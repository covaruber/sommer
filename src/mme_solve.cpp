#include "mme_sp2_support.h"
#include "mme_direct_covariance.h"

// ===========================================================================
// Known-covariance mixed-model equations (mmes(solveOnly=TRUE)).
//
// Iteration-on-data PCG for
//   C = W'R^{-1}W + blockdiag_k(Lambda_k (x) A_k^{-1}),  W = [X Z_1 ... Z_m],
// where C is never assembled. R is block diagonal (one block per residual
// block), with block precisions shared by blocks with identical coordinates.
// ===========================================================================

namespace {

struct SommerCsc {
  int nrow = 0;
  int ncol = 0;
  Rcpp::IntegerVector p;
  Rcpp::IntegerVector i;
  Rcpp::NumericVector x;

  explicit SommerCsc(const Rcpp::S4 & m){
    if(!m.is("dgCMatrix")){
      Rcpp::stop("mme_pcg_solve() expects dgCMatrix inputs.");
    }
    Rcpp::IntegerVector dim = m.slot("Dim");
    nrow = dim[0];
    ncol = dim[1];
    p = m.slot("p");
    i = m.slot("i");
    x = m.slot("x");
  }
};

struct SommerPrecEntry {
  int group;
  int pos;
  int local;
  int row;
  double value;
};

}

// [[Rcpp::export]]
Rcpp::List mme_pcg_solve(const Rcpp::S4 & design,
                         const Rcpp::NumericVector & y,
                         const int nFixed,
                         const Rcpp::IntegerVector & termStart,
                         const Rcpp::IntegerVector & termLevels,
                         const Rcpp::List & AiI,
                         const Rcpp::List & lambdaI,
                         const Rcpp::List & residualCovStruct,
                         const arma::vec & residualPar,
                         const Rcpp::IntegerVector & residualBlock,
                         const Rcpp::IntegerVector & residualIndex,
                         const Rcpp::NumericVector & weights,
                         const double tol,
                         const int maxIter,
                         const int denseFixedMax,
                         const int maxResidualBlock,
                         const Rcpp::NumericVector & start,
                         const bool verbose){

  typedef std::chrono::steady_clock Clock;
  const Clock::time_point setupBegin = Clock::now();
  int residualRefresh = 0;
  if(residualCovStruct.hasAttribute("pcgResidualRefresh")){
    residualRefresh = Rcpp::as<int>(residualCovStruct.attr("pcgResidualRefresh"));
    if(residualRefresh < 0){
      Rcpp::stop("PCG residual refresh interval must be non-negative.");
    }
  }

  const SommerCsc D(design);
  const int n = D.nrow;
  const int neq = D.ncol;
  const int * Dp = D.p.begin();
  const int * Di = D.i.begin();
  const double * Dx = D.x.begin();

  if(y.size() != n || residualBlock.size() != n || residualIndex.size() != n){
    Rcpp::stop("Response, residual blocks and residual coordinates must have one entry per observation.");
  }
  if(weights.size() != 0 && weights.size() != n){
    Rcpp::stop("Weights must be empty or have one entry per observation.");
  }
  if(nFixed < 0 || nFixed > neq){
    Rcpp::stop("Invalid number of fixed-effect columns.");
  }
  if(!std::isfinite(tol) || tol <= 0.0){
    Rcpp::stop("The PCG tolerance must be positive and finite.");
  }
  if(start.size() != 0 && start.size() != neq){
    Rcpp::stop("Starting values must be empty or have one entry per equation.");
  }
  for(int r = 0; r < n; ++r){
    if(!std::isfinite(y[r])){
      Rcpp::stop("The response contains non-finite values.");
    }
  }

  // ---- Random terms ------------------------------------------------------
  const int nTerms = termStart.size();
  if(termLevels.size() != nTerms || AiI.size() != nTerms || lambdaI.size() != nTerms){
    Rcpp::stop("Random-term descriptions have inconsistent lengths.");
  }
  std::vector<SommerCsc> Ai;
  Ai.reserve(static_cast<std::size_t>(nTerms));
  std::vector<arma::mat> Lambda(static_cast<std::size_t>(nTerms));
  for(int k = 0; k < nTerms; ++k){
    Ai.emplace_back(Rcpp::as<Rcpp::S4>(AiI[k]));
    Lambda[static_cast<std::size_t>(k)] = Rcpp::as<arma::mat>(lambdaI[k]);
    const arma::mat & L = Lambda[static_cast<std::size_t>(k)];
    const int nl = termLevels[k];
    if(L.n_rows != L.n_cols || L.n_rows < 1 || !L.is_finite()){
      Rcpp::stop("Random-term precision matrices must be finite and square.");
    }
    if(Ai.back().nrow != nl || Ai.back().ncol != nl){
      Rcpp::stop("A relationship precision matrix does not match its number of levels.");
    }
    if(termStart[k] < nFixed ||
       static_cast<long long>(termStart[k]) + static_cast<long long>(L.n_rows) * nl > neq){
      Rcpp::stop("Random-term coefficient ranges exceed the design matrix.");
    }
  }

  // ---- Row-major copy of the design for W x --------------------------------
  const std::size_t nnzD = static_cast<std::size_t>(Dp[neq]);
  std::vector<std::size_t> rowPtr(static_cast<std::size_t>(n) + 1, 0);
  for(std::size_t e = 0; e < nnzD; ++e){
    ++rowPtr[static_cast<std::size_t>(Di[e]) + 1];
  }
  for(int r = 0; r < n; ++r){
    rowPtr[static_cast<std::size_t>(r) + 1] += rowPtr[static_cast<std::size_t>(r)];
  }
  std::vector<int> rowCol(nnzD);
  std::vector<double> rowVal(nnzD);
  {
    std::vector<std::size_t> fill(rowPtr.begin(), rowPtr.end() - 1);
    for(int c = 0; c < neq; ++c){
      for(int e = Dp[c]; e < Dp[c+1]; ++e){
        const std::size_t pos = fill[static_cast<std::size_t>(Di[e])]++;
        rowCol[pos] = c;
        rowVal[pos] = Dx[e];
      }
    }
  }

  // ---- Residual covariance blocks ------------------------------------------
  if(residualPar.n_elem < 1 || !residualPar.is_finite()){
    Rcpp::stop("Residual covariance parameters must be finite and include log_sigma2.");
  }
  Rcpp::List rfactors = residualCovStruct["factors"];
  const int nRF = rfactors.size();
  const int residualDim = Rcpp::as<int>(residualCovStruct["dim"]);
  std::vector<std::vector<int>> digit(static_cast<std::size_t>(nRF),
                                      std::vector<int>(static_cast<std::size_t>(residualDim)));
  std::vector<int> owner(static_cast<std::size_t>(nRF), -1);
  std::vector<std::vector<arma::mat>> pieces(static_cast<std::size_t>(nRF));
  {
    long long trailing = residualDim;
    for(int f = 0; f < nRF; ++f){
      Rcpp::List fl = Rcpp::as<Rcpp::List>(rfactors[f]);
      const int dim = Rcpp::as<int>(fl["dim"]);
      trailing /= dim;
      for(int local = 0; local < residualDim; ++local){
        digit[static_cast<std::size_t>(f)][static_cast<std::size_t>(local)] =
          static_cast<int>((local / trailing) % dim);
      }
      if(fl.containsElementNamed("section_factor")){
        owner[static_cast<std::size_t>(f)] = Rcpp::as<int>(fl["section_factor"]) - 1;
      }
      const int start1 = fl.containsElementNamed("par_start") ? Rcpp::as<int>(fl["par_start"]) : 1;
      const int end1 = fl.containsElementNamed("par_end") ? Rcpp::as<int>(fl["par_end"]) : 0;
      arma::vec localPar;
      if(end1 >= start1){
        localPar = residualPar.subvec(static_cast<arma::uword>(start1-1),
                                      static_cast<arma::uword>(end1-1));
      }
      const int nSections = fl.containsElementNamed("section_levels") ? Rf_length(fl["section_levels"]) : 0;
      std::vector<arma::mat> & slot = pieces[static_cast<std::size_t>(f)];
      if(nSections == 0){
        slot.assign(1, directEvalFactor(fl, localPar));
      }else{
        if(owner[static_cast<std::size_t>(f)] < 0){
          Rcpp::stop("A sectioned residual factor has no section owner.");
        }
        slot.resize(static_cast<std::size_t>(nSections));
        for(int t = 0; t < nSections; ++t){
          const int listed = directSectionListed(fl, t);
          slot[static_cast<std::size_t>(t)] = listed < 0
            ? arma::eye<arma::mat>(dim, dim)
            : directEvalFactor(fl, directSectionSlice(fl, localPar, listed));
        }
      }
    }
  }
  const double residualScale = std::exp(residualPar(0));

  int nBlocks = 0;
  for(int r = 0; r < n; ++r){
    if(residualBlock[r] < 1 || residualIndex[r] < 1 || residualIndex[r] > residualDim){
      Rcpp::stop("Residual block/local indices are out of range.");
    }
    nBlocks = std::max(nBlocks, residualBlock[r]);
  }
  std::vector<std::size_t> blockPtr(static_cast<std::size_t>(nBlocks) + 1, 0);
  for(int r = 0; r < n; ++r){
    ++blockPtr[static_cast<std::size_t>(residualBlock[r])];
  }
  for(int b = 0; b < nBlocks; ++b){
    blockPtr[static_cast<std::size_t>(b) + 1] += blockPtr[static_cast<std::size_t>(b)];
  }
  std::vector<int> blockRows(static_cast<std::size_t>(n));
  {
    std::vector<std::size_t> fill(blockPtr.begin(), blockPtr.end() - 1);
    for(int r = 0; r < n; ++r){
      blockRows[fill[static_cast<std::size_t>(residualBlock[r] - 1)]++] = r;
    }
  }

  std::map<std::vector<int>, int> patternOf;
  std::vector<arma::mat> patternInv;
  std::vector<int> blockPattern(static_cast<std::size_t>(nBlocks), -1);
  int largestBlock = 0;
  std::vector<int> key;
  for(int b = 0; b < nBlocks; ++b){
    const std::size_t first = blockPtr[static_cast<std::size_t>(b)];
    const int m = static_cast<int>(blockPtr[static_cast<std::size_t>(b) + 1] - first);
    if(m == 0){
      continue;
    }
    if(m > maxResidualBlock){
      Rcpp::stop("A residual block has " + std::to_string(m) +
                 " records, above the solveOnly limit of " + std::to_string(maxResidualBlock) +
                 "; the known-variance solver needs block-diagonal residuals with moderate blocks.");
    }
    largestBlock = std::max(largestBlock, m);
    key.resize(static_cast<std::size_t>(m));
    for(int a = 0; a < m; ++a){
      key[static_cast<std::size_t>(a)] = residualIndex[blockRows[first + static_cast<std::size_t>(a)]] - 1;
    }
    std::map<std::vector<int>, int>::iterator found = patternOf.find(key);
    if(found != patternOf.end()){
      blockPattern[static_cast<std::size_t>(b)] = found->second;
      continue;
    }
    std::vector<const arma::mat *> piece(static_cast<std::size_t>(nRF));
    for(int f = 0; f < nRF; ++f){
      const std::size_t ff = static_cast<std::size_t>(f);
      std::size_t section = 0;
      if(pieces[ff].size() > 1 || owner[ff] >= 0){
        section = static_cast<std::size_t>(
          digit[static_cast<std::size_t>(owner[ff])][static_cast<std::size_t>(key[0])]);
      }
      piece[ff] = &pieces[ff][section];
    }
    arma::mat R0(static_cast<arma::uword>(m), static_cast<arma::uword>(m));
    for(int a1 = 0; a1 < m; ++a1){
      for(int a2 = 0; a2 < m; ++a2){
        double value = residualScale;
        for(int f = 0; f < nRF; ++f){
          const std::size_t ff = static_cast<std::size_t>(f);
          value *= (*piece[ff])(digit[ff][static_cast<std::size_t>(key[static_cast<std::size_t>(a1)])],
                                digit[ff][static_cast<std::size_t>(key[static_cast<std::size_t>(a2)])]);
        }
        R0(a1, a2) = value;
      }
    }
    arma::mat inverse;
    if(!arma::inv_sympd(inverse, 0.5 * (R0 + R0.t()))){
      Rcpp::stop("A residual covariance block is not positive definite at the supplied parameters.");
    }
    const int id = static_cast<int>(patternInv.size());
    patternInv.push_back(inverse);
    patternOf.emplace(key, id);
    blockPattern[static_cast<std::size_t>(b)] = id;
  }

  std::vector<double> sw(static_cast<std::size_t>(n), 1.0);
  if(weights.size() == n){
    for(int r = 0; r < n; ++r){
      if(!std::isfinite(weights[r]) || weights[r] <= 0.0){
        Rcpp::stop("Weights must be positive and finite.");
      }
      sw[static_cast<std::size_t>(r)] = std::sqrt(weights[r]);
    }
  }

  // s = R^{-1} t, with R^{-1} = S P S blockwise (S = diag(sqrt(weights))).
  auto applyRinv = [&](const double * t, double * s){
    #pragma omp parallel
    {
      arma::vec tb;
      #pragma omp for schedule(dynamic, 512)
      for(int b = 0; b < nBlocks; ++b){
        const std::size_t first = blockPtr[static_cast<std::size_t>(b)];
        const int m = static_cast<int>(blockPtr[static_cast<std::size_t>(b) + 1] - first);
        if(m == 0){
          continue;
        }
        const arma::mat & P = patternInv[static_cast<std::size_t>(blockPattern[static_cast<std::size_t>(b)])];
        if(m == 1){
          const std::size_t r = static_cast<std::size_t>(blockRows[first]);
          s[r] = sw[r] * P(0,0) * sw[r] * t[r];
          continue;
        }
        tb.set_size(static_cast<arma::uword>(m));
        for(int a = 0; a < m; ++a){
          const std::size_t r = static_cast<std::size_t>(blockRows[first + static_cast<std::size_t>(a)]);
          tb(static_cast<arma::uword>(a)) = sw[r] * t[r];
        }
        const arma::vec sb = P * tb;
        for(int a = 0; a < m; ++a){
          const std::size_t r = static_cast<std::size_t>(blockRows[first + static_cast<std::size_t>(a)]);
          s[r] = sw[r] * sb(static_cast<arma::uword>(a));
        }
      }
    }
  };

  auto applyDesign = [&](const double * v, double * t){
    #pragma omp parallel for schedule(static)
    for(int r = 0; r < n; ++r){
      double acc = 0.0;
      for(std::size_t e = rowPtr[static_cast<std::size_t>(r)]; e < rowPtr[static_cast<std::size_t>(r) + 1]; ++e){
        acc += rowVal[e] * v[rowCol[e]];
      }
      t[r] = acc;
    }
  };

  auto applyDesignT = [&](const double * s, double * out){
    #pragma omp parallel for schedule(dynamic, 1024)
    for(int c = 0; c < neq; ++c){
      double acc = 0.0;
      for(int e = Dp[c]; e < Dp[c+1]; ++e){
        acc += Dx[e] * s[Di[e]];
      }
      out[c] = acc;
    }
  };

  // out += vec(A_k U Lambda_k) for every random term, U = coefficients as levels x q.
  auto addPriors = [&](const double * v, double * out){
    for(int k = 0; k < nTerms; ++k){
      const std::size_t kk = static_cast<std::size_t>(k);
      const std::size_t off = static_cast<std::size_t>(termStart[k]);
      const int nl = termLevels[k];
      const int q = static_cast<int>(Lambda[kk].n_rows);
      const double * L = Lambda[kk].memptr();
      const int * Ap = Ai[kk].p.begin();
      const int * Aidx = Ai[kk].i.begin();
      const double * Ax = Ai[kk].x.begin();
      #pragma omp parallel
      {
        std::vector<double> au(static_cast<std::size_t>(q));
        #pragma omp for schedule(dynamic, 1024)
        for(int l = 0; l < nl; ++l){
          std::fill(au.begin(), au.end(), 0.0);
          for(int e = Ap[l]; e < Ap[l+1]; ++e){
            const std::size_t row = static_cast<std::size_t>(Aidx[e]);
            const double value = Ax[e];
            for(int a = 0; a < q; ++a){
              au[static_cast<std::size_t>(a)] +=
                value * v[off + static_cast<std::size_t>(a) * static_cast<std::size_t>(nl) + row];
            }
          }
          for(int a = 0; a < q; ++a){
            double acc = 0.0;
            for(int b2 = 0; b2 < q; ++b2){
              acc += au[static_cast<std::size_t>(b2)] * L[b2 + a * q];
            }
            out[off + static_cast<std::size_t>(a) * static_cast<std::size_t>(nl) + static_cast<std::size_t>(l)] += acc;
          }
        }
      }
    }
  };

  std::vector<double> tbuf(static_cast<std::size_t>(n));
  std::vector<double> sbuf(static_cast<std::size_t>(n));
  std::uint64_t operatorApplications = 0;
  auto applyC = [&](const std::vector<double> & v, std::vector<double> & out){
    ++operatorApplications;
    applyDesign(v.data(), tbuf.data());
    applyRinv(tbuf.data(), sbuf.data());
    applyDesignT(sbuf.data(), out.data());
    addPriors(v.data(), out.data());
  };

  // ---- Block-Jacobi preconditioner -----------------------------------------
  // Groups: all fixed effects (dense) or one per fixed column, and one q x q
  // block per random-term level coupling its covariance coordinates.
  std::vector<int> groupOf(static_cast<std::size_t>(neq), -1);
  std::vector<int> posOf(static_cast<std::size_t>(neq), 0);
  std::vector<int> groupSize;
  std::vector<std::size_t> groupColPtr(1, 0);
  std::vector<int> groupCols;
  groupCols.reserve(static_cast<std::size_t>(neq));
  auto openGroup = [&](){
    groupSize.push_back(0);
  };
  auto addToGroup = [&](const int col){
    const int g = static_cast<int>(groupSize.size()) - 1;
    groupOf[static_cast<std::size_t>(col)] = g;
    posOf[static_cast<std::size_t>(col)] = groupSize.back()++;
    groupCols.push_back(col);
  };
  auto closeGroup = [&](){
    groupColPtr.push_back(groupCols.size());
  };
  if(nFixed > 0 && nFixed <= denseFixedMax){
    openGroup();
    for(int c = 0; c < nFixed; ++c){ addToGroup(c); }
    closeGroup();
  }
  for(int k = 0; k < nTerms; ++k){
    const int nl = termLevels[k];
    const int q = static_cast<int>(Lambda[static_cast<std::size_t>(k)].n_rows);
    for(int l = 0; l < nl; ++l){
      openGroup();
      for(int a = 0; a < q; ++a){
        addToGroup(termStart[k] + a * nl + l);
      }
      closeGroup();
    }
  }
  for(int c = 0; c < neq; ++c){
    if(groupOf[static_cast<std::size_t>(c)] < 0){
      openGroup();
      addToGroup(c);
      closeGroup();
    }
  }
  const int nGroups = static_cast<int>(groupSize.size());
  std::vector<std::size_t> groupStore(static_cast<std::size_t>(nGroups) + 1, 0);
  for(int g = 0; g < nGroups; ++g){
    const std::size_t s = static_cast<std::size_t>(groupSize[static_cast<std::size_t>(g)]);
    groupStore[static_cast<std::size_t>(g) + 1] = groupStore[static_cast<std::size_t>(g)] + s * s;
  }
  std::vector<double> store(groupStore.back(), 0.0);

  {
    std::vector<SommerPrecEntry> entries;
    for(int b = 0; b < nBlocks; ++b){
      const std::size_t first = blockPtr[static_cast<std::size_t>(b)];
      const int m = static_cast<int>(blockPtr[static_cast<std::size_t>(b) + 1] - first);
      if(m == 0){
        continue;
      }
      const arma::mat & P = patternInv[static_cast<std::size_t>(blockPattern[static_cast<std::size_t>(b)])];
      entries.clear();
      for(int a = 0; a < m; ++a){
        const int r = blockRows[first + static_cast<std::size_t>(a)];
        for(std::size_t e = rowPtr[static_cast<std::size_t>(r)]; e < rowPtr[static_cast<std::size_t>(r) + 1]; ++e){
          const int c = rowCol[e];
          entries.push_back({groupOf[static_cast<std::size_t>(c)], posOf[static_cast<std::size_t>(c)], a, r, rowVal[e]});
        }
      }
      std::sort(entries.begin(), entries.end(),
                [](const SommerPrecEntry & lhs, const SommerPrecEntry & rhs){ return lhs.group < rhs.group; });
      std::size_t runStart = 0;
      while(runStart < entries.size()){
        std::size_t runEnd = runStart + 1;
        while(runEnd < entries.size() && entries[runEnd].group == entries[runStart].group){ ++runEnd; }
        const int g = entries[runStart].group;
        const std::size_t s = static_cast<std::size_t>(groupSize[static_cast<std::size_t>(g)]);
        double * B = &store[groupStore[static_cast<std::size_t>(g)]];
        for(std::size_t e1 = runStart; e1 < runEnd; ++e1){
          const SommerPrecEntry & u1 = entries[e1];
          for(std::size_t e2 = runStart; e2 < runEnd; ++e2){
            const SommerPrecEntry & u2 = entries[e2];
            const double rinv = sw[static_cast<std::size_t>(u1.row)] *
              P(static_cast<arma::uword>(u1.local), static_cast<arma::uword>(u2.local)) *
              sw[static_cast<std::size_t>(u2.row)];
            B[static_cast<std::size_t>(u1.pos) + static_cast<std::size_t>(u2.pos) * s] +=
              u1.value * rinv * u2.value;
          }
        }
        runStart = runEnd;
      }
    }
  }

  for(int k = 0; k < nTerms; ++k){
    const std::size_t kk = static_cast<std::size_t>(k);
    const int nl = termLevels[k];
    const int q = static_cast<int>(Lambda[kk].n_rows);
    const int * Ap = Ai[kk].p.begin();
    const int * Aidx = Ai[kk].i.begin();
    const double * Ax = Ai[kk].x.begin();
    for(int l = 0; l < nl; ++l){
      double aii = 0.0;
      for(int e = Ap[l]; e < Ap[l+1]; ++e){
        if(Aidx[e] == l){ aii += Ax[e]; }
      }
      const int g = groupOf[static_cast<std::size_t>(termStart[k] + l)];
      double * B = &store[groupStore[static_cast<std::size_t>(g)]];
      for(int a = 0; a < q; ++a){
        for(int b2 = 0; b2 < q; ++b2){
          B[a + b2 * q] += Lambda[kk](static_cast<arma::uword>(a), static_cast<arma::uword>(b2)) * aii;
        }
      }
    }
  }

  int zeroPivots = 0;
  #pragma omp parallel for schedule(dynamic, 1024) reduction(+:zeroPivots)
  for(int g = 0; g < nGroups; ++g){
    const arma::uword s = static_cast<arma::uword>(groupSize[static_cast<std::size_t>(g)]);
    double * B = &store[groupStore[static_cast<std::size_t>(g)]];
    if(s == 1){
      if(B[0] > 0.0 && std::isfinite(B[0])){
        B[0] = 1.0 / B[0];
      }else{
        B[0] = 1.0;
        ++zeroPivots;
      }
      continue;
    }
    arma::mat M(B, s, s, false, true);
    arma::mat inverse;
    if(!arma::inv_sympd(inverse, 0.5 * (M + M.t()))){
      inverse = arma::pinv(0.5 * (M + M.t()));
      ++zeroPivots;
    }
    M = inverse;
  }

  auto applyPrec = [&](const std::vector<double> & r, std::vector<double> & z){
    #pragma omp parallel for schedule(dynamic, 1024)
    for(int g = 0; g < nGroups; ++g){
      const std::size_t s = static_cast<std::size_t>(groupSize[static_cast<std::size_t>(g)]);
      const int * cols = &groupCols[groupColPtr[static_cast<std::size_t>(g)]];
      const double * B = &store[groupStore[static_cast<std::size_t>(g)]];
      for(std::size_t a = 0; a < s; ++a){
        double acc = 0.0;
        for(std::size_t b2 = 0; b2 < s; ++b2){
          acc += B[a + b2 * s] * r[static_cast<std::size_t>(cols[b2])];
        }
        z[static_cast<std::size_t>(cols[a])] = acc;
      }
    }
  };

  auto dot = [&](const std::vector<double> & a, const std::vector<double> & b){
    double sum = 0.0;
    #pragma omp parallel for schedule(static) reduction(+:sum)
    for(int j = 0; j < neq; ++j){
      sum += a[static_cast<std::size_t>(j)] * b[static_cast<std::size_t>(j)];
    }
    return sum;
  };

  // ---- PCG -----------------------------------------------------------------
  const std::size_t N = static_cast<std::size_t>(neq);
  std::vector<double> rhs(N);
  {
    std::vector<double> yv(y.begin(), y.end());
    applyRinv(yv.data(), sbuf.data());
    applyDesignT(sbuf.data(), rhs.data());
  }
  const double rhsNorm = std::sqrt(dot(rhs, rhs));

  std::vector<double> x(N, 0.0);
  if(start.size() == neq){
    std::copy(start.begin(), start.end(), x.begin());
  }
  if(rhsNorm == 0.0){ std::fill(x.begin(), x.end(), 0.0); }
  std::vector<double> r(N), z(N), pdir(N), q(N);
  applyC(x, q);
  for(std::size_t j = 0; j < N; ++j){ r[j] = rhs[j] - q[j]; }

  const int iterLimit = maxIter > 0 ? maxIter : std::max(100, std::min(2 * neq, 10000));
  std::vector<double> history;
  history.reserve(static_cast<std::size_t>(std::min(iterLimit, 10000)) + 1);
  double relres = rhsNorm > 0.0 ? std::sqrt(dot(r, r)) / rhsNorm : 0.0;
  history.push_back(relres);
  bool converged = rhsNorm == 0.0 || relres < tol;
  int iterations = 0;
  int residualReplacements = 0;

  const double setupSeconds = std::chrono::duration<double>(Clock::now() - setupBegin).count();
  const Clock::time_point solveBegin = Clock::now();

  if(verbose){
    Rcpp::Rcout << "MME solve-only path: iteration-on-data PCG (block-Jacobi)" << arma::endl;
  }

  if(!converged){
    applyPrec(r, z);
    pdir = z;
    double rz = dot(r, z);
    if(!(rz > 0.0) || !std::isfinite(rz)){
      Rcpp::stop("PCG breakdown: non-positive preconditioned residual.");
    }
    while(iterations < iterLimit){
      ++iterations;
      applyC(pdir, q);
      const double pq = dot(pdir, q);
      if(!(pq > 0.0) || !std::isfinite(pq)){
        Rcpp::stop("PCG breakdown: the coefficient matrix is not positive definite at the supplied parameters.");
      }
      const double alpha = rz / pq;
      #pragma omp parallel for schedule(static)
      for(int j = 0; j < neq; ++j){
        x[static_cast<std::size_t>(j)] += alpha * pdir[static_cast<std::size_t>(j)];
        r[static_cast<std::size_t>(j)] -= alpha * q[static_cast<std::size_t>(j)];
      }
      relres = std::sqrt(dot(r, r)) / rhsNorm;
      const bool refresh = relres < tol ||
        (residualRefresh > 0 && iterations % residualRefresh == 0);
      if(refresh){
        applyC(x, q);
        for(std::size_t j = 0; j < N; ++j){ r[j] = rhs[j] - q[j]; }
        relres = std::sqrt(dot(r, r)) / rhsNorm;
        ++residualReplacements;
      }
      history.push_back(relres);
      if(verbose && iterations % 50 == 0){
        Rcpp::Rcout << "  PCG iteration " << iterations << ": relative residual " << relres << std::endl;
      }
      if(relres < tol){
        converged = true;
        break;
      }
      if(iterations % 25 == 0){
        Rcpp::checkUserInterrupt();
      }
      applyPrec(r, z);
      const double rzNew = dot(r, z);
      if(!(rzNew > 0.0) || !std::isfinite(rzNew)){
        Rcpp::stop("PCG breakdown: non-positive preconditioned residual update.");
      }
      const double beta = refresh ? 0.0 : rzNew / rz;
      rz = rzNew;
      #pragma omp parallel for schedule(static)
      for(int j = 0; j < neq; ++j){
        pdir[static_cast<std::size_t>(j)] = z[static_cast<std::size_t>(j)] + beta * pdir[static_cast<std::size_t>(j)];
      }
    }
  }
  applyC(x, q);
  for(std::size_t j = 0; j < N; ++j){ r[j] = rhs[j] - q[j]; }
  relres = rhsNorm > 0.0 ? std::sqrt(dot(r, r)) / rhsNorm : 0.0;
  if(!std::isfinite(relres)){
    Rcpp::stop("PCG produced a non-finite true residual.");
  }
  converged = relres <= tol;
  if(!history.empty()){ history.back() = relres; }
  const double solveSeconds = std::chrono::duration<double>(Clock::now() - solveBegin).count();

  std::vector<double> fitted(static_cast<std::size_t>(n));
  applyDesign(x.data(), fitted.data());

  if(verbose){
    Rcpp::Rcout << "PCG " << (converged ? "converged" : "did NOT converge")
                << " in " << iterations << " iterations (relative residual " << relres
                << "; " << neq << " equations, " << patternInv.size()
                << " residual block patterns)." << std::endl;
  }

  return Rcpp::List::create(
    Rcpp::Named("solution") = Rcpp::wrap(x),
    Rcpp::Named("fitted") = Rcpp::wrap(fitted),
    Rcpp::Named("iterations") = iterations,
    Rcpp::Named("relres") = relres,
    Rcpp::Named("trueRelativeResidual") = relres,
    Rcpp::Named("operatorApplications") = static_cast<double>(operatorApplications),
    Rcpp::Named("residualReplacements") = residualReplacements,
    Rcpp::Named("warmStarted") = start.size() == neq,
    Rcpp::Named("converged") = converged,
    Rcpp::Named("history") = Rcpp::wrap(history),
    Rcpp::Named("rhsNorm") = rhsNorm,
    Rcpp::Named("residualBlocks") = nBlocks,
    Rcpp::Named("residualPatterns") = static_cast<int>(patternInv.size()),
    Rcpp::Named("largestResidualBlock") = largestBlock,
    Rcpp::Named("preconditionerGroups") = nGroups,
    Rcpp::Named("singularPreconditionerBlocks") = zeroPivots,
    Rcpp::Named("setupSeconds") = setupSeconds,
    Rcpp::Named("solveSeconds") = solveSeconds
  );
}
