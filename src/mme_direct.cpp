#include "mme_sp2_support.h"
#include "mme_progress.h"

// ====================================================================
// Direct-inversion (observation-space) REML/ML engine.
//
// ai_mme_sp2() solves the Henderson mixed-model equations in
// coefficient space: efficient when there are many more records than
// coefficients (e.g. pedigree BLUP). ai_reml_direct_sp2() instead
// inverts the n x n phenotypic covariance V directly, which is
// efficient when there are many more coefficients than records (e.g.
// marker/SNP-BLUP models with p >> n). It consumes the SAME covStruct
// CovarianceFactor-v2 descriptor contract produced by vsm(), so the
// same covariance-model parameterizations (identity/diag/us/ar1/fa/
// rr/custom ownm(), etc.) are available in both engines.
//
// Derivatives of each small covariance-shaping factor are obtained by
// central finite differences instead of duplicating every native
// analytic derivative formula from ai_mme_sp2(): factors are always
// small (their dimension is the number of levels/traits/lags in a
// Kronecker term, never the number of records or marker coefficients),
// so this is computationally negligible and keeps this engine compact
// and independently auditable.
// ====================================================================

#include "mme_direct_covariance.h"

// [[Rcpp::export]]
arma::mat evaluate_covstruct_cpp(const Rcpp::List & covStruct, const arma::vec & par){
  return directEvaluateDescriptor(covStruct, par);
}

#include "mme_direct_em.h"

// [[Rcpp::export]]
Rcpp::List ai_reml_direct_sp2(const arma::sp_mat & X,
                              const Rcpp::List & ZI,
                              const arma::vec & Zind,
                              const Rcpp::List & AiI,
                              const arma::sp_mat & y0,
                              const arma::sp_mat & H,
                              const bool & useH,
                              const arma::uvec & residualBlockI,
                              const arma::uvec & residualIndexI,
                              int nIters, double tolParConvLL,
                              double tolParConvNorm, double tolParInv,
                              const Rcpp::List & covStructI,
                              const arma::vec & weightEmInf,
                              const arma::vec & weightInf,
                              const bool & verbose,
                              const int & computePev = 0,
                              const bool & reml = true,
                              const bool & responsePrepared = false,
                              const double & preparedMean = 0.0,
                              const double & preparedSd = 1.0,
                              const bool & preparedIntercept = false
){

  if(computePev != 0 && computePev != 2){
    Rcpp::stop("computePev must be 0 (no PEV) or 2 (full PEV, small models only).");
  }

  const int nRRe = covStructI.size();
  if(nRRe < 1){
    Rcpp::stop("At least one covariance descriptor (the residual structure) is required.");
  }
  const int nRe = nRRe - 1;
  const int residualStruct = nRe;
  const int nX = X.n_cols;
  const int nR = y0.n_rows;

  if(y0.n_cols != 1){
    Rcpp::stop("ai_reml_direct_sp2() currently supports a single response column; use the long-format vsm(usm(trait), ...) convention for multi-trait models, exactly as with the Henderson engine.");
  }
  if(residualBlockI.n_elem != static_cast<arma::uword>(nR) ||
     residualIndexI.n_elem != static_cast<arma::uword>(nR)){
    Rcpp::stop("Residual block/index vectors must have one entry per observation.");
  }

  double vary = responsePrepared ? preparedSd * preparedSd : arma::var(arma::vec(arma::mat(y0).col(0)));
  if(!std::isfinite(vary) || vary <= 0.0){
    Rcpp::stop("Response variance must be positive and finite.");
  }
  const double stdy = responsePrepared ? preparedSd : std::sqrt(vary);

  arma::vec cShift(X.n_cols, arma::fill::zeros);
  bool hasIntercept = false;
  if(responsePrepared){
    hasIntercept = preparedIntercept;
    if(hasIntercept && X.n_cols > 0){
      cShift(0) = 1.0;
    }
  } else {
    hasIntercept = find_intercept_vector(X, cShift);
  }

  double muy = 0.0;
  if(responsePrepared){
    muy = preparedMean;
  } else if(hasIntercept){
    muy = arma::mean(arma::vec(arma::mat(y0).col(0)));
  }

  arma::vec y = arma::vec(arma::mat(y0).col(0));
  if(!responsePrepared){
    if(hasIntercept){
      y = (y - muy) / stdy;
    } else {
      y = y / stdy;
    }
  }

  bool intercept = hasIntercept;
  const arma::mat Xd = arma::mat(X);
  arma::uword rankX = arma::rank(Xd);

  // ------------------------------------------------------------
  // Parse covariance descriptors (identical contract to ai_mme_sp2).
  // ------------------------------------------------------------
  std::vector<Rcpp::List> covDescriptor(nRRe);
  arma::field<arma::vec> covPar(nRRe);
  arma::field<arma::vec> covFree(nRRe);
  arma::field<arma::vec> covLower(nRRe);
  arma::field<arma::vec> covUpper(nRRe);
  arma::field<arma::mat> theta(nRRe);

  for(int i = 0; i < nRRe; ++i){
    if(Rf_isNull(covStructI[i])){
      Rcpp::stop("NULL covariance descriptor supplied to ai_reml_direct_sp2().");
    }
    Rcpp::List cs = Rcpp::as<Rcpp::List>(covStructI[i]);
    if(!cs.containsElementNamed("type") || Rcpp::as<std::string>(cs["type"]) != "kron"){
      Rcpp::stop("ai_reml_direct_sp2() accepts only the generic type='kron' covariance descriptor.");
    }
    if(!cs.containsElementNamed("descriptor_version") || Rcpp::as<int>(cs["descriptor_version"]) < 2){
      Rcpp::stop("ai_reml_direct_sp2() requires CovarianceFactor descriptor_version >= 2.");
    }
    covDescriptor[i] = cs;
    if(i == residualStruct){
      const Rcpp::List factors = cs["factors"];
      for(int axis = 0; axis < factors.size(); ++axis){
        if(sommer::zeroSpecificFactor(Rcpp::as<Rcpp::List>(factors[axis]))){
          Rcpp::stop("Zero-specific reduced rank is supported for random effects only; residual covariance must be positive definite.");
        }
      }
    }
    covPar(i) = Rcpp::as<arma::vec>(cs["par"]);
    Rcpp::LogicalVector freeR = cs["free"];
    if(freeR.size() != static_cast<int>(covPar(i).n_elem)){
      Rcpp::stop("covStruct$free and covStruct$par have inconsistent lengths.");
    }
    covFree(i).set_size(covPar(i).n_elem);
    for(arma::uword k = 0; k < covPar(i).n_elem; ++k){
      covFree(i)(k) = freeR[static_cast<int>(k)] ? 1.0 : 0.0;
    }
    covLower(i) = sommer::covarianceVarianceLowerBounds(
      Rcpp::List::create(cs), tolParInv, std::log(vary));
    covUpper(i) = arma::vec(covPar(i).n_elem, arma::fill::value(std::numeric_limits<double>::infinity()));

    covPar(i)(0) -= std::log(vary);
    if(covFree(i)(0) > 0.5) {
      covLower(i)(0) = std::log(std::max(1.0e-8, tolParInv)) - std::log(vary);
    }

    theta(i) = (i == residualStruct && directResidualCompact(cs))
      ? arma::mat()
      : directEvaluateDescriptor(cs, covPar(i));
  }

  arma::vec nVc(nRRe);
  for(int i = 0; i < nRRe; ++i){ nVc(i) = static_cast<double>(covPar(i).n_elem); }
  const int nVcTotal = static_cast<int>(arma::accu(nVc));
  arma::vec nVcEnd = nVc;
  for(int i = 0; i < nRRe; ++i){
    arma::uvec toSum = arma::regspace<arma::uvec>(0, 1, i);
    nVcEnd(i) = arma::accu(nVc(toSum));
  }
  arma::vec nVcStart = nVcEnd - nVc + 1;

  // ------------------------------------------------------------
  // Build and cache Z_a K Z_b' blocks per random structure, once.
  // Constant across REML iterations: only the small Sigma_i(a,b)
  // scalars change, never these n x n record-space blocks.
  // ------------------------------------------------------------
  std::vector<arma::field<arma::mat>> Bcache(nRe);
  std::vector<std::vector<arma::mat>> Zdense(nRe);
  std::vector<arma::uword> qOf(nRe);
  std::vector<arma::uvec> levelStart(nRe), levelEnd(nRe); // per-block coefficient offsets within u

  int lastOffset = nX;
  arma::field<arma::mat> partitions(nRe);

  for(int i = 0; i < nRe; ++i){
    arma::uvec zidx = arma::find(Zind == (i+1));
    const arma::uword q = zidx.n_elem;
    qOf[i] = q;
    if(q != theta(i).n_rows){
      Rcpp::stop("Number of Z design blocks does not match the covariance descriptor dimension for a random structure.");
    }

    Zdense[i].resize(q);
    arma::vec starts(q), ends(q);
    for(arma::uword a = 0; a < q; ++a){
      Zdense[i][a] = arma::mat(Rcpp::as<arma::sp_mat>(ZI[static_cast<int>(zidx(a))]));
      const int levels = Zdense[i][a].n_cols;
      starts(a) = lastOffset + 1;
      ends(a) = lastOffset + levels;
      lastOffset += levels;
    }
    partitions(i) = arma::join_rows(starts, ends);

    arma::sp_mat AiSp = convertSparse(AiI(i));
    const bool aiIsIdentity = isIdentity_spmat(AiSp);
    arma::mat Kdense;
    if(!aiIsIdentity){
      arma::mat Adense = arma::mat(AiSp);
      bool ok = eigenSpdInverse(Adense, Kdense);
      if(!ok){
        arma::mat bend = nearPDcpp(Adense, 100, 1.0e-6, 1.0e-7);
        ok = eigenSpdInverse(bend, Kdense);
        if(!ok){
          Rcpp::stop("Unable to invert a random-effect relationship (precision) matrix for the direct-inversion engine.");
        }
      }
    }

    std::vector<arma::mat> ZK(q);
    for(arma::uword a = 0; a < q; ++a){
      ZK[a] = aiIsIdentity ? Zdense[i][a] : (Zdense[i][a] * Kdense);
    }

    Bcache[i].set_size(q, q);
    for(arma::uword a = 0; a < q; ++a){
      for(arma::uword b = 0; b < q; ++b){
        Bcache[i](a,b) = ZK[a] * Zdense[i][b].t();
      }
    }
  }
  const int Nu = lastOffset - nX;
  const int nEffects = nX + Nu;

  // ------------------------------------------------------------
  // Residual block/local-index cache (same generic mapping as
  // ai_mme_sp2, but here it directly assembles record-space R).
  // ------------------------------------------------------------
  const bool residualCompact = directResidualCompact(covDescriptor[residualStruct]);
  const int residualDim = Rcpp::as<int>(covDescriptor[residualStruct]["dim"]);
  if(residualDim < 1){
    Rcpp::stop("Residual covariance dimension must be positive.");
  }
  int nBlocksR = 0;
  for(int rr = 0; rr < nR; ++rr){
    if(residualBlockI(rr) < 1 || residualIndexI(rr) < 1 || static_cast<int>(residualIndexI(rr)) > residualDim){
      Rcpp::stop("Residual block/local indices are out of range.");
    }
    nBlocksR = std::max(nBlocksR, static_cast<int>(residualBlockI(rr)));
  }
  std::vector<std::vector<std::pair<int,int>>> blockRows(nBlocksR);
  for(int rr = 0; rr < nR; ++rr){
    blockRows[residualBlockI(rr)-1].push_back({rr, static_cast<int>(residualIndexI(rr))-1});
  }

  arma::mat HsInv;
  if(useH){
    arma::mat Hd = arma::mat(H);
    arma::mat Hs = arma::chol(Hd, "upper");
    HsInv = arma::inv(arma::trimatu(Hs));
  }

  auto buildR = [&](const arma::mat & thetaR) -> arma::mat {
    arma::mat R0(nR, nR, arma::fill::zeros);
    for(int b = 0; b < nBlocksR; ++b){
      const std::vector<std::pair<int,int>> & rows = blockRows[b];
      for(std::size_t p1 = 0; p1 < rows.size(); ++p1){
        for(std::size_t p2 = 0; p2 < rows.size(); ++p2){
          R0(rows[p1].first, rows[p2].first) = thetaR(rows[p1].second, rows[p2].second);
        }
      }
    }
    if(useH){ return HsInv.t() * R0 * HsInv; }
    return R0;
  };

  Rcpp::List residualFactors = covDescriptor[residualStruct]["factors"];
  const int nResidualFactors = residualFactors.size();
  std::vector<std::vector<arma::uword>> residualDigit(
    static_cast<std::size_t>(nResidualFactors),
    std::vector<arma::uword>(static_cast<std::size_t>(residualDim)));
  std::vector<int> residualOwner(static_cast<std::size_t>(nResidualFactors), -1);
  {
    arma::uword trailing = static_cast<arma::uword>(residualDim);
    for(int f = 0; f < nResidualFactors; ++f){
      Rcpp::List fl = Rcpp::as<Rcpp::List>(residualFactors[f]);
      const arma::uword dim = static_cast<arma::uword>(Rcpp::as<int>(fl["dim"]));
      trailing /= dim;
      for(int local = 0; local < residualDim; ++local){
        residualDigit[static_cast<std::size_t>(f)][static_cast<std::size_t>(local)] =
          (static_cast<arma::uword>(local) / trailing) % dim;
      }
      if(fl.containsElementNamed("section_factor")){
        residualOwner[static_cast<std::size_t>(f)] = Rcpp::as<int>(fl["section_factor"]) - 1;
      }
    }
  }

  // Record-space R (derivativeK < 0) or dR/dpar_k from the small factor matrices.
  auto buildRCompact = [&](const arma::vec & par, const int derivativeK) -> arma::mat {
    std::vector<std::vector<arma::mat>> pieces(static_cast<std::size_t>(nResidualFactors));
    int derivativeFactor = -1;
    int derivativeSection = -1;
    const int k1 = derivativeK + 1;
    for(int f = 0; f < nResidualFactors; ++f){
      Rcpp::List fl = Rcpp::as<Rcpp::List>(residualFactors[f]);
      const int start1 = fl.containsElementNamed("par_start") ? Rcpp::as<int>(fl["par_start"]) : 1;
      const int end1 = fl.containsElementNamed("par_end") ? Rcpp::as<int>(fl["par_end"]) : 0;
      arma::vec localPar;
      if(end1 >= start1){
        localPar = par.subvec(static_cast<arma::uword>(start1-1), static_cast<arma::uword>(end1-1));
      }
      const bool mine = derivativeK > 0 && end1 >= start1 && k1 >= start1 && k1 <= end1;
      const int nSections = fl.containsElementNamed("section_levels") ? Rf_length(fl["section_levels"]) : 0;
      std::vector<arma::mat> & slot = pieces[static_cast<std::size_t>(f)];
      if(nSections == 0){
        slot.assign(1, mine
          ? directFactorD1(fl, localPar, static_cast<arma::uword>(k1-start1))
          : directEvalFactor(fl, localPar));
        if(mine){ derivativeFactor = f; }
        continue;
      }
      const int npar = Rcpp::as<int>(fl["section_npar"]);
      const int q = Rcpp::as<int>(fl["dim"]);
      slot.resize(static_cast<std::size_t>(nSections));
      for(int t = 0; t < nSections; ++t){
        const int listed = directSectionListed(fl, t);
        if(listed < 0){
          slot[static_cast<std::size_t>(t)] = arma::eye<arma::mat>(q, q);
          continue;
        }
        const arma::vec slice = directSectionSlice(fl, localPar, listed);
        if(mine && (k1 - start1) / npar == listed){
          slot[static_cast<std::size_t>(t)] =
            directFactorD1(fl, slice, static_cast<arma::uword>((k1 - start1) % npar));
          derivativeFactor = f;
          derivativeSection = t;
        }else{
          slot[static_cast<std::size_t>(t)] = directEvalFactor(fl, slice);
        }
      }
    }
    if(derivativeK > 0 && derivativeFactor < 0){
      Rcpp::stop("Residual covariance derivative parameter does not belong to any factor.");
    }

    const double scale = std::exp(par(0));
    arma::mat R0(nR, nR, arma::fill::zeros);
    for(int b = 0; b < nBlocksR; ++b){
      const std::vector<std::pair<int,int>> & rows = blockRows[b];
      if(rows.empty()){ continue; }
      const std::size_t first = static_cast<std::size_t>(rows[0].second);
      std::vector<const arma::mat *> blockPiece(static_cast<std::size_t>(nResidualFactors));
      bool zeroBlock = false;
      for(int f = 0; f < nResidualFactors; ++f){
        const std::size_t ff = static_cast<std::size_t>(f);
        std::size_t section = 0;
        if(pieces[ff].size() > 1 || residualOwner[ff] >= 0){
          section = residualDigit[static_cast<std::size_t>(residualOwner[ff])][first];
          if(f == derivativeFactor && derivativeSection >= 0 &&
             static_cast<int>(section) != derivativeSection){
            zeroBlock = true;
          }
        }
        blockPiece[ff] = &pieces[ff][section];
      }
      if(zeroBlock){ continue; }
      for(std::size_t p1 = 0; p1 < rows.size(); ++p1){
        const std::size_t l1 = static_cast<std::size_t>(rows[p1].second);
        for(std::size_t p2 = 0; p2 < rows.size(); ++p2){
          const std::size_t l2 = static_cast<std::size_t>(rows[p2].second);
          double value = scale;
          for(int f = 0; f < nResidualFactors && value != 0.0; ++f){
            const std::size_t ff = static_cast<std::size_t>(f);
            value *= (*blockPiece[ff])(residualDigit[ff][l1], residualDigit[ff][l2]);
          }
          R0(rows[p1].first, rows[p2].first) = value;
        }
      }
    }
    if(useH){ return arma::mat(HsInv.t() * R0 * HsInv); }
    return R0;
  };

  auto buildResidual = [&](const arma::vec & par, const int derivativeK) -> arma::mat {
    if(residualCompact){
      return buildRCompact(par, derivativeK);
    }
    return derivativeK < 0
      ? buildR(directEvaluateDescriptor(covDescriptor[residualStruct], par))
      : buildR(directDescriptorD1(covDescriptor[residualStruct], par,
                                  static_cast<arma::uword>(derivativeK)));
  };

  // ------------------------------------------------------------
  // Main covariance assembly: V(par) and its derivatives.
  // ------------------------------------------------------------
  auto buildV = [&](const arma::field<arma::vec> & par) -> arma::mat {
    arma::mat V = buildResidual(par(residualStruct), -1);
    for(int i = 0; i < nRe; ++i){
      arma::mat Sigma = directEvaluateDescriptor(covDescriptor[i], par(i));
      const arma::uword q = qOf[i];
      for(arma::uword a = 0; a < q; ++a){
        for(arma::uword b = 0; b < q; ++b){
          if(Sigma(a,b) != 0.0){ V += Sigma(a,b) * Bcache[i](a,b); }
        }
      }
    }
    return V;
  };

  auto buildD = [&](const arma::field<arma::vec> & par, const int iStruct, const arma::uword k) -> arma::mat {
    if(iStruct == residualStruct){
      return buildResidual(par(residualStruct), static_cast<int>(k));
    }
    arma::mat dSigma = directDescriptorD1(covDescriptor[iStruct], par(iStruct), k);
    const arma::uword q = qOf[iStruct];
    arma::mat D(nR, nR, arma::fill::zeros);
    for(arma::uword a = 0; a < q; ++a){
      for(arma::uword b = 0; b < q; ++b){
        if(dSigma(a,b) != 0.0){ D += dSigma(a,b) * Bcache[iStruct](a,b); }
      }
    }
    return D;
  };

  auto structureIsPD = [&](int iStruct, const arma::vec & par) -> bool {
    if(!par.is_finite()){ return false; }
    // Test shape separately from scale so a random variance can reach its
    // original-response-scale floor without failing a standardized PD cutoff.
    arma::vec shape = par;
    shape(0) = 0.0;
    if(iStruct == residualStruct && residualCompact){
      return directDescriptorPD(covDescriptor[iStruct], shape,
        std::numeric_limits<double>::min());
    }
    arma::mat Th = directEvaluateDescriptor(covDescriptor[iStruct], shape);
    if(!Th.is_finite()){ return false; }
    arma::vec eigval;
    bool ok = arma::eig_sym(eigval, arma::symmatu(Th));
    if(iStruct != residualStruct){
      const Rcpp::List factors = covDescriptor[iStruct]["factors"];
      bool reducedRank = false;
      for(int axis = 0; axis < factors.size(); ++axis){
        reducedRank = reducedRank ||
          sommer::zeroSpecificFactor(Rcpp::as<Rcpp::List>(factors[axis]));
      }
      if(reducedRank){
        const double scale = std::max(1.0, arma::norm(Th, "inf"));
        return ok && eigval.is_finite() && eigval.n_elem > 0 &&
          eigval.min() >= -1.0e-10 * scale && eigval.max() > 0.0;
      }
    }
    return ok && eigval.is_finite() && eigval.n_elem > 0 && eigval.min() > 0.0;
  };

  // Constant across iterations: only free (non-fixed) parameters are
  // updated, and their bounds never change.
  arma::vec freeUnlisted(nVcTotal), lowerUnlisted(nVcTotal), upperUnlisted(nVcTotal);
  for(int i = 0; i < nRRe; ++i){
    freeUnlisted.subvec(nVcStart(i)-1, nVcEnd(i)-1) = covFree(i);
    lowerUnlisted.subvec(nVcStart(i)-1, nVcEnd(i)-1) = covLower(i);
    upperUnlisted.subvec(nVcStart(i)-1, nVcEnd(i)-1) = covUpper(i);
  }
  arma::uvec freeIdx = arma::find(freeUnlisted > 0.5);
  const VarianceParameterMap vcMap =
    readVarianceParameterMap(covStructI, static_cast<arma::uword>(nVcTotal));
  arma::mat tangent(nVcTotal, freeIdx.n_elem, arma::fill::zeros);
  for(arma::uword j = 0; j < freeIdx.n_elem; ++j){ tangent(freeIdx(j), j) = 1.0; }
  if(vcMap.active){
    tangent = vcMap.T;
    for(int k = 0; k < nVcTotal; ++k){
      if(freeUnlisted(k) < 0.5 && arma::any(tangent.row(k) != 0.0)){
        Rcpp::stop("The direct-engine constraint map contains a free direction through a fixed parameter.");
      }
    }
    const arma::uvec used = arma::find(arma::sum(arma::abs(tangent), 0).t() > 0.0);
    tangent = tangent.cols(used);
  }
  auto feasibleTangent = [&](const arma::mat & basis, const arma::vec & point,
                             const arma::vec & direction) -> arma::mat {
    std::vector<arma::uword> active;
    for(int k = 0; k < nVcTotal; ++k){
      const double slack = 1.0e-10 * (1.0 + std::abs(point(k)));
      if(((point(k) <= lowerUnlisted(k) + slack && direction(k) < 0.0) ||
          (point(k) >= upperUnlisted(k) - slack && direction(k) > 0.0)) &&
         arma::any(basis.row(k) != 0.0)){
        active.push_back(k);
      }
    }
    if(active.empty()){ return basis; }
    return basis * arma::null(basis.rows(arma::uvec(active)));
  };

  // Standardizing y by stdy shifts log|V| and log|Q| by constants that
  // cancel in yPy but not in the reported likelihood; add the Jacobian
  // correction back so llik/AIC/BIC are reported on the original scale.
  const double llikScaleCorrection =
    -(reml ? (static_cast<double>(nR) - static_cast<double>(rankX)) : static_cast<double>(nR)) * std::log(stdy);

  arma::mat monitor(nVcTotal, nIters, arma::fill::zeros);
  arma::vec llik(nIters, arma::fill::zeros);
  arma::mat InfMatInvWorking(nVcTotal, nVcTotal, arma::fill::eye);
  bool convergence = false;
  int lastIter = 0;
  arma::vec beta(nX, arma::fill::zeros);
  arma::mat VarBeta(nX, nX, arma::fill::zeros);
  arma::vec Py(nR, arma::fill::zeros);
  sommer::MMEProgress progress(nRRe);
  const double scoreTolerance = std::sqrt(std::max(tolParConvLL, tolParConvNorm));
  const arma::vec stepCaps = sommer::covarianceWorkingStepCaps(covStructI);
  double projectedScoreNorm = std::numeric_limits<double>::infinity();
  int lineSearchHalvings = 0;
  std::string stopReason = "iteration-limit";
  const sommer::control::Transport controls(covStructI);
  sommer::control::Proposal proposal;
  proposal.adaptiveRecovery = controls.logical("adaptiveRecovery");
  const bool adaptive = proposal.adaptiveRecovery;
  proposal.blend = controls.blend("Unknown direct-engine EM blending method.");
  const std::string blendMethod = sommer::control::blendName(proposal.blend);
  const bool parameterExpanded = proposal.blend == sommer::control::Blend::PXEM;
  sommer::AIRecoveryPolicy recoveryPolicy;
  bool previousAI = false;
  int previousHalvings = 0;
  int emPreparations = 0, closedSteps = 0, generalizedSteps = 0;
  int qEvaluations = 0, aiUpdateSolves = 0, likelihoodEvaluations = 0;
  int emFallbacks = 0, pxFallbacks = 0, acceptedRecoveries = 0;
  int regularizedAttempts = 0, regularizedAccepted = 0;
  std::vector<int> recoveryIterations;
  std::vector<std::string> recoveryReasons;
  std::vector<double> regularizedRidges;
  std::vector<int> regularizedIterations;
  arma::vec weightHistory(nIters, arma::fill::zeros);
  arma::vec qHistory(nIters, arma::fill::value(NA_REAL));
  arma::vec pxHistory(nIters, arma::fill::value(NA_REAL));
  arma::ivec mstepHistory(nIters, arma::fill::value(NA_INTEGER));
  arma::ivec mstepEvaluations(nIters, arma::fill::zeros);
  arma::ivec pxConvergenceHistory(nIters, arma::fill::value(NA_INTEGER));
  arma::ivec pxEvaluationHistory(nIters, arma::fill::zeros);
  arma::mat expansionHistory(nRe, nIters, arma::fill::zeros);
  std::vector<arma::uword> emWidths(nRe);
  for(int term = 0; term < nRe; ++term) { emWidths[term] = Zdense[term][0].n_cols; }
  arma::mat pendingVi, pendingP;
  double pendingLikelihood = NA_REAL;
  bool pending = false;

  auto evaluatePoint = [&](const arma::field<arma::vec> & point, arma::mat & Vi,
                           arma::mat & P, double & likelihood) -> bool {
    ++likelihoodEvaluations;
    arma::mat V = arma::symmatu(buildV(point));
    if(!V.is_finite() || !arma::inv_sympd(Vi, V)){ return false; }
    const arma::mat VX = Vi * Xd;
    const arma::mat Q = Xd.t() * VX;
    arma::mat Qinv;
    if(!arma::inv_sympd(Qinv, Q)){ return false; }
    P = arma::symmatu(Vi - VX * Qinv * VX.t());
    double logDetV, signV, logDetQ = 0.0, signQ = 1.0;
    arma::log_det(logDetV, signV, V);
    if(reml){ arma::log_det(logDetQ, signQ, Q); }
    const double degrees = reml ? nR - rankX : nR;
    likelihood = -0.5 * (logDetV + logDetQ + arma::as_scalar(y.t() * P * y) +
      degrees * std::log(2.0 * arma::datum::pi)) + llikScaleCorrection;
    return signV > 0.0 && signQ > 0.0 && std::isfinite(likelihood);
  };

  for(int iIter = 0; iIter < nIters; ++iIter){
    Rcpp::checkUserInterrupt();
    lastIter = iIter;
    arma::mat Vi, P;
    double likelihood;
    if(pending) {
      Vi = std::move(pendingVi);
      P = std::move(pendingP);
      likelihood = pendingLikelihood;
      pending = false;
    } else if(!evaluatePoint(covPar, Vi, P, likelihood)){
      Rcpp::stop("Direct-engine V or X'V^{-1}X is numerically singular; check covariance structures and fixed-effect collinearity.");
    }
    Py = P * y;
    llik(iIter) = likelihood;
    arma::vec thetaUnlisted(nVcTotal);
    for(int i = 0; i < nRRe; ++i){
      thetaUnlisted.subvec(nVcStart(i)-1, nVcEnd(i)-1) = covPar(i);
      lowerUnlisted(nVcStart(i)-1) = std::min(covPar(i)(0), covLower(i)(0));
    }
    monitor.col(iIter) = thetaUnlisted;

    // Score and average information.
    arma::vec score(nVcTotal, arma::fill::zeros);
    arma::mat AI(nVcTotal, nVcTotal, arma::fill::zeros);
    std::vector<arma::vec> vList(nVcTotal), PvList(nVcTotal);

    int g = 0;
    for(int i = 0; i < nRRe; ++i){
      for(arma::uword k = 0; k < covPar(i).n_elem; ++k){
        arma::mat Dk = buildD(covPar, i, k);
        arma::vec vk = Dk * Py;
        const double traceTerm = reml ? arma::accu(P % Dk) : arma::accu(Vi % Dk);
        score(g) = 0.5 * (arma::dot(Py, vk) - traceTerm);
        vList[g] = vk;
        PvList[g] = P * vk;
        g++;
      }
    }
    double weight = weightEmInf(iIter);
    if(weight < 0.0 || weight > 1.0 || !std::isfinite(weight)) {
      Rcpp::stop("Direct-engine blend weights must be finite and between 0 and 1.");
    }
    if(weight < 1.0) {
      for(int a = 0; a < nVcTotal; ++a){
        for(int b = a; b < nVcTotal; ++b){
          AI(a,b) = 0.5 * arma::dot(vList[a], PvList[b]);
          AI(b,a) = AI(a,b);
        }
      }
    }
    if(!score.is_finite() || !AI.is_finite()){
      Rcpp::stop("Non-finite direct-engine covariance score or average information.");
    }
    const arma::mat scoreTangent = feasibleTangent(tangent, thetaUnlisted,
      arma::vec(tangent * (tangent.t() * score)));
    const arma::vec projectedScore = scoreTangent.t() * score;
    projectedScoreNorm = 0.0;
    for(arma::uword j = 0; j < projectedScore.n_elem; ++j){
      arma::vec v(nR, arma::fill::zeros), pv(nR, arma::fill::zeros);
      for(int k = 0; k < nVcTotal; ++k) {
        v += scoreTangent(k, j) * vList[k];
        pv += scoreTangent(k, j) * PvList[k];
      }
      const double curvature = 0.5 * arma::dot(v, pv);
      projectedScoreNorm = std::max(projectedScoreNorm, std::abs(projectedScore(j)) /
        std::sqrt(std::max(std::numeric_limits<double>::min(), std::abs(curvature))));
    }
    std::string recoveryReason;
    bool regularized = false;
    bool ordinaryFallback = false;
    double progressTrust = NA_REAL, progressLine = NA_REAL;
    if(adaptive && iIter < nIters - 1) {
      const char * reason = recoveryPolicy.observe(previousAI, previousHalvings,
        iIter > 0 ? likelihood - llik(iIter - 1) : 0.0, tolParConvLL,
        projectedScoreNorm <= scoreTolerance);
      if(reason) { recoveryReason = reason; weight = 1.0; recoveryPolicy.reset(); }
    }
    auto printProgress = [&]() {
      weightHistory(iIter) = weight;
      if(!verbose) { return; }
      if(iIter == 0) {
        Rcpp::Rcout << "Direct path: dense phenotypic V Cholesky ("
                    << (residualCompact ? "compact" : "explicit") << " residual)"
                    << (adaptive ? ", AI-first EM/PX-GEM recovery" : "") << arma::endl;
        progress.header();
      }
      progress.row(iIter + 1, likelihood, iIter > 0 ? llik(iIter - 1) : NA_REAL,
        regularized && weight == 0.0 ? "AI-reg" :
          sommer::MMEProgress::update(weight, parameterExpanded && !ordinaryFallback),
        regularized && weight == 0.0 ? "info-solve" :
          sommer::MMEProgress::policy(adaptive, parameterExpanded, weight, recoveryReason),
        progressTrust, progressLine);
    };
    const bool likelihoodSmall = iIter > 0 &&
      std::abs(likelihood - llik(iIter-1)) < tolParConvLL;
    if(scoreTangent.n_cols == 0 || (likelihoodSmall && projectedScoreNorm <= scoreTolerance)){
      convergence = true;
      stopReason = "converged";
      printProgress();
      break;
    }

    arma::vec delta(nVcTotal, arma::fill::zeros);
    arma::vec cappedDelta(nVcTotal, arma::fill::zeros);
    arma::mat feasible = scoreTangent;
    if(adaptive && weight < 1.0 && iIter < nIters - 1) {
      ++aiUpdateSolves;
      const sommer::FeasibleAIDirection proposal = sommer::feasibleAIDirection(
        AI, score, arma::mat(), tangent, thetaUnlisted, lowerUnlisted, stepCaps);
      regularizedAttempts += proposal.regularizedAttempts;
      if(proposal.reliable) {
        delta = proposal.direction;
        cappedDelta = proposal.cappedDirection;
        regularized = proposal.ridge > 0.0;
        if(regularized) {
          regularizedRidges.push_back(proposal.ridge);
          regularizedIterations.push_back(iIter + 1);
        }
      } else {
        weight = 1.0;
        recoveryReason = "information-solve";
      }
    }
    for(int activeSet = 0; !adaptive && weight < 1.0 && activeSet <= nVcTotal && feasible.n_cols > 0; ++activeSet){
      const arma::mat information = feasible.t() * AI * feasible;
      const arma::vec gradient = feasible.t() * score;
      arma::vec scaling = arma::sqrt(arma::clamp(arma::abs(information.diag()),
        std::numeric_limits<double>::min(), std::numeric_limits<double>::max()));
      const arma::mat normalized = information / (scaling * scaling.t());
      const arma::vec rhs = gradient / scaling;
      arma::vec direction;
      ++aiUpdateSolves;
      bool solved = arma::solve(direction, normalized, rhs,
        arma::solve_opts::likely_sympd + arma::solve_opts::no_approx);
      if(!solved){ direction = arma::pinv(normalized) * rhs; }
      delta = feasible * (direction / scaling);
      if(!delta.is_finite() || arma::dot(delta, score) <= 0.0){
        delta = feasible * (gradient / arma::square(scaling));
      }
      if(!delta.is_finite()){ Rcpp::stop("Non-finite direct-engine covariance direction."); }
      const arma::mat next = feasibleTangent(feasible, thetaUnlisted, delta);
      if(next.n_cols == feasible.n_cols){ break; }
      feasible = next;
      delta.zeros();
    }
    sommer::DirectEMProposal em;
    bool momentsReady = false;
    auto prepareEM = [&]() {
      if(momentsReady) { return; }
      double workspaceEntries = 24.0 * static_cast<double>(nR) * nR +
        8.0 * static_cast<double>(nVcTotal) * nVcTotal;
      for(int term = 0; term < nRe; ++term) {
        workspaceEntries += 12.0 * static_cast<double>(qOf[term]) * qOf[term];
      }
      if(covStructI.hasAttribute("denseMemoryMB")) {
        const double budget = Rcpp::as<double>(covStructI.attr("denseMemoryMB"));
        if(workspaceEntries * sizeof(double) > budget * 1024.0 * 1024.0) {
          Rcpp::stop("Exact direct EM moments exceed denseMemoryMB. Increase the engine memory budget; no approximate moments or diagonal AI substitution are used.");
        }
      }
      em = sommer::directEMProposal(covDescriptor, covPar, freeUnlisted,
        lowerUnlisted, upperUnlisted, tangent, vcMap.active,
        Bcache, emWidths, Xd, y, Vi, P, reml,
        [&](const arma::vec & par) { return buildResidual(par, -1); },
        parameterExpanded);
      momentsReady = true;
      ++emPreparations;
      if(em.closed) { ++closedSteps; } else { ++generalizedSteps; }
      qEvaluations += em.evaluations;
      qHistory(iIter) = em.qGain;
      pxHistory(iIter) = em.pxGain;
      mstepHistory(iIter) = em.converged;
      mstepEvaluations(iIter) = em.ordinaryEvaluations;
      if(parameterExpanded && nRe > 0) {
        pxConvergenceHistory(iIter) = em.pxConverged;
        pxEvaluationHistory(iIter) = em.pxEvaluations;
      }
      if(nRe > 0) { expansionHistory.col(iIter) = em.expansion; }
    };
    const arma::vec aiDirection = delta;
    if(weight > 0.0 && iIter < nIters - 1) {
      prepareEM();
      delta = (1.0 - weight) * aiDirection + weight * (em.point - thetaUnlisted);
    }
    const double relativeStep = arma::norm(delta, 2) / (1.0 + arma::norm(thetaUnlisted, 2));
    if(projectedScoreNorm <= scoreTolerance && relativeStep < tolParConvNorm){
      convergence = true;
      stopReason = "converged";
      printProgress();
      break;
    }
    // Never return a proposal whose likelihood and score were not evaluated.
    if(iIter == nIters - 1){ printProgress(); break; }

    bool accepted = false;
    const int halvingsBefore = lineSearchHalvings;
    const int phases = adaptive && weight == 0.0 &&
      projectedScoreNorm > scoreTolerance ? 2 : 1;
    for(int phase = 0; phase < phases && !accepted; ++phase) {
      if(phase == 1) {
        weight = 1.0;
        recoveryReason = regularized ? "information-solve" : "line-search";
        prepareEM();
        delta = em.point - thetaUnlisted;
        recoveryPolicy.reset();
      }
      const bool mixed = weight > 0.0 && weight < 1.0;
      const int attempts = 1 + (mixed ? 1 : 0) +
        (parameterExpanded && weight > 0.0 && nRe > 0 ? 1 : 0);
      for(int attempt = 0; attempt < attempts && !accepted; ++attempt) {
        if(attempt > 0) {
          weight = 1.0;
          if(mixed && attempt == 1) {
            delta = em.point - thetaUnlisted;
            ++emFallbacks;
          } else {
            delta = em.ordinary - thetaUnlisted;
            ++pxFallbacks;
            ordinaryFallback = true;
          }
        }
        const bool fullEM = weight == 1.0 || attempt > 0;
        double alpha = adaptive && fullEM ? 1.0 : weightInf(iIter);
        for(int k = 0; k < nVcTotal; ++k) {
          const double capped = adaptive && weight == 0.0 ? cappedDelta(k) : delta(k);
          if(!fullEM && capped != 0.0) {
            alpha = std::min(alpha, stepCaps(k) / std::abs(capped));
          }
          if(delta(k) < 0.0 && std::isfinite(lowerUnlisted(k))) {
            alpha = std::min(alpha, std::max(0.0, (thetaUnlisted(k) - lowerUnlisted(k)) / -delta(k)));
          } else if(delta(k) > 0.0 && std::isfinite(upperUnlisted(k))) {
            alpha = std::min(alpha, std::max(0.0, (upperUnlisted(k) - thetaUnlisted(k)) / delta(k)));
          }
        }
        progressTrust = alpha;
        for(int search = 0; search < 30 && alpha > 0.0; ++search) {
          const arma::vec candidate = thetaUnlisted + alpha * delta;
          arma::field<arma::vec> trial(nRRe);
          bool valid = candidate.is_finite() && arma::all(candidate >= lowerUnlisted) &&
            arma::all(candidate <= upperUnlisted);
          for(int i = 0; i < nRRe && valid; ++i) {
            trial(i) = candidate.subvec(nVcStart(i)-1, nVcEnd(i)-1);
            valid = structureIsPD(i, trial(i));
          }
          arma::mat trialVi, trialP;
          double trialLikelihood;
          const double minimumGain = adaptive && (fullEM || regularized) ? 1.0e-10 : 0.0;
          const bool evaluated = valid && evaluatePoint(trial, trialVi, trialP, trialLikelihood);
          const bool stationaryEM = fullEM && projectedScoreNorm <= scoreTolerance &&
            evaluated && trialLikelihood >= likelihood -
              1.0e-10 * (1.0 + std::abs(likelihood));
          if(evaluated &&
             ((trialLikelihood > likelihood + minimumGain &&
               (fullEM || trialLikelihood >= likelihood +
                 1.0e-4 * alpha * arma::dot(score, delta))) || stationaryEM)) {
            covPar = trial;
            pendingVi = std::move(trialVi);
            pendingP = std::move(trialP);
            pendingLikelihood = trialLikelihood;
            pending = true;
            accepted = true;
            progressLine = alpha / progressTrust;
            if(attempt > 0) { weight = 1.0; }
            break;
          }
          alpha *= 0.5;
          ++lineSearchHalvings;
        }
      }
    }
    if(adaptive && weight == 1.0) {
      recoveryIterations.push_back(iIter + 1);
      recoveryReasons.push_back(recoveryReason);
      if(accepted) { ++acceptedRecoveries; }
    }
    if(accepted && regularized && weight == 0.0) { ++regularizedAccepted; }
    previousAI = accepted && weight == 0.0;
    previousHalvings = lineSearchHalvings - halvingsBefore;
    printProgress();
    if(!accepted){
      stopReason = "line-search-stalled";
      Rcpp::warning("Direct-engine optimization stalled at the last accepted parameters; convergence was not certified.");
      break;
    }
  }
  for(int i = 0; i < nRRe; ++i){
    if(!(i == residualStruct && residualCompact)){
      theta(i) = directEvaluateDescriptor(covDescriptor[i], covPar(i));
    }
  }

  // ------------------------------------------------------------
  // Final BLUE/BLUP/PEV computation at the converged (or last) point.
  // ------------------------------------------------------------
  arma::mat V = buildV(covPar);
  V = arma::symmatu(V);
  arma::mat Vi;
  bool okV = arma::inv_sympd(Vi, V);
  arma::mat Djit = arma::eye<arma::mat>(nR, nR);
  for(int jit = 0; !okV && jit < 4; ++jit){
    V = V + Djit * (tolParInv * std::pow(10.0, jit));
    okV = arma::inv_sympd(Vi, V);
  }
  if(!okV){ Rcpp::stop("Final V is numerically singular."); }

  arma::mat VX = Vi * Xd;
  arma::mat Q = Xd.t() * VX;
  arma::mat Qinv;
  bool okQ = arma::inv_sympd(Qinv, Q);
  arma::mat Dxj = arma::eye<arma::mat>(nX, nX);
  for(int jit = 0; !okQ && jit < 4; ++jit){
    Q = Q + Dxj * (tolParInv * std::pow(10.0, jit));
    okQ = arma::inv_sympd(Qinv, Q);
  }
  if(!okQ){ Rcpp::stop("Final X'V^{-1}X is numerically singular."); }

  arma::mat P = arma::symmatu(Vi - VX * Qinv * VX.t());
  Py = P * y;
  beta = Qinv * (Xd.t() * (Vi * y));
  VarBeta = Qinv;

  // Final inference uses clean AI, never the EM/GEM proposal or its curvature.
  {
    arma::mat AIFinal(nVcTotal, nVcTotal, arma::fill::zeros);
    std::vector<arma::vec> vListFinal(nVcTotal), PvListFinal(nVcTotal);
    int g = 0;
    for(int i = 0; i < nRRe; ++i){
      for(arma::uword k = 0; k < covPar(i).n_elem; ++k){
        arma::mat Dk = buildD(covPar, i, k);
        arma::vec vk = Dk * Py;
        vListFinal[g] = vk;
        PvListFinal[g] = P * vk;
        g++;
      }
    }
    for(int a = 0; a < nVcTotal; ++a){
      for(int b = a; b < nVcTotal; ++b){
        AIFinal(a,b) = 0.5 * arma::dot(vListFinal[a], PvListFinal[b]);
        AIFinal(b,a) = AIFinal(a,b);
      }
    }
    InfMatInvWorking.zeros();
    if(tangent.n_cols > 0){
      // Log-scale information vanishes at a variance boundary. Normalize
      // before inversion so that unrelated uncertainty is not truncated.
      const arma::mat information = tangent.t() * AIFinal * tangent;
      const arma::vec scaling = arma::sqrt(arma::clamp(arma::abs(information.diag()),
        std::numeric_limits<double>::min(), std::numeric_limits<double>::max()));
      const arma::mat normalized = information / (scaling * scaling.t());
      InfMatInvWorking = tangent *
        (arma::pinv(normalized) / (scaling * scaling.t())) * tangent.t();
    }
  }

  arma::field<arma::mat> uList(nRe), uPevList(nRe);
  arma::vec u(Nu, arma::fill::zeros);

  for(int i = 0; i < nRe; ++i){
    arma::mat Sigma = theta(i);
    const arma::uword q = qOf[i];
    arma::sp_mat AiSp = convertSparse(AiI(i));
    const bool aiIsIdentity = isIdentity_spmat(AiSp);
    arma::mat Kdense;
    if(!aiIsIdentity){
      arma::mat Adense = arma::mat(AiSp);
      bool ok = eigenSpdInverse(Adense, Kdense);
      if(!ok){
        arma::mat bend = nearPDcpp(Adense, 100, 1.0e-6, 1.0e-7);
        ok = eigenSpdInverse(bend, Kdense);
      }
    }
    const int levels = Zdense[i][0].n_cols;
    arma::mat uMat(levels, q, arma::fill::zeros);
    std::vector<arma::vec> ZtPy(q);
    for(arma::uword b = 0; b < q; ++b){ ZtPy[b] = Zdense[i][b].t() * Py; }
    for(arma::uword a = 0; a < q; ++a){
      arma::vec acc(levels, arma::fill::zeros);
      for(arma::uword b = 0; b < q; ++b){
        acc += Sigma(a,b) * ZtPy[b];
      }
      uMat.col(a) = aiIsIdentity ? acc : (Kdense * acc);
    }
    uList(i) = uMat;
    for(arma::uword a = 0; a < q; ++a){
      u.subvec(static_cast<arma::uword>(partitions(i)(a,0)-1-nX), static_cast<arma::uword>(partitions(i)(a,1)-1-nX)) = uMat.col(a);
    }
    uPevList(i) = arma::mat();
  }

  if(computePev == 2){
    if(Nu > 5000){
      Rcpp::stop("computePev=2 (full PEV) requires Nu <= 5000 coefficients in this engine; use computePev=0 for larger marker/coefficient models.");
    }
    arma::mat Zfull(nR, Nu, arma::fill::zeros);
    arma::mat Gfull(Nu, Nu, arma::fill::zeros);
    for(int i = 0; i < nRe; ++i){
      arma::mat Sigma = theta(i);
      const arma::uword q = qOf[i];
      arma::sp_mat AiSp = convertSparse(AiI(i));
      const bool aiIsIdentity = isIdentity_spmat(AiSp);
      arma::mat Kdense;
      if(!aiIsIdentity){
        arma::mat Adense = arma::mat(AiSp);
        bool ok = eigenSpdInverse(Adense, Kdense);
        if(!ok){
          arma::mat bend = nearPDcpp(Adense, 100, 1.0e-6, 1.0e-7);
          ok = eigenSpdInverse(bend, Kdense);
        }
      }else{
        Kdense = arma::eye<arma::mat>(Zdense[i][0].n_cols, Zdense[i][0].n_cols);
      }
      for(arma::uword a = 0; a < q; ++a){
        const int startA = static_cast<int>(partitions(i)(a,0)) - 1 - nX;
        const int endA = static_cast<int>(partitions(i)(a,1)) - 1 - nX;
        Zfull.cols(startA, endA) = Zdense[i][a];
        for(arma::uword b = 0; b < q; ++b){
          const int startB = static_cast<int>(partitions(i)(b,0)) - 1 - nX;
          const int endB = static_cast<int>(partitions(i)(b,1)) - 1 - nX;
          Gfull.submat(startA, startB, endA, endB) = Sigma(a,b) * Kdense;
        }
      }
    }
    arma::mat GZt = Gfull * Zfull.t();
    arma::mat VarU = GZt * P * GZt.t();
    arma::mat PevFull = Gfull - VarU;
    for(int i = 0; i < nRe; ++i){
      const arma::uword q = qOf[i];
      const int levels = Zdense[i][0].n_cols;
      arma::mat eMat(levels, q, arma::fill::zeros);
      for(arma::uword a = 0; a < q; ++a){
        const int startA = static_cast<int>(partitions(i)(a,0)) - 1 - nX;
        const int endA = static_cast<int>(partitions(i)(a,1)) - 1 - nX;
        eMat.col(a) = arma::diagvec(PevFull.submat(startA, startA, endA, endA)) * vary;
      }
      uPevList(i) = eMat;
    }
  }

  // ------------------------------------------------------------
  // Rescale to the original response scale and report natural-scale
  // covariance parameters + delta-method SE, identical convention to
  // ai_mme_sp2().
  // ------------------------------------------------------------
  for(int i = 0; i < nRRe; ++i){ theta(i) = theta(i) * vary; }
  beta = beta * stdy;
  if(hasIntercept && muy != 0.0){
    for(arma::uword j = 0; j < beta.n_elem; ++j){
      beta(j) += muy * cShift(j);
    }
  }
  VarBeta = VarBeta * vary;
  u = u * stdy;
  for(int i = 0; i < nRe; ++i){ uList(i) = uList(i) * stdy; }

  arma::field<arma::vec> covParOut(nRRe);
  arma::field<arma::vec> covParWorking(nRRe);
  arma::vec finalJacobian(nVcTotal, arma::fill::ones);
  arma::mat monitorOut = monitor;
  for(int i = 0; i < nRRe; ++i){
    arma::vec localJac;
    covParOut(i) = directReportParameters(covDescriptor[i], covPar(i), vary, localJac);
    covParWorking(i) = covPar(i);
    covParWorking(i)(0) += std::log(vary);
    finalJacobian.subvec(nVcStart(i)-1, nVcEnd(i)-1) = localJac;
    for(int iter = 0; iter <= lastIter; ++iter){
      arma::vec localWork = monitor.submat(nVcStart(i)-1, iter, nVcEnd(i)-1, iter);
      arma::vec dummyJac;
      monitorOut.submat(nVcStart(i)-1, iter, nVcEnd(i)-1, iter) =
        directReportParameters(covDescriptor[i], localWork, vary, dummyJac);
    }
  }
  arma::mat theta_se = arma::diagmat(finalJacobian) * InfMatInvWorking * arma::diagmat(finalJacobian);

  const double AIC = (-2.0 * llik(lastIter)) + (2.0 * nX);
  const double BIC = (-2.0 * llik(lastIter)) + (std::log(static_cast<double>(nR)) * nX);

  arma::vec bu(nEffects, arma::fill::zeros);
  bu.subvec(0, nX-1) = beta;
  if(Nu > 0){ bu.subvec(nX, nEffects-1) = u; }

  return Rcpp::List::create(
    Rcpp::Named("llik") = llik.subvec(0, lastIter),
    Rcpp::Named("b") = beta,
    Rcpp::Named("u") = u,
    Rcpp::Named("bu") = bu,
    Rcpp::Named("VarBeta") = VarBeta,
    Rcpp::Named("theta") = theta,
    Rcpp::Named("covPar") = covParOut,
    Rcpp::Named("covParWorking") = covParWorking,
    Rcpp::Named("theta_se") = theta_se,
    Rcpp::Named("monitor") = monitorOut.cols(0, lastIter),
    Rcpp::Named("uList") = uList,
    Rcpp::Named("uPevList") = uPevList,
    Rcpp::Named("AIC") = AIC,
    Rcpp::Named("BIC") = BIC,
    Rcpp::Named("convergence") = convergence,
    Rcpp::Named("partitions") = partitions,
    Rcpp::Named("CiMode") = 0,
    Rcpp::Named("Ci") = arma::sp_mat(),
    Rcpp::Named("engineDiagnostics") = Rcpp::List::create(
      Rcpp::Named("direct") = Rcpp::List::create(
        Rcpp::Named("optimizer") = adaptive ? "ai-first" :
          arma::all(weightEmInf == 0.0) ? "ai" :
          arma::all(weightEmInf == 1.0) ? blendMethod : "blend",
        Rcpp::Named("stopReason") = stopReason,
        Rcpp::Named("finalProjectedScoreNorm") = projectedScoreNorm,
        Rcpp::Named("scoreTolerance") = scoreTolerance,
        Rcpp::Named("lineSearchHalvings") = lineSearchHalvings,
        Rcpp::Named("aiUpdateSolves") = aiUpdateSolves,
        Rcpp::Named("likelihoodEvaluations") = likelihoodEvaluations),
      Rcpp::Named("blend") = Rcpp::List::create(
        Rcpp::Named("method") = blendMethod,
        Rcpp::Named("weights") = weightHistory.head(lastIter + 1),
        Rcpp::Named("qGains") = qHistory.head(lastIter + 1),
        Rcpp::Named("closedFormSteps") = closedSteps,
        Rcpp::Named("generalizedSteps") = generalizedSteps,
        Rcpp::Named("eStepEvaluations") = emPreparations,
        Rcpp::Named("qEvaluations") = qEvaluations,
        Rcpp::Named("numericalMstepConverged") = mstepHistory.head(lastIter + 1),
        Rcpp::Named("numericalMstepEvaluations") = mstepEvaluations.head(lastIter + 1),
        Rcpp::Named("pxQGains") = pxHistory.head(lastIter + 1),
        Rcpp::Named("pxMstepConverged") = pxConvergenceHistory.head(lastIter + 1),
        Rcpp::Named("pxMstepEvaluations") = pxEvaluationHistory.head(lastIter + 1),
        Rcpp::Named("pxLogExpansion") = expansionHistory.cols(0, lastIter),
        Rcpp::Named("emFallbacks") = emFallbacks,
        Rcpp::Named("pxEMFallbacks") = pxFallbacks),
      Rcpp::Named("recovery") = Rcpp::List::create(
        Rcpp::Named("automatic") = adaptive,
        Rcpp::Named("method") = blendMethod,
        Rcpp::Named("likelihood") = reml ? "REML" : "ML",
        Rcpp::Named("steps") = static_cast<int>(recoveryIterations.size()),
        Rcpp::Named("accepted") = acceptedRecoveries,
        Rcpp::Named("iterations") = recoveryIterations,
        Rcpp::Named("reasons") = recoveryReasons,
        Rcpp::Named("regularizedAttempts") = regularizedAttempts,
        Rcpp::Named("stopReason") = stopReason,
        Rcpp::Named("finalProjectedScoreNorm") = projectedScoreNorm,
        Rcpp::Named("scoreTolerance") = scoreTolerance,
        Rcpp::Named("regularizedAcceptedSteps") = regularizedAccepted,
        Rcpp::Named("regularizedIterations") = regularizedIterations,
        Rcpp::Named("regularizedRidges") = regularizedRidges)),
    Rcpp::Named("residualCovarianceCompact") = residualCompact
  );
}
