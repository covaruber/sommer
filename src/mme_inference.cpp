#include "mme_sp2_support.h"
#include "mme_progress.h"

arma::mat predict_mmes_vcov_cpp(Rcpp::List model, const arma::sp_mat & Dmat);

// [[Rcpp::export]]
Rcpp::List post_mme_Cinverse_cpp(Rcpp::List model, const int mode = 1){

  if(mode < 0 || mode > 2){
    Rcpp::stop("mode must be 0, 1, or 2.");
  }

  if(model.containsElementNamed("retainCoefficientMatrix") &&
     !Rcpp::as<bool>(model["retainCoefficientMatrix"])){
    Rcpp::stop("Coefficient matrix C was not retained; refit with retainCoefficientMatrix=TRUE for this method.");
  }
  if(!model.containsElementNamed("C") || Rf_isNull(model["C"])){
    Rcpp::stop("The fitted model does not contain C. Refit with ai_mme_sp2() returning C.");
  }
  if(!model.containsElementNamed("Cscale")){
    Rcpp::stop("The fitted model does not contain Cscale.");
  }
  if(!model.containsElementNamed("partitions")){
    Rcpp::stop("The fitted model does not contain partitions.");
  }

  arma::sp_mat C = Rcpp::as<arma::sp_mat>(model["C"]);
  const double Cscale = Rcpp::as<double>(model["Cscale"]);
  Rcpp::List partitions = model["partitions"];

  const int nEffects = static_cast<int>(C.n_rows);

  if(C.n_rows != C.n_cols){
    Rcpp::stop("Stored C must be square.");
  }

  if(mode == 0){
    arma::sp_mat emptyCi;
    Rcpp::List uPevList(partitions.size());
    for(int i = 0; i < partitions.size(); ++i){
      uPevList[i] = arma::mat();
    }
    model["Ci"] = emptyCi;
    model["uPevList"] = uPevList;
    model["CiComputed"] = false;
    model["CiMode"] = 0;
    return model;
  }

  if(model.containsElementNamed("coefficientMap") &&
     !Rf_isNull(model["coefficientMap"])){
    const arma::sp_mat mapping = Rcpp::as<arma::sp_mat>(model["coefficientMap"]);
    const arma::uword originalEffects = mapping.n_rows;
    arma::sp_mat Ci(originalEffects, originalEffects);
    if(mode == 2){
      // Invert in engine coordinates before expanding zero-constrained effects.
      const Rcpp::List engineModel = Rcpp::List::create(
        Rcpp::Named("C")=model["C"], Rcpp::Named("Cscale")=model["Cscale"]);
      const arma::sp_mat identity = arma::speye<arma::sp_mat>(C.n_rows, C.n_rows);
      const arma::sp_mat engineInverse(predict_mmes_vcov_cpp(engineModel, identity));
      Ci = mapping * engineInverse * mapping.t();
    }else{
      std::vector<bool> usedRows(originalEffects, false);
      for(arma::sp_mat::const_iterator entry = mapping.begin(); entry != mapping.end(); ++entry){
        if(*entry != 0.0) usedRows[entry.row()] = true;
      }
      std::vector<arma::uword> activeRows;
      for(arma::uword rowIndex = 0; rowIndex < originalEffects; ++rowIndex){
        if(usedRows[rowIndex]) activeRows.push_back(rowIndex);
      }
      for(arma::uword first = 0; first < activeRows.size(); first += 32){
        const arma::uword count = std::min<arma::uword>(32, activeRows.size() - first);
        arma::sp_mat selection(count, originalEffects);
        for(arma::uword row = 0; row < count; ++row){ selection(row, activeRows[first + row]) = 1.0; }
        const arma::mat block = predict_mmes_vcov_cpp(model, selection);
        for(arma::uword row = 0; row < count; ++row){
          Ci(activeRows[first + row], activeRows[first + row]) = block(row, row);
        }
      }
    }
    Rcpp::List pevs(partitions.size());
    for(int term = 0; term < partitions.size(); ++term){
      const arma::mat partition = Rcpp::as<arma::mat>(partitions[term]);
      const arma::uword width = partition(0, 1) - partition(0, 0) + 1;
      arma::mat pev(width, partition.n_rows);
      for(arma::uword column = 0; column < partition.n_rows; ++column){
        for(arma::uword row = 0; row < width; ++row){
          const arma::uword index = partition(column, 0) - 1 + row;
          pev(row, column) = Ci(index, index);
        }
      }
      pevs[term] = pev;
    }
    model["Ci"] = Ci;
    model["uPevList"] = pevs;
    model["CiComputed"] = mode == 2;
    model["CiMode"] = mode;
    return model;
  }

  typedef Eigen::SparseMatrix<double, Eigen::ColMajor, int> EigenSpMat;
  typedef Eigen::Triplet<double, int> EigenTriplet;
  typedef Eigen::SimplicialLDLT<
    EigenSpMat,
    Eigen::Lower,
    SommerSparseOrdering
  > EigenLDLT;

  EigenSpMat Ce(nEffects, nEffects);
  std::vector<EigenTriplet> triplets;
  triplets.reserve(static_cast<std::size_t>(C.n_nonzero));

  for(arma::sp_mat::const_iterator it = C.begin(); it != C.end(); ++it){
    triplets.emplace_back(
      static_cast<int>(it.row()),
      static_cast<int>(it.col()),
      (*it)
    );
  }

  Ce.setFromTriplets(triplets.begin(), triplets.end());
  Ce.makeCompressed();

  EigenLDLT Cfactor;
  Cfactor.compute(Ce);

  if(Cfactor.info() != Eigen::Success){
    Rcpp::stop("Sparse LDLT factorisation of stored C failed.");
  }

  Rcpp::List uPevList(partitions.size());
  arma::sp_mat Ci;

  if(mode == 1){

    struct SelectedInverseSubset {
      std::vector< std::vector<int> > rows;
      std::vector< std::vector<double> > values;
      std::vector<int> originalToPermuted;
    };

    const Eigen::VectorXd D = Cfactor.vectorD();
    const int n = static_cast<int>(D.size());

    EigenSpMat Lmat = Cfactor.matrixL();
    Lmat.makeCompressed();

    SelectedInverseSubset subset;
    subset.rows.resize(n);
    subset.values.resize(n);
    subset.originalToPermuted.assign(n, -1);

    const auto & perm = Cfactor.permutationP();

    for(int i = 0; i < n; ++i){
      const int p = perm.indices()(i);
      if(p < 0 || p >= n){
        Rcpp::stop("Invalid LDLT permutation.");
      }
      subset.originalToPermuted[i] = p;
    }

    for(int col = 0; col < n; ++col){
      subset.rows[col].push_back(col);

      for(EigenSpMat::InnerIterator it(Lmat, col); it; ++it){
        if(it.row() > col){
          subset.rows[col].push_back(it.row());
        }
      }

      std::sort(subset.rows[col].begin(), subset.rows[col].end());

      subset.rows[col].erase(
        std::unique(subset.rows[col].begin(), subset.rows[col].end()),
        subset.rows[col].end()
      );

      subset.values[col].assign(subset.rows[col].size(), 0.0);
    }

    auto getPermuted = [&](int a, int b, double & value) -> bool {
      const int col = std::min(a,b);
      const int row = std::max(a,b);
      const std::vector<int> & rr = subset.rows[col];

      auto pos = std::lower_bound(rr.begin(), rr.end(), row);

      if(pos == rr.end() || (*pos) != row){
        return false;
      }

      const std::size_t idx =
        static_cast<std::size_t>(pos - rr.begin());

      value = subset.values[col][idx];
      return true;
    };

    auto setPermuted = [&](int row, int col, double value) {
      if(row < col){
        std::swap(row,col);
      }

      const std::vector<int> & rr = subset.rows[col];
      auto pos = std::lower_bound(rr.begin(), rr.end(), row);

      if(pos == rr.end() || (*pos) != row){
        Rcpp::stop("Internal Takahashi pattern error.");
      }

      const std::size_t idx =
        static_cast<std::size_t>(pos - rr.begin());

      subset.values[col][idx] = value;
    };

    for(int i = n - 1; i >= 0; --i){

      if(!std::isfinite(D(i)) || D(i) <= 0.0){
        Rcpp::stop("Non-positive/non-finite LDLT pivot during Takahashi inversion.");
      }

      std::vector<int> neighbours;
      std::vector<double> lvalues;

      for(EigenSpMat::InnerIterator it(Lmat, i); it; ++it){
        if(it.row() > i){
          neighbours.push_back(it.row());
          lvalues.push_back(it.value());
        }
      }

      for(std::size_t jj = 0; jj < neighbours.size(); ++jj){
        const int j = neighbours[jj];
        double sum = 0.0;

        for(std::size_t kk = 0; kk < neighbours.size(); ++kk){
          double zkj = 0.0;

          if(!getPermuted(neighbours[kk], j, zkj)){
            Rcpp::stop("Takahashi recurrence requested an unavailable entry.");
          }

          sum += lvalues[kk] * zkj;
        }

        setPermuted(j, i, -sum);
      }

      double diagCorrection = 0.0;

      for(std::size_t kk = 0; kk < neighbours.size(); ++kk){
        double zki = 0.0;

        if(!getPermuted(neighbours[kk], i, zki)){
          Rcpp::stop("Unable to retrieve Takahashi off-diagonal entry.");
        }

        diagCorrection += lvalues[kk] * zki;
      }

      setPermuted(i, i, (1.0 / D(i)) - diagCorrection);
    }

    auto getOriginal = [&](int originalRow, int originalCol, double & value) -> bool {
      const int a = subset.originalToPermuted[originalRow];
      const int b = subset.originalToPermuted[originalCol];
      return getPermuted(a, b, value);
    };

    for(int i = 0; i < partitions.size(); ++i){

      arma::mat p = Rcpp::as<arma::mat>(partitions[i]);

      if(p.n_rows == 0 || p.n_cols < 2){
        Rcpp::stop("Invalid random-effect partition matrix.");
      }

      const arma::uword blockSize =
        static_cast<arma::uword>(p(0,1) - p(0,0) + 1);

      arma::mat eMat(blockSize, p.n_rows);

      for(arma::uword j = 0; j < p.n_rows; ++j){

        const arma::uword first =
          static_cast<arma::uword>(p(j,0) - 1);

        const arma::uword last =
          static_cast<arma::uword>(p(j,1) - 1);

        if(last < first || last >= C.n_rows){
          Rcpp::stop("Random-effect partition is outside stored C.");
        }

        if((last-first+1) != blockSize){
          Rcpp::stop("Inconsistent random-effect block sizes.");
        }

        for(arma::uword k = first; k <= last; ++k){

          double cii = 0.0;

          if(!getOriginal(
               static_cast<int>(k),
               static_cast<int>(k),
               cii
             )){
            Rcpp::stop("Required diagonal inverse entry missing from Takahashi subset.");
          }

          eMat(k-first, j) = cii * Cscale;
        }
      }

      uPevList[i] = eMat;
    }

    Ci.reset();
  }

  if(mode == 2){

    Eigen::MatrixXd identity =
      Eigen::MatrixXd::Identity(
        static_cast<Eigen::Index>(nEffects),
        static_cast<Eigen::Index>(nEffects)
      );

    Eigen::MatrixXd CiEig = Cfactor.solve(identity);

    if(Cfactor.info() != Eigen::Success){
      Rcpp::stop("Sparse LDLT solve failed while computing full C inverse.");
    }

    arma::mat CiDense(
      static_cast<arma::uword>(nEffects),
      static_cast<arma::uword>(nEffects)
    );

    std::copy(
      CiEig.data(),
      CiEig.data() + CiEig.size(),
      CiDense.memptr()
    );

    CiDense = 0.5 * (CiDense + CiDense.t());
    CiDense *= Cscale;
    Ci = arma::sp_mat(CiDense);

    for(int i = 0; i < partitions.size(); ++i){

      arma::mat p = Rcpp::as<arma::mat>(partitions[i]);

      const arma::uword blockSize =
        static_cast<arma::uword>(p(0,1) - p(0,0) + 1);

      arma::mat eMat(blockSize, p.n_rows);

      for(arma::uword j = 0; j < p.n_rows; ++j){

        const arma::uword first =
          static_cast<arma::uword>(p(j,0) - 1);

        const arma::uword last =
          static_cast<arma::uword>(p(j,1) - 1);

        eMat.col(j) =
          arma::diagvec(
            Ci.submat(first, first, last, last)
          );
      }

      uPevList[i] = eMat;
    }
  }

  model["Ci"] = Ci;
  model["uPevList"] = uPevList;
  model["CiComputed"] = (mode == 2);
  model["CiMode"] = mode;

  return model;
}

class SommerPredictionFactor {
  typedef Eigen::SparseMatrix<double, Eigen::ColMajor, int> EigenSpMat;
  typedef Eigen::Triplet<double, int> EigenTriplet;
  typedef Eigen::SimplicialLDLT<
    EigenSpMat, Eigen::Lower, SommerSparseOrdering
  > EigenLDLT;

  EigenLDLT factor;
  arma::sp_mat mapping;
  bool mapped;
  double scale;
  int nEffects;

public:
  explicit SommerPredictionFactor(const Rcpp::List & model) : mapped(false) {
  if(model.containsElementNamed("retainCoefficientMatrix") &&
     !Rcpp::as<bool>(model["retainCoefficientMatrix"])){
    Rcpp::stop("Coefficient matrix C was not retained; refit with retainCoefficientMatrix=TRUE for this method.");
  }
  if(!model.containsElementNamed("C") || Rf_isNull(model["C"])){
    Rcpp::stop("The fitted model does not contain C. Refit with ai_mme_sp2() returning C.");
  }
  if(!model.containsElementNamed("Cscale")){
    Rcpp::stop("The fitted model does not contain Cscale.");
  }

  arma::sp_mat C = Rcpp::as<arma::sp_mat>(model["C"]);
  scale = Rcpp::as<double>(model["Cscale"]);
  if(model.containsElementNamed("coefficientMap") &&
     !Rf_isNull(model["coefficientMap"])){
    mapping = Rcpp::as<arma::sp_mat>(model["coefficientMap"]);
    mapped = true;
  }

  if(C.n_rows != C.n_cols){
    Rcpp::stop("Stored C must be square.");
  }

  if(C.n_rows > static_cast<arma::uword>(std::numeric_limits<int>::max())){
    Rcpp::stop("Stored C is too large for sparse prediction solves.");
  }
  nEffects = static_cast<int>(C.n_rows);
  if(mapped && mapping.n_cols != C.n_rows){
    Rcpp::stop("Coefficient map does not conform to stored C.");
  }

  EigenSpMat Ce(nEffects, nEffects);
  std::vector<EigenTriplet> triplets;
  triplets.reserve(static_cast<std::size_t>(C.n_nonzero));

  for(arma::sp_mat::const_iterator it = C.begin(); it != C.end(); ++it){
    triplets.emplace_back(
      static_cast<int>(it.row()),
      static_cast<int>(it.col()),
      (*it)
    );
  }

  Ce.setFromTriplets(triplets.begin(), triplets.end());
  Ce.makeCompressed();

  factor.compute(Ce);

  if(factor.info() != Eigen::Success){
    Rcpp::stop("Sparse LDLT factorisation of stored C failed.");
  }
  }

  arma::mat covariance(const arma::sp_mat & Dmat) const {
  arma::sp_mat design = Dmat;
  if(mapped){
    if(design.n_cols != mapping.n_rows){
      Rcpp::stop("D must have one column per original mixed-model coefficient.");
    }
    design = design * mapping;
  }
  if(design.n_cols != static_cast<arma::uword>(nEffects)){
    Rcpp::stop("D must have one column per mixed-model effect (nrow/ncol of stored C).");
  }

  const int k = static_cast<int>(design.n_rows);

  Eigen::MatrixXd rhs = Eigen::MatrixXd::Zero(nEffects, k);
  for(arma::sp_mat::const_iterator it = design.begin(); it != design.end(); ++it){
    rhs(static_cast<Eigen::Index>(it.col()), static_cast<Eigen::Index>(it.row())) = (*it);
  }

  Eigen::MatrixXd X = factor.solve(rhs);

  if(factor.info() != Eigen::Success){
    Rcpp::stop("Sparse LDLT solve failed while computing prediction variances.");
  }

  const arma::mat Xarma(X.data(), static_cast<arma::uword>(nEffects), static_cast<arma::uword>(k));

  arma::mat vcov = (design * Xarma) * scale;
  vcov = 0.5 * (vcov + vcov.t());

  return vcov;
  }
};

// [[Rcpp::export]]
arma::mat predict_mmes_vcov_cpp(Rcpp::List model, const arma::sp_mat & Dmat){
  return SommerPredictionFactor(model).covariance(Dmat);
}

// Operation-local factorization; never retained in a fitted/serialized model.
// [[Rcpp::export]]
SEXP mmes_vcov_factor_cpp(Rcpp::List model){
  Rcpp::XPtr<SommerPredictionFactor> pointer(new SommerPredictionFactor(model), true);
  pointer.attr("class") = "sommer_prediction_factor";
  return pointer;
}

// [[Rcpp::export]]
arma::mat mmes_vcov_solve_cpp(SEXP pointer, const arma::sp_mat & Dmat){
  if(TYPEOF(pointer) != EXTPTRSXP ||
     !Rf_inherits(pointer, "sommer_prediction_factor") ||
     R_ExternalPtrAddr(pointer) == NULL){
    Rcpp::stop("Prediction factor is invalid; create a new operation-local solver.");
  }
  Rcpp::XPtr<SommerPredictionFactor> factor(pointer);
  return factor->covariance(Dmat);
}
