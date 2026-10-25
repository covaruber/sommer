#include "mme_sp2_support.h"
#include "mme_progress.h"

namespace {

struct RelationshipMarkers {
  arma::uvec columns;
  arma::vec frequencies;

  RelationshipMarkers(const arma::mat & X, double minMAF)
    : columns(X.n_cols), frequencies(X.n_cols) {
    if(X.n_rows < 2 || X.n_cols == 0) {
      Rcpp::stop("X must have at least two individuals and one marker.");
    }
    if(!std::isfinite(minMAF) || minMAF < 0 || minMAF > 0.5) {
      Rcpp::stop("min.MAF must be finite and between zero and 0.5.");
    }
    arma::uword count = 0;
    for(arma::uword j = 0; j < X.n_cols; ++j) {
      if(j % 1024 == 0) { Rcpp::checkUserInterrupt(); }
      const double * x = X.colptr(j);
      double total = 0;
      bool variable = false;
      for(arma::uword i = 0; i < X.n_rows; ++i) {
        if(!std::isfinite(x[i]) || x[i] < -1 || x[i] > 1) {
          Rcpp::stop("X must contain finite genotypes in [-1, 1] after imputation.");
        }
        total += x[i] + 1;
        variable = variable || x[i] != x[0];
      }
      const double p = total / (2.0 * X.n_rows);
      if(variable && std::min(p, 1 - p) > minMAF) {
        columns[count] = j;
        frequencies[count++] = p;
      }
    }
    if(count == 0) {
      Rcpp::stop("No variable markers remain after min.MAF filtering.");
    }
    columns.resize(count);
    frequencies.resize(count);
  }
};

arma::mat relationshipKernel(const arma::mat & X,
                             const RelationshipMarkers & markers,
                             bool additive, bool method) {
  arma::mat W(X.n_rows, markers.columns.n_elem);
  double denominator = 0;
  for(arma::uword k = 0; k < W.n_cols; ++k) {
    if(k % 1024 == 0) { Rcpp::checkUserInterrupt(); }
    const double p = markers.frequencies[k];
    const double pq2 = 2 * p * (1 - p);
    const double * x = X.colptr(markers.columns[k]);
    double * w = W.colptr(k);
    if(additive) {
      const double center = method ? 2 * p : arma::mean(X.col(markers.columns[k]));
      const double offset = method ? 1 : 0;
      for(arma::uword i = 0; i < W.n_rows; ++i) {
        w[i] = x[i] + offset - center;
      }
      denominator += pq2;
    }else if(method) {
      // Orthogonal dominance coding also applies to fractional genotypes.
      for(arma::uword i = 0; i < W.n_rows; ++i) {
        w[i] = (2 * p - 1) * x[i] - x[i] * x[i] + pq2;
      }
      denominator += pq2 * pq2;
    }else{
      for(arma::uword i = 0; i < W.n_rows; ++i) {
        w[i] = 1 - std::abs(x[i]) - pq2;
      }
      denominator += pq2 * (1 - pq2);
    }
  }
  arma::mat result = W * W.t();
  if(additive && !method) {
    denominator = arma::mean(result.diag());
  }
  if(!std::isfinite(denominator) || denominator <= 0) {
    Rcpp::stop("The relationship matrix normalization must be finite and positive.");
  }
  result /= denominator;
  return result;
}

void validateRelationshipMatrix(const arma::mat & matrix, const char * name) {
  if(matrix.n_rows == 0 || matrix.n_rows != matrix.n_cols || !matrix.is_finite()) {
    Rcpp::stop("%s must be a finite, non-empty square matrix.", name);
  }
  if(!matrix.is_symmetric() && !matrix.is_symmetric(1e-10)) {
    Rcpp::stop("%s must be symmetric.", name);
  }
}

arma::mat relationshipCholesky(arma::mat & matrix, double ridge, const char * name) {
  if(!matrix.is_finite()) {
    Rcpp::stop("%s must contain only finite values.", name);
  }
  arma::mat factor;
  if(!arma::chol(factor, matrix)) {
    matrix.diag() += ridge;
    if(!matrix.is_finite() || !arma::chol(factor, matrix)) {
      Rcpp::stop("%s is not positive definite even after adding tolparinv to its diagonal.", name);
    }
  }
  return factor;
}

arma::mat relationshipSolve(const arma::mat & factor, const arma::mat & rhs) {
  arma::mat intermediate, result;
  const auto options = arma::solve_opts::fast + arma::solve_opts::no_approx;
  if(!arma::solve(intermediate, arma::trimatl(factor.t()), rhs, options) ||
     !arma::solve(result, arma::trimatu(factor), intermediate, options) ||
     !result.is_finite()) {
    Rcpp::stop("The relationship matrix Cholesky solve failed.");
  }
  return result;
}

arma::mat relationshipInverse(const arma::mat & factor) {
  if(factor.n_rows > static_cast<arma::uword>(std::numeric_limits<arma::blas_int>::max())) {
    Rcpp::stop("The relationship matrix exceeds the LAPACK integer limit.");
  }
  arma::mat result = factor;
  arma::blas_int n = static_cast<arma::blas_int>(factor.n_rows);
  char upper = 'U';
  arma::blas_int info = 0;
  // POTRI reuses the Cholesky factor without a second factorization.
  arma::lapack::potri(&upper, &n, result.memptr(), &n, &info);
  if(info != 0 || !result.is_finite()) {
    Rcpp::stop("The relationship matrix Cholesky inverse failed.");
  }
  result = arma::symmatu(result);
  return result;
}

}

// [[Rcpp::export]]
arma::mat amat(const arma::mat & Xo, const bool & vanraden, double minMAF) {
  const RelationshipMarkers markers(Xo, minMAF);
  return relationshipKernel(Xo, markers, true, vanraden);
}

// [[Rcpp::export]]
arma::mat dmat(const arma::mat & Xo, const bool & nishio, double minMAF) {
  const RelationshipMarkers markers(Xo, minMAF);
  return relationshipKernel(Xo, markers, false, nishio);
}

// [[Rcpp::export]]
arma::mat emat(const arma::mat & X1, const arma::mat & X2) {
  if(X1.n_rows != X2.n_rows || X1.n_cols != X2.n_cols) {
    Rcpp::stop("Epistatic relationship matrices must have identical dimensions.");
  }
  return X1 % X2;
}

// [[Rcpp::export]]
arma::mat emat_genotypes(const arma::mat & Xo, bool nishio,
                         const std::string & type, double minMAF) {
  if(type != "A#A" && type != "A#D" && type != "D#D") {
    Rcpp::stop("type must be one of 'A#A', 'A#D' or 'D#D'.");
  }
  const RelationshipMarkers markers(Xo, minMAF);
  arma::mat E = relationshipKernel(Xo, markers, type != "D#D", type != "D#D" || nishio);
  if(type == "A#D") {
    E %= relationshipKernel(Xo, markers, false, nishio);
  }else{
    E = arma::square(E);
  }
  return E;
}

// [[Rcpp::export]]
arma::mat hmat(const arma::mat & A, const arma::mat & G22,
               const arma::vec & index, double tolparinv,
               double tau, double omega) {
  validateRelationshipMatrix(A, "A");
  validateRelationshipMatrix(G22, "G22");
  if(index.n_elem != A.n_rows || !index.is_finite() ||
     arma::any((index != 0) % (index != 1))) {
    Rcpp::stop("index must contain one zero or one per individual in A.");
  }
  if(!std::isfinite(tolparinv) || tolparinv <= 0 ||
     !std::isfinite(tau) || !std::isfinite(omega)) {
    Rcpp::stop("tolparinv must be finite and positive; tau and omega must be finite.");
  }
  const arma::uvec index1 = arma::find(index == 1);
  const arma::uvec index2 = arma::find(index == 0);
  if(index2.n_elem != G22.n_rows) {
    Rcpp::stop("G22 dimensions must match the genotyped individuals in index.");
  }
  const arma::mat A22 = A.submat(index2, index2);
  arma::mat regularizedA22 = A22;
  const arma::mat factorA = relationshipCholesky(regularizedA22, tolparinv, "A22");
  arma::mat H22 = G22;
  const arma::mat factorG = relationshipCholesky(H22, tolparinv, "G22");
  if(omega == 1 && tau > 0) {
    // Preserve G's ridge retry, but avoid inverting G and then its inverse.
    H22 /= tau;
  }else{
    arma::mat precision = tau * relationshipInverse(factorG) +
      (1 - omega) * relationshipInverse(factorA);
    const arma::mat factorH = relationshipCholesky(precision, tolparinv, "H22 precision");
    H22 = relationshipInverse(factorH);
  }
  if(!H22.is_finite()) {
    Rcpp::stop("The hybrid genotyped relationship block is not finite.");
  }
  arma::mat H = A;
  if(!index1.is_empty()) {
    const arma::mat B = relationshipSolve(factorA, A.submat(index2, index1));
    const arma::mat cross = B.t() * (H22 - A22);
    arma::mat correction = cross * B;
    correction = arma::symmatu(correction);
    H.submat(index1, index1) += correction;
    H.submat(index1, index2) += cross;
    H.submat(index2, index1) += cross.t();
  }
  H.submat(index2, index2) = H22;
  if(!H.is_finite()) {
    Rcpp::stop("The hybrid relationship matrix is not finite.");
  }
  return H;
}

// [[Rcpp::export]]
arma::cube scorecalc(const arma::mat & Mimv,
                     const arma::mat & Ymv, // Y is provided as multitrait
                     const arma::mat & Zmv, // Z is provided as univariate
                     const arma::mat & Xmv, // X is provided as univariate
                     const arma::mat & Vinv, // multivariate inverse of V
                     int nt, double minMAF
) {

  double tolparinv = 0.00001;

  //
  arma::rowvec pf = mean(Mimv+1)/2; // allele frequency of p
  arma::rowvec qf = 1 - pf; // allele frequency of q
  double MAF = min(arma::join_rows(pf,qf)); // calculate MAF

  // start calculation
  double n = Ymv.n_rows;
  arma::mat ZMimv = Zmv * Mimv;
  arma::mat XZMimv = join_rows( Xmv, ZMimv);
  double p = XZMimv.n_cols;
  double v1 = 1;
  double v2 = n - p;
  // create Wi = (XZ Vi ZX)-1
  arma::mat Winv;
  arma::mat W = XZMimv.t() * (Vinv * XZMimv);
  arma::sp_mat D = arma::speye<arma::sp_mat>(W.n_cols,W.n_cols);
  arma::inv_sympd(Winv,W); // try to invert normally
  if(Winv.n_rows == 0){ // if fails try to invert with diag(1e-3)
    W = W + (D*tolparinv);
    arma::inv_sympd(Winv,W);
  }

  // Initialize result cube
  arma::cube result(Mimv.n_cols, nt, 3, arma::fill::zeros);

  if (Winv.n_rows > 0 && MAF > minMAF) {
    // Main calculations
    arma::mat XZMimvVy = XZMimv.t() * (Vinv * Ymv); // XZM' Vi y
    arma::colvec b = Winv * XZMimvVy;               // (XZ'VinvXZ)^-1 XZ'VinvY
    arma::colvec e = Ymv - (XZMimv * b);            // Residuals: Y - Xb
    arma::mat mVar = (e.t() * (Vinv * e)) / v2;     // Residual variance
    double mVarAsDouble = mVar(0, 0);
    arma::mat bVar = Winv * mVarAsDouble;           // Beta variance

    // Extract the right fixed effect
    arma::mat bn(b.n_rows, 1); // Empty vector
    for (int i = 0; i < b.n_rows; ++i) {
      bn(i) = i; // Fill it with their own position
    }
    arma::uvec ps = find(bn > (Xmv.n_cols - 1)); // Index for good markers
    arma::colvec bMarker = b(ps);               // Beta for marker
    arma::mat bMarkerVar = bVar(ps, ps);        // Variance for beta
    arma::vec SEMarker = arma::sqrt(diagvec(bMarkerVar));
    arma::vec fStat = arma::pow(bMarker / SEMarker, 2); // F statistic
    arma::vec x = v2 / (v2 + v1 * fStat);               // Probability scores

    // Fill the result cube
    result.slice(0) = arma::reshape(x, Mimv.n_cols, nt);         // Scores
    result.slice(1) = arma::reshape(bMarker, Mimv.n_cols, nt);  // Beta coefficients
    result.slice(2) = arma::reshape(SEMarker, Mimv.n_cols, nt); // Standard errors
  }

  return result;
}

// [[Rcpp::export]]
arma::cube gwasForLoop(const arma::mat & M, // marker matrix
                       const arma::mat & Y, // Y is provided as multitrait
                       const arma::mat & Z, // Z is provided as univariate
                       const arma::mat & X, // X is provided as univariate
                       const arma::mat & Vinv, // multivariate inverse of V
                       double minMAF,
                       bool display_progress=true
) {
  int nt = Y.n_cols;
  arma::mat Dnt = arma::eye<arma::mat>(nt,nt) ; // diagonal of nt dimensions
  // multivariate versions
  arma::mat Ymv = arma::vectorise(Y.t(),0); // multivariate Y
  // arma::mat Ymv = Ymvt.t();
  arma::mat Zmv = arma::kron(Z,Dnt); // kronecker
  arma::mat Xmv = arma::kron(X,Dnt); // kronecker

  // start calculation
  int n_marker = M.n_cols;
  arma::vec dummy(nt, arma::fill::ones);
  arma::uvec pos = arma::find(dummy > 0); // index for good markers
  arma::cube results(n_marker, nt, 3, arma::fill::zeros); // 3D array: markers x traits x (scores & bMarkers)

  // start for loop for each marker
  Progress p(n_marker, display_progress);
  for (int i = 0; i < n_marker; ++i) {
    if (Progress::check_abort() ){
      return arma::cube(0, 0, 0); // return an empty cube on abort
    }
    p.increment();
    arma::mat Mi = M.col(i); // extract marker i
    arma::mat Mimv = arma::kron(Mi,Dnt); // kronecker for multivariate
    arma::cube prov = scorecalc(Mimv, Ymv, Zmv, Xmv, Vinv, nt, minMAF);

    // Assign slices
    results.slice(0).row(i) = prov.slice(0).t(); // Scores
    results.slice(1).row(i) = prov.slice(1).t(); // bMarkers
    results.slice(2).row(i) = prov.slice(2).t(); // SEs
  }

  return results;
}
