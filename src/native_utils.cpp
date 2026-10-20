#include "mme_sp2_support.h"
#include "mme_progress.h"

// [[Rcpp::export]]
const std::string currentDateTime() {
  time_t     now = time(0);
  struct tm  tstruct;
  char       buf[80];
  tstruct = *localtime(&now);
  // Visit http://en.cppreference.com/w/cpp/chrono/c/strftime
  // for more information about date/time format
  strftime(buf, sizeof(buf), "%Y-%m-%d.%X", &tstruct);

  return buf;
}

// [[Rcpp::export]]
arma::vec seqCpp(const int & a,
                 const int & b){
  int c = b-a+1,i,counter;
  arma::vec d(c);
  counter = a;
  for(i=0; i < c; i++){
    d[i] = counter;
    counter++;
  }
  return d;
}

// [[Rcpp::export]]
arma::vec mat_to_vecCpp(const arma::mat & x,
                        const arma::mat & x2){
  // x is the matrix to be passed to a vector form in the output (out)
  // x2 is a mtrix of constraints to indicate wheter the value to be passed should be pass intefer (>0) or not passed (=0)
  int ncol = x.n_cols;
  arma::uvec nent2 = find(x2 > 0); int nent3 = nent2.n_elem;
  Rcpp::NumericVector out(nent3);
  // std::vector<bool> out2(nent3, true); // create position vector
  int counter = 0;
  int i, j;
  for (i = 0; i < ncol; i++){
    for (j = 0; j < ncol; j++){
      if (i > j){}else{
        // only extract the variance component if it was planned to be estimated
        if(x2(i,j) > 0){
          out[counter] = x(i,j);
          counter++;
        }
      }
    }
  }
  return out;
}

// [[Rcpp::export]]
arma::mat vec_to_matCpp(const arma::vec & x,
                        const arma::mat & x2){
  // x is the vector to be passed to a matrix form in the output (out)
  // x2 is a matrix of constraints to indicate wheter the value to be passed should be pass intefer (>0) or not passed (=0)
  int ncol = x2.n_cols;
  arma::uvec nent2 = find(x2 > 0);
  arma::mat out(ncol,ncol);
  //
  int counter = 0;
  int i, j;
  for (j = 0; j < ncol; j++){
    for (i = 0; i < ncol; i++){
      if (i > j){}else{
        // only extract the variance component if it was planned to be estimated
        if(x2(i,j) > 0){
          out(i,j) = x(counter);
          counter++;
        }
      }
    }
  }
  return out;
}

// [[Rcpp::export]]
arma::cube vec_to_cubeCpp(const arma::vec & x,
                          const Rcpp::List & g){

  int nge = g.size();
  arma::mat uuu = Rcpp::as<arma::mat>(g[0]);
  int ncols = uuu.n_cols;
  arma::cube Ge(ncols,ncols,nge); // copy GeI and we will replace the correct ones
  int i, j, k;
  int counter = 0;
  for(k = 0; k < nge; k++){ // FOR EACH RANDOM EFFECT access the matrix
    arma::mat x2 = Rcpp::as<arma::mat>(g[k]);
    int ncol = x2.n_cols;
    arma::mat x3(ncol,ncol);
    for (i = 0; i < ncol; i++){
      for (j = 0; j < ncol; j++){
        if (i > j){}else{//only upper triangular
          // only extract the variance component if it was planned to be estimated
          if(x2(i,j) > 0){
            x3(i,j) = x(counter);
            counter++;
          }else{x3(i,j) = 0;}
        }
      }
    }
    x3 = arma::symmatu(x3);
    Ge.slice(k) = x3;
  }
  return Ge;
}

// [[Rcpp::export]]
arma::vec varCols(const arma::mat & x){
  int nrow = x.n_rows, ncol = x.n_cols;
  Rcpp::NumericVector out(ncol);

  for (int j = 0; j < ncol; j++) {

    double mean = 0;
    double M2 = 0;
    int n;
    double delta, xx;

    for (int i = 0; i < nrow; i++) {
      n = i+1;
      xx = x(i,j);
      delta = xx - mean;
      mean += delta/n;
      M2 = M2 + delta*(xx-mean);
    }
    out(j) = M2/(n-1);
  }
  return out;
}

// [[Rcpp::export]]
arma::mat scaleCpp(const arma::mat & x) { // scale a matrix

  arma::vec sds = 1 / sqrt(varCols(x));
  arma::mat sdsD = diagmat(sds);

  int n = x.n_rows;
  arma::mat D = arma::eye<arma::mat>(n,n) ;
  arma::vec  v = arma::ones<arma::vec>(n);
  arma::mat vtv = ( v * v.t() ) / n;
  arma::mat xs = ((D - vtv) * x) * sdsD ;

  return xs;
}

// [[Rcpp::export]]
arma::mat makeFull(const arma::mat & X) {
  arma::mat U;
  arma::vec s;
  arma::mat V;
  arma::svd(U,s,V,X);

  int ncols0 = X.n_cols;
  arma::uvec s2 = find( s > 1e-8 );
  int ncols = s2.max() + 1;
  arma::uvec indices(ncols);
  for (int i = 0; i < ncols; ++i) {
    indices[i] = 1;
  }
  arma::mat Xf = U.cols(find(indices == 1));
  if(ncols == ncols0){
    arma::mat Xf = X;
  }
  return Xf;
}

// [[Rcpp::export]]
bool isIdentity_mat(const arma::mat x){
  int N = x.n_rows;
  for (int row = 0; row < N; row++){
    for (int col = 0; col < N; col++){
      if (row == col && x(row,col) != 1)
        return false;
      else if (row != col && x(row,col) != 0)
        return false;
    }
  }
  return true;
}

// [[Rcpp::export]]
bool isIdentity_spmat(const arma::sp_mat x){
  int N = x.n_rows;
  for (int row = 0; row < N; row++){
    for (int col = 0; col < N; col++){
      if (row == col && x(row,col) != 1)
        return false;
      else if (row != col && x(row,col) != 0)
        return false;
    }
  }
  return true;
}

// [[Rcpp::export]]
bool isDiagonal_mat(const arma::mat x){
  int N = x.n_rows;
  for (int row = 0; row < N; row++){
    for (int col = 0; col < N; col++){
      if (row != col && x(row,col) != 0)
        return false;
    }
  }
  return true;
}

// [[Rcpp::export]]
bool isDiagonal_spmat(const arma::sp_mat x){
  int N = x.n_rows;
  for (int row = 0; row < N; row++){
    for (int col = 0; col < N; col++){
      if (row != col && x(row,col) != 0)
        return false;
    }
  }
  return true;
}

// [[Rcpp::export]]
arma::sp_mat convertSparse(Rcpp::S4 mat) {
  // https://gallery.rcpp.org/articles/armadillo-sparse-matrix/
  // mat is an S4 (R) sparse matrix to be converted to an Armadillo sp_mat
  // obtain dim, i, p. x from S4 object
  Rcpp::IntegerVector dims = mat.slot("Dim");
  arma::urowvec i = Rcpp::as<arma::urowvec>(mat.slot("i"));
  arma::urowvec p = Rcpp::as<arma::urowvec>(mat.slot("p"));
  arma::vec x     = Rcpp::as<arma::vec>(mat.slot("x"));

  int nrow = dims[0], ncol = dims[1];

  // use Armadillo sparse matrix constructor
  arma::sp_mat res(i, p, x, nrow, ncol);
  // dsCMatrix stores one triangle only; downstream code expects both.
  if(mat.is("dsCMatrix")){
    arma::sp_mat d = arma::diagmat(res);
    res = res + res.t() - d;
  }
  return(res);
}

// [[Rcpp::export]]
arma::vec mat_to_vecCpp2(const arma::mat & x,
                         const arma::mat & x2){
  // x is the matrix to be passed to a vector form in the output (out)
  // x2 is a mtrix of constraints to indicate wheter the value to be passed should be pass intefer (>0) or not passed (=0)
  int ncol = x.n_cols;
  arma::uvec nent2 = find(x2 > 0); int nent3 = nent2.n_elem;
  Rcpp::NumericVector out(nent3);
  // std::vector<bool> out2(nent3, true); // create position vector
  int counter = 0;
  int i, j;
  for (j = 0; j < ncol; j++){
    for (i = 0; i < ncol; i++){
      if (i > j){}else{
        // only extract the variance component if it was planned to be estimated
        if(x2(i,j) > 0){
          out[counter] = x(i,j);
          counter++;
        }
      }
    }
  }
  return out;
}

// [[Rcpp::export]]
arma::mat nearPDcpp(const arma::mat X0,
                    const int & maxit,
                    const double & eig_tol,
                    const double & conv_tol){

  // Strict positive-definite version of nearPDcpp().
  // It first projects onto the positive-semidefinite cone and then
  // adds only the diagonal shift needed to guarantee strict PD.

  if(X0.n_rows != X0.n_cols){
    Rcpp::stop("nearPDcpp requires a square matrix.");
  }

  if(X0.n_rows == 0){
    return X0;
  }

  if(!X0.is_finite()){
    Rcpp::stop("nearPDcpp received a matrix containing non-finite values.");
  }

  const arma::uword n = X0.n_rows;

  arma::mat X = arma::symmatu(X0);

  double matrixScale = arma::norm(X, "inf");
  if(!std::isfinite(matrixScale) || matrixScale < 1.0){
    matrixScale = 1.0;
  }

  const double machineFloor =
    100.0 * std::numeric_limits<double>::epsilon() * matrixScale;

  // ai_mme_sp2 currently requires lambda_min > 1e-10.
  // Keep a small safety margin above that threshold.
  const double strictPdFloor =
    std::max(1.0e-9, machineFloor);

  // Scalar covariance structure.
  if(n == 1){
    if(X(0,0) < strictPdFloor){
      X(0,0) = strictPdFloor;
    }
    return X;
  }

  arma::mat D_S(n, n, arma::fill::zeros);
  arma::mat Y = X;

  const double safeEigTol = std::max(0.0, eig_tol);
  const double safeConvTol = std::max(0.0, conv_tol);

  for(int iter = 0; iter < maxit; ++iter){

    Y = X;

    arma::mat R =
      arma::symmatu(
        Y - D_S
      );

    arma::vec eigval;
    arma::mat eigvec;

    const bool eigOK =
      arma::eig_sym(
        eigval,
        eigvec,
        R
      );

    if(!eigOK || eigval.n_elem != n || !eigval.is_finite()){
      Rcpp::stop("nearPDcpp eigendecomposition failed.");
    }

    double eigScale = arma::abs(eigval).max();
    if(!std::isfinite(eigScale) || eigScale < 1.0){
      eigScale = 1.0;
    }

    const double psdTol =
      safeEigTol * eigScale;

    arma::vec projectedEig = eigval;

    for(arma::uword j = 0; j < projectedEig.n_elem; ++j){
      if(projectedEig(j) <= psdTol){
        projectedEig(j) = 0.0;
      }
    }

    X =
      eigvec
      *
      arma::diagmat(projectedEig)
      *
      eigvec.t();

    X = arma::symmatu(X);

    D_S =
      X - R;

    const double denom =
      std::max(
        arma::norm(Y, "inf"),
        1.0
      );

    const double conv =
      arma::norm(Y - X, "inf")
      /
      denom;

    if(std::isfinite(conv) && conv <= safeConvTol){
      break;
    }
  }

  // Final strict-PD safeguard.
  X = arma::symmatu(X);

  arma::vec finalEig;
  const bool finalEigOK =
    arma::eig_sym(
      finalEig,
      X
    );

  if(!finalEigOK || finalEig.n_elem != n || !finalEig.is_finite()){
    Rcpp::stop("nearPDcpp final eigendecomposition failed.");
  }

  const double minEig =
    finalEig.min();

  if(minEig < strictPdFloor){

    const double shift =
      (strictPdFloor - minEig)
      +
      machineFloor;

    X.diag() +=
      shift;

    X = arma::symmatu(X);
  }

  // Defensive verification using Cholesky.
  arma::mat cholFactor;
  bool cholOK =
    arma::chol(
      cholFactor,
      X,
      "lower"
    );

  if(!cholOK){

    double jitter =
      std::max(
        strictPdFloor,
        machineFloor
      );

    for(int attempt = 0; attempt < 6 && !cholOK; ++attempt){

      X.diag() +=
        jitter;

      X = arma::symmatu(X);

      cholOK =
        arma::chol(
          cholFactor,
          X,
          "lower"
        );

      jitter *=
        10.0;
    }
  }

  if(!cholOK){
    Rcpp::stop(
      "nearPDcpp was unable to obtain a strictly positive-definite matrix."
    );
  }

  return X;
}
