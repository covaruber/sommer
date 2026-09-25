// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// RcppArmadillo pulls in Rcpp.h for us
#define ARMA_DONT_PRINT_ERRORS

#include <RcppArmadillo.h>
#include "stdlib.h"

#include <progress.hpp>

// Eigen sparse LDLT used by ai_mme_sp()
#include <Eigen/SparseCore>
#include <Eigen/SparseCholesky>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/Eigenvalues>
#include <cctype>

// Nested-dissection sparse ordering for ai_mme_sp2()'s LDLT factorisations,
// enabled only when configure detected a usable METIS installation; falls
// back to Eigen's built-in AMD ordering otherwise.
#ifdef SOMMER_HAVE_METIS
#include <Eigen/MetisSupport>
typedef Eigen::MetisOrdering<int> SommerSparseOrdering;
#else
typedef Eigen::AMDOrdering<int> SommerSparseOrdering;
#endif

// Supernodal (BLAS-3) sparse Cholesky backend for ai_mme_sp2()'s
// solver="cholmod" option, provided by R's Matrix package. Matrix.h declares
// the M_cholmod_*() wrappers; stubs.c defines them via lazy R_GetCCallable()
// lookups into the already-loaded Matrix package, so no extra link flags
// are required beyond LinkingTo: Matrix.
#include <Matrix/Matrix.h>
#include <Matrix/stubs.c>

// Standard C++ headers used by the new implementation
#include <vector>
#include <limits>
#include <cmath>
#include <algorithm>
#include <atomic>
#include <cstring>
#include <functional>

#ifdef _OPENMP
#include <omp.h>
#endif

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

// simple example of creating two matrices and
// returning the result of an operatioon on them
//
// via the exports attribute we tell Rcpp to make this function
// available from R
//

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
arma::mat amat(const arma::mat & Xo, const bool & vanraden, double minMAF) {
  
  // remove min.MAF
  arma::rowvec pfreq = mean(Xo+1,0)/2; // frequency of p
  arma::mat pqfreq = arma::join_cols(pfreq,1-pfreq); // frequencies of p and q
  arma::rowvec MAF = min(pqfreq,0); // minor allele freqs
  arma::uvec indexMAF = find(MAF > minMAF); // index for good markers > minMAF
  arma::mat Xo2 = Xo.cols(indexMAF); // new X only with polymorphic markers
  
  // remove monomorphic markers
  arma::rowvec xVar = var(Xo2,0); // column variance
  arma::uvec index = find(xVar > 0); // index for good markers
  arma::mat X = Xo2.cols(index); // new X only with polymorphic markers
  
  // initialize A
  int p = X.n_cols;// number of markers
  int n = X.n_rows;
  arma::mat A(n,n);
  
  if(vanraden == true){ //  regular vanRaden 
    
    arma::rowvec ms012 = mean( X+1, 0 ); // means of columns
    arma::rowvec freq = ms012/2;
    double v = 2 * mean(freq % (1 - freq));
    
    arma::mat one(n, 1, arma::fill::ones);
    arma::mat freqmat = one * freq;
    arma::mat W = (X + 1) - (2 * freqmat);
    //
    arma::mat K = W * W.t();
    A = K/v/p;
    
  }else{ // Endelman (currently we have a bug here)
    
    // IN R: M <- scale(X, center = TRUE, scale = FALSE)
    arma::rowvec ms = mean( X, 0 ); // means of columns
    arma::mat M = X.each_row() - ms;
    // IN R: tcrossprod(M)
    arma::mat K = M * M.t();
    // IN R: K/mean(diag(K))   mean(K.diag())
    double v = mean(diagvec(K));
    A = K/v;
    
  }
  
  return A;
}

// [[Rcpp::export]]
arma::mat dmat(const arma::mat & Xo, const bool & nishio, double minMAF) {
  
  // remove min.MAF
  arma::rowvec pfreq = mean(Xo+1,0)/2; // frequency of p
  arma::mat pqfreq = arma::join_cols(pfreq,1-pfreq); // frequencies of p and q
  arma::rowvec MAF = min(pqfreq,0); // minor allele freqs
  arma::uvec indexMAF = find(MAF > minMAF); // index for good markers > minMAF
  arma::mat Xo2 = Xo.cols(indexMAF); // new X only with polymorphic markers
  
  // remove monomorphic markers
  arma::rowvec xVar = var(Xo2,0); // column variance
  arma::uvec index = find(xVar > 0); // index for good markers
  arma::mat X = Xo2.cols(index); // new X only with polymorphic markers
  
  // initialize A
  // int p = X.n_cols;
  int n = X.n_rows;
  arma::mat D(n,n);
  
  arma::mat Xd = 1 - abs(X);
  
  if(nishio == true){ //  Nishio ans Satoh. (2014)
    
    // IN R: M <- scale(Xd, center = TRUE, scale = FALSE)
    arma::rowvec ms = mean( Xd, 0 ); // means of columns
    arma::mat M = Xd.each_row() - ms; // centered Xd matrix
    // IN R: bAlleleFrequency <- colMeans(X+1)/2; 0-1-2
    arma::rowvec bAlleleFrequency = mean( X+1, 0 )/2; // means of columns
    // IN R: varHW <- sum((2 * bAlleleFrequency * (1 - bAlleleFrequency))^2)
    double varHW = arma::accu(arma::square(2 * bAlleleFrequency % (1 - bAlleleFrequency)));
    // IN R: tcrossprod(M)
    arma::mat K = M * M.t();
    //
    D = K/varHW;
    
  }else{ // Su et al. (2012)
    
    // IN R: M <- scale(X, center = TRUE, scale = FALSE)
    // arma::rowvec ms = mean( Xd, 0 ); // means of columns
    // arma::mat M = Xd.each_row() - ms; // centered Xd matrix
    // IN R: p <- colSums(X+1)/(2*n) # from marker marix in 0,1,2 format
    arma::rowvec p = sum( X+1, 0 )/(2*n); // means of columns
    arma::rowvec q = 1-p;
    // IN R: varHW <- sum(2*p*q * (1-(2*p*q)) )
    arma::rowvec p2q = 2*(p%q);
    double varHW = arma::accu( p2q % (1-p2q) );
    // IN R: Xdpq <- apply(Xd, 1, function(x){ x - (2 * p * q)})
    arma::mat M = Xd.each_row() - p2q;
    // IN R: tcrossprod(M)
    arma::mat K = M * M.t();
    D = K/varHW;
    
  }
  
  return D;
}

// [[Rcpp::export]]
arma::mat emat(const arma::mat & X1, const arma::mat & X2) {
  
  arma::mat E = X1 % X2;
  
  return E;
}

// [[Rcpp::export]]
arma::mat hmat(const arma::mat & A, const arma::mat & G22,
               const arma::vec & index, double tolparinv,
               double tau, double omega) {
  
  arma::uvec index1 = find(index == true); // index for good markers
  arma::uvec index2 = find(index == false); // index for good markers
  // A11 <- A[index, index]
  arma::mat A11 = A.submat(index1,index1);
  // A12 <- A[index, !index]
  arma::mat A12 = A.submat(index1,index2);
  // A21 <- A[!index, index]
  arma::mat A21 = A.submat(index2,index1);
  // A22 <- A[!index, !index]
  arma::mat A22 = A.submat(index2,index2);
  // A22inv = solve(A22)
  arma::mat A22inv(A22.n_cols,A22.n_cols);
  arma::inv_sympd(A22inv,A22); // try to invert normally
  arma::sp_mat Ia = arma::speye<arma::sp_mat>(A22.n_cols,A22.n_cols);
  if(A22inv.n_rows == 0){ // if fails try to invert with diag(1e-3)
    arma::mat A22b = A22 + (Ia*tolparinv);
    arma::inv_sympd(A22inv,A22b);
  }
  // G22inv = try(solve(G22), silent = TRUE)
  arma::mat G22inv(G22.n_cols,G22.n_cols);
  arma::inv_sympd(G22inv,G22); // try to invert normally
  arma::sp_mat Ig = arma::speye<arma::sp_mat>(G22.n_cols,G22.n_cols);
  if(G22inv.n_rows == 0){ // if fails try to invert with diag(1e-3)
    arma::mat G22b = G22 + (Ig*tolparinv);
    arma::inv_sympd(G22inv,G22b);
  }
  //   H22 = solve((tau * G22inv + (1 - omega) * A22inv))
  arma::mat H22p = (tau * G22inv) + ((1 - omega) * A22inv); //constant by matrix product
  arma::mat H22inv(H22p.n_cols,H22p.n_cols);
  arma::inv_sympd(H22inv,H22p); // try to invert normally
  arma::sp_mat Ih = arma::speye<arma::sp_mat>(H22p.n_cols,H22p.n_cols);
  if(H22inv.n_rows == 0){ // if fails try to invert with diag(1e-3)
    arma::mat H22pb = H22p + (Ih*tolparinv);
    arma::inv_sympd(H22inv,H22pb);
  }
  //   H11 = A12 %*% A22inv %*% (H22 - A22) %*% A22inv %*% A21
  arma::mat H11 = A12 * A22inv * (H22inv - A22) * A22inv * A21;
  //   H12 = A12 %*% A22inv %*% (H22 - A22)
  arma::mat H12 = A12 * A22inv * (H22inv - A22);
  //   H21 = (H22 - A22) %*% A22inv %*% A21
  arma::mat H21 = (H22inv - A22) * A22inv * A21;
  //   H22 = (H22 - A22)
  arma::mat H22 = (H22inv - A22);
  //   H = A + cbind(rbind(H11, H21), rbind(H12, H22))
  arma::mat H = A + join_rows(join_cols(H11, H21), join_cols(H12, H22));
  
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

// [[Rcpp::depends(RcppProgress)]]
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

// [[Rcpp::export]]
Rcpp::List newton_di_sp(const arma::sp_mat & Y, const Rcpp::List & X,
                        const Rcpp::List & Gx,
                        const Rcpp::List & Z, const Rcpp::List & K,
                        const Rcpp::List & R, 
                        const Rcpp::List & Ge, const Rcpp::List & GeI, // theta and thetaC
                        const arma::sp_mat & W, const bool & isInvW,
                        int iters, double tolpar, double tolparinv,
                        const bool & ai, const bool & pev,
                        const bool & verbose,const bool & retscaled,
                        const arma::vec & stepweight, // const arma::vec & emupdate,
                        const arma::vec & emweight, const Rcpp::List & thetaConstOri,
                        const arma::vec & thetaIndex) {
  
  time_t before = time(0);
  localtime(&before);
  
  int n_fixed = X.size(); // define nre=number of fixed effects
  int n_random = Z.size(); // define nre=number of random effects
  int n_rcov = R.size(); // define nre=number of residual effects
  int n_re = n_random + n_rcov; // define nre=number of total random effects z+r
  int n_traits = Y.n_cols; // define n_traits=number of traits
  int no = Y.n_rows; // define n_traits=number of traits
  arma::vec n_levels(n_re, arma::fill::ones); // to store the number of columns each Z and R matrix has
  arma::mat diagTrait = arma::eye(n_traits,n_traits);
  // ****************************************************
  // define ZKZ' and R
  // ****************************************************
  // calculate and concatenate ZKZ' and R
  arma::cube ZKZtR(no,no,n_re);
  
  for (int i = 0; i < n_re; ++i) { // for each random effect
    int irw = i - n_random;
    if(i < n_random && n_random > 0){ // if random effect (not residual)
      
      arma::sp_mat zp = Rcpp::as<arma::sp_mat>(Z[i]); // transform as sparse
      n_levels(i) = zp.n_cols; // store the number of columns or levels for this random effect
      
      bool dcheck = isIdentity_mat(arma::mat(Rcpp::as<arma::sp_mat>(K[i])));
      if(dcheck == true){ // if K[i] is diagonal
        if(zp.n_rows == zp.n_cols){//is a square matrix
          bool dcheck2 = isIdentity_spmat(zp);
          if(dcheck2 == true){ // if Z[i] is diagonal
            ZKZtR.slice(i) = arma::mat(Rcpp::as<arma::sp_mat>(K[i]));
          }else{ZKZtR.slice(i) = zp * zp.t(); }
        }else{ // is a rectangular matrix
          ZKZtR.slice(i) = zp * zp.t();
        }
      }else{ // if K[i] is not diagonal
        if(zp.n_rows == zp.n_cols){//is a square matrix
          bool dcheck2 = isIdentity_spmat(zp);
          if(dcheck2 == true){ // if Z[i] is diagonal
            ZKZtR.slice(i) = arma::mat(Rcpp::as<arma::sp_mat>(K[i]));
          }else{ZKZtR.slice(i) = zp * arma::mat(Rcpp::as<arma::sp_mat>(K[i])) * zp.t(); }
        }else{
          ZKZtR.slice(i) = zp * arma::mat(Rcpp::as<arma::sp_mat>(K[i])) * zp.t();
        }
      }
      
    }else{//if is an rcov term
      // bool dcheck3 = isIdentity_mat(W);
      double dcheck3 = accu(W) - W.n_cols;
      if(dcheck3 == 0){ // if W is diagonal no need to multiply
        ZKZtR.slice(i) = Rcpp::as<arma::sp_mat>(R[irw]);
      }else{ // if W (weights) is not diagonal then multiply Wis R Wis
        // arma::vec ws = W.diag();// 1 / sqrt(diagvec(W));
        if(isInvW == true){ // user has provided a squared and inverted W already
          ZKZtR.slice(i) = W * Rcpp::as<arma::sp_mat>(R[irw]) * W;
        }else{ // user has provided only W
          arma::mat Wis = inv(chol(arma::mat(W)));
          ZKZtR.slice(i) = Wis * Rcpp::as<arma::sp_mat>(R[irw]) * Wis.t();
        }
        // arma::vec ws2 = 1/sqrt(ws);// arma::mat Wis = diagmat(ws2); // W inverse squared  // ZKZtR.slice(i) = Wis * Rcpp::as<arma::sp_mat>(R[irw]) * Wis;
      }
      n_levels(i) = ZKZtR.slice(i).n_cols; // store the number of columns or levels for this random effect
    }
  }
  // ****************************************************
  // build multivariate versions of X and Y
  // ****************************************************
  arma::vec Ym = vectorise(arma::mat(Y)); // multivariate Y in original scale
  int nom = Ym.n_rows; // number of observations on the vector-form of multivariate Y
  
  arma::mat Xm;
  for (int i = 0; i < n_fixed; ++i) { // for each fixed effect
    if(i==0){ // build multivariate X for 1st fixed effect
      Xm = kron(Rcpp::as<arma::mat>(Gx[i]), arma::mat(Rcpp::as<arma::sp_mat>(X[i])) );
    }else{ // build multivariate X for 2nd to nth fixed effect and column bind them
      Xm = arma::join_horiz( Xm , kron(Rcpp::as<arma::mat>(Gx[i]), arma::mat(Rcpp::as<arma::sp_mat>(X[i])) ) );
    }
  }
  arma::mat Ys = scaleCpp(arma::mat(Y)); // scale Y using the scaleCpp function made
  arma::vec Ysm = vectorise(Ys); // multivariate Y in scaled form
  // ****************************************************
  // initial VC
  // ****************************************************
  arma::mat base_var = cov(arma::mat(Y)); // matrix of original variance-covariance in responses
  arma::mat sc_var = cov(Ys); // matrix of scaled variance-covariance in responses
  int rankX = Xm.n_rows - rank(Xm); // n - p.x
  // VC matrix with dimensions n_traits x n_traits (sigma)
  // we need one for each random effect (n_re)
  arma::cube sigma(n_traits,n_traits,n_re);
  arma::cube sigma_scaled(n_traits,n_traits,n_re);
  
  arma::field<arma::vec> sigma_ut(n_re); // undefined LIST to store the VC in a vector-form with length n_re (#of random effects)
  arma::field<arma::vec> constraintsL(n_re); // undefined LIST to store the constraints in a vector-form with length n_re (#of random effects)
  arma::field<arma::vec> n_levels_multi_traitL(n_re); // undefined LIST to store the n_levels in a vector form
  int no_vc = 0; // to add and find out how many VC exist in total
  for (int i = 0; i < n_re; ++i) { // for each random effect fill the cube
    sigma.slice(i) = Rcpp::as<arma::mat>(Ge[i]); // take Ge for a random effect (initial VC values) and save them in a slice
    arma::vec oo = mat_to_vecCpp(sigma.slice(i),GeI[i]) ; // extract upper triangular from that slice in a vector form, pass the constraints as 2nd argument
    sigma_ut[i] = oo; // oo is sigma2 in vector form and stored in the list sigma_ut
    constraintsL[i] = mat_to_vecCpp(GeI[i],GeI[i]) ; // who are diagonal and non-diagonal VCs, pass constraints in list form
    n_levels_multi_traitL[i] = constraintsL[i] ;
    no_vc = no_vc + oo.n_elem; // keep adding the #of VC
  }
  // sigma_ut_un will have all VC for all random effects in a single vector
  arma::vec sigma_ut_un; // vector to unlist the LIST of VC for all random effects
  arma::vec constraints; // vector to unlist constraints
  arma::vec n_levels_multi_trait; // vector to unlist constraints
  for(int i=0; i < n_re ; i++){ // for each random effect unlist
    sigma_ut_un = join_cols(sigma_ut_un,sigma_ut[i]); // column bind vectors so we end up with a very long vector with all VC
    constraints = join_cols(constraints,constraintsL[i]); // column bind vectors so we end up with a very long vector with all constraints
    arma::vec provX = n_levels_multi_traitL[i];
    arma::vec hpos = provX;
    for(int h=0; h < provX.n_cols ; h++){ // for each random effect unlist
      hpos(h) = n_levels(i);
    }
    n_levels_multi_trait = join_cols(n_levels_multi_trait,(provX/provX) % hpos);
  }
  arma::vec sigmaF_ut_un = sigma_ut_un; // make a copy for fixed-value vc's when we use constraints
  arma::vec coef_ut_un = sigma_ut_un; // make a 2nd copy of the same vector for stabilization
  arma::vec coef_ut_un_explode = sigma_ut_un; // make a 3rd copy of the same vector for checking issues with vc going too early outside the parameter space
  
  int  kk = sigma_ut_un.n_elem; // how many VCs are in the model?
  arma::vec llstore(iters); // container for LL
  arma::vec pos(sigma_ut_un.n_elem, arma::fill::zeros); // create an index vector with as many 0's as VCs
  
  // ****************************************************
  // dummy matrices for multivariate derivatives
  // ****************************************************
  int tot = n_re*n_traits*n_traits; // maximum number of variance components
  arma::vec re_mapper(tot); // mapper to know which VC belongs to each random effect
  arma::cube deriv_dummy(n_traits,n_traits,tot);
  int counter3 = 0;
  for(int i=0; i < n_re; i++){ // for each random effect
    arma::mat prov = Rcpp::as<arma::mat>(GeI[i]);
    int ncol = prov.n_cols; // traits
    
    for (int k = 0; k < ncol; k++){ // go through GeI(i) and make the dummy derivatives where there's a value > 0
      for (int j = 0; j < ncol; j++){
        if (k > j){}else{
          // only extract the variance component if it was planned to be estimated
          if(prov(k,j) > 0){
            arma::mat prov4(ncol,ncol,arma::fill::zeros);
            prov4(k,j)=1;
            prov4 = arma::symmatu(prov4);
            deriv_dummy.slice(counter3) = prov4;
            re_mapper[counter3] = i;
            counter3++;
          }
        }
      }
    }
    
  }
  // ****************************************************
  // ****************************************************
  // ##### iterative algorithm starts
  // ****************************************************
  // ****************************************************
  // Rcpp::List PdViList(kk); // list to store the multivariate derivatives * P or PVi=P*dZKZ'/ds
  
  arma::vec v(nom, arma::fill::ones); // generate enough ones for an identity matrix of dimensions nt x nt
  arma::mat Vi(nom,nom); // V or phenotypic variance matrix
  arma::mat P(nom,nom); // to fill the projection matrix
  arma::sp_mat D = arma::speye<arma::sp_mat>(nom,nom);
  arma::vec seqrankX = seqCpp(0,rankX-1); // will be used to keep only the eigen values for indices 1 to rankX
  arma::vec seqkk = seqCpp(0,kk-1);
  arma::vec popo = arma::vec(rankX, arma::fill::zeros);
  for(int i=0; i < rankX; i++){popo(i) = 1;}
  arma::mat Inf(kk,kk,arma::fill::zeros); // to store second derivatives (information matrix)
  arma::mat InfEM(kk,kk,arma::fill::zeros); // to store second derivatives (information matrix)
  arma::mat InfJoin(kk,kk,arma::fill::zeros); // to store second derivatives (information matrix)
  arma::mat InfJoin_inv(kk,kk,arma::fill::zeros); // to store second derivatives (information matrix)
  
  arma::mat Infw(kk,kk,arma::fill::zeros); // weights for AI information matrix
  arma::mat InfEMw(kk,kk,arma::fill::zeros); // weights for EM information matrix
  
  arma::vec score(kk); // vector to store first derivatives, the product Y'PViPY - tr(PVi) = dL/ds
  arma::mat Inf_inv; // to store the inverse of the information matrix
  arma::vec eigval2; // will be used for the decomposition of P, within the algorithm
  arma::mat eigvec2; // will be used for the decomposition of P
  arma::mat sigma_store(sigma_ut_un.n_elem,iters); // to store variance comp through the different iterations
  arma::mat sigma_perc_change(sigma_ut_un.n_elem,iters); // to store percent change of variance components
  arma::mat llik_store(1,iters); // to store llik through the different iterations
  
  arma::mat beta, fitted, residuals; // empty matrices for ..
  Rcpp::List VarU(n_random); // list object for the BLUP variances
  Rcpp::List PevU(n_random); // list object for the BLUP PEVs
  Rcpp::List U(n_random); // list object for the BLUPs
  
  arma::vec vdD(n_traits,arma::fill::ones);
  arma::mat dD = arma::diagmat(vdD);
  arma::mat sigma_cov;
  arma::mat tXVXi; // var-cov fixed effects
  arma::mat V(nom,nom); // V or phenotypic variance matrix
  
  bool convergence = false;
  bool last_iteration = false;
  int cycle, cycle2, ikk;
  double ldet, llik, llik0, delta_llik, checkP, seconds; // to store likelihoods and determinants
  // ###############
  // LOOP for cycles
  // ###############
  for(cycle=0; cycle < iters; cycle++){ // for each cycle
    
    for (int i = 0; i < n_re; ++i) {  // for each random effect in the formula
      sigma_ut[i] = mat_to_vecCpp(sigma.slice(i),Rcpp::as<arma::mat>(GeI[i])) ; // extract upper triangular in a vector form
    } // sigma_ut is a LIST
    arma::vec sigmatwo; // create a vector for variance components
    for(int i=0; i < n_re ; i++){ // for each random effect
      sigmatwo = join_cols(sigmatwo,sigma_ut[i]); // column bind to make a vector of vectors
    } // sigmatwo now has all VCs in a vector
    
    // multivariate ZKZ' and V
    
    int i;
    for(i=0; i < n_re; i++){ // loop for filling the multivariate ZGZ' and V
      // listGs.slice(i) = prov;
      if(last_iteration == true){ // if is the last iteration multivariate ZKZ is opposite
        if(i == 0){
          V = arma::kron(ZKZtR.slice(i),sigma.slice(i));
        }else{V = V + arma::kron(ZKZtR.slice(i),sigma.slice(i));}
      }else{
        if(i == 0){
          V = arma::kron(sigma.slice(i),ZKZtR.slice(i));
        }else{V = V + arma::kron(sigma.slice(i),ZKZtR.slice(i));}
      }
    }
    // invert V and P (projection matrix)
    
    arma::inv_sympd(Vi,V); // try to invert normally
    if(Vi.n_rows == 0){ // if fails try to invert with diag(1e-3)
      V = V + (D*tolparinv);
      arma::inv_sympd(Vi,V);
      if(Vi.n_rows == 0){// if fails try to invert with diag(1e-2)
        V = V + (D*(tolparinv*10));
        arma::inv_sympd(Vi,V);
        if(Vi.n_rows == 0){ // if fails try to invert with diag(1e-1)
          V = V + (D*(tolparinv*100));
          arma::inv_sympd(Vi,V);
          if(Vi.n_rows == 0){ // finally, if fails try to invert with diag(1e-3)
            // Rcpp::Rcout << "System is singular (V). Stopping the job. Try a bigger number of tolParInv." << arma::endl;
            Rcpp::stop("System is singular (V). Aborting the job. Try a bigger number of tolParInv.");
            // return 0;
          }
        }
      }
    }
    // if last iteration let's make Xm in the opposite direction
    if(last_iteration == true){
      for (int i = 0; i < n_fixed; ++i) {
        if(i==0){
          Xm = kron(arma::mat(Rcpp::as<arma::sp_mat>(X[i])), Rcpp::as<arma::mat>(Gx[i]) );
        }else{
          Xm = arma::join_horiz( Xm , arma::kron( arma::mat(Rcpp::as<arma::sp_mat>(X[i])), Rcpp::as<arma::mat>(Gx[i]) ) );
        }
      }
      Ym = vectorise(Y.t());
      Ysm = vectorise(Ys.t());
    }
    arma::mat VX = Vi * Xm; // VX
    arma::mat tXVX = Xm.t() * VX; // X'VX
    
    arma::mat tXVXVX; // X'VXVX
    tXVXVX = arma::solve(tXVX, VX.t()); // X'VXVX
    arma::solve(tXVXVX,tXVX,VX.t());
    if(tXVXVX.n_rows == 0){ // if fails try to invert with diag(1e-6)
      arma::solve(tXVXVX,tXVX + (D*(tolparinv)),VX.t());
      if(tXVXVX.n_rows == 0){// if fails try to invert with diag(1e-5)
        arma::solve(tXVXVX,tXVX + (D*(tolparinv*10)),VX.t());
        if(tXVXVX.n_rows == 0){ // if fails try to invert with diag(1e-4)
          arma::solve(tXVXVX,tXVX + (D*(tolparinv*100)),VX.t());
          if(tXVXVX.n_rows == 0){ // finally stop
            // Rcpp::Rcout << "System is singular (tXVXVX). Aborting the job. Try a bigger number of tolParInv." << arma::endl;
            Rcpp::stop("System is singular (tXVXVX). Aborting the job. Try a bigger number of tolParInv.");
            // return 0;
          }
        }
      }
    }
    
    // projection matrix
    P = Vi - (VX*tXVXVX); // V - V(XVX)-V
    
    if(last_iteration == false){
      
      arma::vec rss = Ysm.t() * (P * Ysm); // yPy = scalar RSS
      
      double rankXorss = arma::as_scalar(rankX/rss); // (n-p)/y'Py
      double rssorankX = arma::as_scalar(rss/rankX); // y'Py/(n-p)
      
      sigmatwo = sigmatwo * rssorankX;
      
      // weight the projection matrix to provide stability
      coef_ut_un(arma::find(pos == 0)) =  sigmatwo(arma::find(pos == 0)); // VC1[which(pos==0)] = VC2[which(pos==0)]
      coef_ut_un(arma::find(pos == 1)) = log(sigmatwo(arma::find(pos == 1))); // VC1[which(pos==1)] = log(VC2[which(pos==1)])
      
      // calculate the log-likelihood
      P = P * rankXorss; // P * [(n-p)/y'Py]
      rss = rankX; // yPy = n-p
      arma::eig_sym(eigval2, eigvec2, P); // VlV
      eigval2 = sort(eigval2,"descend"); // sort eigen vectors
      eigval2 = eigval2(arma::find(popo == 1));//(find(seqrankX < rankX)); // only take the values from 1 to
      checkP = eigval2.min();
      if(checkP < 0){ // if any eigen value is < 0 recalculate P
        P = P + (D * (tolpar - eigval2.min())) ;
        eigval2 = eigval2 + tolpar - eigval2.min();
      }
      ldet = accu(log(eigval2)); // sum(log(lambda))
      llik = ldet/2 - (arma::as_scalar(rss)/2); // llik = [sum(log(lambda))/2] - [(n-p)/2]
      
      if(cycle == 0){llik0 = llik;}
      delta_llik = llik - llik0;
      llik0 = llik;
      
      // use the stabilization
      arma::vec var_components(kk, arma::fill::ones); // VC = rep(0,nVC)
      double check00 = accu(pos); // accu is like sum() in R
      if(check00 > 0){  // if there's 1's in the pos vector
        arma::uvec ind = find(pos == 1); // which are 1's
        var_components(ind) = sigmatwo(ind); // var_components[which(pos==1)] = sigmatwo[which(pos==1)]
      }
      
      // calculate first derivatives (dL/ds = score)
      
      arma::cube PdViList(nom,nom,kk); // list to store the multivariate derivatives * P or PVi=P*dZKZ'/ds
      for(int i=0; i < kk; i++){
        int re = re_mapper(i);
        arma::mat zkzp = ZKZtR.slice(re); // it repeats the same ZKZtR if is a vc for the same random effect
        arma::mat PdVi = P * kron(deriv_dummy.slice(i),zkzp); // multivariate dVi = dZKZ'/ds
        if(ai && cycle > 2){
          score[i] = - (0.5 * arma::as_scalar(trace(PdVi))) + (0.5 * arma::as_scalar((Ysm.t() * PdVi * P * Ysm)));
        }else{
          score[i] = arma::as_scalar(Ysm.t() * PdVi * P * Ysm) - accu(diagvec(PdVi));
        }
        PdViList.slice(i) = PdVi;
      }
      // theta(k) * dL/ds  ..... are scalar values
      score = score % var_components; // to be used later for updating the variance components
      // if all goes well var_components is just ones
      
      // calculate second derivatives (AverageInformation)
      // Fisher's Information tr(PVi * PVi) .... A*=Vi=dV/ds .... [Vi Vj'] si sj ; TT is the list of derivatives for all random effects - trait combos
      
      // if(emupdate(cycle) == 0){ // if user wants an EM update (1st derivatives) . It works but it didn't speed up the algorithm when using EM. This leads to don't have information matrix and therefore SE for variance components.
      for (int i = 0; i < kk; i++){
        for (int j = 0; j < kk; j++){
          if (i > j){}else{//only upper triangular
            if(ai && cycle > 2){ // if average information
              Inf(i,j) = 0.5 * arma::as_scalar(Ysm.t() * PdViList.slice(i) * P * PdViList.slice(j) * P * (P * Ysm)); // j is .t() ?
            }else{ // if newton raphson
              Inf(i,j) = accu(PdViList.slice(i) % PdViList.slice(j).t()) * arma::as_scalar(var_components(i)) * arma::as_scalar(var_components(j));
            }
          }
        }
      }
      Inf = arma::symmatu(Inf); // copy lower in upper triangular
      Inf_inv = arma::pinv(Inf, 1.490116e-08); // Inverse of Fishers or information matrix
      
      if(Inf_inv.n_rows == 0){ // if fails
        // Rcpp::Rcout << "System is singular (Inf_inv). Aborting the job. Try a bigger number of tolParInv." << arma::endl;
        // return 0;
        Rcpp::stop("System is singular (Inf_inv). Aborting the job. Try a bigger number of tolParInv.");
      }
      // }
      
      // vector to store the update = F- * sigma(k) * dL/ds
      arma::vec delta(kk);
      
      // if(emupdate(cycle) == 1){ // if user wants an EM update (1st derivatives)
      InfEM.diag() = (coef_ut_un % coef_ut_un) / n_levels_multi_trait;  // I.em inverse
      InfEM = arma::pinv( InfEM ,  1.490116e-08 ); // I.em
      arma::vec emw(kk); // vectors for weights
      arma::vec aiw(kk);
      for(ikk=0; ikk < kk; ikk++){
        emw(ikk)= emweight(cycle);
        aiw(ikk)= 1 - emweight(cycle);
      }
      Infw.diag() = aiw;  // put weights in diagonal fill::value is still not available in this version
      InfEMw.diag() = emw; //
      InfJoin = (Inf*Infw)+(InfEM*InfEMw); // joint information matrix
      InfJoin_inv = arma::pinv(InfJoin, 1.490116e-08); // inverse the joint information matrix
      delta = InfJoin_inv * score; //update for variance components where: delta = Information.inv * dL/ds
      // delta = (coef_ut_un % score % coef_ut_un)/n_levels; // previous way I was calculating the deltas
      // }else{ // if user wants an information*score update
      //   delta = Inf_inv * score; //update for variance components where: delta = Information.inv * dL/ds
      // }
      
      // ^^^^^^^^^^^^^^^^^^
      // ^^^^^^^^^^^^^^^^^^
      // parameter restrain
      // GeI values
      // 0 not estimated
      // 1 positive
      // 2 unconstrained
      // 3 fixed
      arma::vec coef_ut_unC = coef_ut_un + (stepweight(cycle) * delta); // provisional new variance components
      arma::uvec restrain = find(constraints == 1 && coef_ut_unC < 0); // which vcs are negative and should be positive
      arma::vec cc = coef_ut_unC(restrain); // extract the ones that suppose to be positive
      // arma::vec cc2 = cc(find(cc < 0)); // identify var comp < 0 (1's)
      if(cc.n_elem > 0){ // we have to restrain
        // rest0 = '(';  rest1=cc.n_elem; rest2 = 'restrained)';
        arma::uvec no_restrain = find((constraints == 1 && coef_ut_unC > 0) || (constraints > 1)); // indices of columns that are OK to use (no restrain)
        arma::mat Inf_norestrain = Inf.submat(no_restrain,no_restrain); // subset of Information matrix
        arma::mat Inf_norestrain_inv; // define the inverse of the information matrix
        arma::inv(Inf_norestrain_inv, Inf_norestrain); // Inverse of Fishers (subset of Inf)
        if(Inf_norestrain_inv.n_rows == 0){ // if fails
          // Rcpp::Rcout << "System is singular (Inf_norestrain_inv). Stopping the job. Try a bigger number of tolParInv." << arma::endl;
          // return 0;
          Rcpp::stop("System is singular (Inf_norestrain_inv). Aborting the job. Try a bigger number of tolParInv.");
        }
        arma::vec scorenorestrain = score(no_restrain); // subset of scores (1st derivatives)
        arma::vec coef_ut_un_norestrain = coef_ut_un(no_restrain); // subset of vc
        arma::vec deltanorestrain; //  define the delta for no restrained
        
        //
        arma::mat InfEM_norestrain = InfEM.submat(no_restrain,no_restrain); // subset of Information matrix
        arma::mat Infw_norestrain = Infw.submat(no_restrain,no_restrain); // subset of Information matrix
        arma::mat InfEMw_norestrain = InfEMw.submat(no_restrain,no_restrain); // subset of Information matrix
        arma::mat InfJoin_norestrain = InfJoin.submat(no_restrain,no_restrain); // subset of Information matrix
        arma::mat InfJoin_inv_norestrain = InfJoin_inv.submat(no_restrain,no_restrain); // subset of Information matrix
        // if(emupdate(cycle) == 1){ // if user wants an EM update (1st derivatives)
        InfJoin_norestrain = (Inf_norestrain*Infw_norestrain)+(InfEM_norestrain*InfEMw_norestrain); // joint information matrix
        InfJoin_inv_norestrain = arma::pinv(InfJoin_norestrain, 1.490116e-08); // inverse the joint information matrix
        deltanorestrain = InfJoin_inv_norestrain * scorenorestrain; //update for variance components where: delta = Information.inv * dL/ds
        // deltanorestrain = (coef_ut_un_norestrain % scorenorestrain % coef_ut_un_norestrain)/n_levels;
        // }else{ // if user wants an information*score update
        //   deltanorestrain = Inf_norestrain_inv * scorenorestrain; //update variance components
        // }
        delta(no_restrain) = deltanorestrain;
        delta(restrain) = delta(restrain)*0;
        
      }//else just keep going
      // end of parameter restrain
      // ^^^^^^^^^^^^^^^^^^
      // ^^^^^^^^^^^^^^^^^^
      coef_ut_un = coef_ut_un + (stepweight(cycle) * delta);
      //
      // constraint the parameters that should be positive and are going negative
      if(cc.n_elem > 0){
        coef_ut_un(restrain) = coef_ut_un(restrain)*0; // the ones that still go below zero and shouldn't let's fix them
      }
      // weight the projection matrix to provide stability
      sigmatwo(arma::find(pos == 0)) =  coef_ut_un(arma::find(pos == 0)); // index of pos
      sigmatwo(arma::find(pos == 1)) = exp(coef_ut_un(arma::find(pos == 1)));
      // the fixed paramters are forced to be the original value
      sigmatwo(find(constraints == 3)) = sigmaF_ut_un(find(constraints == 3));
      // bring back sigma as a list
      sigma = vec_to_cubeCpp(sigmatwo, GeI);
      // check if likelihood has reached it's maximum and stop if so
      llstore(cycle) = llik;
      // get current time
      time_t now = time(0);
      tm *ltm = localtime(&now);
      // keep track of time difference between iterations
      seconds = difftime(now,before);
      // update before time to be the now to be used in the next iteration
      time_t before = time(0);
      localtime(&before);
      // store paramaters
      sigma_store.col(cycle) = sigmatwo;
      if(cycle > 0){
        sigma_perc_change.col(cycle) = ((coef_ut_un/sigma_store.col(cycle-1))-1) * 100; // percent change
      }
      llik_store(cycle) = llik;
      // return output to the console
      if(verbose == true){ //  arma::cout
        if(cycle == 0){Rcpp::Rcout << "iteration   " << " LogLik   " << "  wall    " << "cpu(sec)   " << "restrained" << arma::endl;}
        Rcpp::Rcout << "    " << cycle+1 << "      " <<  llik << "   " << ltm->tm_hour << ":" << ltm->tm_min << ":" << ltm->tm_sec << "      " << seconds << "           " << cc.n_elem << arma::endl;
      }
      // define the end of the algorithm          std::setprecision(5) <<
      if(((cycle > 2) && (delta_llik < tolpar)) || cycle == iters-1 ){ // tolpar*10
        cycle2 = cycle;
        if((cycle > 2) && (delta_llik < tolpar)){convergence = true;}
        last_iteration = true;
        cycle = iters-2;
        // if user wants to get scaled results we don't bring back to original scale
        if(retscaled == false){
          for (int i = 0; i < n_re; ++i) {
            sigma_scaled.slice(i) = sigma.slice(i); // just make a copy before we scale to normal units
            sigma.slice(i) = (sigma.slice(i)%base_var)/sc_var ;
          }
        }
        
        // Fisher inverse
        arma::mat  FI = Inf/2;
        arma::vec myone(pos.n_elem,arma::fill::ones);
        arma::vec sp = ((sigmatwo - myone) % pos) + myone;
        arma::mat FI_c = FI / (sp * sp.t());
        sigma_cov = pinv(FI_c);
        if(sigma_cov.n_rows == 0){ // if fails
          // Rcpp::Rcout << "System is singular (sigma_cov). Aborting the job." << arma::endl;
          // return 0;
          Rcpp::stop("System is singular (sigma_cov). Aborting the job. Try a bigger number of tolParInv.");
        }
      }
      
    }else{// if we are in the last iteration now we calculate u, PEV, B, XB
      
      arma::inv(tXVXi,tXVX);
      if(tXVXi.n_rows == 0){ // if fails try to invert with diag(1e-6)
        arma::inv(tXVXi,tXVX+(D*(tolparinv)));
        if(tXVXi.n_rows == 0){// if fails try to invert with diag(1e-5)
          arma::inv(tXVXi,tXVX+(D*(tolparinv*10)));
          if(tXVXi.n_rows == 0){
            // Rcpp::Rcout << "System is singular (tXVXi). Aborting the job. Try a bigger number of tolParInv." << arma::endl;
            // return 0;
            Rcpp::stop("System is singular (tXVXi). Aborting the job. Try a bigger number of tolParInv.");
          }
        }
      }
      // arma::vec Ym_rw = vectorise(Y.t());
      if(retscaled == true){// if we have to return scaled results we use Yms
        beta = tXVXi * ((Xm.t() * Vi) * Ysm);
      }else{ // we return in normal scale
        beta = tXVXi * ((Xm.t() * Vi) * Ym);
      }
      
      // beta.reshape(X.n_cols,n_traits);
      fitted = Xm * beta;
      residuals = Ym - fitted;
      // residuals.reshape(no,n_traits);
      // arma::vec residuals_rw = vectorise(residuals.t());
      arma::mat Vie = Vi * residuals;
      if(n_random > 0){
        for(i=0; i < n_random; i++){
          arma::mat Zprov = arma::mat(Rcpp::as<arma::sp_mat>(Z[i]));
          
          // arma::mat Zprov2 = Rcpp::as<arma::mat>(Z[i]);
          arma::mat Ki = arma::mat(Rcpp::as<arma::sp_mat>(K[i]));
          // arma::mat Ki2 = (Zprov.t() * Zprov) * 0; //
          // Ki2.diag() = arma::ones<arma::vec>(Ki2.n_cols);
          // arma::mat Ki2 = arma::mat(arma::speye( Zprov.n_cols, Zprov.n_cols ));
          arma::mat VarK;
          arma::mat ZKfv;
          // double VarKscalar = arma::as_scalar(sigma.slice(i));
          
          // IMPORTANT
          // for rrBLUP models we had to allow a K matrix to be a 1 x 1 matrix so dimensions do not match with Z
          if(Ki.n_cols == Zprov.n_cols){ // if a regular random effect
            // Rcpp::Rcout << "regular" << arma::endl;
            VarK = arma::kron(arma::mat(Rcpp::as<arma::sp_mat>(K[i])),sigma.slice(i)); // Gu * var.u
            ZKfv = VarK * arma::kron(Zprov.t(),dD); // G Z'
          }else{ // if huge matrix from models like rrBLUP we need to create a diagonal to calculate VarK and BLUPs
            // Rcpp::Rcout << "rrBLUP" << arma::endl;
            // VarK = arma::kron(Ki2,sigma.slice(i)); // Gu * var.u
            ZKfv = arma::kron(Zprov.t(),dD*sigma.slice(i)); // G Z'
          }
          
          U(i) = ZKfv * Vie; // BLUP = Z' G Vi (Y - Xb)
          if(pev==true){
            
            if(Ki.n_cols == Zprov.n_cols){ // if a regular random effect
              VarU(i) = ZKfv * (P * ZKfv.t()); // var(u) = Z' G [Vi - (VX*tXVXVX)] G Z'
              PevU(i) = VarK - Rcpp::as<arma::mat>(VarU(i)); // PEV = G - var(u)
            }else{ // 
              VarU(i) = ZKfv * (P * ZKfv.t()); // var(u) = Z' G [Vi - (VX*tXVXVX)] G Z'
              // TO BE FIXED
              // not sure how to get the PEV without constructing VarK due to high-memory requirements in rrBLUP models with potentially millions of SNPs
              PevU(i) = Rcpp::as<arma::mat>(VarU(i)); // PEV = G - var(u)
            }
            
          }
        }
      }
    }
    
  }
  // ****************************************************
  // end of algorithm
  // ****************************************************
  arma::vec dd,ee;
  for (int i = 0; i < n_re; ++i) {
    dd = join_cols(dd,mat_to_vecCpp(base_var,Rcpp::as<arma::mat>(GeI[i]))) ; // extract upper triangular in a vector form
    ee = join_cols(ee,mat_to_vecCpp(sc_var,Rcpp::as<arma::mat>(GeI[i]))) ; // extract upper triangular in a vector form
  }
  
  arma::mat FISH = (sigma_cov % (dd*dd.t())) / (ee*ee.t()); // bring back to original scale
  // recalculate V and P with original sigma values
  double AIC = (-2 * llik) + (2 * Xm.n_cols);
  double ny = Ym.n_elem;
  double BIC = (-2 * llik) + (log(ny) * Xm.n_cols);
  // monitor
  sigma_store.each_col() %= dd;
  sigma_store.each_col() /= ee;
  // arma::mat monitor = sigma_store.cols(0, cycle2); // join_cols(llik_store,sigma_store);
  // arma::uvec indices(cycle2,arma::fill::ones);
  // arma::mat monitor2 = monitor.cols(find(indices == 1));
  // arma::mat monitor2 = monitor.cols(0, cycle2);
  arma::mat sigma_perc_change2;
  if(iters > 1){
    sigma_perc_change2 = sigma_perc_change.cols(1, cycle2); // indicate first and last column to subset to return at the end
  }else{
    sigma_perc_change2 = sigma_perc_change; // indicate first and last column to subset to return at the end
  }
  
  // code to put variance-covariance components in original structure defined in thetaC and also blups
  // this piece is what makes the direct inversion to return the same output than the mme algorithm from Jensen
  arma::field<arma::mat> newTheta(thetaConstOri.size());
  arma::field<arma::mat> uList(thetaConstOri.size()-1), uPevList(thetaConstOri.size()-1), PEVs(thetaConstOri.size()-1); // no residual effects stored
  arma::mat u; // join blups in a single matrix format
  arma::field<arma::mat> partitions(thetaConstOri.size()-1); // store indices for each random effect
  arma::vec end, start;
  int counter4 = 0;
  int value = 0;
  for (int i = 0; i < thetaConstOri.size(); ++i) { // for each random effect
    arma::uvec effsToUse = find(thetaIndex == (i+1) ); // which thetas we should use
    arma::mat thetaConstOriIth = thetaConstOri(i); // get effect i
    arma::mat newThetaIth(thetaConstOriIth.n_rows,thetaConstOriIth.n_rows); // to store new thetas
    arma::mat newThetaIthProv; // to store the sigma ith slice provisionally
    arma::mat blupTable, pevTable, partitionsTable, pevFull; // 
    for (int j = 0; j < thetaConstOriIth.n_cols; ++j) { // for each row in theta
      for (int k = 0; k < thetaConstOriIth.n_rows; ++k) { // for each col in theta
        if(thetaConstOriIth(j,k) != 0){ // if was estimated fill the matrix objects
          newThetaIthProv = sigma(counter4);
          if(j==k){ // variance components (both random or residual)
            newThetaIth(j,k) = newThetaIthProv(0,0);
            if( i < (thetaConstOri.size()-1) ){ // if ~random not residual
              // Rcpp::Rcout << (thetaConstOri.size()-1) << arma::endl;
              arma::mat provPev;
              if(pev==true){
                provPev = Rcpp::as<arma::mat>(PevU(counter4));
                arma::mat zerosMatFill(pevFull.n_rows,provPev.n_cols,arma::fill::zeros);
                pevFull = arma::join_cols( arma::join_rows(pevFull, zerosMatFill),arma::join_rows(zerosMatFill.t(),provPev) );
              }
              if( blupTable.n_cols <= j){ // not yet populated
                blupTable = arma::join_rows(blupTable, Rcpp::as<arma::mat>(U[counter4]) );
                if(pev==true){
                  pevTable = arma::join_rows(pevTable, provPev.diag());
                }
              }else{ // the column is already populated
                blupTable.col(j) = blupTable.col(j) + Rcpp::as<arma::mat>(U[counter4]);
                if(pev==true){
                  pevTable.col(j) = pevTable.col(j) + provPev.diag();
                }
              }
              
            }
          }else{ // covariance components
            newThetaIth(j,k) = newThetaIthProv(0,0);
            newThetaIth(k,j) = newThetaIthProv(0,0);
            if( i < (thetaConstOri.size()-1) ){ // if random effect
              arma::vec indexCovJ(blupTable.n_rows, arma::fill::ones);
              indexCovJ=indexCovJ*j;
              arma::vec indexCovK(blupTable.n_rows, arma::fill::ones);
              indexCovK=indexCovK*k;
              arma::vec indexCovJK = arma::join_cols(indexCovJ,indexCovK);
              // find which cov components belong where
              arma::uvec isJ = arma::find(indexCovJK == j);
              arma::uvec isK = arma::find(indexCovJK == k);
              arma::mat provBlup = Rcpp::as<arma::mat>(U[counter4]);
              arma::mat provPev;
              if(pev==true){
                provPev = Rcpp::as<arma::mat>(PevU(counter4));
                provPev = provPev.diag();
              }
              // add cov blup and cov pev 
              if(blupTable.n_cols <= j){ // not yet populated
                blupTable = arma::join_rows(blupTable, provBlup.rows(isJ) );
                if(pev==true){
                  pevTable = arma::join_rows(pevTable, provPev(isJ) );
                }
              }else{ // the column is already populated
                blupTable.col(j) = blupTable.col(j) + provBlup.rows(isJ);
                if(pev==true){
                  pevTable.col(j) = pevTable.col(j) + provPev(isJ);
                }
              }
              if( blupTable.n_cols <= k){ // not yet populated
                blupTable = arma::join_rows(blupTable, provBlup.rows(isK) );
                if(pev==true){
                  pevTable = arma::join_rows(pevTable, provPev(isK) );
                }
              }else{
                blupTable.col(k) = blupTable.col(k) + provBlup.rows(isK);
                if(pev==true){
                  pevTable.col(k) = pevTable.col(k) + provPev(isK);
                }
              }
            }
          }
          counter4++;
        }
      }
    }
    newTheta(i) = newThetaIth;
    if( i < (thetaConstOri.size()-1) ){
      uList(i) = blupTable;
      if(pev==true){
        uPevList(i) = pevTable;
        PEVs(i) = pevFull;
      }
      int nrbt= blupTable.n_rows;
      for (int l = 0; l < blupTable.n_cols; ++l) {
        if( (i==0) && (l==0) ){
          start = Xm.n_cols + (l*nrbt) + 1; // index of where the random effect starts
          end = Xm.n_cols + (nrbt * (l+1)) ; // index of where the random effect ends
        }else{
          start = value + 1;
          end= start + nrbt - 1;
        }
        partitionsTable = arma::join_cols(partitionsTable, arma::join_rows(start,end));
        value = partitionsTable.max();
      }
      u = join_cols(u, arma::vectorise(blupTable) ); // join blups in a single matrix
      partitions(i) = partitionsTable;
      
    }
  }
  arma::mat bu = join_cols(beta, u );
  arma::mat Ci(bu.n_rows,bu.n_rows, arma::fill::zeros);
  Ci.submat(0, 0, beta.n_rows-1, beta.n_rows-1 ) = tXVXi;
  // Rcpp::Rcout << "good2" << arma::endl;
  if(pev==true){
    for (int i = 0; i < partitions.size(); ++i) {//for each major random effect
      Ci.submat(partitions(i).min()-1, partitions(i).min()-1, partitions(i).max()-1, partitions(i).max()-1 ) = PEVs(i);
    }
  }
  
  // ****************************************************
  // return the results
  // ****************************************************
  
  return Rcpp::List::create(
    Rcpp::Named("llik") =  llik_store.cols(0, cycle2) ,
    Rcpp::Named("b") = beta,
    Rcpp::Named("u") = u,
    Rcpp::Named("bu") = bu,
    Rcpp::Named("Ci") = Ci,
    Rcpp::Named("theta") = newTheta,
    Rcpp::Named("theta_se") = FISH, // inverse of fisher's information
    // Rcpp::Named("theta_scaled") = sigma_scaled,
    Rcpp::Named("InfMat") = Inf, // dL2
    Rcpp::Named("monitor") = sigma_store.cols(0, cycle2),
    Rcpp::Named("uList") = uList,
    Rcpp::Named("uPevList") = uPevList,
    Rcpp::Named("AIC") = AIC,
    Rcpp::Named("BIC") = BIC,
    Rcpp::Named("convergence") = convergence,
    Rcpp::Named("partitions") = partitions,
    Rcpp::Named("fitted") = fitted,
    Rcpp::Named("residuals") = residuals,
    Rcpp::Named("dL") = score, // dL
    Rcpp::Named("percChange") = sigma_perc_change2
  // Rcpp::Named("Vi") = Vi,
  // Rcpp::Named("P") = P,
  // Rcpp::Named("u_var") = VarU,
  
  );
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

// [[Rcpp::export]]
Rcpp::List post_mme_Cinverse_cpp(Rcpp::List model, const int mode = 1){

  if(mode < 0 || mode > 2){
    Rcpp::stop("mode must be 0, 1, or 2.");
  }

  if(!model.containsElementNamed("C")){
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

// [[Rcpp::export]]
arma::mat predict_mmes_vcov_cpp(Rcpp::List model, const arma::sp_mat & Dmat){

  // Exact prediction-variance approach that avoids computeCi altogether:
  // Var(D %*% bu) = D C^{-1} D' is obtained by solving C X = D' for the
  // handful of rows actually requested (k = nrow(Dmat)), instead of either
  // the Takahashi selected-inverse subset (computeCi=1, incomplete for
  // arbitrary linear combinations) or a full n x n inverse (computeCi=2,
  // wasteful when k << nEffects). Cost is k sparse triangular solves.
  if(!model.containsElementNamed("C")){
    Rcpp::stop("The fitted model does not contain C. Refit with ai_mme_sp2() returning C.");
  }
  if(!model.containsElementNamed("Cscale")){
    Rcpp::stop("The fitted model does not contain Cscale.");
  }

  arma::sp_mat C = Rcpp::as<arma::sp_mat>(model["C"]);
  const double Cscale = Rcpp::as<double>(model["Cscale"]);

  if(C.n_rows != C.n_cols){
    Rcpp::stop("Stored C must be square.");
  }

  const int nEffects = static_cast<int>(C.n_rows);

  if(static_cast<int>(Dmat.n_cols) != nEffects){
    Rcpp::stop("D must have one column per mixed-model effect (nrow/ncol of stored C).");
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

  const int k = static_cast<int>(Dmat.n_rows);

  Eigen::MatrixXd rhs = Eigen::MatrixXd::Zero(nEffects, k);
  for(arma::sp_mat::const_iterator it = Dmat.begin(); it != Dmat.end(); ++it){
    rhs(static_cast<Eigen::Index>(it.col()), static_cast<Eigen::Index>(it.row())) = (*it);
  }

  Eigen::MatrixXd X = Cfactor.solve(rhs);

  if(Cfactor.info() != Eigen::Success){
    Rcpp::stop("Sparse LDLT solve failed while computing prediction variances.");
  }

  const arma::mat Xarma(X.data(), static_cast<arma::uword>(nEffects), static_cast<arma::uword>(k));

  arma::mat vcov = (Dmat * Xarma) * Cscale;
  vcov = 0.5 * (vcov + vcov.t());

  return vcov;
}


// [[Rcpp::export]]
Rcpp::List ai_mme_sp(const arma::sp_mat & X, const Rcpp::List & ZI,  const arma::vec & Zind,
                     const Rcpp::List & AiI, const arma::sp_mat & y0,
                     const Rcpp::List & SI, const Rcpp::List & partitionsS,
                     const arma::sp_mat & H, const bool & useH,
                     int nIters, double tolParConvLL, double tolParConvNorm,
                     double tolParInv, const Rcpp::List & thetaI,
                     const Rcpp::List & thetaCI, const arma::mat & thetaF,
                     const arma::vec & addScaleParam, const arma::vec & weightEmInf,
                     const arma::vec & weightInf, const bool & verbose
){
  
  time_t before = time(0);
  localtime(&before);
  // define element sizes
  int nSs = SI.size(); // number of residual inverse matrices
  int nZs = ZI.size(); // number of random effects
  int nRe;
  if(nZs > 0){
    nRe = Zind.max(); // number of actual random effects specified in random
  }else{nRe =0;}
  int nZsFake = 1; // a fake value in case there's no random effects we avoid a bad allocation error
  int nReFake = 1; // a fake value in case there's no random effects we avoid a bad allocation error
  int nRRe = thetaI.size(); // number of random + residual effects
  int nX = X.n_cols;// number of fixed effects
  int nR = y0.n_rows; // number of records
  // find variance and mean for the response
  double vary2 = arma::mean(arma::var(arma::square(y0)));
  double vary = arma::mean(arma::var(y0));
  double stdy = arma::mean(arma::stddev(y0));
  double muy = arma::mean(arma::mean(y0));
  arma::sp_mat y = arma::sp_mat(scaleCpp(arma::mat(y0)));
  bool intercept = false;
  // Rcpp::Rcout << intercept << arma::endl;
  if(arma::accu(X.col(0)) == X.n_rows){ // there's an intercept
    intercept = true;
  }
  // create a list to store the symmetric version of thetaC
  arma::field<arma::mat> theta(nRRe), thetaC(nRRe);
  for (int i = 0; i < nRRe; ++i) { // create a copy of thetas
    theta[i]=Rcpp::as<arma::mat>(thetaI[i])/vary ; //
    thetaC[i]=Rcpp::as<arma::mat>(thetaCI[i]) ; //
  }
  // move Z to sparse arma objects
  int nZsAl; // integer to define the allocation of Z
  if(nZs > 0){
    nZsAl = nZs; // if there's random effects the nZs to allocate is equal to Z.size
  }else{
    nZsAl = nZsFake; // otherwise at least we allocate 1 element to avoid the program to crash
  }
  arma::field<arma::sp_mat> Z(nZsAl); // allocate size of Z
  if(nZs > 0){ // if there's random effects
    for (int i = 0; i < nZs; ++i) { // for each Z
      Z(i)=convertSparse(ZI(i)); // convert the matrix to sparse and store in the field
    }
  }
  // delete ZI;
  // move S inverse to sparse arma objects
  arma::field<arma::sp_mat> Si(nSs); // allocate size of Si
  for (int i = 0; i < nSs; ++i) {
    Si(i)=convertSparse(SI(i)); // convert the matrix to sparse and store in the field
  }
  // move Ai to sparse arma objects
  int nReAl;
  if(nZs > 0){
    nReAl = nRe;
  }else{
    nReAl = nReFake;
  }
  arma::field<arma::sp_mat> Ai(nReAl); // allocate size of Ai field
  if(nZs > 0){
    for (int i = 0; i < nRe; ++i) {
      Ai(i)=convertSparse(AiI(i)); // convert the matrix to sparse and store in the field
    }
  }
  // delete AiI;
  // calculate log determinants of Ai's
  arma::rowvec logDetA(nReAl);
  if(nZs > 0){ // of there's random effects
    for (int i = 0; i < nRe; ++i) { // for each random effect
      typedef Eigen::SparseMatrix<double, Eigen::ColMajor, int> RelationshipSpMat;
      typedef Eigen::Triplet<double, int> RelationshipTriplet;
      typedef Eigen::SimplicialLDLT<
        RelationshipSpMat,
        Eigen::Lower,
        Eigen::AMDOrdering<int>
      > RelationshipLDLT;

      if(Ai(i).n_rows > static_cast<arma::uword>(std::numeric_limits<int>::max())){
        Rcpp::stop("Relationship inverse is too large for Eigen's sparse index type.");
      }

      RelationshipSpMat relationshipPrecision(
        static_cast<int>(Ai(i).n_rows),
        static_cast<int>(Ai(i).n_cols)
      );
      std::vector<RelationshipTriplet> relationshipEntries;
      relationshipEntries.reserve(static_cast<std::size_t>(Ai(i).n_nonzero));

      for(arma::sp_mat::const_iterator entry = Ai(i).begin();
          entry != Ai(i).end();
          ++entry){
        relationshipEntries.emplace_back(
          static_cast<int>(entry.row()),
          static_cast<int>(entry.col()),
          *entry
        );
      }

      relationshipPrecision.setFromTriplets(
        relationshipEntries.begin(),
        relationshipEntries.end()
      );
      relationshipPrecision.makeCompressed();

      RelationshipLDLT relationshipFactor;
      relationshipFactor.compute(relationshipPrecision);
      if(relationshipFactor.info() != Eigen::Success){
        Rcpp::stop("Sparse LDLT factorisation of a relationship inverse failed.");
      }

      const Eigen::VectorXd relationshipPivots = relationshipFactor.vectorD();
      double relationshipLogDet = 0.0;
      for(Eigen::Index pivotIndex = 0;
          pivotIndex < relationshipPivots.size();
          ++pivotIndex){
        const double pivot = relationshipPivots(pivotIndex);
        if(!std::isfinite(pivot) || pivot <= 0.0){
          Rcpp::stop("Relationship inverse has a non-positive LDLT pivot.");
        }
        relationshipLogDet += std::log(pivot);
      }

      logDetA(i) = -relationshipLogDet;
    }
  }
  
  // define partitions (only used if random effects exist)
  int last = X.n_cols;
  arma::field<arma::mat> partitions(nReAl); // store indices of the random effects
  arma::vec zsAva;
  int Nu = 0;
  if(nZs > 0){ //if there's random effects (Z matrices) check where each starts and ends
    zsAva = unique(Zind);
    for (int i = 0; i < nRe; ++i) { // for each effect
      arma::uvec indexZind = find(Zind == (i+1) ); // which Z matrices to use , +1 because of the way indeces are used in C++
      int nIndexZind = indexZind.size(); //  number of Z matrices to use
      arma::vec Nus(nIndexZind); // vector to store number of columns in each Z matrix
      // for each matrix in this random effect
      for (int j = 0; j < nIndexZind; ++j) {
        int jj = indexZind(j); // thake the jj matrix
        arma::sp_mat Zprov = Z(jj); // put it in a provisional object
        Nus(j)=Zprov.n_cols; // calculate the number of columns
      }
      arma::vec end = Nus; // define ends and starts
      for (int k = 0; k < nIndexZind; ++k) { // for each effect
        arma::uvec toSum = arma::regspace<arma::uvec>(0,  1,  k); // equivalent to seq()
        end(k)=arma::accu(Nus(toSum));
      }
      arma::vec ones(nIndexZind, arma::fill::ones);
      arma::vec lastM(nIndexZind, arma::fill::value(last));
      arma::vec start = end - Nus + ones;
      start = start + lastM;// adjust start by adding # of fixed effects
      end = end + lastM;//adjust end by adding # of fixed effects
      partitions(i) = arma::join_rows(start,end);
      last = end.max();
      Nu = Nu + accu(Nus);
    }
  }// end of if statement when random effects exist
  
  // define the number of variance components to estimate per random effect structure
  arma::vec nVc(nRRe);
  for (int i = 0; i < nRRe; ++i) {
    arma::mat thetaCprov = thetaC[i];
    arma::uvec nVcProv = find(thetaCprov > 0);
    nVc(i) = nVcProv.size();
  }
  int nVcTotal = accu(nVc); // total number of variance components
  // assign a start and an end index to each covariance structure using the #of VC
  arma::vec nVcEnd = nVc;
  for (int i = 0; i < nRRe; ++i) {
    arma::uvec toSum = arma::regspace<arma::uvec>(0,  1,  i); // equivalent to seq()
    nVcEnd(i)=arma::accu(nVc(toSum));
  }
  arma::vec nVcStart = nVcEnd - nVc + 1;
  // move constraints to vector form binding the columns
  arma::vec thetaCUnlisted;
  for (int i = 0; i < nRRe; ++i) {
    thetaCUnlisted = join_cols(thetaCUnlisted,mat_to_vecCpp2(thetaC[i],thetaC[i]));
  }
  // removing complex structures how many effects are really there
  arma::vec nUsTotal(nReAl);
  if(nZs > 0){ //
    for (int i = 0; i < nRe; ++i) {
      arma::mat partitionsProv = partitions(i);
      nUsTotal(i) = partitionsProv(0,1) - partitionsProv(0,0) + 1;
    }
  }
  // define objects to store theta and llik across iterations
  arma::mat monitor(nVcTotal,nIters); // matrix to store variance components
  arma::mat percChange(nVcTotal,nIters); // matrix to store variance components
  arma::rowvec llik(nIters); // store log likellihood values
  
  int nEffects = Nu+nX;
  int nEffectsPlusY = nEffects + 1;
  arma::mat Mchol; // (nEffectsPlusY,nEffectsPlusY)
  arma::umat PM_mat;
  arma::mat MWuchol;
  arma::umat PMWu_mat;
  
  arma::sp_mat M(nEffectsPlusY,nEffectsPlusY), W(nR,nEffects), Wy(nR,nEffectsPlusY), C(nEffects,nEffects), Ci(nEffects,nEffects);
  arma::vec u(Nu), b(nX), bu(nEffects);
  arma::mat buWu(nEffects,nVcTotal);
  arma::mat avInf(nVcTotal,nVcTotal);
  arma::mat emInf(nVcTotal,nVcTotal);
  arma::mat InfMat(nVcTotal,nVcTotal);
  arma::mat InfMatInv(nVcTotal,nVcTotal);
  bool convergence = false;
  double seconds;
  arma::sp_mat XWjxZWj(nEffects,nVcTotal), WiWj(nVcTotal,nVcTotal);
  arma::mat Mchol_XZ;//, Wu2;
  arma::sp_mat I = arma::speye<arma::sp_mat>(nEffectsPlusY,nEffectsPlusY);
  arma::mat I2 = arma::eye(nEffects+nVcTotal,nEffects+nVcTotal);
  arma::vec delta(nVcTotal), delta_minus1(nVcTotal);
  // objects for constraints
  arma::mat percDelta(nVcTotal,nIters,arma::fill::zeros); // store % change of the delta with respect to the previous iteration
  arma::mat normMonitor(3,nIters); // store in each iteration the 3 stopping criteria of Madsen and Jensen
  arma::mat toBoundary(nIters,nVcTotal, arma::fill::zeros ); // store which values have been set to the boundary value
  arma::vec sumToBoundary(nVcTotal, arma::fill::zeros ); // to apply sum across iterations and if a VC goes to the boundary 3 times it is fixed to the boundary
  arma::sp_mat Ri(nR,nR); // matrix to store R inverse
  arma::sp_mat Hs(H.n_cols,H.n_cols); // square of H matrix
  if(useH == true){ // do cholesky decomposition of H if user wants to use weights
    Rcpp::Rcout << "Using the weights matrix " << arma::endl;
    Hs = arma::sp_mat(chol(arma::mat(H)));
  }
  arma::vec dLuOut;//(nVcTotal); // we will join cols
  ////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////
  // START ITERATIVE ALGORITHM
  ////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////
  
  
  for (int iIter = 0; iIter < nIters; ++iIter) {
    
    // ###########################
    // # 1) absorption of M into y to obtain y'Py and logDetC
    // # PAPER FORMULA from Jensen and Madsen 1997, Gilmour et al., 1995
    // # expand coefficient matrix (C) to have the response variable
    // # M = W' Ri W # with W = [X Z y]
    // #
    // #     [X'RiX  X'RiZ     X'Riy ]
    // # M = [Z'RiX  Z'RiZ+Gi  Z'Riy ]
    // #     [y'RiX  y'RiZ     y'Riy ]
    // #
    // # where Gi = Ai*(s2e/s2u) = (A*s2u)*s2e = kronecker(Ai,solve(s2u))
    // #
    // # lambda = solve(theta) # inverse of var-covar matrices
    // # MChol = chol(M)
    // # yPy = MChol[n,n] # where n is the last element of the matrix
    // # logDetC = 2 * E log(diag(MChol))
    // ###########################
    arma::vec thetaResidualsVec = mat_to_vecCpp2(theta(nRRe-1),thetaC[(nRRe-1)]);
    if(thetaResidualsVec.n_elem != Si.size()){
      Rcpp::stop("The number of residual parameters does not match the residual covariance bases.");
    }
    arma::sp_mat Rmat(nR,nR);
    for (int i = 0; i < nSs; ++i) {
      Rmat += Si(i) * thetaResidualsVec(i);
    }
    arma::mat Ridense;
    bool okR = arma::inv_sympd(Ridense, arma::mat(Rmat));
    if(okR == false){
      arma::mat Rpd = nearPDcpp(arma::symmatu(arma::mat(Rmat)), 100, 1e-06, 1e-07);
      Rmat = arma::sp_mat(Rpd);
      okR = arma::inv_sympd(Ridense, Rpd);
      if(okR == false){
        Rcpp::stop("Inversion of the residual covariance matrix failed.");
      }
    }
    Ri = arma::sp_mat(Ridense);
    // adjust R inverse if user provides weights
    if(useH == true){
      Ri = Hs *  Ri * Hs.t();
    }
    // Rcpp::Rcout << ltm->tm_hour << ":" << ltm->tm_min << ":" << ltm->tm_sec << "      " << seconds << "           " << arma::endl;
    /////////////////////////////////
    // form the mixed model equations
    if(iIter == 0){ // only form W and Wy once in the first iteration
      W = X;
      // if random effects exist
      if(nZs > 0){
        for (int i = 0; i < nZs; ++i) {
          W = arma::join_rows( W, Z(i) );
        }
      }
      Wy = arma::join_rows(W,y);
    }
    
    M = Wy.t() * Ri * Wy;
    
    arma::field<arma::sp_mat> lambda(nReAl); // to store theta inverses
    arma::field<arma::sp_mat> GI(nReAl); // to store kron(thetainv,Ainv) 
    if(nZs > 0){
      for (int i = 0; i < nRe; ++i) {
        arma::mat bend = nearPDcpp(arma::symmatu(theta(i)), 100, 1e-06, 1e-07);
        // lambda(i) = arma::sp_mat( dddd );
        // arma::mat bend = arma::eye(theta(i).n_rows,theta(i).n_rows) * 1e-6;
        lambda(i) = arma::sp_mat( arma::pinv( bend ) );
        GI(i) = kron(lambda(i), Ai(i) );
        arma::mat partitionsP = partitions(i);
        int ff = partitionsP(0,0) - 1;
        int ll = partitionsP(partitionsP.n_rows-1,1) - 1;
        M.submat( ff, ff, ll, ll ) = M.submat( ff, ff, ll, ll ) + GI(i);
      }
    }
    
    // Rcpp::Rcout << "all good" << arma::endl;
    // arma::chol(Mchol, PM_mat, arma::mat(M), "upper", "matrix");
    // Rcpp::Rcout <<  Mchol.n_cols << arma::endl;
    // Rcpp::Rcout <<  Mchol.n_rows << arma::endl;
    
    bool okChol = arma::chol(Mchol, arma::mat(M));
    if(okChol == false){
      if(verbose == true){
        Rcpp::Rcout << "Making M positive definite " << arma::endl;
      }
      arma::mat Mp = arma::symmatu(arma::mat(M));
      Mp = nearPDcpp(Mp, 100, 1e-06, 1e-07);
      Mchol = arma::chol(Mp) ;
      if(verbose == true){
        Rcpp::Rcout << "Cholesky of M succeeded " << arma::endl;
      }
    }
    arma::vec yPy =arma::square(Mchol.submat( Mchol.n_rows-1, Mchol.n_cols-1, Mchol.n_rows-1,  Mchol.n_cols-1 ));
    Mchol_XZ = Mchol.submat( 0,0, Mchol.n_rows-2,  Mchol.n_cols-2 ); // M without y portion (last row and column of M)
    arma::vec My = Mchol.submat( 0,Mchol.n_rows-1, Mchol.n_cols-2,  Mchol.n_cols-1 );
    double logDetC = 2 * accu(log(Mchol_XZ.diag()));
    // ###########################
    // # 1.1) calculate the log-likelihood
    // # PAPER FORMULA (Lee and Van der Werf, 2006)    #
    // # LL = -0.5 [((Nr-Nb-Nu-...)*ln(s2e)) - ln|C| + ln|Au| + ... + (Nu*ln(s2u)) + ... + y'Py ]
    // # PAPER FORMULA (Jensen and Madsen, 1997)
    // # LL = -0.5 [ln|C| + ln|R| + (ln|A.u| +ln|theta.u|) + ... + y'Py ]
    // # where | | is the determinant of a matrix
    // #       A.u is the pure relationship matrix for the uth random effect
    // #       theta.u is the vc matrix for the uth random effect
    // ###########################
    
    double llikp=0;
    if(nZs > 0){
      for (int i = 0; i < nRe; ++i) {
        double val;
        double sign;
        bool ok1 = log_det(val, sign, theta(i));  // form 2
        if(ok1 == false){ Rcpp::Rcout << "log determinant failed " << arma::endl;};
        llikp = llikp + (nUsTotal(i)*val*sign) + (logDetA(i)*theta(i).n_rows);
      }
    }
    double val;
    double sign;
    bool ok2 = log_det(val, sign, arma::mat(Rmat));
    if(ok2 == false){ Rcpp::Rcout << "log determinant of R failed " << arma::endl;};
    double logDetR = val * sign;
    llik(iIter) = (- 0.5) * ( llikp + logDetC + logDetR + arma::as_scalar(yPy) );
    // ###########################
    // # 2) backsubstitute to get b and u (CORRECT)
    // # use the results from the absorption to obtain BLUE & BLUPs
    // # b = backsolve(MChol[,rest],MChol[,last])
    // ###########################
    
    // arma::spsolve_factoriser SF;
    // bool status = SF.factorise(arma::sp_mat(Mchol_XZ));
    // if(status == false) { Rcpp::Rcout << "factorisation failed" << arma::endl; }
    // double rcond_value = SF.rcond();
    // bool solution1_ok = SF.solve(bu,My);
    
    arma::spsolve(bu, arma::sp_mat(Mchol_XZ) , My, "lapack" );  // use LAPACK  solver
    
    arma::uvec bInd = arma::regspace<arma::uvec>(0,  1,  (nX-1)); // which ones are BLUEs, equivalent to seq()
    b = bu(bInd); // move BLUEs to a different vector
    if(nZs > 0){ // move BLUPs to a different vector
      arma::uvec uInd = arma::regspace<arma::uvec>((nX),  1,  (nX+Nu-1)); // equivalent to seq()
      u = bu(uInd);
    }
    // ###########################
    // # 3) calculate Wu (working variates)
    // # PAPER FORMULA (Notes on Estimation of Genetic Parameters from Van der Werf)
    // # wu = Zu/s2u; we = e/s2e
    // # PAPER FORMULA (Jensen and Madsen, 1997)
    // # U = [u1 | u2 | ... | ui]
    // # US = U * lambda
    // # Wu.ii = Zui*USi # for variance component
    // # Wu.ij = Zui*USj + Zuj*USi # for covariance component
    // # Wr.j = Rj * Rinv * e  # for residual variance component
    // ###########################
    
    arma::sp_mat Wu;
    arma::field<arma::sp_mat> uSinv(nReAl); // field of sparse matrices to store products u * thetainv
    
    if(nZs > 0){ // if random effects exist
      for(int iR = 0; iR < nRe; ++iR){ // for each random effect u
        arma::mat partitionsP = partitions(iR); // access the partition for the iR random effect
        arma::sp_mat U(partitionsP(0,1) - partitionsP(0,0) + 1, partitionsP.n_rows); // all BLUPs all the ith random effect
        for(int iRow = 0; iRow < partitionsP.n_rows; ++iRow){ // for each partition row
          arma::uvec usedPartition = arma::regspace<arma::uvec>((partitionsP(iRow,0)-1),  1, (partitionsP(iRow,1)-1)  ); // equivalent to seq()
          U.col(iRow) = bu(usedPartition); // move vector to a matrix: u -> [u1|u2|...]
        }
        // # [a || m] [s2a || sam] = [s2a a + sam m  || sam a + s2m m]
        // #          [sam || s2m]
        arma::sp_mat uSinvProv = U * lambda(iR); // 
        uSinv(iR) = uSinvProv;
        arma::mat thetaCprov = thetaC[iR];
        // // for the ij var comp we calculate the Wu
        arma::uvec useZind = find(Zind == iR+1); // which Z matrices we should use for this random effect
        arma::vec nVcForIr = mat_to_vecCpp2(thetaC[iR],thetaC[iR]);
        arma::sp_mat WuiR(nR, nVcForIr.size()); // store Wu for random effect iR (we need to all, vc and cov)
        arma::sp_mat ZuiR(nR, useZind.size()); // store Wu for random effect iR (we need only for vc)
        int counterWu = 0;
        int counterZu = 0;
        //
        for(int iRow = 0; iRow < lambda(iR).n_rows; ++iRow){
          for(int iCol = 0; iCol < lambda(iR).n_cols; ++iCol){
            if(thetaCprov(iRow,iCol) > 0){ // if vc has to be estimated
              if(iRow == iCol){ // variance component
                // Wu
                WuiR.col(counterWu) = Z(useZind(iRow)) * uSinvProv.col(iRow);
                counterWu ++;
                // Zu
                ZuiR.col(counterZu) = Z(useZind(iRow)) * U.col(iRow) ;
                counterZu ++;
              }else{ // covariance component
                // Wu
                WuiR.col(counterWu) = ( Z(useZind(iCol)) * uSinvProv.col(iRow) )  -  ( Z(useZind(iRow)) * uSinvProv.col(iCol) )  ;
                counterWu ++;
              }
            } // end of thetaCprov(iRow,iCol) > 0
          } // end of icol loop
        } // end of irow loop
        Wu = arma::join_rows(Wu,WuiR);
      } // end of loop for each random effect
    } // end of condition when random effects exist
    // calculate residuals
    arma::vec e = y - (arma::sp_mat(W.submat(0,0,W.n_rows-1,W.n_cols-1)) * bu);
    // Working variates for residual VCs
    for(int iS = 0; iS < Si.size(); ++iS){
      Wu = arma::join_rows(Wu , Si(iS) * Ri * arma::sp_mat(e) );
    }
    // ###########################
    // # 4) absorption of M into Wu (2 VAR, 1 COV) to obtain Wu' P Wu  which is the AI matrix
    // # we had to change the avInf to avInf/sigmas
    // # PAPER FORMULA (Smith, 1995) Differentiation of the Cholesky Algorithm
    // # avInf.ij = ((chol(M))[n,n])^2 # the square of the last diagonal element of the cholesky factorization
    // # where:
    // #        [X'RiX  X'RiZ     X'Riwj ]
    // # M.Wu = [Z'RiX  Z'RiZ+Gi  Z'Riwj ]
    // #        [wk'RiX wk'RiZ    wk'Riwj]
    // # where:
    // # wi: working variate i
    // # wk: working variate k
    // # Ri: is R inverse
    // # and the the part corresponding to X and Z is the coefficient matrix C
    // # AI = (M.Wu.chol)^2
    // ###########################
    
    XWjxZWj = W.t() * Ri * Wu ;// [X'Riwj Z'Riwj]' # C12 upper right
    WiWj = Wu.t() * Ri * Wu ;//  wk'Riwj # C22 lower right
    // // solve method !!
    // // similar to arma::spsolve(bu, arma::sp_mat(Mchol_XZ) , My, "lapack" ); but My = XZRiy and XWjxZWj = XZRi.Wu
    arma::spsolve(buWu, arma::sp_mat(M.submat( 0,0, M.n_rows-2,  M.n_cols-2 )), arma::mat(XWjxZWj), "lapack" );  // use LAPACK  solver
    avInf = WiWj - (buWu.t()*XWjxZWj); // WuWu' - bu.Wu'*W.Wu
    
    // cholesky method!! requires scaling of the response to work
    // arma::mat MWu = arma::join_cols(
    //   arma::join_rows(arma::mat(M.submat( 0,0, M.n_rows-2,  M.n_cols-2 )),arma::mat(XWjxZWj0) ),
    //   arma::join_rows(arma::mat(XWjxZWj0.t()), arma::mat(WiWj) )
    // );
    // arma::chol(MWuchol, PMWu_mat, MWu, "lower", "matrix");
    // MWuchol = PMWu_mat.t() * MWuchol;
    // avInf = MWuchol.submat( MWuchol.n_cols-Wu.n_cols, MWuchol.n_cols-Wu.n_cols, MWuchol.n_cols-1, MWuchol.n_cols-1);
    // avInf = avInf * avInf.t();
    
    // ##########################
    // # 5) get 1st derivatives (dL/ds2i) from MME-version
    // # PAPER FORMULA (Lee and Van der Werf, 2006)
    // # dL/ds2u = -0.5 [(Nu/s2u) - (tr(AiCuu)/s4u) -  (e/s2e)'(Zu/s2u)]
    // # dL/ds2e = -0.5 [((Nr-Nb)/s2e) - [(Nu - (tr(AiCuu)/s2u))*(1/s2e)] - ... - (e/s2e)'(e/s2e)]
    // #
    // # PAPER FORMULA (Jensen and Madsen, 1997)
    // # dL/ds2u = (q.i * lambda) - (lambda * (T + S) * lambda)  Eq. 18
    // # dL/ds2e = tr(Rij*Ri) - tr(Ci*W'*Ri*Rij*Ri*W) - (e'*Ri*Rij*Ri*e)
    // ###########################
    
    // get the inverse of the coefficient matrix
    arma::vec v(Mchol.n_cols-1, arma::fill::ones);//option 1
    arma::mat D = diagmat(v); // option 1
    arma::mat Cichol = arma::solve( trimatu(Mchol.submat(0,0,Mchol.n_rows-2,Mchol.n_cols-2) ), D);  // indicate that A is triangular; option 1
    // multiply by it's transpose
    arma::sp_mat Cip = arma::sp_mat(Cichol);
    Ci = Cip * Cip.t() ;
    arma::field<arma::mat> emInfList(nRRe);
    arma::vec dLu;//(nVcTotal); // we will join cols // dLu(Wu.n_cols);//(nVcTotal); // we will join cols
    if(nZs > 0){ // if random effects exist
      for(int iR = 0; iR < nRe; ++iR){ // for each random effect u
        arma::mat thetaCprov = thetaC[iR];
        arma::sp_mat traces(lambda(iR).n_rows,lambda(iR).n_cols);
        for(int iRow = 0; iRow < lambda(iR).n_rows; ++iRow){
          for(int iCol = 0; iCol < lambda(iR).n_cols; ++iCol){
            if(thetaCprov(iRow,iCol) > 0){ // if vc has to be estimated
              arma::mat partitionsP = partitions(iR);
              // X.submat( first_row, first_col, last_row, last_col )
              double trAiCuu = arma::trace(  Ai(iR) * Ci.submat(partitionsP(iRow,0)-1, partitionsP(iCol,0)-1, partitionsP(iRow,1)-1, partitionsP(iCol,1)-1 )  );
              traces(iRow,iCol) = trAiCuu;
            }else{
              traces(iRow,iCol) = 0;
            }
          }// end of loop for iCol
        }// end of loop for iRow
        traces = arma::symmatu(traces); // copy upper in lower triangular
        //  first derivatives = dL/ds2u = (q.i * lambda) - (lambda * (T + S) * lambda)    where S=UAiU and we use U.lambda
        arma::sp_mat dLuProv = (arma::as_scalar(nUsTotal(iR)) * lambda(iR) ) - ( uSinv(iR).t() * Ai(iR) * uSinv(iR) ) - ( lambda(iR) * traces * lambda(iR) );
        // althernative EM update
        // current(theta)   -   update(delta)  but we need to decompose the update(delta) = Iem * vech(dLu/ds2u) , Iem is then of dimensions equal to vech(dLu/ds2u)
        // theta[[iR]] - (theta[[iR]]%*%dLuProv%*%theta[[iR]])/Nus[iR]    Eq.34
        arma::vec thetaUnlisted = mat_to_vecCpp2( theta(iR),thetaCprov);
        arma::mat thetaUnlistedMat = diagmat(thetaUnlisted);
        arma::mat emInfInvProvExt = ( thetaUnlistedMat * thetaUnlistedMat.t() ) / arma::as_scalar(nUsTotal(iR));
        emInfList(iR) = arma::pinv(emInfInvProvExt, tolParInv);
        dLu = join_cols( dLu, mat_to_vecCpp2(arma::mat(dLuProv),thetaCprov) );
      }// end of loop for each random effect
    }// end of condition if random effects exist
    arma::vec dLe(Si.size());
    arma::sp_mat eProv = arma::sp_mat(e);
    for(int iS = 0; iS < Si.size(); ++iS){ // Rij <- S[[iS]]%*%Ri
      arma::sp_mat Sprov = Si(iS);
      dLe(iS) = ( arma::trace( Sprov *Ri) - arma::trace( Ci * W.t() * Ri * Sprov * Ri * W ) ) - arma::as_scalar( eProv.t() * Ri * Sprov * Ri * eProv );
    }
    arma::vec thetaRUnlisted = mat_to_vecCpp2(theta(nRRe-1),thetaC[nRRe-1]);
    arma::mat thetaRUnlistedMat = diagmat(thetaRUnlisted);
    arma::mat emInfInvRProvExt = ( thetaRUnlistedMat * thetaRUnlistedMat.t() ) / arma::as_scalar(nR);
    emInfList(nRRe-1) = arma::pinv(emInfInvRProvExt, tolParInv); // we invert because this is equivalen to the inverse of the information and we need the information
    
    dLu = join_cols( dLu, dLe );// join the random and residual first derivatives in a single vector
    
    for(int i = 0; i < nRRe; ++i){
      emInf.submat(nVcStart(i)-1, nVcStart(i)-1, nVcEnd(i)-1, nVcEnd(i)-1 ) = emInfList(i);
    }
    
    // ###########################
    // # 6) update the variance paramters using the Newton method
    // # PAPER FORMULA (Lee and Van der Werf, 2006)
    // # theta.n+1 = theta.n + (AInfi * dL/ds2)
    // ###########################
    
    arma::vec thetaUnlisted, thetaCUnlisted;
    for(int i = 0; i < nRRe; ++i){
      thetaUnlisted = arma::join_cols(thetaUnlisted, mat_to_vecCpp2( theta(i),thetaC(i)) );
      thetaCUnlisted = arma::join_cols(thetaCUnlisted, mat_to_vecCpp2(thetaC(i),thetaC(i)) );
    }
    // Rcpp::Rcout << "thetaCUnlisted" << thetaCUnlisted << arma::endl;
    // create the 'weight' EM information matrix (TO BE USED LATER WITHIN THE OPTIMIZATION)
    arma::vec v2(nVcTotal, arma::fill::ones) ;
    arma::mat weightEmInfMat = diagmat(v2) * arma::as_scalar(weightEmInf(iIter));
    arma::mat weightAiInfMat = diagmat(v2) * (1 - arma::as_scalar(weightEmInf(iIter)));
    // Joint information matrix and update
    //                  AVERAGE INFORMATION                         +       EXPECTATION MAXIMIZATION
    InfMat = (weightAiInfMat * avInf) + (weightEmInfMat * emInf);
    InfMatInv = arma::pinv(InfMat, tolParInv); // inverse of the information matrix
    delta = (InfMatInv * arma::as_scalar(weightInf(iIter))) * dLu; // delta = I- * dLu/dLx
    // new values for variance components theta.i+1 = theta.i + delta
    arma::vec  expectedNewTheta = thetaUnlisted - delta;
    
    // #######################
    // # 7) APPLY CONSTRAINTS to VC
    // # suggestions from Madsen and Jensen (1997) and Gilmour (2019)
    // #######################
    // A) apply constraints for fixed and positive
    for(int i = 0; i < thetaCUnlisted.size() ; ++i){
      if(thetaCUnlisted(i) == 1){
        if(expectedNewTheta(i) < 1e-10){
          // Rcpp::Rcout << "Restraining to small value" << arma::endl;
          expectedNewTheta(i)=1e-10;
          toBoundary(iIter,i)=1; // toBoundary(nIters,nVcTotal)
          // change to fixed the ones that constantly (3 times) go to the boundaries
          for(int j = 0; j < toBoundary.n_cols ; ++j){
            sumToBoundary(j) =arma::accu(toBoundary.col(j)) ;
          }
          // to force
          arma::uvec toBeForced = arma::find(sumToBoundary >= 2);
          if(toBeForced.n_elem > 0){
            thetaCUnlisted(toBeForced) = thetaCUnlisted(toBeForced) - thetaCUnlisted(toBeForced) + 3;
          }
        }// end of if(expectedNewTheta(i) < 1e-10)
      }// end of positive constraints
      // any vc outide the search space should come back
      if(expectedNewTheta(i) > 1){ // since we scale the modell we don't allow to explain more than 90%
        // Rcpp::Rcout << "Restraining to small value" << arma::endl;
        expectedNewTheta(i)=1;
        toBoundary(iIter,i)=1; // toBoundary(nIters,nVcTotal)
      }// end of if(expectedNewTheta(i) > 1)
      if(thetaCUnlisted(i) == 3){
        arma::vec thetaUnlistedPlusAddScaleParam;
        if(iIter == 0){
          thetaUnlistedPlusAddScaleParam = arma::join_cols(thetaUnlisted,addScaleParam); // expectedNewTheta
        }else{
          thetaUnlistedPlusAddScaleParam = arma::join_cols(monitor.col((iIter-1)),addScaleParam);
        }
        // theta.i            = scaleParameter.selected      *  Theta
        expectedNewTheta(i) = arma::as_scalar(thetaF.row(i) *  thetaUnlistedPlusAddScaleParam);
      }
    }
    // B)  if there's constrained or to boundary VC we need to partition InfMat and make a different update
    arma::vec thetaCUnlisted0 = thetaCUnlisted;
    thetaCUnlisted0(arma::find(thetaCUnlisted0 < 3)) = thetaCUnlisted0(arma::find(thetaCUnlisted0 < 3)) * 0; // we make everything zero except fully constrained
    arma::uvec constrained = arma::find((thetaCUnlisted0+sumToBoundary) > 0);
    arma::uvec unconstrained = arma::find((thetaCUnlisted0+sumToBoundary) <= 0);//arma::find(thetaCUnlisted != 3);
    arma::mat InfMat_uu, InfMat_ff, InfMat_uf, InfMatInv_uu, InfMatInv_ff;
    arma::vec dLu_uu, dLu_ff, delta_uu, delta_ff;
    if(constrained.n_elem > 0){   // Rcpp::Rcout << "Updates using constrained Information matrix" << arma::endl;
      InfMat_uu = InfMat(unconstrained,unconstrained);
      InfMat_ff = InfMat(constrained,constrained);
      InfMat_uf = InfMat(unconstrained,constrained);
      dLu_uu = dLu(unconstrained);
      dLu_ff = dLu(constrained);
      InfMatInv_ff = pinv(InfMat_ff,tolParInv);
      InfMatInv_uu = pinv(InfMat_uu,tolParInv);
      delta_ff = InfMatInv_ff * dLu_ff;
      delta_uu = InfMatInv_uu * (dLu_uu - ( InfMat_uf * (InfMatInv_ff*delta_ff) ) );
      delta(unconstrained) = delta_uu;
      expectedNewTheta(unconstrained) = thetaUnlisted(unconstrained) - delta(unconstrained);
    }
    // C) quantify delta changes
    if(iIter == 0){
      delta_minus1 = delta;
    }else{
      percDelta.col(iIter) =delta/delta_minus1; 
    }
    // D) if not positive-definite change to EM update
    arma::vec pdCheck;
    for(int i = 0; i < nRRe; ++i){ // caculate eigen values
      arma::uvec toFill = arma::regspace<arma::uvec>(nVcStart(i)-1,  1,  nVcEnd(i)-1); // equivalent to seq()
      arma::mat thetaProvNew = vec_to_matCpp(expectedNewTheta(toFill), thetaC[i] );
      arma::mat thetaProvNewS = arma::symmatu(thetaProvNew);
      arma::vec eigval;
      arma::mat eigvec;
      arma::eig_sym(eigval, eigvec, thetaProvNewS);  // find 5 eigenvalues/eigenvectors
      pdCheck = arma::join_cols(pdCheck,eigval);
    }
    arma::uvec eigenFind = find(pdCheck < 0);
    if(eigenFind.n_elem > 0){
      if(verbose == true){ //
        Rcpp::Rcout << "Updated VC is not positive definite, changing to EM step" << arma::endl;
      }
      InfMat = (0.5 * avInf) + (0.5 * emInf);
      if(constrained.n_elem > 0){
        Rcpp::Rcout << "Update using constraints" << arma::endl;
        InfMat_uu = InfMat(unconstrained,unconstrained);
        InfMat_ff = InfMat(constrained,constrained);
        dLu_uu = dLu(unconstrained);
        dLu_ff = dLu(constrained);
        InfMatInv_uu = pinv(InfMat_uu,tolParInv);
        delta_uu = InfMatInv_uu * dLu_uu;
        delta(unconstrained) = delta_uu;
        // delta(constrained) = delta(constrained) * 0;
        expectedNewTheta(unconstrained) = thetaUnlisted(unconstrained) - delta(unconstrained);
      }else{
        InfMatInv = arma::pinv(InfMat, tolParInv); // inverse of the information matrix
        expectedNewTheta = thetaUnlisted - (InfMatInv * dLu);
      }
    }
    // #######################
    // # 8) Bring back theta in vector form to matrix form and save for monitor
    // #######################
    monitor.col(iIter) = expectedNewTheta;
    for(int i = 0; i < nRRe; ++i){
      arma::uvec toFill = arma::regspace<arma::uvec>(nVcStart(i)-1,  1,  nVcEnd(i)-1); // equivalent to seq()
      arma::mat thetaProvNew = vec_to_matCpp(expectedNewTheta(toFill), thetaC[i] );
      arma::mat thetaProvNewPD = nearPDcpp(arma::symmatu(thetaProvNew), 100, 1e-06, 1e-07); // maxit=100, eig_tol = 1e-06, conv_tol = 1e-07
      // arma::mat thetaProvNewPD = thetaProvNew + (arma::eye(thetaProvNew.n_rows,thetaProvNew.n_rows)*tolParInv); // nearPDcpp(arma::symmatu(thetaProvNew), 100, 1e-06, 1e-07);
      theta(i) = thetaProvNewPD; // arma::symmatu(thetaProvNew);
      //
      arma::mat thetaCProvNew = vec_to_matCpp(thetaCUnlisted(toFill), thetaC[i] );
      thetaC(i) = thetaCProvNew;
    }
    // #######################
    // # 9) Stopping criteria
    // #######################
    // get current time
    time_t now = time(0);
    tm *ltm = localtime(&now);
    // keep track of time difference between iterations
    seconds = difftime(now,before);
    // update before time to be the now to be used in the next iteration
    time_t before = time(0);
    localtime(&before);
    
    // E) stopping criteria (norms)
    normMonitor(0,iIter) = arma::norm(delta(unconstrained), 1) ; // stopping criteria 1
    normMonitor(1,iIter) = arma::norm(dLu(unconstrained),1) ; // stopping criteria 1
    arma::vec nVcTotalv(1);
    nVcTotalv(0) = nVcTotal;
    arma::vec stopCriteria3 = (InfMatInv.diag()/arma::as_scalar(arma::sqrt(nVcTotalv)) ) % dLu;
    normMonitor(2,iIter) = arma::norm(stopCriteria3(unconstrained),1) ; // stopping criteria 1
    
    arma::uvec restrained = arma::find(toBoundary.row(iIter) > 0);
    if(verbose == true){ //
      if(iIter == 0){Rcpp::Rcout << "iteration   " << " LogLik   " << "  wall    " << "cpu(sec)   " << "restrained" << arma::endl;}
      Rcpp::Rcout << "    " << iIter+1 << "      " <<  llik(iIter) << "   " << ltm->tm_hour << ":" << ltm->tm_min << ":" << ltm->tm_sec << "      " << seconds << "           " <<  restrained.n_elem << arma::endl;
    }
    if(iIter > 0){
      double delta_llik = llik(iIter) - llik(iIter-1);
      // if( (  (delta_llik < tolParConv)) || (iIter == nIters)  ){ // || changePer < .001
      if( (normMonitor(2,iIter) <  tolParConvNorm) || (delta_llik < tolParConvLL)  || (iIter == nIters) ) {
        // if(delta_llik < tolParConv){
        if( (normMonitor(2,iIter) < tolParConvNorm) || (delta_llik < tolParConvLL) ){
          convergence = true;
        }
        monitor = monitor.cols(0,iIter);
        normMonitor = normMonitor.cols(0,iIter);
        percDelta = percDelta.cols(0,iIter);
        llik = llik.cols(0,iIter);
        break; //
      }
    }
    dLuOut = dLu;
  }// end of iterative optimization
  
  
  double AIC = (-2 * llik((llik.n_cols-1))) + (2 * nX);
  double BIC = (-2 * llik((llik.n_cols-1))) + (log(nR) * nX);
  // move constraints to vector form binding the columns
  arma::vec thetaCUnlistedFinal;
  for (int i = 0; i < nRRe; ++i) {
    thetaCUnlistedFinal = join_cols(thetaCUnlistedFinal,mat_to_vecCpp2(thetaC[i],thetaC[i]));
  }
  // bring back to original scale the variance components
  for (int i = 0; i < nRRe; ++i) {
    theta(i) = theta(i)*vary;
  }
  bu = bu*stdy;
  b = b*stdy;
  // Rcpp::Rcout << intercept << arma::endl;
  if(intercept==true){ // if true mu is required in the first position only
    b(0,0) = b(0,0) + muy;
    bu(0,0) = bu(0,0) + muy;
  }else{ // if false mu is required all over b
    bu.submat(0, 0, (b.n_rows-1), 0) = bu.submat(0, 0, (b.n_rows-1), 0) + muy;
    b = b + muy;
  }
  Ci = Ci*vary;
  // InfMat=InfMat/stdy;//*(1/(vary2/2));
  InfMatInv=InfMatInv*vary;//*(vary2/2);
  Mchol_XZ=Mchol_XZ/stdy;
  monitor=monitor*vary;
  // dLuOut=dLuOut/vary;
  // move the effects from vector to a field with matrices
  arma::field<arma::mat> uList(nRe), uPevList(nRe); // store indices of the random effects
  for (int i = 0; i < nRe; ++i) {
    arma::mat partitionsP = partitions(i); // partition of random effect i
    arma::mat uMat( arma::as_scalar( partitionsP(0,1) - partitionsP(0,0) + 1 ), partitionsP.n_rows );
    arma::mat eMat( arma::as_scalar( partitionsP(0,1) - partitionsP(0,0) + 1 ), partitionsP.n_rows );
    for (int j = 0; j < partitionsP.n_rows; ++j) {
      uMat.col(j) = bu.submat( arma::as_scalar(partitionsP(j,0)-1), 0, arma::as_scalar(partitionsP(j,1)-1), 0 );
      eMat.col(j) = arma::diagvec(  Ci.submat( arma::as_scalar(partitionsP(j,0)-1), arma::as_scalar(partitionsP(j,0)-1), arma::as_scalar(partitionsP(j,1)-1), arma::as_scalar(partitionsP(j,1)-1) ) );
    }
    uList(i) = uMat;
    uPevList(i) = eMat;
  }
  // return results in a list form
  return Rcpp::List::create(
    Rcpp::Named("llik") = llik,
    // Rcpp::Named("M") = M,
    Rcpp::Named("W") = W,
    Rcpp::Named("b") = b,
    Rcpp::Named("u") = u,
    Rcpp::Named("bu") = bu,
    Rcpp::Named("Ci") = Ci,
    Rcpp::Named("theta") = theta,
    Rcpp::Named("theta_se") = InfMatInv,
    Rcpp::Named("InfMat") = InfMat, //InfMat,
    Rcpp::Named("monitor") = monitor,
    // Rcpp::Named("constraints") = thetaCUnlistedFinal,
    Rcpp::Named("uList") = uList,
    Rcpp::Named("uPevList") = uPevList,
    Rcpp::Named("AIC") = AIC,
    Rcpp::Named("BIC") = BIC,
    Rcpp::Named("convergence") = convergence,
    Rcpp::Named("partitions") = partitions,
    Rcpp::Named("percDelta") = percDelta,
    Rcpp::Named("normMonitor") = normMonitor,
    Rcpp::Named("toBoundary") = toBoundary,
    Rcpp::Named("dLu") = dLuOut,
    Rcpp::Named("Cchol") = Mchol_XZ
    // Rcpp::Named("PMWu_mat") = PMWu_mat, 
    // Rcpp::Named("MWuchol") = MWuchol
  );
  
}

// [[Rcpp::export]]
Rcpp::List MNR(const arma::mat & Y, const Rcpp::List & X,
               const Rcpp::List & Gx,
               const Rcpp::List & Z, const Rcpp::List & K,
               const Rcpp::List & R, const Rcpp::List & Ge,
               const Rcpp::List & GeI, const arma::mat & W, const bool & isInvW,
               int iters, double tolpar, double tolparinv,
               const bool & ai, const bool & pev,
               const bool & verbose,const bool & retscaled,
               const arma::vec & stepweight, // const arma::vec & emupdate,
               const arma::vec & emweight) {
  
  time_t before = time(0);
  localtime(&before);
  
  int n_fixed = X.size(); // define nre=number of fixed effects
  int n_random = Z.size(); // define nre=number of random effects
  int n_rcov = R.size(); // define nre=number of residual effects
  int n_re = n_random + n_rcov; // define nre=number of total random effects z+r
  int n_traits = Y.n_cols; // define n_traits=number of traits
  int no = Y.n_rows; // define n_traits=number of traits
  arma::vec n_levels(n_re, arma::fill::ones); // to store the number of columns each Z and R matrix has
  // ****************************************************
  // define ZKZ' and R
  // ****************************************************
  // calculate and concatenate ZKZ' and R
  arma::cube ZKZtR(no,no,n_re);
  
  for (int i = 0; i < n_re; ++i) { // for each random effect
    int irw = i - n_random;
    if(i < n_random && n_random > 0){ // if random effect (not residual)
      
      arma::sp_mat zp = Rcpp::as<arma::sp_mat>(Z[i]); // transform as sparse
      n_levels(i) = zp.n_cols; // store the number of columns or levels for this random effect
      bool dcheck = isIdentity_mat(Rcpp::as<arma::mat>(K[i]));
      if(dcheck == true){ // if K[i] is diagonal
        if(zp.n_rows == zp.n_cols){//is a square matrix
          bool dcheck2 = isIdentity_spmat(zp);
          if(dcheck2 == true){ // if Z[i] is diagonal
            ZKZtR.slice(i) = Rcpp::as<arma::mat>(K[i]);
          }else{ZKZtR.slice(i) = zp * zp.t(); }
        }else{ // is a rectangular matrix
          ZKZtR.slice(i) = zp * zp.t();
        }
      }else{ // if K[i] is not diagonal
        if(zp.n_rows == zp.n_cols){//is a square matrix
          bool dcheck2 = isIdentity_spmat(zp);
          if(dcheck2 == true){ // if Z[i] is diagonal
            ZKZtR.slice(i) = Rcpp::as<arma::mat>(K[i]);
          }else{ZKZtR.slice(i) = zp * Rcpp::as<arma::mat>(K[i]) * zp.t(); }
        }else{
          ZKZtR.slice(i) = zp * Rcpp::as<arma::mat>(K[i]) * zp.t();
        }
      }
      
    }else{//if is an rcov term
      // bool dcheck3 = isIdentity_mat(W);
      double dcheck3 = accu(W) - W.n_cols;
      if(dcheck3 == 0){ // if W is diagonal no need to multiply
        ZKZtR.slice(i) = Rcpp::as<arma::sp_mat>(R[irw]);
      }else{ // if W (weights) is not diagonal then multiply Wis R Wis
        // arma::vec ws = W.diag();// 1 / sqrt(diagvec(W));
        if(isInvW == true){ // user has provided a squared and inverted W already
          ZKZtR.slice(i) = W * Rcpp::as<arma::sp_mat>(R[irw]) * W;
        }else{ // user has provided only W
          arma::mat Wis = inv(chol(W));
          ZKZtR.slice(i) = Wis * Rcpp::as<arma::sp_mat>(R[irw]) * Wis.t();
        }
        // arma::vec ws2 = 1/sqrt(ws);// arma::mat Wis = diagmat(ws2); // W inverse squared  // ZKZtR.slice(i) = Wis * Rcpp::as<arma::sp_mat>(R[irw]) * Wis;
      }
      n_levels(i) = ZKZtR.slice(i).n_cols; // store the number of columns or levels for this random effect
    }
  }
  // ****************************************************
  // build multivariate versions of X and Y
  // ****************************************************
  arma::vec Ym = vectorise(Y); // multivariate Y in original scale
  int nom = Ym.n_rows; // number of observations on the vector-form of multivariate Y
  
  arma::mat Xm;
  for (int i = 0; i < n_fixed; ++i) { // for each fixed effect
    if(i==0){ // build multivariate X for 1st fixed effect
      Xm = kron(Rcpp::as<arma::mat>(Gx[i]), Rcpp::as<arma::mat>(X[i]));
    }else{ // build multivariate X for 2nd to nth fixed effect and column bind them
      Xm = arma::join_horiz( Xm , kron(Rcpp::as<arma::mat>(Gx[i]), Rcpp::as<arma::mat>(X[i])) );
    }
  }
  
  arma::mat Ys = scaleCpp(Y); // scale Y using the scaleCpp function made
  arma::vec Ysm = vectorise(Ys); // multivariate Y in scaled form
  // ****************************************************
  // initial VC
  // ****************************************************
  arma::mat base_var = cov(Y); // matrix of original variance-covariance in responses
  arma::mat sc_var = cov(Ys); // matrix of scaled variance-covariance in responses
  int rankX = Xm.n_rows - rank(Xm); // n - p.x
  // VC matrix with dimensions n_traits x n_traits (sigma)
  // we need one for each random effect (n_re)
  arma::cube sigma(n_traits,n_traits,n_re);
  arma::cube sigma_scaled(n_traits,n_traits,n_re);
  
  arma::field<arma::vec> sigma_ut(n_re); // undefined LIST to store the VC in a vector-form with length n_re (#of random effects)
  arma::field<arma::vec> constraintsL(n_re); // undefined LIST to store the constraints in a vector-form with length n_re (#of random effects)
  arma::field<arma::vec> n_levels_multi_traitL(n_re); // undefined LIST to store the n_levels in a vector form
  for (int i = 0; i < n_re; ++i) { // for each random effect fill the cube
    sigma.slice(i) = Rcpp::as<arma::mat>(Ge[i]); // take Ge for a random effect (initial VC values) and save them in a slice
    arma::vec oo = mat_to_vecCpp(sigma.slice(i),GeI[i]) ; // extract upper triangular from that slice in a vector form, pass the constraints as 2nd argument
    sigma_ut[i] = oo; // oo is sigma2 in vector form and stored in the list sigma_ut
    constraintsL[i] = mat_to_vecCpp(GeI[i],GeI[i]) ; // who are diagonal and non-diagonal VCs, pass constraints in list form
    n_levels_multi_traitL[i] = constraintsL[i] ;
  }
  // sigma_ut_un will have all VC for all random effects in a single vector
  arma::vec sigma_ut_un; // vector to unlist the LIST of VC for all random effects
  arma::vec constraints; // vector to unlist constraints
  arma::vec n_levels_multi_trait; // vector to unlist constraints
  for(int i=0; i < n_re ; i++){ // for each random effect unlist
    sigma_ut_un = join_cols(sigma_ut_un,sigma_ut[i]); // column bind vectors so we end up with a very long vector with all VC
    constraints = join_cols(constraints,constraintsL[i]); // column bind vectors so we end up with a very long vector with all constraints
    arma::vec provX = n_levels_multi_traitL[i];
    arma::vec hpos = provX;
    for(int h=0; h < provX.n_cols ; h++){ // for each random effect unlist
      hpos(h) = n_levels(i);
    }
    n_levels_multi_trait = join_cols(n_levels_multi_trait,(provX/provX) % hpos);
  }
  arma::vec sigmaF_ut_un = sigma_ut_un; // make a copy for fixed-value vc's when we use constraints
  arma::vec coef_ut_un = sigma_ut_un; // make a 2nd copy of the same vector for stabilization
  
  int  kk = sigma_ut_un.n_elem; // how many VCs are in the model?
  arma::vec llstore(iters); // container for LL
  arma::vec pos(sigma_ut_un.n_elem, arma::fill::zeros); // create an index vector with as many 0's as VCs
  
  // ****************************************************
  // dummy matrices for multivariate derivatives
  // ****************************************************
  int tot = n_re*n_traits*n_traits; // maximum number of variance components
  arma::vec re_mapper(tot); // mapper to know which VC belongs to each random effect
  arma::cube deriv_dummy(n_traits,n_traits,tot);
  int counter3 = 0;
  for(int i=0; i < n_re; i++){ // for each random effect
    arma::mat prov = Rcpp::as<arma::mat>(GeI[i]);
    int ncol = prov.n_cols; // traits
    
    for (int k = 0; k < ncol; k++){ // go through GeI(i) and make the dummy derivatives where there's a value > 0
      for (int j = 0; j < ncol; j++){
        if (k > j){}else{
          // only extract the variance component if it was planned to be estimated
          if(prov(k,j) > 0){
            arma::mat prov4(ncol,ncol,arma::fill::zeros);
            prov4(k,j)=1;
            prov4 = arma::symmatu(prov4);
            deriv_dummy.slice(counter3) = prov4;
            re_mapper[counter3] = i;
            counter3++;
          }
        }
      }
    }
    
  }
  // ****************************************************
  // ****************************************************
  // ##### iterative algorithm starts
  // ****************************************************
  // ****************************************************
  // Rcpp::List PdViList(kk); // list to store the multivariate derivatives * P or PVi=P*dZKZ'/ds
  
  arma::mat Vi(nom,nom); // V or phenotypic variance matrix
  arma::mat P(nom,nom); // to fill the projection matrix
  arma::sp_mat D = arma::speye<arma::sp_mat>(nom,nom);
  arma::vec popo = arma::vec(rankX, arma::fill::zeros);
  for(int i=0; i < rankX; i++){popo(i) = 1;}
  arma::mat Inf(kk,kk,arma::fill::zeros); // to store second derivatives (information matrix)
  arma::mat InfEM(kk,kk,arma::fill::zeros); // to store second derivatives (information matrix)
  arma::mat InfJoin(kk,kk,arma::fill::zeros); // to store second derivatives (information matrix)
  arma::mat InfJoin_inv(kk,kk,arma::fill::zeros); // to store second derivatives (information matrix)
  
  arma::mat Infw(kk,kk,arma::fill::zeros); // weights for AI information matrix
  arma::mat InfEMw(kk,kk,arma::fill::zeros); // weights for EM information matrix
  
  arma::vec score(kk); // vector to store first derivatives, the product Y'PViPY - tr(PVi) = dL/ds
  arma::vec eigval2; // will be used for the decomposition of P, within the algorithm
  arma::mat sigma_store(sigma_ut_un.n_elem,iters); // to store variance comp through the different iterations
  arma::mat sigma_perc_change(sigma_ut_un.n_elem,iters); // to store percent change of variance components
  arma::mat llik_store(1,iters); // to store llik through the different iterations
  
  arma::mat beta, fitted, residuals; // empty matrices for ..
  Rcpp::List VarU(n_random); // list object for the BLUP variances
  Rcpp::List PevU(n_random); // list object for the BLUP PEVs
  Rcpp::List U(n_random); // list object for the BLUPs
  
  arma::vec vdD(n_traits,arma::fill::ones);
  arma::mat dD = arma::diagmat(vdD);
  arma::mat sigma_cov;
  arma::mat tXVXi; // var-cov fixed effects
  
  bool convergence = false;
  bool last_iteration = false;
  int cycle, cycle2, ikk;
  double ldet, llik, llik0, delta_llik, checkP, seconds; // to store likelihoods and determinants
  // ###############
  // LOOP for cycles
  // ###############
  for(cycle=0; cycle < iters; cycle++){ // for each cycle
    
    for (int i = 0; i < n_re; ++i) {  // for each random effect in the formula
      sigma_ut[i] = mat_to_vecCpp(sigma.slice(i),Rcpp::as<arma::mat>(GeI[i])) ; // extract upper triangular in a vector form
    } // sigma_ut is a LIST
    arma::vec sigmatwo; // create a vector for variance components
    for(int i=0; i < n_re ; i++){ // for each random effect
      sigmatwo = join_cols(sigmatwo,sigma_ut[i]); // column bind to make a vector of vectors
    } // sigmatwo now has all VCs in a vector
    
    // multivariate ZKZ' and V
    arma::mat V(nom,nom); // V or phenotypic variance matrix
    int i;
    for(i=0; i < n_re; i++){ // loop for filling the multivariate ZGZ' and V
      // listGs.slice(i) = prov;
      if(last_iteration == true){ // if is the last iteration multivariate ZKZ is opposite
        if(i == 0){
          V = arma::kron(ZKZtR.slice(i),sigma.slice(i));
        }else{V = V + arma::kron(ZKZtR.slice(i),sigma.slice(i));}
      }else{
        if(i == 0){
          V = arma::kron(sigma.slice(i),ZKZtR.slice(i));
        }else{V = V + arma::kron(sigma.slice(i),ZKZtR.slice(i));}
      }
    }
    // invert V and P (projection matrix)
    
    arma::inv_sympd(Vi,V); // try to invert normally
    if(Vi.n_rows == 0){ // if fails try to invert with diag(1e-3)
      V = V + (D*tolparinv);
      arma::inv_sympd(Vi,V);
      if(Vi.n_rows == 0){// if fails try to invert with diag(1e-2)
        V = V + (D*(tolparinv*10));
        arma::inv_sympd(Vi,V);
        if(Vi.n_rows == 0){ // if fails try to invert with diag(1e-1)
          V = V + (D*(tolparinv*100));
          arma::inv_sympd(Vi,V);
          if(Vi.n_rows == 0){ // finally, if fails try to invert with diag(1e-3)
            // Rcpp::Rcout << "System is singular (V). Stopping the job. Try a bigger number of tolParInv." << arma::endl;
            Rcpp::stop("System is singular (V). Aborting the job. Try a bigger number of tolParInv.");
            // return 0;
          }
        }
      }
    }
    // if last iteration let's make Xm in the opposite direction
    if(last_iteration == true){
      for (int i = 0; i < n_fixed; ++i) {
        if(i==0){
          Xm = kron(Rcpp::as<arma::mat>(X[i]), Rcpp::as<arma::mat>(Gx[i]) );
        }else{
          Xm = arma::join_horiz( Xm , arma::kron( Rcpp::as<arma::mat>(X[i]), Rcpp::as<arma::mat>(Gx[i]) ) );
        }
      }
      Ym = vectorise(Y.t());
      Ysm = vectorise(Ys.t());
    }
    arma::mat VX = Vi * Xm; // VX
    arma::mat tXVX = Xm.t() * VX; // X'VX
    
    arma::mat tXVXVX; // X'VXVX
    arma::solve(tXVXVX,tXVX,VX.t()); // X'VXVX (was computed twice; the status-form call already provides tXVXVX)
    if(tXVXVX.n_rows == 0){ // if fails try to invert with diag(1e-6)
      arma::solve(tXVXVX,tXVX + (D*(tolparinv)),VX.t());
      if(tXVXVX.n_rows == 0){// if fails try to invert with diag(1e-5)
        arma::solve(tXVXVX,tXVX + (D*(tolparinv*10)),VX.t());
        if(tXVXVX.n_rows == 0){ // if fails try to invert with diag(1e-4)
          arma::solve(tXVXVX,tXVX + (D*(tolparinv*100)),VX.t());
          if(tXVXVX.n_rows == 0){ // finally stop
            // Rcpp::Rcout << "System is singular (tXVXVX). Aborting the job. Try a bigger number of tolParInv." << arma::endl;
            Rcpp::stop("System is singular (tXVXVX). Aborting the job. Try a bigger number of tolParInv.");
            // return 0;
          }
        }
      }
    }
    
    // projection matrix
    P = Vi - (VX*tXVXVX); // V - V(XVX)-V
    
    if(last_iteration == false){
      
      arma::vec rss = Ysm.t() * (P * Ysm); // yPy = scalar RSS
      
      double rankXorss = arma::as_scalar(rankX/rss); // (n-p)/y'Py
      double rssorankX = arma::as_scalar(rss/rankX); // y'Py/(n-p)
      
      sigmatwo = sigmatwo * rssorankX;
      
      // weight the projection matrix to provide stability
      coef_ut_un(arma::find(pos == 0)) =  sigmatwo(arma::find(pos == 0)); // VC1[which(pos==0)] = VC2[which(pos==0)]
      coef_ut_un(arma::find(pos == 1)) = log(sigmatwo(arma::find(pos == 1))); // VC1[which(pos==1)] = log(VC2[which(pos==1)])
      
      // calculate the log-likelihood
      P = P * rankXorss; // P * [(n-p)/y'Py]
      rss = rankX; // yPy = n-p
      arma::eig_sym(eigval2, P); // VlV; eigenvectors were never used, values-only is cheaper
      eigval2 = sort(eigval2,"descend"); // sort eigen vectors
      eigval2 = eigval2(arma::find(popo == 1));//(find(seqrankX < rankX)); // only take the values from 1 to
      checkP = eigval2.min();
      if(checkP < 0){ // if any eigen value is < 0 recalculate P
        P = P + (D * (tolpar - eigval2.min())) ;
        eigval2 = eigval2 + tolpar - eigval2.min();
      }
      ldet = accu(log(eigval2)); // sum(log(lambda))
      llik = ldet/2 - (arma::as_scalar(rss)/2); // llik = [sum(log(lambda))/2] - [(n-p)/2]
      
      if(cycle == 0){llik0 = llik;}
      delta_llik = llik - llik0;
      llik0 = llik;
      
      // use the stabilization
      arma::vec var_components(kk, arma::fill::ones); // VC = rep(0,nVC)
      double check00 = accu(pos); // accu is like sum() in R
      if(check00 > 0){  // if there's 1's in the pos vector
        arma::uvec ind = find(pos == 1); // which are 1's
        var_components(ind) = sigmatwo(ind); // var_components[which(pos==1)] = sigmatwo[which(pos==1)]
      }
      
      // calculate first derivatives (dL/ds = score)
      arma::vec Py = P * Ysm; // constant across the vc/AI loops below, hoisted out of them
      
      arma::cube PdViList(nom,nom,kk); // list to store the multivariate derivatives * P or PVi=P*dZKZ'/ds
      for(int i=0; i < kk; i++){
        int re = re_mapper(i);
        arma::mat zkzp = ZKZtR.slice(re); // it repeats the same ZKZtR if is a vc for the same random effect
        arma::mat PdVi = P * kron(deriv_dummy.slice(i),zkzp); // multivariate dVi = dZKZ'/ds
        if(ai && cycle > 2){
          score[i] = - (0.5 * arma::as_scalar(trace(PdVi))) + (0.5 * arma::as_scalar((Ysm.t() * PdVi * Py)));
        }else{
          score[i] = arma::as_scalar(Ysm.t() * PdVi * Py) - accu(diagvec(PdVi));
        }
        PdViList.slice(i) = PdVi;
      }
      // theta(k) * dL/ds  ..... are scalar values
      score = score % var_components; // to be used later for updating the variance components
      // if all goes well var_components is just ones
      
      // calculate second derivatives (AverageInformation)
      // Fisher's Information tr(PVi * PVi) .... A*=Vi=dV/ds .... [Vi Vj'] si sj ; TT is the list of derivatives for all random effects - trait combos
      
      // if(emupdate(cycle) == 0){ // if user wants an EM update (1st derivatives) . It works but it didn't speed up the algorithm when using EM. This leads to don't have information matrix and therefore SE for variance components.
      for (int i = 0; i < kk; i++){
        for (int j = 0; j < kk; j++){
          if (i > j){}else{//only upper triangular
            if(ai && cycle > 2){ // if average information
              Inf(i,j) = 0.5 * arma::as_scalar(Ysm.t() * PdViList.slice(i) * P * PdViList.slice(j) * Py); // j is .t() ?
            }else{ // if newton raphson
              Inf(i,j) = accu(PdViList.slice(i) % PdViList.slice(j).t()) * arma::as_scalar(var_components(i)) * arma::as_scalar(var_components(j));
            }
          }
        }
      }
      Inf = arma::symmatu(Inf); // copy lower in upper triangular
      // Note: Inf_inv (pinv of Inf) used to be computed here purely for a
      // singularity check whose result was never used afterward (the actual
      // update below uses InfJoin_inv) - removed as dead computation.
      
      // vector to store the update = F- * sigma(k) * dL/ds
      arma::vec delta(kk);
      
      // if(emupdate(cycle) == 1){ // if user wants an EM update (1st derivatives)
      InfEM.diag() = (coef_ut_un % coef_ut_un) / n_levels_multi_trait;  // I.em inverse
      InfEM = arma::pinv( InfEM ,  1.490116e-08 ); // I.em
      arma::vec emw(kk); // vectors for weights
      arma::vec aiw(kk);
      for(ikk=0; ikk < kk; ikk++){
        emw(ikk)= emweight(cycle);
        aiw(ikk)= 1 - emweight(cycle);
      }
      Infw.diag() = aiw;  // put weights in diagonal fill::value is still not available in this version
      InfEMw.diag() = emw; //
      InfJoin = (Inf*Infw)+(InfEM*InfEMw); // joint information matrix
      InfJoin_inv = arma::pinv(InfJoin, 1.490116e-08); // inverse the joint information matrix
      delta = InfJoin_inv * score; //update for variance components where: delta = Information.inv * dL/ds
      // delta = (coef_ut_un % score % coef_ut_un)/n_levels; // previous way I was calculating the deltas
      // }else{ // if user wants an information*score update
      //   delta = Inf_inv * score; //update for variance components where: delta = Information.inv * dL/ds
      // }
      
      // ^^^^^^^^^^^^^^^^^^
      // ^^^^^^^^^^^^^^^^^^
      // parameter restrain
      // GeI values
      // 0 not estimated
      // 1 positive
      // 2 unconstrained
      // 3 fixed
      arma::vec coef_ut_unC = coef_ut_un + (stepweight(cycle) * delta); // provisional new variance components
      arma::uvec restrain = find(constraints == 1 && coef_ut_unC < 0); // which vcs are negative and should be positive
      arma::vec cc = coef_ut_unC(restrain); // extract the ones that suppose to be positive
      // arma::vec cc2 = cc(find(cc < 0)); // identify var comp < 0 (1's)
      if(cc.n_elem > 0){ // we have to restrain
        // rest0 = '(';  rest1=cc.n_elem; rest2 = 'restrained)';
        arma::uvec no_restrain = find((constraints == 1 && coef_ut_unC > 0) || (constraints > 1)); // indices of columns that are OK to use (no restrain)
        arma::mat Inf_norestrain = Inf.submat(no_restrain,no_restrain); // subset of Information matrix
        // Note: Inf_norestrain_inv (plain inverse, used only for a singularity
        // check whose result was never used afterward) removed as dead
        // computation - the actual update uses InfJoin_inv_norestrain below.
        arma::vec scorenorestrain = score(no_restrain); // subset of scores (1st derivatives)
        arma::vec deltanorestrain; //  define the delta for no restrained
        
        //
        arma::mat InfEM_norestrain = InfEM.submat(no_restrain,no_restrain); // subset of Information matrix
        arma::mat Infw_norestrain = Infw.submat(no_restrain,no_restrain); // subset of Information matrix
        arma::mat InfEMw_norestrain = InfEMw.submat(no_restrain,no_restrain); // subset of Information matrix
        arma::mat InfJoin_norestrain = InfJoin.submat(no_restrain,no_restrain); // subset of Information matrix
        arma::mat InfJoin_inv_norestrain = InfJoin_inv.submat(no_restrain,no_restrain); // subset of Information matrix
        // if(emupdate(cycle) == 1){ // if user wants an EM update (1st derivatives)
        InfJoin_norestrain = (Inf_norestrain*Infw_norestrain)+(InfEM_norestrain*InfEMw_norestrain); // joint information matrix
        InfJoin_inv_norestrain = arma::pinv(InfJoin_norestrain, 1.490116e-08); // inverse the joint information matrix
        deltanorestrain = InfJoin_inv_norestrain * scorenorestrain; //update for variance components where: delta = Information.inv * dL/ds
        // deltanorestrain = (coef_ut_un_norestrain % scorenorestrain % coef_ut_un_norestrain)/n_levels;
        // }else{ // if user wants an information*score update
        //   deltanorestrain = Inf_norestrain_inv * scorenorestrain; //update variance components
        // }
        delta(no_restrain) = deltanorestrain;
        delta(restrain) = delta(restrain)*0;
        
      }//else just keep going
      // end of parameter restrain
      // ^^^^^^^^^^^^^^^^^^
      // ^^^^^^^^^^^^^^^^^^
      coef_ut_un = coef_ut_un + (stepweight(cycle) * delta);
      //
      // constraint the parameters that should be positive and are going negative
      if(cc.n_elem > 0){
        coef_ut_un(restrain) = coef_ut_un(restrain)*0; // the ones that still go below zero and shouldn't let's fix them
      }
      // weight the projection matrix to provide stability
      sigmatwo(arma::find(pos == 0)) =  coef_ut_un(arma::find(pos == 0)); // index of pos
      sigmatwo(arma::find(pos == 1)) = exp(coef_ut_un(arma::find(pos == 1)));
      // the fixed paramters are forced to be the original value
      sigmatwo(find(constraints == 3)) = sigmaF_ut_un(find(constraints == 3));
      // bring back sigma as a list
      sigma = vec_to_cubeCpp(sigmatwo, GeI);
      // check if likelihood has reached it's maximum and stop if so
      llstore(cycle) = llik;
      // get current time
      time_t now = time(0);
      tm *ltm = localtime(&now);
      // keep track of time difference between iterations
      seconds = difftime(now,before);
      // update before time to be the now to be used in the next iteration
      time_t before = time(0);
      localtime(&before);
      // store paramaters
      sigma_store.col(cycle) = sigmatwo;
      if(cycle > 0){
        sigma_perc_change.col(cycle) = ((coef_ut_un/sigma_store.col(cycle-1))-1) * 100; // percent change
      }
      llik_store(cycle) = llik;
      // return output to the console
      if(verbose == true){ //  arma::cout
        if(cycle == 0){Rcpp::Rcout << "iteration   " << " LogLik   " << "  wall    " << "cpu(sec)   " << "restrained" << arma::endl;}
        Rcpp::Rcout << "    " << cycle+1 << "      " <<  llik << "   " << ltm->tm_hour << ":" << ltm->tm_min << ":" << ltm->tm_sec << "      " << seconds << "           " << cc.n_elem << arma::endl;
      }
      // define the end of the algorithm          std::setprecision(5) <<
      if(((cycle > 2) && (delta_llik < tolpar)) || cycle == iters-1 ){ // tolpar*10
        cycle2 = cycle;
        if((cycle > 2) && (delta_llik < tolpar)){convergence = true;}
        last_iteration = true;
        cycle = iters-2;
        // if user wants to get scaled results we don't bring back to original scale
        if(retscaled == false){
          for (int i = 0; i < n_re; ++i) {
            sigma_scaled.slice(i) = sigma.slice(i); // just make a copy before we scale to normal units
            sigma.slice(i) = (sigma.slice(i)%base_var)/sc_var ;
          }
        }
        
        // Fisher inverse
        arma::mat  FI = Inf/2;
        arma::vec myone(pos.n_elem,arma::fill::ones);
        arma::vec sp = ((sigmatwo - myone) % pos) + myone;
        arma::mat FI_c = FI / (sp * sp.t());
        sigma_cov = pinv(FI_c);
        if(sigma_cov.n_rows == 0){ // if fails
          // Rcpp::Rcout << "System is singular (sigma_cov). Aborting the job." << arma::endl;
          // return 0;
          Rcpp::stop("System is singular (sigma_cov). Aborting the job. Try a bigger number of tolParInv.");
        }
      }
      
    }else{// if we are in the last iteration now we calculate u, PEV, B, XB
      
      arma::inv(tXVXi,tXVX);
      if(tXVXi.n_rows == 0){ // if fails try to invert with diag(1e-6)
        arma::inv(tXVXi,tXVX+(D*(tolparinv)));
        if(tXVXi.n_rows == 0){// if fails try to invert with diag(1e-5)
          arma::inv(tXVXi,tXVX+(D*(tolparinv*10)));
          if(tXVXi.n_rows == 0){
            // Rcpp::Rcout << "System is singular (tXVXi). Aborting the job. Try a bigger number of tolParInv." << arma::endl;
            // return 0;
            Rcpp::stop("System is singular (tXVXi). Aborting the job. Try a bigger number of tolParInv.");
          }
        }
      }
      // arma::vec Ym_rw = vectorise(Y.t());
      if(retscaled == true){// if we have to return scaled results we use Yms
        beta = tXVXi * ((Xm.t() * Vi) * Ysm);
      }else{ // we return in normal scale
        beta = tXVXi * ((Xm.t() * Vi) * Ym);
      }
      
      // beta.reshape(X.n_cols,n_traits);
      fitted = Xm * beta;
      residuals = Ym - fitted;
      // residuals.reshape(no,n_traits);
      // arma::vec residuals_rw = vectorise(residuals.t());
      arma::mat Vie = Vi * residuals;
      if(n_random > 0){
        for(i=0; i < n_random; i++){
          arma::mat Zprov = arma::mat(Rcpp::as<arma::sp_mat>(Z[i]));
          
          // arma::mat Zprov2 = Rcpp::as<arma::mat>(Z[i]);
          arma::mat Ki = Rcpp::as<arma::mat>(K[i]);
          // arma::mat Ki2 = (Zprov.t() * Zprov) * 0; //
          // Ki2.diag() = arma::ones<arma::vec>(Ki2.n_cols);
          // arma::mat Ki2 = arma::mat(arma::speye( Zprov.n_cols, Zprov.n_cols ));
          arma::mat VarK;
          arma::mat ZKfv;
          // double VarKscalar = arma::as_scalar(sigma.slice(i));
          
          // IMPORTANT
          // for rrBLUP models we had to allow a K matrix to be a 1 x 1 matrix so dimensions do not match with Z
          if(Ki.n_cols == Zprov.n_cols){ // if a regular random effect
            // Rcpp::Rcout << "regular" << arma::endl;
            VarK = arma::kron(Rcpp::as<arma::mat>(K[i]),sigma.slice(i)); // Gu * var.u
            ZKfv = VarK * arma::kron(Zprov.t(),dD); // G Z'
          }else{ // if huge matrix from models like rrBLUP we need to create a diagonal to calculate VarK and BLUPs
            // Rcpp::Rcout << "rrBLUP" << arma::endl;
            // VarK = arma::kron(Ki2,sigma.slice(i)); // Gu * var.u
            ZKfv = arma::kron(Zprov.t(),dD*sigma.slice(i)); // G Z'
          }
          
          U(i) = ZKfv * Vie; // BLUP = Z' G Vi (Y - Xb)
          if(pev==true){
            
            if(Ki.n_cols == Zprov.n_cols){ // if a regular random effect
              VarU(i) = ZKfv * (P * ZKfv.t()); // var(u) = Z' G [Vi - (VX*tXVXVX)] G Z'
              PevU(i) = VarK - Rcpp::as<arma::mat>(VarU(i)); // PEV = G - var(u)
            }else{
              VarU(i) = ZKfv * (P * ZKfv.t()); // var(u) = Z' G [Vi - (VX*tXVXVX)] G Z'
              // TO BE FIXED
              // not sure how to get the PEV without constructing VarK due to high-memory requirements in rrBLUP models with potentially millions of SNPs
              PevU(i) = Rcpp::as<arma::mat>(VarU(i)); // PEV = G - var(u)
            }
            
          }
        }
      }
    }
    
  }
  // ****************************************************
  // end of algorithm
  // ****************************************************
  arma::vec dd,ee;
  for (int i = 0; i < n_re; ++i) {
    dd = join_cols(dd,mat_to_vecCpp(base_var,Rcpp::as<arma::mat>(GeI[i]))) ; // extract upper triangular in a vector form
    ee = join_cols(ee,mat_to_vecCpp(sc_var,Rcpp::as<arma::mat>(GeI[i]))) ; // extract upper triangular in a vector form
  }
  arma::mat FISH = (sigma_cov % (dd*dd.t())) / (ee*ee.t()); // bring back to original scale
  // recalculate V and P with original sigma values
  double AIC = (-2 * llik) + (2 * Xm.n_cols);
  double ny = Ym.n_elem;
  double BIC = (-2 * llik) + (log(ny) * Xm.n_cols);
  // monitor
  sigma_store.each_col() %= dd;
  sigma_store.each_col() /= ee;
  arma::mat monitor = join_cols(llik_store,sigma_store);
  // arma::uvec indices(cycle2,arma::fill::ones);
  // arma::mat monitor2 = monitor.cols(find(indices == 1));
  arma::mat monitor2 = monitor.cols(0, cycle2);
  arma::mat sigma_perc_change2;
  if(iters > 1){
    sigma_perc_change2 = sigma_perc_change.cols(1, cycle2); // indicate first and last column to subset to return at the end
  }else{
    sigma_perc_change2 = sigma_perc_change; // indicate first and last column to subset to return at the end
  }
  
  // ****************************************************
  // return the results
  // ****************************************************
  
  return Rcpp::List::create(
    Rcpp::Named("Vi") = Vi,
    Rcpp::Named("P") = P,
    Rcpp::Named("sigma") = sigma,
    Rcpp::Named("sigma_scaled") = sigma_scaled,
    Rcpp::Named("sigmaSE") = FISH,
    Rcpp::Named("Beta") = beta,
    Rcpp::Named("VarBeta") = tXVXi,
    Rcpp::Named("U") = U,
    Rcpp::Named("VarU") = VarU,
    Rcpp::Named("PevU") = PevU,
    Rcpp::Named("fitted") = fitted,
    Rcpp::Named("residuals") = residuals,
    Rcpp::Named("AIC") = AIC,
    Rcpp::Named("BIC") = BIC,
    Rcpp::Named("convergence") = convergence,
    Rcpp::Named("monitor") = monitor2,
    Rcpp::Named("percChange") = sigma_perc_change2,
    Rcpp::Named("dL") = score,
    Rcpp::Named("dL2") = Inf
  );
}

// Reusable arma::sp_mat <-> Eigen::SparseMatrix<double> bridges for the
// pieces of ai_mme_sp2() that are being migrated off Armadillo incrementally.
static inline Eigen::SparseMatrix<double> armaSparseToEigenGlobal(
    const arma::sp_mat & src){
  std::vector<Eigen::Triplet<double>> triplets;
  triplets.reserve(src.n_nonzero);
  for(arma::sp_mat::const_iterator it = src.begin(); it != src.end(); ++it){
    triplets.emplace_back(
      static_cast<int>(it.row()),
      static_cast<int>(it.col()),
      *it
    );
  }
  Eigen::SparseMatrix<double> out(
    static_cast<int>(src.n_rows),
    static_cast<int>(src.n_cols)
  );
  out.setFromTriplets(triplets.begin(), triplets.end());
  out.makeCompressed();
  return out;
}

static inline arma::sp_mat eigenSparseToArmaGlobal(
    const Eigen::SparseMatrix<double> & src){
  const arma::uword nnz = static_cast<arma::uword>(src.nonZeros());
  arma::umat locations(2, nnz);
  arma::vec values(nnz);
  arma::uword k = 0;
  for(int col = 0; col < src.outerSize(); ++col){
    for(Eigen::SparseMatrix<double>::InnerIterator it(src, col); it; ++it){
      locations(0,k) = static_cast<arma::uword>(it.row());
      locations(1,k) = static_cast<arma::uword>(it.col());
      values(k) = it.value();
      ++k;
    }
  }
  return arma::sp_mat(
    locations,
    values,
    static_cast<arma::uword>(src.rows()),
    static_cast<arma::uword>(src.cols())
  );
}

// Sparse LDLT log-determinant of a sparse SPD precision/relationship matrix,
// avoiding the O(n^2) densification that arma::log_det() would require.
static inline double sparseSpdLogDet(const arma::sp_mat & A){
  const Eigen::SparseMatrix<double> Ae = armaSparseToEigenGlobal(A);
  Eigen::SimplicialLDLT<
    Eigen::SparseMatrix<double>,
    Eigen::Lower,
    SommerSparseOrdering
  > factor;
  factor.compute(Ae);
  if(factor.info() != Eigen::Success){
    Rcpp::stop("Sparse LDLT factorisation failed while computing a log-determinant.");
  }
  const Eigen::VectorXd pivots = factor.vectorD();
  double logDet = 0.0;
  for(Eigen::Index k = 0; k < pivots.size(); ++k){
    const double pivot = pivots(k);
    if(!std::isfinite(pivot) || pivot <= 0.0){
      Rcpp::stop("Non-positive LDLT pivot while computing a sparse log-determinant.");
    }
    logDet += std::log(pivot);
  }
  return logDet;
}

// Dense SPD inverse via Eigen::LLT, used in place of arma::inv_sympd() for
// the small per-random-effect covariance matrices in ai_mme_sp2().
static inline bool eigenSpdInverse(const arma::mat & A, arma::mat & out){
  const arma::uword n = A.n_rows;
  // Armadillo and Eigen dense matrices are both column-major, so the
  // conversion is a straight buffer copy rather than an element loop.
  const Eigen::Map<const Eigen::MatrixXd> Ae(A.memptr(), n, n);
  Eigen::LLT<Eigen::MatrixXd> llt(Ae);
  if(llt.info() != Eigen::Success){
    return false;
  }
  const Eigen::MatrixXd inv = llt.solve(Eigen::MatrixXd::Identity(n, n));
  out.set_size(n, n);
  std::copy(inv.data(), inv.data() + inv.size(), out.memptr());
  return true;
}

// [[Rcpp::export]]
Rcpp::List ai_mme_sp2(const arma::sp_mat & X, const Rcpp::List & ZI,
                     const arma::vec & Zind, const Rcpp::List & AiI,
                     const arma::sp_mat & y0,
                     const arma::sp_mat & H, const bool & useH,
                     const arma::uvec & residualBlockI,
                     const arma::uvec & residualIndexI,
                     int nIters, double tolParConvLL, double tolParConvNorm,
                     double tolParInv, const Rcpp::List & covStructI,
                     const arma::vec & weightEmInf,
                     const arma::vec & weightInf, const bool & verbose,
                     const int & computeCi = 0,
                     const std::string & solver = "ldlt",
                     const double & pcgTol = 1.0e-8,
                     const int & pcgMaxIters = 0,
                     const int & pcgTraceProbes = 8,
                     const int & pcgLanczosSteps = 20,
                     const bool & reml = true
){

  if(computeCi < 0 || computeCi > 2){
    Rcpp::stop(
      "computeCi must be 0, 1, or 2: "
      "0 = no final C-inverse work; "
      "1 = Takahashi sparse inverse subset; "
      "2 = full C inverse."
    );
  }

  // Linear-system backend used for C x = rhs solves.
  // In solver="pcg" mode the MME matrix C is NEVER factorised by LDLT:
  // log|C| is estimated by deterministic-probe stochastic Lanczos quadrature
  // and C^{-1} trace terms by Hutchinson probes solved with PCG.
  std::string solverName = solver;
  std::transform(solverName.begin(), solverName.end(), solverName.begin(),
                 [](unsigned char c){ return static_cast<char>(std::tolower(c)); });
  if(solverName != "ldlt" && solverName != "pcg" && solverName != "cholmod"){
    Rcpp::stop("solver must be 'ldlt', 'pcg', or 'cholmod'.");
  }
  if(!std::isfinite(pcgTol) || pcgTol <= 0.0){
    Rcpp::stop("pcgTol must be positive and finite.");
  }
  if(pcgMaxIters < 0){
    Rcpp::stop("pcgMaxIters must be >= 0; 0 selects an automatic limit.");
  }
  if(pcgTraceProbes < 1){
    Rcpp::stop("pcgTraceProbes must be >= 1.");
  }
  if(pcgLanczosSteps < 2){
    Rcpp::stop("pcgLanczosSteps must be >= 2.");
  }
  if(solverName == "pcg" && computeCi == 1){
    Rcpp::stop("computeCi=1 requires LDLT/Takahashi. Use computeCi=0 for a genuinely factorisation-free PCG fit, or computeCi=2 for an explicit PCG full inverse (small systems only).");
  }
  if(!reml && solverName != "ldlt" && solverName != "cholmod"){
    Rcpp::stop("reml=FALSE (maximum likelihood) currently requires solver='ldlt' or solver='cholmod'.");
  }

  if(verbose){
#ifdef _OPENMP
    Rcpp::Rcout
      << "OpenMP available: up to "
      << omp_get_max_threads()
      << " threads."
      << arma::endl;
#else
    Rcpp::Rcout
      << "OpenMP unavailable: this build runs solver kernels serially."
      << arma::endl;
#endif
  }

  time_t before = time(0);
  localtime(&before);

  const int nZs = ZI.size();
  const int nRRe = covStructI.size();
  if(nRRe < 1){
    Rcpp::stop("At least one covariance descriptor (the residual structure) is required.");
  }
  const int nRe = nRRe - 1;
  const int nZsFake = 1;
  const int nReFake = 1;
  const int nX = X.n_cols;
  const int nR = y0.n_rows;

  if(residualBlockI.n_elem != static_cast<arma::uword>(nR) ||
     residualIndexI.n_elem != static_cast<arma::uword>(nR)){
    Rcpp::stop("Residual block/index vectors must have one entry per observation.");
  }

  double vary = arma::mean(arma::var(y0));
  double stdy = arma::mean(arma::stddev(y0));
  double muy = arma::mean(arma::mean(y0));
  if(!std::isfinite(vary) || vary <= 0.0){
    Rcpp::stop("Response variance must be positive and finite.");
  }

  arma::mat responseScaled = arma::mat(y0);
  const arma::rowvec responseMeans = arma::mean(responseScaled, 0);
  const arma::rowvec responseSds = arma::stddev(responseScaled, 0, 0);
  responseScaled.each_row() -= responseMeans;
  responseScaled.each_row() /= responseSds;
  arma::sp_mat y = arma::sp_mat(responseScaled);
  bool intercept = false;
  if(X.n_cols > 0 && arma::accu(X.col(0)) == X.n_rows){
    intercept = true;
  }

  // ------------------------------------------------------------------
  // Generic arbitrary-Kronecker covariance descriptor interface.
  //
  // Every random/residual covariance structure is represented as
  //
  //   Sigma = exp(log_sigma2) * K1 (x) K2 (x) ... (x) Km
  //
  // with unconstrained working parameters. Each factor is a universal
  // CovarianceFactor descriptor carrying evaluator, derivative, reporting,
  // and trust-cap specifications. Statistical model names are not consumed
  // by the REML/MME optimizer.
  // ------------------------------------------------------------------
  arma::field<arma::mat> theta(nRRe), thetaC(nRRe);
  arma::field<arma::vec> covPar(nRRe);
  arma::field<arma::vec> covConstraint(nRRe);
  arma::field<arma::vec> covScale(nRRe);
  arma::field<arma::vec> covLower(nRRe);
  arma::field<arma::vec> covUpper(nRRe);

  std::vector<std::string> covType(
    static_cast<std::size_t>(nRRe),
    "kron"
  );

  std::vector<Rcpp::List> covDescriptor;
  covDescriptor.reserve(static_cast<std::size_t>(nRRe));

  // ================================================================
  // Generic covariance engine: native backend
  // ================================================================
  // This layer owns the optimized built-in covariance primitives.  The
  // REML/MME solver below never dispatches on statistical model names; it
  // only calls evaluateFactor() / derivativeFactor().
  auto evalNativeFactor =
    [&](const Rcpp::List & f,
        const arma::vec & localPar,
        const std::string & op) -> arma::mat {

      const arma::uword q =
        static_cast<arma::uword>(
          Rcpp::as<int>(f["dim"])
        );

      if(op == "identity"){
        return arma::eye<arma::mat>(q,q);
      }

      if(op == "diag"){
        if(localPar.n_elem + 1 != q){
          Rcpp::stop("Malformed diagonal covariance factor.");
        }
        arma::vec d(q, arma::fill::ones);
        for(arma::uword k = 0; k < localPar.n_elem; ++k){
          d(k+1) = std::exp(localPar(k));
        }
        return arma::diagmat(d);
      }

      if(op == "ar1"){
        if(localPar.n_elem != 1){
          Rcpp::stop("Malformed AR1 covariance factor.");
        }
        const double rho = std::tanh(localPar(0));
        arma::mat K(q,q,arma::fill::zeros);
        for(arma::uword i = 0; i < q; ++i){
          for(arma::uword j = 0; j < q; ++j){
            const arma::uword d = (i > j ? i-j : j-i);
            K(i,j) = std::pow(rho, static_cast<double>(d));
          }
        }
        return K;
      }

      if(op == "us"){
        Rcpp::IntegerVector rr = f["us_row"];
        Rcpp::IntegerVector cc = f["us_col"];
        Rcpp::LogicalVector dd = f["us_diag"];

        if(localPar.n_elem != static_cast<arma::uword>(rr.size()) ||
           rr.size() != cc.size() || rr.size() != dd.size()){
          Rcpp::stop("Malformed unstructured covariance factor.");
        }

        arma::mat L(q,q,arma::fill::zeros);
        L(0,0) = 1.0;

        for(arma::uword k = 0; k < localPar.n_elem; ++k){
          const arma::uword i =
            static_cast<arma::uword>(rr[static_cast<int>(k)] - 1);
          const arma::uword j =
            static_cast<arma::uword>(cc[static_cast<int>(k)] - 1);
          if(dd[static_cast<int>(k)]){
            L(i,j) = std::exp(localPar(k));
          }else{
            L(i,j) = localPar(k);
          }
        }

        return L * L.t();
      }


      if(op == "cor_uniform"){
        if(q < 2 || localPar.n_elem != 1){
          Rcpp::stop("Malformed compound-symmetry/uniform-correlation factor.");
        }

        const double lo =
          -1.0
          /
          static_cast<double>(q - 1);

        const double eta = localPar(0);

        const double s =
          eta >= 0.0
          ?
          1.0 / (1.0 + std::exp(-eta))
          :
          std::exp(eta) / (1.0 + std::exp(eta));

        const double rho =
          lo
          +
          (1.0 - lo) * s;

        arma::mat K(
          q,
          q,
          arma::fill::value(rho)
        );

        K.diag().ones();

        return K;
      }

      if(op == "corh"){
        if(q < 2 || localPar.n_elem != q){
          Rcpp::stop("Malformed heterogeneous uniform-correlation factor.");
        }

        const double lo =
          -1.0
          /
          static_cast<double>(q - 1);

        const double eta = localPar(0);

        const double s =
          eta >= 0.0
          ?
          1.0 / (1.0 + std::exp(-eta))
          :
          std::exp(eta) / (1.0 + std::exp(eta));

        const double rho =
          lo
          +
          (1.0 - lo) * s;

        arma::mat C(
          q,
          q,
          arma::fill::value(rho)
        );

        C.diag().ones();

        arma::vec variances(
          q,
          arma::fill::ones
        );

        for(arma::uword k = 1; k < q; ++k){
          variances(k) =
            std::exp(
              localPar(k)
            );
        }

        arma::vec sd =
          arma::sqrt(
            variances
          );

        return
          arma::diagmat(sd)
          *
          C
          *
          arma::diagmat(sd);
      }

      if(op == "arp"){
        const int order =
          Rcpp::as<int>(
            f["order"]
          );

        if(order < 1 ||
           localPar.n_elem != static_cast<arma::uword>(order) ||
           q <= static_cast<arma::uword>(order)){
          Rcpp::stop("Malformed AR(p) covariance factor.");
        }

        // Reflection/PACF coefficients -> stable AR coefficients.
        arma::vec pacf(
          static_cast<arma::uword>(order),
          arma::fill::zeros
        );

        for(int j = 0; j < order; ++j){
          pacf(static_cast<arma::uword>(j)) =
            std::tanh(
              localPar(static_cast<arma::uword>(j))
            );
        }

        arma::vec phi;

        for(int m = 1; m <= order; ++m){

          arma::vec next(
            static_cast<arma::uword>(m),
            arma::fill::zeros
          );

          next(static_cast<arma::uword>(m-1)) =
            pacf(static_cast<arma::uword>(m-1));

          if(m > 1){

            for(int j = 0; j < m-1; ++j){

              next(static_cast<arma::uword>(j)) =
                phi(static_cast<arma::uword>(j))
                -
                pacf(static_cast<arma::uword>(m-1))
                *
                phi(static_cast<arma::uword>(m-2-j));
            }
          }

          phi = next;
        }

        // Solve Yule-Walker equations for rho_1,...,rho_p with rho_0=1.
        arma::mat A(
          static_cast<arma::uword>(order),
          static_cast<arma::uword>(order),
          arma::fill::zeros
        );

        arma::vec b(
          static_cast<arma::uword>(order),
          arma::fill::zeros
        );

        for(int kk = 1; kk <= order; ++kk){

          A(
            static_cast<arma::uword>(kk-1),
            static_cast<arma::uword>(kk-1)
          ) += 1.0;

          for(int jj = 1; jj <= order; ++jj){

            const int d =
              std::abs(
                kk - jj
              );

            const double pj =
              phi(
                static_cast<arma::uword>(jj-1)
              );

            if(d == 0){

              b(
                static_cast<arma::uword>(kk-1)
              ) += pj;

            }else{

              A(
                static_cast<arma::uword>(kk-1),
                static_cast<arma::uword>(d-1)
              ) -= pj;
            }
          }
        }

        arma::vec rInitial;

        bool ok =
          arma::solve(
            rInitial,
            A,
            b
          );

        if(!ok || !rInitial.is_finite()){
          Rcpp::stop("Unable to solve Yule-Walker equations for AR(p) factor.");
        }

        arma::vec rho(
          q,
          arma::fill::zeros
        );

        rho(0) = 1.0;

        for(int h = 1; h <= order; ++h){
          rho(static_cast<arma::uword>(h)) =
            rInitial(static_cast<arma::uword>(h-1));
        }

        for(arma::uword h = static_cast<arma::uword>(order+1);
            h < q;
            ++h){

          double value = 0.0;

          for(int jj = 1; jj <= order; ++jj){

            value +=
              phi(static_cast<arma::uword>(jj-1))
              *
              rho(
                h
                -
                static_cast<arma::uword>(jj)
              );
          }

          rho(h) = value;
        }

        arma::mat K(
          q,
          q,
          arma::fill::zeros
        );

        for(arma::uword i = 0; i < q; ++i){
          for(arma::uword j = 0; j < q; ++j){
            const arma::uword d =
              i > j
              ?
              i-j
              :
              j-i;

            K(i,j) =
              rho(d);
          }
        }

        return K;
      }

      if(op == "ma"){
        const int order =
          Rcpp::as<int>(
            f["order"]
          );

        if(order < 1 ||
           localPar.n_elem != static_cast<arma::uword>(order) ||
           q <= static_cast<arma::uword>(order)){
          Rcpp::stop("Malformed MA(q) covariance factor.");
        }

        arma::vec coef(
          static_cast<arma::uword>(order + 1),
          arma::fill::zeros
        );

        coef(0) = 1.0;

        for(int j = 1; j <= order; ++j){
          coef(static_cast<arma::uword>(j)) =
            localPar(static_cast<arma::uword>(j-1));
        }

        arma::vec rho(
          q,
          arma::fill::zeros
        );

        double gamma0 = 0.0;

        for(int j = 0; j <= order; ++j){
          const double c =
            coef(static_cast<arma::uword>(j));
          gamma0 += c*c;
        }

        if(!std::isfinite(gamma0) || gamma0 <= 0.0){
          Rcpp::stop("Invalid MA covariance normalization.");
        }

        rho(0) = 1.0;

        for(int h = 1; h <= order; ++h){

          double gamma = 0.0;

          for(int j = 0; j <= order-h; ++j){

            gamma +=
              coef(static_cast<arma::uword>(j))
              *
              coef(static_cast<arma::uword>(j+h));
          }

          rho(static_cast<arma::uword>(h)) =
            gamma / gamma0;
        }

        arma::mat K(
          q,
          q,
          arma::fill::zeros
        );

        for(arma::uword i = 0; i < q; ++i){
          for(arma::uword j = 0; j < q; ++j){

            const arma::uword d =
              i > j
              ?
              i-j
              :
              j-i;

            K(i,j) =
              d <= static_cast<arma::uword>(order)
              ?
              rho(d)
              :
              0.0;
          }
        }

        return K;
      }

      if(op == "corg"){
        Rcpp::IntegerVector rr =
          f["corg_row"];

        Rcpp::IntegerVector cc =
          f["corg_col"];

        if(localPar.n_elem != static_cast<arma::uword>(rr.size()) ||
           rr.size() != cc.size()){
          Rcpp::stop("Malformed general-correlation factor.");
        }

        arma::mat A(
          q,
          q,
          arma::fill::eye
        );

        for(arma::uword k = 0; k < localPar.n_elem; ++k){

          const arma::uword i =
            static_cast<arma::uword>(
              rr[static_cast<int>(k)] - 1
            );

          const arma::uword j =
            static_cast<arma::uword>(
              cc[static_cast<int>(k)] - 1
            );

          A(i,j) =
            localPar(k);
        }

        arma::mat S =
          A * A.t();

        arma::vec sd =
          arma::sqrt(
            S.diag()
          );

        arma::mat denom = sd * sd.t();
        arma::mat K = S % arma::pow(denom, -1.0);

        K.diag().ones();

        return
          0.5
          *
          (
            K
            +
            K.t()
          );
      }


      if(op == "fa"){
        const int order =
          Rcpp::as<int>(
            f["order"]
          );

        const int nload =
          Rcpp::as<int>(
            f["fa_nload"]
          );

        Rcpp::IntegerVector rr =
          f["fa_row"];

        Rcpp::IntegerVector cc =
          f["fa_col"];

        Rcpp::LogicalVector dd =
          f["fa_diag"];

        if(order < 1 ||
           nload < 1 ||
           rr.size() != nload ||
           cc.size() != nload ||
           dd.size() != nload ||
           localPar.n_elem !=
             static_cast<arma::uword>(
               nload
               +
               static_cast<int>(q)
               -
               1
             )){
          Rcpp::stop("Malformed factor-analytic covariance factor.");
        }

        arma::mat L(
          q,
          static_cast<arma::uword>(order),
          arma::fill::zeros
        );

        double referenceLogLoading = 0.0;

        for(int a = 0; a < nload; ++a){
          if(rr[a] == 1 && cc[a] == 1 && dd[a]){
            referenceLogLoading =
              localPar(static_cast<arma::uword>(a));
            break;
          }
        }

        const double twiceReferenceLogLoading =
          2.0 * referenceLogLoading;

        const double logScale =
          twiceReferenceLogLoading > 0.0
          ?
          twiceReferenceLogLoading
          +
          std::log1p(std::exp(-twiceReferenceLogLoading))
          :
          std::log1p(std::exp(twiceReferenceLogLoading));

        const double inverseReferenceSd =
          std::exp(-0.5 * logScale);

        for(int a = 0; a < nload; ++a){

          const arma::uword i =
            static_cast<arma::uword>(
              rr[a] - 1
            );

          const arma::uword j =
            static_cast<arma::uword>(
              cc[a] - 1
            );

          const double value =
            dd[a]
            ?
            std::exp(
              localPar(
                static_cast<arma::uword>(a)
              )
              -
              0.5 * logScale
            )
            :
            localPar(
              static_cast<arma::uword>(a)
            )
            *
            inverseReferenceSd;

          L(i,j) =
            value;
        }

        arma::vec psi(
          q,
          arma::fill::zeros
        );

        psi(0) =
          std::exp(-logScale);

        for(arma::uword i = 1; i < q; ++i){

          psi(i) =
            std::exp(
              localPar(
                static_cast<arma::uword>(nload)
                +
                i
                -
                1
              )
              -
              logScale
            );
        }

        return
          L * L.t()
          +
          arma::diagmat(
            psi
          );
      }

      if(op == "ante"){
        const int ncoef =
          Rcpp::as<int>(
            f["ante_ncoef"]
          );

        Rcpp::IntegerVector rr =
          f["ante_row"];

        Rcpp::IntegerVector cc =
          f["ante_col"];

        if(ncoef < 1 ||
           rr.size() != ncoef ||
           cc.size() != ncoef ||
           localPar.n_elem !=
             static_cast<arma::uword>(
               ncoef
               +
               static_cast<int>(q)
               -
               1
             )){
          Rcpp::stop("Malformed antedependence covariance factor.");
        }

        arma::mat T(
          q,
          q,
          arma::fill::eye
        );

        for(int a = 0; a < ncoef; ++a){

          const arma::uword i =
            static_cast<arma::uword>(
              rr[a] - 1
            );

          const arma::uword j =
            static_cast<arma::uword>(
              cc[a] - 1
            );

          T(i,j) =
            -localPar(
              static_cast<arma::uword>(a)
            );
        }

        arma::vec innovation(
          q,
          arma::fill::ones
        );

        for(arma::uword i = 1; i < q; ++i){

          innovation(i) =
            std::exp(
              localPar(
                static_cast<arma::uword>(ncoef)
                +
                i
                -
                1
              )
            );
        }

        arma::mat Ti =
          arma::inv(
            arma::trimatl(T)
          );

        arma::mat M =
          Ti
          *
          arma::diagmat(
            innovation
          )
          *
          Ti.t();

        const double scale =
          M(0,0);

        return
          M / scale;
      }


      Rcpp::stop("Unsupported native covariance evaluator opcode: " + op);
      return arma::mat();
    };



  // Universal factor evaluator. Backends are deliberately small and stable:
  //   native : optimized built-in primitive selected by evaluator$op
  //   fixed  : a supplied fixed covariance matrix
  //   R      : arbitrary callback fun(par) -> covariance matrix
  // Adding a future covariance model therefore does not require changing
  // ai_mme_sp2(): it can be expressed through the R backend immediately and
  // promoted to a native opcode later only if performance justifies it.
  auto evalFactor =
    [&](const Rcpp::List & f,
        const arma::vec & localPar) -> arma::mat {

      if(!f.containsElementNamed("evaluator")){
        Rcpp::stop("CovarianceFactor is missing evaluator specification.");
      }

      Rcpp::List spec =
        Rcpp::as<Rcpp::List>(f["evaluator"]);

      const std::string backend =
        Rcpp::as<std::string>(spec["backend"]);

      const arma::uword q =
        static_cast<arma::uword>(Rcpp::as<int>(f["dim"]));

      arma::mat K;

      if(backend == "native"){
        if(!spec.containsElementNamed("op")){
          Rcpp::stop("Native CovarianceFactor evaluator is missing op.");
        }
        const std::string op =
          Rcpp::as<std::string>(spec["op"]);
        K = evalNativeFactor(f, localPar, op);

      }else if(backend == "fixed"){
        if(!spec.containsElementNamed("matrix")){
          Rcpp::stop("Fixed CovarianceFactor evaluator is missing matrix.");
        }
        K = Rcpp::as<arma::mat>(spec["matrix"]);

      }else if(backend == "R"){
        if(!spec.containsElementNamed("fun") || Rf_isNull(spec["fun"])){
          Rcpp::stop("R CovarianceFactor evaluator is missing fun(par).");
        }
        Rcpp::Function fun = spec["fun"];
        SEXP ans = fun(Rcpp::wrap(localPar));
        K = Rcpp::as<arma::mat>(ans);

      }else{
        Rcpp::stop("Unknown CovarianceFactor evaluator backend: " + backend);
      }

      if(K.n_rows != q || K.n_cols != q || !K.is_finite()){
        Rcpp::stop("CovarianceFactor evaluator returned an invalid matrix.");
      }

      K = 0.5 * (K + K.t());
      return K;
    };

  // Universal first-derivative evaluator.  The AI algorithm only needs first
  // covariance derivatives.  Native analytic derivatives are used for cheap
  // high-frequency primitives; arbitrary callbacks may provide dfun(par,k);
  // otherwise the covariance engine differentiates the small factor matrix by
  // a central finite difference.  The likelihood itself is never numerically
  // differentiated.
  auto factorD1 =
    [&](const Rcpp::List & f,
        const arma::vec & localPar,
        const arma::uword k) -> arma::mat {

      if(k >= localPar.n_elem){
        Rcpp::stop("Invalid CovarianceFactor derivative parameter index.");
      }
      if(!f.containsElementNamed("derivative")){
        Rcpp::stop("CovarianceFactor is missing derivative specification.");
      }

      Rcpp::List spec =
        Rcpp::as<Rcpp::List>(f["derivative"]);
      const std::string backend =
        Rcpp::as<std::string>(spec["backend"]);
      const arma::uword q =
        static_cast<arma::uword>(Rcpp::as<int>(f["dim"]));

      if(backend == "native"){
        if(!spec.containsElementNamed("op")){
          Rcpp::stop("Native derivative specification is missing op.");
        }
        const std::string op =
          Rcpp::as<std::string>(spec["op"]);

        if(op == "diag"){
          arma::mat D(q,q,arma::fill::zeros);
          D(k+1,k+1) = std::exp(localPar(k));
          return D;
        }

        if(op == "ar1"){
          if(k != 0 || localPar.n_elem != 1){
            Rcpp::stop("Malformed native AR1 derivative request.");
          }
          const double rho = std::tanh(localPar(0));
          const double drho = 1.0 - rho*rho;
          arma::mat D(q,q,arma::fill::zeros);
          for(arma::uword i = 0; i < q; ++i){
            for(arma::uword j = 0; j < q; ++j){
              const arma::uword d = (i > j ? i-j : j-i);
              if(d > 0){
                D(i,j) =
                  static_cast<double>(d)
                  * std::pow(rho, static_cast<double>(d-1))
                  * drho;
              }
            }
          }
          return D;
        }

        if(op == "us"){
          Rcpp::IntegerVector rr = f["us_row"];
          Rcpp::IntegerVector cc = f["us_col"];
          Rcpp::LogicalVector dd = f["us_diag"];
          if(localPar.n_elem != static_cast<arma::uword>(rr.size()) ||
             rr.size() != cc.size() || rr.size() != dd.size()){
            Rcpp::stop("Malformed native US derivative metadata.");
          }
          arma::mat L(q,q,arma::fill::zeros);
          L(0,0) = 1.0;
          for(arma::uword a = 0; a < localPar.n_elem; ++a){
            const arma::uword i =
              static_cast<arma::uword>(rr[static_cast<int>(a)] - 1);
            const arma::uword j =
              static_cast<arma::uword>(cc[static_cast<int>(a)] - 1);
            L(i,j) = dd[static_cast<int>(a)]
              ? std::exp(localPar(a)) : localPar(a);
          }
          arma::mat dL(q,q,arma::fill::zeros);
          const arma::uword i =
            static_cast<arma::uword>(rr[static_cast<int>(k)] - 1);
          const arma::uword j =
            static_cast<arma::uword>(cc[static_cast<int>(k)] - 1);
          dL(i,j) = dd[static_cast<int>(k)]
            ? std::exp(localPar(k)) : 1.0;
          return dL * L.t() + L * dL.t();
        }

        Rcpp::stop("Unsupported native derivative opcode: " + op);
      }

      if(backend == "R"){
        if(!spec.containsElementNamed("fun") || Rf_isNull(spec["fun"])){
          Rcpp::stop("R derivative specification is missing fun(par,k).");
        }
        Rcpp::Function dfun = spec["fun"];
        SEXP ans = dfun(Rcpp::wrap(localPar), static_cast<int>(k + 1));
        arma::mat D = Rcpp::as<arma::mat>(ans);
        if(D.n_rows != q || D.n_cols != q || !D.is_finite()){
          Rcpp::stop("CovarianceFactor derivative callback returned an invalid matrix.");
        }
        return 0.5 * (D + D.t());
      }

      if(backend == "numeric"){
        double relStep = 1.0e-6;
        if(spec.containsElementNamed("rel_step")){
          relStep = Rcpp::as<double>(spec["rel_step"]);
        }
        if(!std::isfinite(relStep) || relStep <= 0.0){
          Rcpp::stop("Numerical derivative rel_step must be positive and finite.");
        }
        const double h = relStep * (1.0 + std::abs(localPar(k)));
        arma::vec plus = localPar;
        arma::vec minus = localPar;
        plus(k) += h;
        minus(k) -= h;
        return (evalFactor(f, plus) - evalFactor(f, minus)) / (2.0*h);
      }

      if(backend == "none"){
        Rcpp::stop("Derivative requested for a CovarianceFactor with no parameters.");
      }

      Rcpp::stop("Unknown CovarianceFactor derivative backend: " + backend);
      return arma::mat();
    };

  auto evaluateDescriptor =
    [&](const Rcpp::List & cs,
        const arma::vec & par) -> arma::mat {

      if(par.n_elem < 1){
        Rcpp::stop("Covariance descriptor must contain log_sigma2.");
      }

      Rcpp::List factors = cs["factors"];
      arma::mat K(1,1,arma::fill::ones);

      for(int fidx = 0; fidx < factors.size(); ++fidx){
        Rcpp::List f = Rcpp::as<Rcpp::List>(factors[fidx]);

        const int start1 =
          f.containsElementNamed("par_start")
          ?
          Rcpp::as<int>(f["par_start"])
          :
          1;

        const int end1 =
          f.containsElementNamed("par_end")
          ?
          Rcpp::as<int>(f["par_end"])
          :
          0;

        arma::vec localPar;
        if(end1 >= start1){
          localPar =
            par.subvec(
              static_cast<arma::uword>(start1 - 1),
              static_cast<arma::uword>(end1 - 1)
            );
        }

        K = arma::kron(K, evalFactor(f, localPar));
      }

      return std::exp(par(0)) * K;
    };

  auto descriptorD1 =
    [&](const Rcpp::List & cs,
        const arma::vec & par,
        const arma::uword k) -> arma::mat {

      if(k == 0){
        return evaluateDescriptor(cs, par);
      }

      Rcpp::List factors = cs["factors"];
      arma::mat K(1,1,arma::fill::ones);
      bool found = false;

      for(int fidx = 0; fidx < factors.size(); ++fidx){
        Rcpp::List f = Rcpp::as<Rcpp::List>(factors[fidx]);

        const int start1 =
          f.containsElementNamed("par_start")
          ?
          Rcpp::as<int>(f["par_start"])
          :
          1;

        const int end1 =
          f.containsElementNamed("par_end")
          ?
          Rcpp::as<int>(f["par_end"])
          :
          0;

        arma::vec localPar;
        if(end1 >= start1){
          localPar =
            par.subvec(
              static_cast<arma::uword>(start1 - 1),
              static_cast<arma::uword>(end1 - 1)
            );
        }

        arma::mat piece;

        const int k1 = static_cast<int>(k) + 1;
        if(end1 >= start1 && k1 >= start1 && k1 <= end1){
          piece =
            factorD1(
              f,
              localPar,
              static_cast<arma::uword>(k1 - start1)
            );
          found = true;
        }else{
          piece = evalFactor(f, localPar);
        }

        K = arma::kron(K, piece);
      }

      if(!found){
        Rcpp::stop("Covariance derivative parameter does not belong to any factor.");
      }

      return std::exp(par(0)) * K;
    };

  for(int i = 0; i < nRRe; ++i){

    if(Rf_isNull(covStructI[i])){
      Rcpp::stop("NULL covariance descriptor supplied to ai_mme_sp2().");
    }

    Rcpp::List cs =
      Rcpp::as<Rcpp::List>(
        covStructI[i]
      );

    if(!cs.containsElementNamed("type") ||
       Rcpp::as<std::string>(cs["type"]) != "kron"){
      Rcpp::stop("ai_mme_sp2() accepts only the generic type='kron' covariance descriptor.");
    }
    if(!cs.containsElementNamed("descriptor_version") ||
       Rcpp::as<int>(cs["descriptor_version"]) < 2){
      Rcpp::stop("ai_mme_sp2() requires CovarianceFactor descriptor_version >= 2.");
    }

    covDescriptor.push_back(cs);

    covPar(i) =
      Rcpp::as<arma::vec>(
        cs["par"]
      );

    Rcpp::LogicalVector freeR =
      cs["free"];

    if(freeR.size() != static_cast<int>(covPar(i).n_elem)){
      Rcpp::stop("covStruct$free and covStruct$par have inconsistent lengths.");
    }

    covConstraint(i).set_size(covPar(i).n_elem);
    for(arma::uword k = 0; k < covPar(i).n_elem; ++k){
      covConstraint(i)(k) =
        freeR[static_cast<int>(k)]
        ?
        2.0
        :
        3.0;
    }

    covScale(i) =
      arma::vec(
        covPar(i).n_elem,
        arma::fill::zeros
      );

    covLower(i) =
      arma::vec(
        covPar(i).n_elem,
        arma::fill::value(
          -std::numeric_limits<double>::infinity()
        )
      );

    covUpper(i) =
      arma::vec(
        covPar(i).n_elem,
        arma::fill::value(
          std::numeric_limits<double>::infinity()
        )
      );

    // Standardise only the product-level variance:
    // log(sigma2_internal) = log(sigma2_original) - log(var(y)).
    covPar(i)(0) -= std::log(vary);
    
    // Give the product-level log(sigma2) working parameter (always
    // covPar(i)(0), see vsm()) the same variance floor legacy type-1
    // parameters receive via applyVarianceBounds()/vcFloor, converted to
    // this working/standardised scale. Without this, a boundary-approaching
    // variance component (true value near/at zero) has no lower bound at
    // all, so an unconstrained Newton/AI step can jump arbitrarily far
    // past the PD-check floor in one iteration; the structure then gets
    // flagged "not PD" and, for single-parameter structures, permanently
    // frozen (see the covPar(iStruct).n_elem==1 fallback below).
    covLower(i)(0) =
      std::log(std::max(1.0e-8, tolParInv)) - std::log(vary);

    theta(i) =
      evaluateDescriptor(
        cs,
        covPar(i)
      );

    thetaC(i) =
      arma::zeros<arma::mat>(
        theta(i).n_rows,
        theta(i).n_cols
      );
  }

  auto covarianceD1 =
    [&](const int iStruct,
        const arma::uword k) -> arma::mat {

      return
        descriptorD1(
          covDescriptor[static_cast<std::size_t>(iStruct)],
          covPar(iStruct),
          k
        );
    };

  std::vector< std::vector<arma::mat> > covarianceDerivativeCache(
    static_cast<std::size_t>(nRRe)
  );

  std::vector< std::vector<bool> > covarianceDerivativeReady(
    static_cast<std::size_t>(nRRe)
  );

  auto cachedCovarianceD1 =
    [&](const int iStruct,
        const arma::uword k) -> const arma::mat & {

      const std::size_t structureOffset =
        static_cast<std::size_t>(iStruct);

      if(!covarianceDerivativeReady[structureOffset][k]){
        covarianceDerivativeCache[structureOffset][k] =
          covarianceD1(iStruct, k);
        covarianceDerivativeReady[structureOffset][k] = true;
      }

      return covarianceDerivativeCache[structureOffset][k];
    };

  auto evaluateStructure =
    [&](const int iStruct,
        const arma::vec & par) -> arma::mat {

      return
        evaluateDescriptor(
          covDescriptor[static_cast<std::size_t>(iStruct)],
          par
        );
    };

  auto tryEvaluateStructure =
    [&](const int iStruct,
        const arma::vec & par,
        arma::mat & value) -> bool {

      if(!par.is_finite()){
        return false;
      }

      try{
        value =
          evaluateStructure(
            iStruct,
            par
          );
      }catch(const std::exception &){
        return false;
      }

      return
        value.n_rows
        ==
        theta(iStruct).n_rows
        &&
        value.n_cols
        ==
        theta(iStruct).n_cols
        &&
        value.is_finite();
    };

  // move Z to sparse arma objects
  int nZsAl; // integer to define the allocation of Z
  if(nZs > 0){
    nZsAl = nZs; // if there's random effects the nZs to allocate is equal to Z.size
  }else{
    nZsAl = nZsFake; // otherwise at least we allocate 1 element to avoid the program to crash
  }
  arma::field<arma::sp_mat> Z(nZsAl); // allocate size of Z
  if(nZs > 0){ // if there's random effects
    for (int i = 0; i < nZs; ++i) { // for each Z
      Z(i)=convertSparse(ZI(i)); // convert the matrix to sparse and store in the field
    }
  }
  // delete ZI;
  // Residual covariance is now built directly from residual block/local-index
  // metadata plus the same generic Kronecker descriptor used for random G.
  // No residual basis-list expansion is required.
  // move Ai to sparse arma objects
  int nReAl;
  if(nZs > 0){
    nReAl = nRe;
  }else{
    nReAl = nReFake;
  }
  arma::field<arma::sp_mat> Ai(nReAl); // allocate size of Ai field
  if(nZs > 0){
    for (int i = 0; i < nRe; ++i) {
      Ai(i)=convertSparse(AiI(i)); // convert the matrix to sparse and store in the field
    }
  }
  // calculate log determinants of Ai's
  arma::rowvec logDetA(nReAl);
  if(nZs > 0){ // of there's random effects
    for (int i = 0; i < nRe; ++i) { // for each random effect
      // Sparse LDLT log-determinant: avoids densifying the (potentially
      // large) relationship/pedigree precision matrix Ai.
      logDetA(i) = -sparseSpdLogDet(Ai(i));
    }
  }
  
  // define partitions (only used if random effects exist)
  int last = X.n_cols;
  arma::field<arma::mat> partitions(nReAl); // store indices of the random effects
  int Nu = 0;
  if(nZs > 0){ //if there's random effects (Z matrices) check where each starts and ends
    for (int i = 0; i < nRe; ++i) { // for each effect
      arma::uvec indexZind = find(Zind == (i+1) ); // which Z matrices to use , +1 because of the way indeces are used in C++
      int nIndexZind = indexZind.size(); //  number of Z matrices to use
      arma::vec Nus(nIndexZind); // vector to store number of columns in each Z matrix
      // for each matrix in this random effect
      for (int j = 0; j < nIndexZind; ++j) {
        int jj = indexZind(j); // thake the jj matrix
        arma::sp_mat Zprov = Z(jj); // put it in a provisional object
        Nus(j)=Zprov.n_cols; // calculate the number of columns
      }
      arma::vec end = Nus; // define ends and starts
      for (int k = 0; k < nIndexZind; ++k) { // for each effect
        arma::uvec toSum = arma::regspace<arma::uvec>(0,  1,  k); // equivalent to seq()
        end(k)=arma::accu(Nus(toSum));
      }
      arma::vec ones(nIndexZind, arma::fill::ones);
      arma::vec lastM(nIndexZind, arma::fill::value(last));
      arma::vec start = end - Nus + ones;
      start = start + lastM;// adjust start by adding # of fixed effects
      end = end + lastM;//adjust end by adding # of fixed effects
      partitions(i) = arma::join_rows(start,end);
      last = end.max();
      Nu = Nu + accu(Nus);
    }
  }// end of if statement when random effects exist
  
  // define the number of optimizer parameters per covariance structure
  arma::vec nVc(nRRe);
  for (int i = 0; i < nRRe; ++i) {
    nVc(i) =
      static_cast<double>(
        covPar(i).n_elem
      );
  }
  int nVcTotal = accu(nVc); // total number of variance components
  // assign a start and an end index to each covariance structure using the #of VC
  arma::vec nVcEnd = nVc;
  for (int i = 0; i < nRRe; ++i) {
    arma::uvec toSum = arma::regspace<arma::uvec>(0,  1,  i); // equivalent to seq()
    nVcEnd(i)=arma::accu(nVc(toSum));
  }
  arma::vec nVcStart = nVcEnd - nVc + 1;
  // removing complex structures how many effects are really there
  arma::vec nUsTotal(nReAl);
  if(nZs > 0){ //
    for (int i = 0; i < nRe; ++i) {
      arma::mat partitionsProv = partitions(i);
      nUsTotal(i) = partitionsProv(0,1) - partitionsProv(0,0) + 1;
    }
  }
  // define objects to store theta and llik across iterations
  arma::mat monitor(nVcTotal,nIters); // matrix to store variance components
  arma::rowvec llik(nIters); // store log likellihood values
  
  int nEffects = Nu+nX;
  // Stage 2 Eigen migration: W and C are the per-iteration hot-path sparse
  // objects (assembled every iteration and factorised every iteration).
  // They are native Eigen sparse matrices so the LDLT/PCG solvers consume
  // them directly, with no per-iteration Armadillo<->Eigen round trip.
  Eigen::SparseMatrix<double> W(nR,nEffects), C(nEffects,nEffects);
  arma::sp_mat Ci(nEffects,nEffects);
  arma::vec u(Nu), b(nX), bu(nEffects);
  arma::mat avInf(nVcTotal,nVcTotal);
  arma::mat emInf(nVcTotal,nVcTotal);
  arma::mat InfMat(nVcTotal,nVcTotal);
  arma::mat InfMatInv(nVcTotal,nVcTotal);
  bool convergence = false;
  double seconds;
  arma::vec delta(nVcTotal), delta_minus1(nVcTotal);
  // objects for constraints
  arma::mat percDelta(nVcTotal,nIters,arma::fill::zeros); // store % change of the delta with respect to the previous iteration
  arma::mat normMonitor(3,nIters); // store in each iteration the 3 stopping criteria of Madsen and Jensen
  arma::mat toBoundary(nIters,nVcTotal, arma::fill::zeros ); // store which values have been set to the boundary value
  arma::vec sumToBoundary(nVcTotal, arma::fill::zeros ); // to apply sum across iterations and if a VC goes to the boundary 3 times it is fixed to the boundary
  arma::sp_mat Hs(H.n_cols,H.n_cols); // square of H matrix
  if(useH == true){ // sparse Cholesky decomposition of H if user wants to use weights
    Rcpp::Rcout << "Using the weights matrix " << arma::endl;
    // Stage 1 Eigen migration: factorise H directly in sparse form instead
    // of densifying it just to call arma::chol().
    const Eigen::SparseMatrix<double> He = armaSparseToEigenGlobal(H);
    Eigen::SimplicialLLT<Eigen::SparseMatrix<double>, Eigen::Upper> HeLLT;
    HeLLT.compute(He);
    if(HeLLT.info() != Eigen::Success){
      Rcpp::stop("Sparse Cholesky factorisation of the weights matrix H failed.");
    }
    Eigen::SparseMatrix<double> HeU = HeLLT.matrixU();
    HeU.makeCompressed();
    Hs = eigenSparseToArmaGlobal(HeU);
  }
  arma::vec dLuOut;//(nVcTotal); // we will join cols

  // ============================================================
  // CHANGE 1: cache objects that never change across iterations
  // ============================================================

  // W = [X Z] is constant.
  {
    std::vector<Eigen::Triplet<double>> Wtriplets;
    const Eigen::SparseMatrix<double> Xeig = armaSparseToEigenGlobal(X);
    Wtriplets.reserve(
      static_cast<std::size_t>(Xeig.nonZeros())
    );
    for(int col = 0; col < Xeig.outerSize(); ++col){
      for(Eigen::SparseMatrix<double>::InnerIterator it(Xeig, col); it; ++it){
        Wtriplets.emplace_back(
          static_cast<int>(it.row()),
          static_cast<int>(it.col()),
          it.value()
        );
      }
    }
    int columnOffset = static_cast<int>(Xeig.cols());
    if(nZs > 0){
      for(int i = 0; i < nZs; ++i){
        const Eigen::SparseMatrix<double> Zeig = armaSparseToEigenGlobal(Z(i));
        Wtriplets.reserve(
          Wtriplets.size() + static_cast<std::size_t>(Zeig.nonZeros())
        );
        for(int col = 0; col < Zeig.outerSize(); ++col){
          for(Eigen::SparseMatrix<double>::InnerIterator it(Zeig, col); it; ++it){
            Wtriplets.emplace_back(
              static_cast<int>(it.row()),
              columnOffset + static_cast<int>(it.col()),
              it.value()
            );
          }
        }
        columnOffset += static_cast<int>(Zeig.cols());
      }
    }
    W.setFromTriplets(Wtriplets.begin(), Wtriplets.end());
    W.makeCompressed();
  }

  std::vector<Eigen::SparseMatrix<double>> AiEigenCache(
    static_cast<std::size_t>(nRe)
  );
  if(solverName == "pcg"){
    for(int iR = 0; iR < nRe; ++iR){
      AiEigenCache[static_cast<std::size_t>(iR)] =
        armaSparseToEigenGlobal(Ai(iR));
      AiEigenCache[static_cast<std::size_t>(iR)].makeCompressed();
    }
  }

  // Dense response vector is constant.
  const arma::vec yDense = arma::vectorise(arma::mat(y));

  // Fixed/random coefficient indexes are constant.
  const arma::uvec bInd =
    arma::regspace<arma::uvec>(
      static_cast<arma::uword>(0),
      static_cast<arma::uword>(1),
      static_cast<arma::uword>(nX - 1)
    );

  arma::uvec uInd;
  if(nZs > 0){
    uInd =
      arma::regspace<arma::uvec>(
        static_cast<arma::uword>(nX),
        static_cast<arma::uword>(1),
        static_cast<arma::uword>(nX + Nu - 1)
      );
  }

  // Cache the Z indexes belonging to every random covariance structure
  // and the coefficient indexes for every partition.
  arma::field<arma::uvec> useZindCache(nReAl);
  std::vector< std::vector<arma::uvec> > partitionIndexCache(
    static_cast<std::size_t>(nRe)
  );
  std::vector< std::vector<arma::uword> > partitionStartCache(
    static_cast<std::size_t>(nRe)
  );
  std::vector< std::vector<arma::uword> > partitionEndCache(
    static_cast<std::size_t>(nRe)
  );

  if(nZs > 0){
    for(int iR = 0; iR < nRe; ++iR){

      useZindCache(iR) =
        arma::find(
          Zind == (iR + 1)
        );

      arma::mat partitionsP =
        partitions(iR);

      partitionIndexCache[static_cast<std::size_t>(iR)].resize(
        static_cast<std::size_t>(partitionsP.n_rows)
      );
      partitionStartCache[static_cast<std::size_t>(iR)].resize(
        static_cast<std::size_t>(partitionsP.n_rows)
      );
      partitionEndCache[static_cast<std::size_t>(iR)].resize(
        static_cast<std::size_t>(partitionsP.n_rows)
      );

      for(arma::uword j = 0; j < partitionsP.n_rows; ++j){

        const arma::uword first =
          static_cast<arma::uword>(
            partitionsP(j,0) - 1
          );

        const arma::uword lastIndex =
          static_cast<arma::uword>(
            partitionsP(j,1) - 1
          );

        partitionStartCache[static_cast<std::size_t>(iR)][static_cast<std::size_t>(j)] =
          first;

        partitionEndCache[static_cast<std::size_t>(iR)][static_cast<std::size_t>(j)] =
          lastIndex;

        partitionIndexCache[static_cast<std::size_t>(iR)][static_cast<std::size_t>(j)] =
          arma::regspace<arma::uvec>(
            first,
            static_cast<arma::uword>(1),
            lastIndex
          );
      }
    }
  }

  // One-time validation that Ai(iR) is dimensionally consistent with every
  // partition block, so the per-iteration G^-1 accumulation loop does not
  // need to repeat these checks on every REML iteration.
  if(nZs > 0){
    for(int iR = 0; iR < nRe; ++iR){
      const std::size_t iCache = static_cast<std::size_t>(iR);
      for(std::size_t j = 0; j < partitionStartCache[iCache].size(); ++j){
        const arma::uword blockSize =
          partitionEndCache[iCache][j] - partitionStartCache[iCache][j] + 1;
        if(blockSize != Ai(iR).n_rows || blockSize != Ai(iR).n_cols){
          Rcpp::stop(
            "Relationship inverse dimensions are inconsistent with a "
            "random-effect MME partition."
          );
        }
      }
    }
  }

  // ============================================================
  // Generic residual block/layout cache
  // ============================================================
  const int residualStructIndex = nRRe - 1;
  const int residualDim =
    static_cast<int>(
      theta(residualStructIndex).n_rows
    );

  if(residualDim < 1){
    Rcpp::stop("Residual covariance dimension must be positive.");
  }

  std::vector<int> residualBlockOfRow(
    static_cast<std::size_t>(nR),
    -1
  );

  std::vector<int> residualLocalOfRow(
    static_cast<std::size_t>(nR),
    -1
  );

  int residualNBlocks = 0;

  for(int rr = 0; rr < nR; ++rr){

    if(residualBlockI(static_cast<arma::uword>(rr)) < 1 ||
       residualIndexI(static_cast<arma::uword>(rr)) < 1){
      Rcpp::stop("Residual block and local indices are 1-based positive integers.");
    }

    const int b =
      static_cast<int>(
        residualBlockI(static_cast<arma::uword>(rr)) - 1
      );

    const int local =
      static_cast<int>(
        residualIndexI(static_cast<arma::uword>(rr)) - 1
      );

    if(local < 0 || local >= residualDim){
      Rcpp::stop("Residual local covariance index exceeds the covariance-product dimension.");
    }

    residualBlockOfRow[static_cast<std::size_t>(rr)] = b;
    residualLocalOfRow[static_cast<std::size_t>(rr)] = local;
    residualNBlocks = std::max(residualNBlocks, b + 1);
  }

  std::vector< std::vector<arma::uword> > residualRowsTmp(
    static_cast<std::size_t>(residualNBlocks)
  );

  for(int rr = 0; rr < nR; ++rr){
    residualRowsTmp[
      static_cast<std::size_t>(
        residualBlockOfRow[static_cast<std::size_t>(rr)]
      )
    ].push_back(
      static_cast<arma::uword>(rr)
    );
  }

  std::vector<arma::uvec> residualKronRows(
    static_cast<std::size_t>(residualNBlocks)
  );

  bool residualKronBlocks = true;

  for(int b = 0; b < residualNBlocks; ++b){

    std::vector< std::pair<int, arma::uword> > ordered;
    ordered.reserve(
      residualRowsTmp[static_cast<std::size_t>(b)].size()
    );

    std::vector<int> seen(
      static_cast<std::size_t>(residualDim),
      0
    );

    for(const arma::uword rr :
        residualRowsTmp[static_cast<std::size_t>(b)]){

      const int local =
        residualLocalOfRow[static_cast<std::size_t>(rr)];

      if(seen[static_cast<std::size_t>(local)] != 0){
        Rcpp::stop("A residual block contains the same covariance-product coordinate more than once.");
      }

      seen[static_cast<std::size_t>(local)] = 1;
      ordered.emplace_back(local, rr);
    }

    std::sort(
      ordered.begin(),
      ordered.end(),
      [](const std::pair<int,arma::uword> & a,
         const std::pair<int,arma::uword> & b){
        return a.first < b.first;
      }
    );

    residualKronRows[static_cast<std::size_t>(b)].set_size(
      ordered.size()
    );

    for(std::size_t j = 0; j < ordered.size(); ++j){
      residualKronRows[static_cast<std::size_t>(b)](
        static_cast<arma::uword>(j)
      ) =
        ordered[j].second;
    }

    if(
        static_cast<int>(ordered.size()) != residualDim
    ){
      residualKronBlocks = false;
    }else{
      for(int local = 0; local < residualDim; ++local){
        if(seen[static_cast<std::size_t>(local)] == 0){
          residualKronBlocks = false;
          break;
        }
      }
    }
  }

  const int residualKronNBlocks =
    residualNBlocks;

  const int residualKronBlockSize =
    residualDim;

  std::vector< std::vector<arma::uword> > residualActiveColumnsTmp(
    static_cast<std::size_t>(residualNBlocks)
  );

  for(int wCol = 0; wCol < W.outerSize(); ++wCol){
    for(Eigen::SparseMatrix<double>::InnerIterator entry(W, wCol); entry; ++entry){
      const int block = residualBlockOfRow[
        static_cast<std::size_t>(entry.row())
      ];
      residualActiveColumnsTmp[static_cast<std::size_t>(block)].push_back(
        static_cast<arma::uword>(entry.col())
      );
    }
  }

  std::vector<arma::uvec> residualActiveColumns(
    static_cast<std::size_t>(residualNBlocks)
  );

  for(int block = 0; block < residualNBlocks; ++block){
    std::vector<arma::uword> & columns =
      residualActiveColumnsTmp[static_cast<std::size_t>(block)];
    std::sort(columns.begin(), columns.end());
    columns.erase(std::unique(columns.begin(), columns.end()), columns.end());
    residualActiveColumns[static_cast<std::size_t>(block)] =
      arma::uvec(columns);
  }

  // Structural diagonality is determined by the covariance factors, not by
  // current parameter values.
  bool residualStructurallyDiagonal = true;

  {
    Rcpp::List rcs =
      covDescriptor[
        static_cast<std::size_t>(
          residualStructIndex
        )
      ];

    Rcpp::List rfactors =
      rcs["factors"];

    for(int fidx = 0; fidx < rfactors.size(); ++fidx){
      Rcpp::List f =
        Rcpp::as<Rcpp::List>(
          rfactors[fidx]
        );

      if(!f.containsElementNamed("structurally_diagonal") ||
         !Rcpp::as<bool>(f["structurally_diagonal"])){
        residualStructurallyDiagonal = false;
        break;
      }
    }
  }

  // Same structural check, per random-effect covariance descriptor. This is
  // a property of the descriptor (which factors it is built from), not of
  // the current parameter values, so it is computed once here rather than
  // every REML iteration. Used to short-circuit the O(q^5) AI-matrix
  // second-derivative loop and the O(q^2) score/trace loop down to their
  // exact closed-form scalar equivalents when a random-effect covariance
  // product is diagonal (e.g. dsm() alone, or a Kronecker product of only
  // diagonal factors).
  std::vector<bool> randomStructurallyDiagonal(
    static_cast<std::size_t>(nRe),
    true
  );

  for(int iR = 0; iR < nRe; ++iR){
    Rcpp::List rcs =
      covDescriptor[static_cast<std::size_t>(iR)];

    Rcpp::List rfactors =
      rcs["factors"];

    bool diagonalHere = true;

    for(int fidx = 0; fidx < rfactors.size(); ++fidx){
      Rcpp::List f =
        Rcpp::as<Rcpp::List>(rfactors[fidx]);

      if(!f.containsElementNamed("structurally_diagonal") ||
         !Rcpp::as<bool>(f["structurally_diagonal"])){
        diagonalHere = false;
        break;
      }
    }

    randomStructurallyDiagonal[static_cast<std::size_t>(iR)] = diagonalHere;
  }

  // Cache whether the optional H Cholesky factor is diagonal.
  bool Hdiag = true;
  arma::vec HdiagSquared(
    nR,
    arma::fill::ones
  );

  if(useH){
    for(arma::sp_mat::const_iterator it = Hs.begin();
        it != Hs.end();
        ++it){

      if(it.row() != it.col() && std::abs(*it) > 0.0){
        Hdiag = false;
      }
    }

    if(Hdiag){
      for(int ii = 0; ii < nR; ++ii){
        const double h = Hs(ii,ii);
        HdiagSquared(ii) = h*h;
      }
    }
  }

  // ------------------------------------------------------------
  // Cache the block design matrices used by the repeated-block/Kronecker
  // residual C-assembly path (including any diagonal-H row weighting).
  // W, Hs, and the block/column membership are all constant across REML
  // iterations, so this "blockW" is identical every iteration; precomputing
  // it once avoids rescanning every active column of W (previously
  // O(sum of active-column nnz over blocks) per iteration) inside the main
  // loop below.
  // ------------------------------------------------------------
  std::vector<arma::mat> residualKronBlockWCache(
    static_cast<std::size_t>(residualNBlocks)
  );

  // Only complete/repeated blocks (residualKronBlocks) guarantee
  // rows.n_elem == residualDim, which the local-coordinate indexing below
  // requires. When blocks are irregular the Rkron path is never taken at
  // runtime, so the cache is simply left empty.
  if(residualKronBlocks){
  for(int block = 0; block < residualNBlocks; ++block){
    const arma::uvec & rows =
      residualKronRows[static_cast<std::size_t>(block)];
    const arma::uvec & columns =
      residualActiveColumns[static_cast<std::size_t>(block)];

    if(columns.n_elem == 0){ continue; }

    arma::mat blockW(
      rows.n_elem,
      columns.n_elem,
      arma::fill::zeros
    );

    for(arma::uword localCol = 0;
        localCol < columns.n_elem;
        ++localCol){
      for(Eigen::SparseMatrix<double>::InnerIterator designEntry(
            W, static_cast<int>(columns(localCol)));
          designEntry; ++designEntry){
        const arma::uword globalRow =
          static_cast<arma::uword>(designEntry.row());
        if(residualBlockOfRow[static_cast<std::size_t>(globalRow)] == block){
          const int localCoordinate = residualLocalOfRow[
            static_cast<std::size_t>(globalRow)
          ];
          blockW(static_cast<arma::uword>(localCoordinate), localCol) =
            designEntry.value();
        }
      }
    }

    if(useH && Hdiag){
      arma::vec hDiagonal(rows.n_elem);
      for(arma::uword localRow = 0; localRow < rows.n_elem; ++localRow){
        hDiagonal(localRow) = Hs(rows(localRow), rows(localRow));
      }
      blockW.each_col() %= hDiagonal;
    }

    residualKronBlockWCache[static_cast<std::size_t>(block)] =
      std::move(blockW);
  }
  }

  // ------------------------------------------------------------
  // Persistent residual (gr,gc,lr,lc) index list.
  //
  // residualKronRows/residualLocalOfRow are purely data-driven (fixed for
  // the whole optimisation), so the set of global row/col pairs feeding
  // Rmat and every residualDerivativeBasis(k) never changes across
  // iterations - only the numeric value read at (lr,lc) does.  Building
  // this list once and reusing it every iteration replaces the former
  // per-iteration index bookkeeping plus single-element arma::sp_mat
  // insertion with one batch sparse construction from cached locations.
  // ------------------------------------------------------------
  arma::umat residualPatternLocations;
  std::vector<int> residualPatternLr;
  std::vector<int> residualPatternLc;

  {
    std::vector<arma::uword> patternRows;
    std::vector<arma::uword> patternCols;

    for(int b = 0; b < residualNBlocks; ++b){
      const arma::uvec & rows =
        residualKronRows[static_cast<std::size_t>(b)];

      for(arma::uword aa = 0; aa < rows.n_elem; ++aa){
        const arma::uword gr = rows(aa);
        const int lr = residualLocalOfRow[static_cast<std::size_t>(gr)];

        for(arma::uword bb = 0; bb < rows.n_elem; ++bb){
          const arma::uword gc = rows(bb);
          const int lc = residualLocalOfRow[static_cast<std::size_t>(gc)];

          patternRows.push_back(gr);
          patternCols.push_back(gc);
          residualPatternLr.push_back(lr);
          residualPatternLc.push_back(lc);
        }
      }
    }

    residualPatternLocations.set_size(2, patternRows.size());
    for(std::size_t k = 0; k < patternRows.size(); ++k){
      residualPatternLocations(0,k) = patternRows[k];
      residualPatternLocations(1,k) = patternCols[k];
    }
  }

  auto buildResidualSparseFromPattern =
    [&](const arma::mat & localMat) -> arma::sp_mat {
      const arma::uword nnz = residualPatternLr.size();
      arma::vec values(nnz);
      for(arma::uword k = 0; k < nnz; ++k){
        values(k) = localMat(
          static_cast<arma::uword>(residualPatternLr[k]),
          static_cast<arma::uword>(residualPatternLc[k])
        );
      }
      return arma::sp_mat(residualPatternLocations, values, nR, nR);
    };


  // CHANGE 2 helper:
  // exact Moore-Penrose inverse of
  // diag(theta^2 / denominator), without forming the diagonal
  // matrix and without an SVD-based pinv().
  // ------------------------------------------------------------

  // ------------------------------------------------------------
  // Information-system helper:
  // solve A x = rhs directly; no inverse is formed.
  // ------------------------------------------------------------
  auto solveInformationSystem =
    [&](arma::vec & solution,
        const arma::mat & A,
        const arma::vec & rhs,
        const std::string & context) -> bool {

      bool ok =
        arma::solve(
          solution,
          A,
          rhs,
          arma::solve_opts::likely_sympd
        );

      if(!ok){
        ok =
          arma::solve(
            solution,
            A,
            rhs
          );
      }

      if(!ok){
        Rcpp::Rcout
          << "Information-system solve failed in "
          << context
          << arma::endl;
      }

      return ok;
    };

  auto buildEmInformationDiagonal =
    [&](const arma::vec & thetaVec,
        const double denominator) -> arma::mat {

      arma::vec infoDiag(
        thetaVec.n_elem,
        arma::fill::zeros
      );

      for(arma::uword j = 0; j < thetaVec.n_elem; ++j){

        const double varianceLike =
          (
            thetaVec(j)
            *
            thetaVec(j)
          )
          /
          denominator;

        // Match pinv(A, tolParInv): singular diagonal values at or
        // below the tolerance contribute zero to the pseudoinverse.
        if(std::isfinite(varianceLike) &&
           varianceLike > tolParInv){

          infoDiag(j) =
            1.0
            /
            varianceLike;
        }
      }

      return
        arma::diagmat(
          infoDiag
        );
    };


  // ============================================================
  // Sparse LDLT / Takahashi selected-inverse infrastructure
  // ============================================================
  typedef Eigen::SparseMatrix<double, Eigen::ColMajor, int> EigenSpMat;
  typedef Eigen::Triplet<double, int> EigenTriplet;
  typedef Eigen::SimplicialLDLT<
    EigenSpMat,
    Eigen::Lower,
    SommerSparseOrdering
  > EigenLDLT;

  typedef Eigen::ConjugateGradient<
    EigenSpMat,
    Eigen::Lower | Eigen::Upper,
    Eigen::DiagonalPreconditioner<double>
  > EigenPCG;

  struct SelectedInverseSubset {
    std::vector< std::vector<int> > rows;
    std::vector< std::vector<double> > values;
    std::vector<int> originalToPermuted;
  };

  auto buildSelectedInverseSubset =
    [&](const EigenLDLT & factor,
        const std::string & context,
        SelectedInverseSubset & out,
        const bool reuseTopology) -> void {

      const Eigen::VectorXd D = factor.vectorD();
      const int n = static_cast<int>(D.size());
      EigenSpMat Lmat = factor.matrixL();
      Lmat.makeCompressed();

      // When reuseTopology is set, `out` already holds the sparsity pattern
      // (rows/originalToPermuted) from a previous call with the same LDLT
      // fill-in pattern, so only the numeric values need to be reset. This
      // avoids a full deep copy of the (potentially large) pattern vectors
      // on every REML iteration.
      if(reuseTopology){
        if(static_cast<int>(out.rows.size()) != n ||
           static_cast<int>(out.originalToPermuted.size()) != n){
          Rcpp::stop("Cached selected-inverse topology is incompatible with the LDLT factor.");
        }
        for(int col = 0; col < n; ++col){
          std::fill(out.values[col].begin(), out.values[col].end(), 0.0);
        }
      }else{
        out.rows.assign(n, std::vector<int>());
        out.values.assign(n, std::vector<double>());
        out.originalToPermuted.assign(n, -1);

        const auto & perm = factor.permutationP();
        // Eigen's convention is P e_i = e_{sigma(i)}, therefore
        // perm.indices()(i) maps ORIGINAL index i -> PERMUTED index.
        for(int i = 0; i < n; ++i){
          const int permutedIndex = perm.indices()(i);
          if(permutedIndex < 0 || permutedIndex >= n){
            Rcpp::stop("Invalid permutation encountered while building selected inverse subset.");
          }
          out.originalToPermuted[i] = permutedIndex;
        }

        for(int col = 0; col < n; ++col){
          out.rows[col].push_back(col);
          for(EigenSpMat::InnerIterator it(Lmat, col); it; ++it){
            const int row = it.row();
            if(row > col){ out.rows[col].push_back(row); }
          }
          std::sort(out.rows[col].begin(), out.rows[col].end());
          out.rows[col].erase(std::unique(out.rows[col].begin(), out.rows[col].end()), out.rows[col].end());
          out.values[col].assign(out.rows[col].size(), 0.0);
        }
      }

      auto getPermuted = [&](int a, int b, double & value) -> bool {
        const int col = std::min(a,b);
        const int row = std::max(a,b);
        const std::vector<int> & rr = out.rows[col];
        auto pos = std::lower_bound(rr.begin(), rr.end(), row);
        if(pos == rr.end() || (*pos) != row){ return false; }
        const std::size_t idx = static_cast<std::size_t>(pos - rr.begin());
        value = out.values[col][idx];
        return true;
      };

      auto setPermuted = [&](int row, int col, double value) {
        if(row < col){ std::swap(row,col); }
        const std::vector<int> & rr = out.rows[col];
        auto pos = std::lower_bound(rr.begin(), rr.end(), row);
        if(pos == rr.end() || (*pos) != row){ Rcpp::stop("Internal selected-inverse pattern error."); }
        const std::size_t idx = static_cast<std::size_t>(pos - rr.begin());
        out.values[col][idx] = value;
      };

      // Takahashi recursions for A = P' L D L' P.
      for(int i = n - 1; i >= 0; --i){
        if(!std::isfinite(D(i)) || D(i) <= 0.0){
          Rcpp::stop("Non-positive LDLT pivot while building selected inverse subset in " + context + ".");
        }

        std::vector<int> neighbours;
        std::vector<double> lvalues;
        for(EigenSpMat::InnerIterator it(Lmat, i); it; ++it){
          if(it.row() > i){
            neighbours.push_back(it.row());
            lvalues.push_back(it.value());
          }
        }

        // Off-diagonal entries first.
        for(std::size_t jj = 0; jj < neighbours.size(); ++jj){
          const int j = neighbours[jj];
          double sum = 0.0;
          for(std::size_t kk = 0; kk < neighbours.size(); ++kk){
            double zkj = 0.0;
            if(!getPermuted(neighbours[kk], j, zkj)){
              Rcpp::stop("Takahashi recurrence requested an entry outside the filled LDLT pattern in " + context + ".");
            }
            sum += lvalues[kk] * zkj;
          }
          setPermuted(j, i, -sum);
        }

        double diagCorrection = 0.0;
        for(std::size_t kk = 0; kk < neighbours.size(); ++kk){
          double zki = 0.0;
          if(!getPermuted(neighbours[kk], i, zki)){
            Rcpp::stop("Unable to retrieve selected inverse off-diagonal during Takahashi recursion in " + context + ".");
          }
          diagCorrection += lvalues[kk] * zki;
        }
        setPermuted(i, i, (1.0 / D(i)) - diagCorrection);
      }
    };

  auto getSelectedInverseOriginal =
    [&](const SelectedInverseSubset & subset,
        int originalRow,
        int originalCol,
        double & value) -> bool {
      if(originalRow < 0 || originalCol < 0 ||
         originalRow >= static_cast<int>(subset.originalToPermuted.size()) ||
         originalCol >= static_cast<int>(subset.originalToPermuted.size())){ return false; }
      const int a = subset.originalToPermuted[originalRow];
      const int b = subset.originalToPermuted[originalCol];
      const int col = std::min(a,b);
      const int row = std::max(a,b);
      const std::vector<int> & rr = subset.rows[col];
      auto pos = std::lower_bound(rr.begin(), rr.end(), row);
      if(pos == rr.end() || (*pos) != row){ return false; }
      const std::size_t idx = static_cast<std::size_t>(pos - rr.begin());
      value = subset.values[col][idx];
      return true;
    };



  // ============================================================
  // CHANGE 5: reusable symbolic factorisation for C
  // ============================================================
  // The numerical values of C change with the variance components,
  // but its sparsity pattern is usually unchanged.  Cache the sparse
  // pattern and reuse Eigen's symbolic analysis whenever possible.
  auto eigenSparsePatternMatches =
    [&](const EigenSpMat & A,
        const std::vector<int> & cachedOuter,
        const std::vector<int> & cachedInner) -> bool {

      if(
          cachedOuter.size()
          !=
          static_cast<std::size_t>(A.outerSize() + 1)
          ||
          cachedInner.size()
          !=
          static_cast<std::size_t>(A.nonZeros())
      ){
        return false;
      }

      const int * outer =
        A.outerIndexPtr();

      const int * inner =
        A.innerIndexPtr();

      for(Eigen::Index j = 0; j < A.outerSize() + 1; ++j){
        if(cachedOuter[static_cast<std::size_t>(j)] != outer[j]){
          return false;
        }
      }

      for(Eigen::Index j = 0; j < A.nonZeros(); ++j){
        if(cachedInner[static_cast<std::size_t>(j)] != inner[j]){
          return false;
        }
      }

      return true;
    };

  auto cacheEigenSparsePattern =
    [&](const EigenSpMat & A,
        std::vector<int> & cachedOuter,
        std::vector<int> & cachedInner){

      cachedOuter.assign(
        A.outerIndexPtr(),
        A.outerIndexPtr() + A.outerSize() + 1
      );

      cachedInner.assign(
        A.innerIndexPtr(),
        A.innerIndexPtr() + A.nonZeros()
      );
    };

  // ============================================================
  // MME linear-solver module
  //
  // Cfactor is the exact direct backend and also supplies the LDLT factors
  // needed by log|C| and Takahashi.  Cpcg is an optional iterative backend
  // for the repeated C x = rhs solves (BLUE/BLUP, sensitivities and rare
  // selected-block fallbacks).  Keeping all solve dispatch here makes future
  // solver backends local to this module rather than the AI-REML code.
  // ============================================================
  EigenLDLT Cfactor;
  EigenPCG Cpcg;
  bool CsymbolicReady = false;
  bool CnumericReady = false;
  bool CpcgReady = false;
  std::vector<int> CouterPattern;
  std::vector<int> CinnerPattern;
  SelectedInverseSubset CselectedTopology;
  bool CselectedTopologyReady = false;

  // ML (reml=FALSE) support: D is the random-effects-only block of C,
  // D = Z'R^{-1}Z + G^{-1} (rows/cols nX..nEffects-1 of C, re-indexed to
  // 0..Nu-1). log|D| and D's own Takahashi selected inverse replace
  // log|C| and CselectedTopology in the likelihood/score computations,
  // dropping the log|X'V^{-1}X| restriction term that makes REML
  // "restricted". Only used when reml==false (validated solver=="ldlt"
  // above); never touched otherwise.
  EigenLDLT Dfactor;
  SelectedInverseSubset DselectedTopology;
  bool DselectedTopologyReady = false;


  // ============================================================
  // CHOLMOD (via R's Matrix package) supernodal backend for solver="cholmod"
  // ============================================================
  // CholmodState owns the cholmod_common/cholmod_factor lifetime so they are
  // released even if Rcpp::stop() unwinds out of this function.
  struct CholmodState {
    cholmod_common common;
    cholmod_factor * factor = nullptr;
    bool started = false;
    void ensureStarted(){
      if(!started){
        M_cholmod_start(&common);
        common.supernodal = CHOLMOD_SUPERNODAL;
        started = true;
      }
    }
    ~CholmodState(){
      if(factor != nullptr){ M_cholmod_free_factor(&factor, &common); }
      if(started){ M_cholmod_finish(&common); }
    }
  };
  CholmodState cholmodState;
  bool CholmodSymbolicReady = false;

  // A cholmod_sparse VIEW over an Eigen sparse matrix's existing compressed-
  // storage arrays (zero-copy): CHOLMOD's analyze/factorize only read A, so
  // a const_cast to CHOLMOD's void* fields is safe.
  auto eigenToCholmodSparseView =
    [](const Eigen::SparseMatrix<double> & A) -> cholmod_sparse {
      cholmod_sparse view;
      std::memset(&view, 0, sizeof(view));
      view.nrow = static_cast<size_t>(A.rows());
      view.ncol = static_cast<size_t>(A.cols());
      view.nzmax = static_cast<size_t>(A.nonZeros());
      view.p = const_cast<int *>(A.outerIndexPtr());
      view.i = const_cast<int *>(A.innerIndexPtr());
      view.x = const_cast<double *>(A.valuePtr());
      // Both matrices are symmetric with both triangles populated; stype=-1
      // tells CHOLMOD to read only the lower triangle and assume symmetry
      // for the rest. stype=0 ("general") would instead make CHOLMOD
      // factorise A*A'.
      view.stype = -1;
      view.itype = CHOLMOD_INT;
      view.xtype = CHOLMOD_REAL;
      view.dtype = CHOLMOD_DOUBLE;
      view.sorted = 1;
      view.packed = 1;
      return view;
    };

  auto solveCVectorCholmod =
    [&](const Eigen::Ref<const Eigen::VectorXd> & rhs,
        const std::string & context) -> Eigen::VectorXd {
      cholmod_dense rhsView;
      std::memset(&rhsView, 0, sizeof(rhsView));
      rhsView.nrow = static_cast<size_t>(rhs.size());
      rhsView.ncol = 1;
      rhsView.nzmax = static_cast<size_t>(rhs.size());
      rhsView.d = static_cast<size_t>(rhs.size());
      rhsView.x = const_cast<double *>(rhs.data());
      rhsView.xtype = CHOLMOD_REAL;
      rhsView.dtype = CHOLMOD_DOUBLE;

      cholmod_dense * solution =
        M_cholmod_solve(CHOLMOD_A, cholmodState.factor, &rhsView, &cholmodState.common);
      if(solution == nullptr){
        Rcpp::stop("CHOLMOD solve failed in " + context + ".");
      }
      Eigen::VectorXd ans =
        Eigen::Map<const Eigen::VectorXd>(static_cast<double *>(solution->x), rhs.size());
      M_cholmod_free_dense(&solution, &cholmodState.common);
      return ans;
    };

  auto solveCMatrixCholmod =
    [&](const Eigen::Ref<const Eigen::MatrixXd> & rhs,
        const std::string & context) -> Eigen::MatrixXd {
      cholmod_dense rhsView;
      std::memset(&rhsView, 0, sizeof(rhsView));
      rhsView.nrow = static_cast<size_t>(rhs.rows());
      rhsView.ncol = static_cast<size_t>(rhs.cols());
      rhsView.nzmax = static_cast<size_t>(rhs.size());
      rhsView.d = static_cast<size_t>(rhs.rows());
      // CHOLMOD dense matrices are column-major, matching Eigen's default
      // storage order, so rhs's buffer can be handed over without copying.
      rhsView.x = const_cast<double *>(rhs.data());
      rhsView.xtype = CHOLMOD_REAL;
      rhsView.dtype = CHOLMOD_DOUBLE;

      cholmod_dense * solution =
        M_cholmod_solve(CHOLMOD_A, cholmodState.factor, &rhsView, &cholmodState.common);
      if(solution == nullptr){
        Rcpp::stop("CHOLMOD multi-RHS solve failed in " + context + ".");
      }
      Eigen::MatrixXd ans =
        Eigen::Map<const Eigen::MatrixXd>(
          static_cast<double *>(solution->x), rhs.rows(), rhs.cols()
        );
      M_cholmod_free_dense(&solution, &cholmodState.common);
      return ans;
    };

  // CHOLMOD backend for D (the random-effects-only block used by ML,
  // reml==false), kept in a separate state/factor from C's - D is a
  // different matrix (Nu x Nu, no fixed-effect rows/cols) with its own
  // sparsity pattern, mirroring how cholmodRState is kept separate from
  // cholmodState for R.
  CholmodState cholmodDState;
  bool CholmodDSymbolicReady = false;
  std::vector<int> DouterPattern;
  std::vector<int> DinnerPattern;

  auto solveDMatrixCholmod =
    [&](const Eigen::Ref<const Eigen::MatrixXd> & rhs,
        const std::string & context) -> Eigen::MatrixXd {
      cholmod_dense rhsView;
      std::memset(&rhsView, 0, sizeof(rhsView));
      rhsView.nrow = static_cast<size_t>(rhs.rows());
      rhsView.ncol = static_cast<size_t>(rhs.cols());
      rhsView.nzmax = static_cast<size_t>(rhs.size());
      rhsView.d = static_cast<size_t>(rhs.rows());
      rhsView.x = const_cast<double *>(rhs.data());
      rhsView.xtype = CHOLMOD_REAL;
      rhsView.dtype = CHOLMOD_DOUBLE;

      cholmod_dense * solution =
        M_cholmod_solve(CHOLMOD_A, cholmodDState.factor, &rhsView, &cholmodDState.common);
      if(solution == nullptr){
        Rcpp::stop("CHOLMOD multi-RHS solve for the random-effects-only matrix D failed in " + context + " (reml=FALSE).");
      }
      Eigen::MatrixXd ans =
        Eigen::Map<const Eigen::MatrixXd>(
          static_cast<double *>(solution->x), rhs.rows(), rhs.cols()
        );
      M_cholmod_free_dense(&solution, &cholmodDState.common);
      return ans;
    };

  auto preparePCG = [&](const EigenSpMat & Ce){
    CpcgReady = false;
    if(solverName != "pcg"){ return; }
    Cpcg.setTolerance(pcgTol);
    const int automaticMax = std::max<int>(1000, static_cast<int>(Ce.rows()));
    Cpcg.setMaxIterations(pcgMaxIters > 0 ? pcgMaxIters : automaticMax);
    Cpcg.compute(Ce);
    if(Cpcg.info() != Eigen::Success){
      Rcpp::stop("PCG setup failed for the MME coefficient matrix C.");
    }
    CpcgReady = true;
  };

  auto solveCVector = [&](const Eigen::Ref<const Eigen::VectorXd> & rhs,
                          const std::string & context) -> Eigen::VectorXd {
    if(solverName == "ldlt"){
      Eigen::VectorXd ans = Cfactor.solve(rhs);
      if(Cfactor.info() != Eigen::Success){
        Rcpp::stop("Sparse LDLT solve failed in " + context + ".");
      }
      return ans;
    }
    if(solverName == "cholmod"){
      return solveCVectorCholmod(rhs, context);
    }
    if(!CpcgReady){
      Rcpp::stop("PCG solve requested before the PCG backend was prepared.");
    }
    Eigen::VectorXd ans = Cpcg.solve(rhs);
    if(Cpcg.info() != Eigen::Success || !ans.allFinite()){
      Rcpp::stop("PCG failed to converge in " + context + ".");
    }
    return ans;
  };

  auto solveCMatrix = [&](const Eigen::Ref<const Eigen::MatrixXd> & rhs,
                          const std::string & context) -> Eigen::MatrixXd {
    Eigen::MatrixXd ans(rhs.rows(), rhs.cols());
    if(solverName == "ldlt"){
      ans = Cfactor.solve(rhs);
      if(Cfactor.info() != Eigen::Success){
        Rcpp::stop("Sparse LDLT multi-RHS solve failed in " + context + ".");
      }
      return ans;
    }
    if(solverName == "cholmod"){
      return solveCMatrixCholmod(rhs, context);
    }
    if(!CpcgReady){
      Rcpp::stop("PCG solve requested before the PCG backend was prepared.");
    }
    // Eigen's iterative solver is run column-by-column so convergence is
    // checked independently for every derivative/right-hand side.
    for(Eigen::Index j = 0; j < rhs.cols(); ++j){
      ans.col(j) = Cpcg.solve(rhs.col(j));
      if(Cpcg.info() != Eigen::Success || !ans.col(j).allFinite()){
        Rcpp::stop("PCG failed to converge for RHS column " +
                   std::to_string(static_cast<long long>(j + 1)) +
                   " in " + context + ".");
      }
    }
    return ans;
  };

  // Computes tr(B * factor^{-1}) using the selected-inverse fast path when
  // available (LDLT only). genericSolve, when non-empty, routes the
  // fallback through it instead of factor.solve() directly - required for
  // a solver="cholmod" factor (its supernodal representation has no
  // Takahashi subset, so `factor`/`subset` are never populated for it).
  auto sparseTraceInverseTimes =
    [&](const arma::sp_mat & B,
        const EigenLDLT & factor,
        const SelectedInverseSubset & subset,
        const std::string & context,
        bool & usedFallback,
        const std::function<Eigen::MatrixXd(const Eigen::Ref<const Eigen::MatrixXd> &, const std::string &)> & genericSolve = nullptr) -> double {
      usedFallback = false;
      double traceValue = 0.0;
      bool allAvailable = true;
      for(arma::sp_mat::const_iterator it = B.begin(); it != B.end(); ++it){
        double zij = 0.0;
        if(!getSelectedInverseOriginal(subset, static_cast<int>(it.col()), static_cast<int>(it.row()), zij)){
          allAvailable = false;
          break;
        }
        traceValue += (*it) * zij;
      }
      if(allAvailable){ return traceValue; }

      usedFallback = true;
      traceValue = 0.0;
      std::vector<arma::uword> activeColumns;
      activeColumns.reserve(static_cast<std::size_t>(B.n_cols));
      for(arma::uword j = 0; j < B.n_cols; ++j){
        if(B.begin_col(j) != B.end_col(j)){ activeColumns.push_back(j); }
      }
      const std::size_t batchSize = 32;
      for(std::size_t batchStart = 0; batchStart < activeColumns.size(); batchStart += batchSize){
        const std::size_t batchEnd = std::min(batchStart + batchSize, activeColumns.size());
        const std::size_t currentBatchSize = batchEnd - batchStart;
        Eigen::MatrixXd rhs = Eigen::MatrixXd::Zero(static_cast<Eigen::Index>(B.n_rows), static_cast<Eigen::Index>(currentBatchSize));
        for(std::size_t k = 0; k < currentBatchSize; ++k){
          const arma::uword sourceCol = activeColumns[batchStart + k];
          for(arma::sp_mat::const_col_iterator it = B.begin_col(sourceCol); it != B.end_col(sourceCol); ++it){
            rhs(static_cast<Eigen::Index>(it.row()), static_cast<Eigen::Index>(k)) = (*it);
          }
        }
        Eigen::MatrixXd solution;
        if(genericSolve){
          solution = genericSolve(rhs, context);
        }else{
          solution = factor.solve(rhs);
          if(factor.info() != Eigen::Success){ Rcpp::stop("Sparse LDLT trace fallback solve failed in " + context + "."); }
        }
        for(std::size_t k = 0; k < currentBatchSize; ++k){
          const arma::uword sourceCol = activeColumns[batchStart + k];
          traceValue += solution(static_cast<Eigen::Index>(sourceCol), static_cast<Eigen::Index>(k));
        }
      }
      return traceValue;
    };

  // ============================================================
  // Factorisation-free PCG trace/log-determinant helpers for C
  // ============================================================
  // Deterministic Rademacher signs give common random numbers across REML
  // iterations.  This is important for the likelihood line search: changes
  // in the approximate log determinant then reflect C, not fresh Monte Carlo
  // noise at every iteration.
  auto pcgProbeSign = [](const unsigned long long row,
                         const unsigned long long probe) -> double {
    unsigned long long x = row + 0x9e3779b97f4a7c15ULL * (probe + 1ULL);
    x ^= x >> 30;
    x *= 0xbf58476d1ce4e5b9ULL;
    x ^= x >> 27;
    x *= 0x94d049bb133111ebULL;
    x ^= x >> 31;
    return (x & 1ULL) ? 1.0 : -1.0;
  };

  Eigen::MatrixXd pcgLanczosInverseGuess;

  auto pcgApproxLogDet = [&](const EigenSpMat & A) -> double {
    const Eigen::Index n = A.rows();
    if(n <= 0){ return 0.0; }
    const int mMax = std::min<int>(pcgLanczosSteps, static_cast<int>(n));
    std::vector<double> probeEstimates(
      static_cast<std::size_t>(pcgTraceProbes),
      0.0
    );
    pcgLanczosInverseGuess.resize(n, pcgTraceProbes);
    std::atomic<bool> invalidProbe(false);

    // Each deterministic Rademacher probe is independent.  Keeping one
    // accumulator per probe makes the final reduction reproducible across
    // OpenMP schedules.
    #pragma omp parallel for if(pcgTraceProbes > 1)
    for(int probe = 0; probe < pcgTraceProbes; ++probe){
      Eigen::VectorXd q(n), qPrev = Eigen::VectorXd::Zero(n);
      for(Eigen::Index i = 0; i < n; ++i){
        q(i) = pcgProbeSign(static_cast<unsigned long long>(i),
                           static_cast<unsigned long long>(probe));
      }
      const double normz = q.norm();
      q /= normz;
      Eigen::MatrixXd lanczosBasis(n, mMax);

      std::vector<double> alpha;
      std::vector<double> beta;
      alpha.reserve(static_cast<std::size_t>(mMax));
      beta.reserve(static_cast<std::size_t>(std::max(0, mMax-1)));
      double betaPrev = 0.0;

      for(int j = 0; j < mMax; ++j){
        lanczosBasis.col(j) = q;
        Eigen::VectorXd w = A * q;
        if(j > 0){ w.noalias() -= betaPrev * qPrev; }
        const double a = q.dot(w);
        w.noalias() -= a * q;
        // A small second local orthogonalisation against qPrev reduces
        // loss of orthogonality without storing the complete Lanczos basis.
        if(j > 0){
          const double corr = qPrev.dot(w);
          w.noalias() -= corr * qPrev;
        }
        const double b = w.norm();
        alpha.push_back(a);
        if(j + 1 >= mMax || b <= 1.0e-14 * std::max(1.0, std::abs(a))){ break; }
        beta.push_back(b);
        qPrev.swap(q);
        q = w / b;
        betaPrev = b;
      }

      const int m = static_cast<int>(alpha.size());
      Eigen::MatrixXd T = Eigen::MatrixXd::Zero(m,m);
      for(int j = 0; j < m; ++j){
        T(j,j) = alpha[static_cast<std::size_t>(j)];
        if(j + 1 < m){
          const double b = beta[static_cast<std::size_t>(j)];
          T(j,j+1) = b;
          T(j+1,j) = b;
        }
      }
      Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(T);
      if(es.info() != Eigen::Success){
        invalidProbe.store(true, std::memory_order_relaxed);
        continue;
      }
      const Eigen::VectorXd eval = es.eigenvalues();
      const Eigen::MatrixXd evec = es.eigenvectors();
      double quad = 0.0;
      Eigen::VectorXd inverseCoefficients = Eigen::VectorXd::Zero(m);
      for(int j = 0; j < m; ++j){
        if(!std::isfinite(eval(j)) || eval(j) <= 0.0){
          invalidProbe.store(true, std::memory_order_relaxed);
          quad = 0.0;
          break;
        }
        const double w0 = evec(0,j);
        quad += w0*w0*std::log(eval(j));
        inverseCoefficients.noalias() +=
          (normz * w0 / eval(j)) * evec.col(j);
      }
      probeEstimates[static_cast<std::size_t>(probe)] = normz*normz*quad;
      pcgLanczosInverseGuess.col(probe).noalias() =
        lanczosBasis.leftCols(m) * inverseCoefficients;
    }

    if(invalidProbe.load(std::memory_order_relaxed)){
      Rcpp::stop("PCG/SLQ failed while estimating log|C|; C may not be positive definite or Lanczos accuracy is insufficient.");
    }

    double total = 0.0;
    for(int probe = 0; probe < pcgTraceProbes; ++probe){
      total += probeEstimates[static_cast<std::size_t>(probe)];
    }
    return total / static_cast<double>(pcgTraceProbes);
  };

  // Z and X=C^{-1}Z are allocated only in PCG mode and reused by every
  // random/residual trace in the current REML iteration.
  Eigen::MatrixXd pcgTraceZ;
  Eigen::MatrixXd pcgTraceX;
  bool reportedOpenMpSlq = false;
  bool reportedOpenMpResidualBlocks = false;

  auto preparePCGTraceProbes = [&](const EigenSpMat & A){
    if(solverName != "pcg"){ return; }
    const Eigen::Index n = A.rows();
    pcgTraceZ.resize(n, pcgTraceProbes);
    pcgTraceX.resize(n, pcgTraceProbes);
    for(int p = 0; p < pcgTraceProbes; ++p){
      for(Eigen::Index i = 0; i < n; ++i){
        pcgTraceZ(i,p) = pcgProbeSign(static_cast<unsigned long long>(i),
                                     static_cast<unsigned long long>(p));
      }
      if(
          pcgLanczosInverseGuess.rows() == n
          &&
          pcgLanczosInverseGuess.cols() == pcgTraceProbes
      ){
        pcgTraceX.col(p) = Cpcg.solveWithGuess(
          pcgTraceZ.col(p),
          pcgLanczosInverseGuess.col(p)
        );
      }else{
        pcgTraceX.col(p) = Cpcg.solve(pcgTraceZ.col(p));
      }
      if(Cpcg.info() != Eigen::Success || !pcgTraceX.col(p).allFinite()){
        Rcpp::stop("PCG failed while preparing Hutchinson trace probes for C^{-1}.");
      }
    }
  };


    if(verbose && !reportedOpenMpSlq){
#ifdef _OPENMP
      Rcpp::Rcout
        << "OpenMP active: parallel SLQ log-determinant probes ("
        << pcgTraceProbes
        << " probes)."
        << arma::endl;
#else
      Rcpp::Rcout
        << "OpenMP unavailable: SLQ log-determinant probes run serially."
        << arma::endl;
#endif
      reportedOpenMpSlq = true;
    }
  auto pcgTraceCInverseTimesSparse = [&](const arma::sp_mat & B) -> double {
    if(solverName != "pcg" || pcgTraceX.cols() != pcgTraceProbes){
      Rcpp::stop("PCG trace probes are not available.");
    }
    double out = 0.0;
    for(int p = 0; p < pcgTraceProbes; ++p){
      double one = 0.0;
      for(arma::sp_mat::const_iterator it = B.begin(); it != B.end(); ++it){
        one += pcgTraceZ(static_cast<Eigen::Index>(it.row()), p)
             * (*it)
             * pcgTraceX(static_cast<Eigen::Index>(it.col()), p);
      }
      out += one;
    }
    return out / static_cast<double>(pcgTraceProbes);
  };


  // ============================================================
  // CHANGE 7: reusable factorisation of the small Kronecker block R0
  // ============================================================
  EigenLDLT RkronFactor;
  bool RkronSymbolicReady = false;
  std::vector<int> RkronOuterPattern;
  std::vector<int> RkronInnerPattern;

  auto factorizeRkronWithCachedPattern =
    [&](const EigenSpMat & Re) -> bool {

      const bool samePattern =
        RkronSymbolicReady
        &&
        eigenSparsePatternMatches(
          Re,
          RkronOuterPattern,
          RkronInnerPattern
        );

      if(!samePattern){

        RkronFactor.analyzePattern(
          Re
        );

        if(RkronFactor.info() != Eigen::Success){
          return false;
        }

        cacheEigenSparsePattern(
          Re,
          RkronOuterPattern,
          RkronInnerPattern
        );

        RkronSymbolicReady =
          true;
      }

      RkronFactor.factorize(
        Re
      );

      if(
          RkronFactor.info() != Eigen::Success
          &&
          samePattern
      ){

        RkronFactor.analyzePattern(
          Re
        );

        if(RkronFactor.info() == Eigen::Success){
          RkronFactor.factorize(
            Re
          );
        }
      }

      return
        RkronFactor.info()
        ==
        Eigen::Success;
    };

  // ============================================================
  // CHANGE 6: reusable symbolic factorisation for generic non-diagonal R
  // ============================================================
  EigenLDLT Rfactor;
  bool RsymbolicReady = false;
  std::vector<int> RouterPattern;
  std::vector<int> RinnerPattern;

  // CHOLMOD backend for R, used when solverName=="cholmod" (Phase 2). R is
  // usually far smaller than C (this path is only reached for irregular
  // residual designs that are neither diagonal nor repeated-block), but the
  // supernodal factor is reused with the same symbolic-pattern-caching
  // discipline as C for consistency.
  CholmodState cholmodRState;
  bool CholmodRSymbolicReady = false;

  auto solveRMatrixCholmod =
    [&](const Eigen::Ref<const Eigen::MatrixXd> & rhs,
        const std::string & context) -> Eigen::MatrixXd {
      cholmod_dense rhsView;
      std::memset(&rhsView, 0, sizeof(rhsView));
      rhsView.nrow = static_cast<size_t>(rhs.rows());
      rhsView.ncol = static_cast<size_t>(rhs.cols());
      rhsView.nzmax = static_cast<size_t>(rhs.size());
      rhsView.d = static_cast<size_t>(rhs.rows());
      rhsView.x = const_cast<double *>(rhs.data());
      rhsView.xtype = CHOLMOD_REAL;
      rhsView.dtype = CHOLMOD_DOUBLE;

      cholmod_dense * solution =
        M_cholmod_solve(CHOLMOD_A, cholmodRState.factor, &rhsView, &cholmodRState.common);
      if(solution == nullptr){
        Rcpp::stop("CHOLMOD solve for the residual covariance matrix R failed in " + context + ".");
      }
      Eigen::MatrixXd ans =
        Eigen::Map<const Eigen::MatrixXd>(
          static_cast<double *>(solution->x), rhs.rows(), rhs.cols()
        );
      M_cholmod_free_dense(&solution, &cholmodRState.common);
      return ans;
    };

  auto factorizeRWithCachedPattern =
    [&](const EigenSpMat & Re) -> bool {

      if(solverName == "cholmod"){
        cholmodRState.ensureStarted();
        cholmod_sparse Rview = eigenToCholmodSparseView(Re);

        const bool samePattern =
          CholmodRSymbolicReady
          &&
          eigenSparsePatternMatches(Re, RouterPattern, RinnerPattern);

        if(!samePattern || cholmodRState.factor == nullptr){
          if(cholmodRState.factor != nullptr){
            M_cholmod_free_factor(&cholmodRState.factor, &cholmodRState.common);
          }
          cholmodRState.factor = M_cholmod_analyze(&Rview, &cholmodRState.common);
          if(cholmodRState.factor == nullptr || cholmodRState.common.status != CHOLMOD_OK){
            return false;
          }
          cacheEigenSparsePattern(Re, RouterPattern, RinnerPattern);
          CholmodRSymbolicReady = true;
        }

        const int factorizeOk =
          M_cholmod_factorize(&Rview, cholmodRState.factor, &cholmodRState.common);
        return factorizeOk != 0 && cholmodRState.common.status == CHOLMOD_OK;
      }

      const bool sameRPattern =
        RsymbolicReady
        &&
        eigenSparsePatternMatches(
          Re,
          RouterPattern,
          RinnerPattern
        );

      if(!sameRPattern){

        Rfactor.analyzePattern(
          Re
        );

        if(Rfactor.info() != Eigen::Success){
          return false;
        }

        cacheEigenSparsePattern(
          Re,
          RouterPattern,
          RinnerPattern
        );

        RsymbolicReady =
          true;
      }

      Rfactor.factorize(
        Re
      );

      if(Rfactor.info() != Eigen::Success && sameRPattern){

        // Defensive retry with a fresh symbolic analysis.
        Rfactor.analyzePattern(
          Re
        );

        if(Rfactor.info() == Eigen::Success){
          Rfactor.factorize(
            Re
          );
        }
      }

      return
        Rfactor.info()
        ==
        Eigen::Success;
    };


  // ============================================================
  // GLOBAL OPTIMIZER SAFEGUARDS FOR TRANSFORMED COVARIANCE COORDINATES
  // ============================================================
  //
  // The generic covariance descriptors use log / tanh / Cholesky working
  // coordinates.  Those mappings guarantee valid covariance structures over
  // very large regions of parameter space, so "is Sigma positive definite?"
  // is no longer an adequate proxy for "is this AI step reasonable?".
  //
  // We therefore use three complementary safeguards:
  //
  //   1) likelihood convergence uses |Delta logLik|, never a signed decrease;
  //   2) every proposed global covariance step is checked at the next pass
  //      through the exact REML likelihood and backtracked geometrically when
  //      it decreases the likelihood;
  //   3) transformed-coordinate trust caps limit one-iteration movement before
  //      the likelihood line search is even attempted.
  //
  // Backtracking is implemented without duplicating the expensive likelihood
  // code: a rejected trial is replaced by a half-step, iIter is decremented,
  // and the same iteration slot is recomputed.  Consequently retries do not
  // consume the user's requested nIters.
  // ============================================================

  bool haveAcceptedLikelihood = false;
  double acceptedLikelihood =
    -std::numeric_limits<double>::infinity();

  bool pendingLineSearch = false;
  arma::vec lineSearchBase;
  arma::vec lineSearchTarget;
  double lineSearchAlpha = 1.0;
  int lineSearchHalvings = 0;

  const int maxLineSearchHalvings = 24;

  // Permit only a negligible numerical drop when deciding whether a trial
  // likelihood is acceptable.
  const double likelihoodAcceptTolerance =
    std::max(
      1.0e-10,
      0.1 * tolParConvLL
    );

  bool lineSearchStalled = false;

  // ------------------------------------------------------------
  // Per-parameter trust caps in WORKING coordinates.
  // ------------------------------------------------------------
  // The covariance model owns these caps through its CovarianceFactor
  // descriptor.  The optimizer no longer contains model-specific branches.
  arma::vec workingStepCap(
    nVcTotal,
    arma::fill::value(2.0)
  );

  {
    arma::uword globalOffset = 0;

    for(int iStruct = 0; iStruct < nRRe; ++iStruct){
      const arma::uword nLocal = covPar(iStruct).n_elem;
      if(nLocal == 0){
        continue;
      }

      // Product-level log(sigma^2).
      workingStepCap(globalOffset) = 1.0;

      Rcpp::List cs =
        covDescriptor[static_cast<std::size_t>(iStruct)];
      Rcpp::List factors = cs["factors"];

      for(int fidx = 0; fidx < factors.size(); ++fidx){
        Rcpp::List f = Rcpp::as<Rcpp::List>(factors[fidx]);
        const int start1 = Rcpp::as<int>(f["par_start"]);
        const int end1 = Rcpp::as<int>(f["par_end"]);
        if(end1 < start1){
          continue;
        }

        arma::vec caps = Rcpp::as<arma::vec>(f["trust_cap"]);
        const int nFactorPar = end1 - start1 + 1;
        if(caps.n_elem != static_cast<arma::uword>(nFactorPar)){
          Rcpp::stop("CovarianceFactor trust_cap has incompatible length.");
        }

        for(int local = 0; local < nFactorPar; ++local){
          const double cap = caps(static_cast<arma::uword>(local));
          if(!std::isfinite(cap) || cap <= 0.0){
            Rcpp::stop("CovarianceFactor trust caps must be positive and finite.");
          }
          workingStepCap(
            globalOffset + static_cast<arma::uword>(start1 - 1 + local)
          ) = cap;
        }
      }

      globalOffset += nLocal;
    }
  }

  ////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////
  // START ITERATIVE ALGORITHM
  ////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////
  
  
    
for (int iIter = 0; iIter < nIters; ++iIter) {

    // ###########################
    // # 1) absorption of M into y to obtain y'Py and logDetC
    // ###########################

    const int residualStruct =
      nRRe - 1;

    const arma::uword nResidualPar =
      covPar(residualStruct).n_elem;

    const arma::mat residualSigma =
      theta(residualStruct);

    for(int iStruct = 0; iStruct < nRRe; ++iStruct){
      const std::size_t structureOffset =
        static_cast<std::size_t>(iStruct);
      const std::size_t parameterCount =
        static_cast<std::size_t>(covPar(iStruct).n_elem);
      covarianceDerivativeCache[structureOffset].resize(parameterCount);
      covarianceDerivativeReady[structureOffset].assign(parameterCount, false);
    }

    if(
        residualSigma.n_rows != static_cast<arma::uword>(residualDim)
        ||
        residualSigma.n_cols != static_cast<arma::uword>(residualDim)
    ){
      Rcpp::stop("Residual covariance dimension changed unexpectedly.");
    }

    arma::field<arma::mat> residualLocalD1(
      nResidualPar
    );

    arma::field<arma::sp_mat> residualDerivativeBasis(
      nResidualPar
    );

    for(arma::uword k = 0; k < nResidualPar; ++k){

      residualLocalD1(k) =
        cachedCovarianceD1(
          residualStruct,
          k
        );

      residualDerivativeBasis(k) =
        buildResidualSparseFromPattern(residualLocalD1(k));
    }

    // ============================================================
    // STRUCTURED R HANDLING
    // ============================================================
    // R is assembled sparsely from the residual covariance bases.
    // If R (and the optional weight transform) are diagonal, all
    // R^{-1} applications are elementwise.  Otherwise R is kept
    // sparse and factorised by Eigen::SimplicialLDLT; no dense
    // inv_sympd(R) is formed.
    // ============================================================

    arma::sp_mat Rmat;

    const bool Rdiag =
      residualStructurallyDiagonal;

    const bool Rkron =
      (
        !Rdiag
        &&
        residualKronBlocks
      );

    arma::vec RdiagInv;
    arma::vec effectiveRiDiag;
    double logDetR = 0.0;
    SelectedInverseSubset Rselected;
    arma::mat RkronInv;

    if(Rdiag){

      arma::vec rd(
        nR,
        arma::fill::zeros
      );

      for(int rr = 0; rr < nR; ++rr){
        const int local =
          residualLocalOfRow[
            static_cast<std::size_t>(rr)
          ];

        rd(static_cast<arma::uword>(rr)) =
          residualSigma(
            static_cast<arma::uword>(local),
            static_cast<arma::uword>(local)
          );
      }

      const double minRd = rd.min();

      if(!rd.is_finite() || minRd <= 0.0){
        Rcpp::stop("Diagonal residual covariance produced a non-positive/non-finite value.");
      }

      RdiagInv = 1.0 / rd;
      logDetR = arma::accu(arma::log(rd));

      if(!useH || Hdiag){
        effectiveRiDiag = RdiagInv;
        if(useH){
          effectiveRiDiag %= HdiagSquared;
        }
      }

    }else if(Rkron){

      arma::sp_mat R0 =
        arma::sp_mat(
          residualSigma
        );

      EigenSpMat R0e(
        residualKronBlockSize,
        residualKronBlockSize
      );

      std::vector<EigenTriplet> tripsR0;
      tripsR0.reserve(
        static_cast<std::size_t>(
          R0.n_nonzero
        )
      );

      for(arma::sp_mat::const_iterator it = R0.begin();
          it != R0.end();
          ++it){

        tripsR0.emplace_back(
          static_cast<int>(it.row()),
          static_cast<int>(it.col()),
          (*it)
        );
      }

      R0e.setFromTriplets(
        tripsR0.begin(),
        tripsR0.end()
      );

      R0e.makeCompressed();

      bool factorOK =
        factorizeRkronWithCachedPattern(
          R0e
        );

      if(!factorOK){

        bool recovered = false;

        for(int k = 0; k < 3 && !recovered; ++k){

          const double jitter =
            tolParInv
            *
            std::pow(
              10.0,
              static_cast<double>(k)
            );

          arma::sp_mat R0try =
            R0
            +
            arma::speye<arma::sp_mat>(
              residualKronBlockSize,
              residualKronBlockSize
            )
            *
            jitter;

          std::vector<EigenTriplet> retryTrips;
          retryTrips.reserve(
            static_cast<std::size_t>(
              R0try.n_nonzero
            )
          );

          EigenSpMat R0tryE(
            residualKronBlockSize,
            residualKronBlockSize
          );

          for(arma::sp_mat::const_iterator it = R0try.begin();
              it != R0try.end();
              ++it){

            retryTrips.emplace_back(
              static_cast<int>(it.row()),
              static_cast<int>(it.col()),
              (*it)
            );
          }

          R0tryE.setFromTriplets(
            retryTrips.begin(),
            retryTrips.end()
          );

          R0tryE.makeCompressed();

          if(
              factorizeRkronWithCachedPattern(
                R0tryE
              )
          ){
            R0 = R0try;
            R0e = R0tryE;
            recovered = true;
          }
        }

        if(!recovered){
          Rcpp::stop(
            "Kronecker residual block factorisation failed."
          );
        }
      }

      const Eigen::VectorXd Dr0 =
        RkronFactor.vectorD();

      double logDetR0 = 0.0;

      for(Eigen::Index j = 0;
          j < Dr0.size();
          ++j){

        const double dj = Dr0(j);

        if(!std::isfinite(dj) || dj <= 0.0){
          Rcpp::stop(
            "Kronecker residual block LDLT produced a non-positive/non-finite pivot."
          );
        }

        logDetR0 += std::log(dj);
      }

      logDetR =
        static_cast<double>(
          residualKronNBlocks
        )
        *
        logDetR0;

      Eigen::MatrixXd Iq =
        Eigen::MatrixXd::Identity(
          residualKronBlockSize,
          residualKronBlockSize
        );

      Eigen::MatrixXd R0invEig =
        RkronFactor.solve(
          Iq
        );

      if(RkronFactor.info() != Eigen::Success){
        Rcpp::stop(
          "Kronecker residual block inverse solve failed."
        );
      }

      RkronInv.set_size(
        static_cast<arma::uword>(
          residualKronBlockSize
        ),
        static_cast<arma::uword>(
          residualKronBlockSize
        )
      );

      std::copy(
        R0invEig.data(),
        R0invEig.data() + R0invEig.size(),
        RkronInv.memptr()
      );

      RkronInv =
        0.5
        *
        (
          RkronInv
          +
          RkronInv.t()
        );

      if(verbose && iIter == 0){
        Rcpp::Rcout
          << "Using exact descriptor-driven repeated-block/Kronecker residual engine: "
          << residualKronNBlocks
          << " blocks x "
          << residualKronBlockSize
          << " covariance coordinates"
          << arma::endl;
      }

    }else{

      // Generic sparse residual path for incomplete/unbalanced blocks.
      Rmat = buildResidualSparseFromPattern(residualSigma);

      if(Rmat.n_rows > static_cast<arma::uword>(std::numeric_limits<int>::max())){
        Rcpp::stop("R is too large for the 32-bit Eigen sparse index type used in ai_mme_sp2().");
      }

      auto armaSparseToEigen =
        [&](const arma::sp_mat & A) -> EigenSpMat {

          EigenSpMat out(
            static_cast<int>(A.n_rows),
            static_cast<int>(A.n_cols)
          );

          std::vector<EigenTriplet> trips;
          trips.reserve(
            static_cast<std::size_t>(
              A.n_nonzero
            )
          );

          for(arma::sp_mat::const_iterator it = A.begin();
              it != A.end();
              ++it){

            trips.emplace_back(
              static_cast<int>(it.row()),
              static_cast<int>(it.col()),
              *it
            );
          }

          out.setFromTriplets(
            trips.begin(),
            trips.end()
          );

          out.makeCompressed();
          return out;
        };

      EigenSpMat Re =
        armaSparseToEigen(
          Rmat
        );

      bool RfactorOK =
        factorizeRWithCachedPattern(
          Re
        );

      if(!RfactorOK){

        bool recovered = false;

        for(int k = 0; k < 3 && !recovered; ++k){

          const double jitter =
            tolParInv
            *
            std::pow(
              10.0,
              static_cast<double>(k)
            );

          arma::sp_mat Rtry =
            Rmat
            +
            arma::speye<arma::sp_mat>(
              nR,
              nR
            )
            *
            jitter;

          Re =
            armaSparseToEigen(
              Rtry
            );

          if(
              factorizeRWithCachedPattern(
                Re
              )
          ){
            Rmat = Rtry;
            recovered = true;
          }
        }

        if(!recovered){
          Rcpp::stop(
            "Sparse LDLT factorisation of the residual covariance matrix R failed."
          );
        }
      }

      if(solverName == "cholmod"){
        logDetR += M_cholmod_factor_ldetA(cholmodRState.factor);
        if(!std::isfinite(logDetR)){
          Rcpp::stop("CHOLMOD produced a non-finite log-determinant for R.");
        }
        // No Takahashi subset for a supernodal factor: Rselected stays
        // empty, and sparseTraceInverseTimes()'s fallback is routed through
        // solveRMatrixCholmod (see the R-trace call site below).
      }else{
        const Eigen::VectorXd Dr =
          Rfactor.vectorD();

        for(Eigen::Index j = 0; j < Dr.size(); ++j){

          if(!std::isfinite(Dr(j)) || Dr(j) <= 0.0){
            Rcpp::stop(
              "Residual sparse LDLT produced a non-positive/non-finite pivot."
            );
          }

          logDetR += std::log(Dr(j));
        }

        buildSelectedInverseSubset(
          Rfactor,
          "R",
          Rselected,
          false
        );
      }
    }

    // Apply the effective residual precision Hs * R^{-1} * Hs'.
    auto applyRiDense = [&](const arma::mat & B) -> arma::mat {
      if(B.n_rows != static_cast<arma::uword>(nR)){
        Rcpp::stop("R^{-1} application received an incompatible RHS.");
      }

      if(effectiveRiDiag.n_elem == static_cast<arma::uword>(nR)){
        arma::mat out = B;
        out.each_col() %= effectiveRiDiag;
        return out;
      }

      arma::mat rhs = B;
      if(useH){ rhs = arma::mat(Hs.t() * rhs); }

      arma::mat solved;
      if(Rdiag){
        solved = rhs;
        solved.each_col() %= RdiagInv;
      }else if(Rkron){

        solved.zeros(
          rhs.n_rows,
          rhs.n_cols
        );

        if(verbose && !reportedOpenMpResidualBlocks){
#ifdef _OPENMP
          Rcpp::Rcout
            << "OpenMP active: parallel repeated residual-block precision applications ("
            << residualKronNBlocks
            << " blocks)."
            << arma::endl;
#else
          Rcpp::Rcout
            << "OpenMP unavailable: repeated residual-block precision applications run serially."
            << arma::endl;
#endif
          reportedOpenMpResidualBlocks = true;
        }

        // Apply the same small R0^{-1} to every independent block.
        #pragma omp parallel for if(residualKronNBlocks > 1)
        for(int b = 0;
            b < residualKronNBlocks;
            ++b){

          const arma::uvec & rows =
            residualKronRows[
              static_cast<std::size_t>(b)
            ];

          const arma::mat solvedBlock =
            RkronInv * rhs.rows(rows);

          for(arma::uword localRow = 0;
              localRow < rows.n_elem;
              ++localRow){
            solved.row(rows(localRow)) = solvedBlock.row(localRow);
          }
        }

      }else{
        Eigen::Map<const Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor> > rhsEig(
          rhs.memptr(), static_cast<Eigen::Index>(rhs.n_rows), static_cast<Eigen::Index>(rhs.n_cols)
        );
        Eigen::MatrixXd solEig;
        if(solverName == "cholmod"){
          solEig = solveRMatrixCholmod(rhsEig, "R^{-1} application");
        }else{
          solEig = Rfactor.solve(rhsEig);
          if(Rfactor.info() != Eigen::Success){
            Rcpp::stop("Sparse residual LDLT solve failed.");
          }
        }
        solved.set_size(rhs.n_rows, rhs.n_cols);
        std::copy(solEig.data(), solEig.data()+solEig.size(), solved.memptr());
      }

      if(useH){ solved = arma::mat(Hs * solved); }
      return solved;
    };

    auto applyRiVec = [&](const arma::vec & v) -> arma::vec {
      arma::mat tmp(v);
      arma::mat ans = applyRiDense(tmp);
      return ans.col(0);
    };

    // ------------------------------------------------------------
    // Build the data part of C and the MME RHS without an explicit Ri.
    // ------------------------------------------------------------
    arma::vec Riy = applyRiVec(yDense);
    double yRiy = arma::dot(yDense, Riy);
    arma::vec rhsMME;

    Eigen::SparseMatrix<double> RiWsp;
    arma::mat RiWdense;
    bool RiWisSparse = (effectiveRiDiag.n_elem == static_cast<arma::uword>(nR));

    if(RiWisSparse){
      const Eigen::Map<const Eigen::VectorXd> effectiveRiDiagEig(
        effectiveRiDiag.memptr(), static_cast<Eigen::Index>(effectiveRiDiag.n_elem)
      );
      RiWsp = effectiveRiDiagEig.asDiagonal() * W;
      C = Eigen::SparseMatrix<double>(W.transpose()) * RiWsp;
      C.makeCompressed();
      const Eigen::Map<const Eigen::VectorXd> RiyEig(
        Riy.memptr(), static_cast<Eigen::Index>(Riy.n_elem)
      );
      const Eigen::VectorXd rhsMMEEig = W.transpose() * RiyEig;
      rhsMME.set_size(rhsMMEEig.size());
      std::copy(rhsMMEEig.data(), rhsMMEEig.data() + rhsMMEEig.size(), rhsMME.memptr());
    }else if(Rkron && (!useH || Hdiag)){
      // A repeated residual block does not couple observations from other
      // blocks.  Restricting each operation to the MME columns active in a
      // block avoids the global nR x nEffects dense RiW temporary.
      C.resize(nEffects, nEffects);
      RiWsp.resize(nR, nEffects);
      std::vector<Eigen::Triplet<double>> Ctriplets;
      std::vector<Eigen::Triplet<double>> RiWspTriplets;

      for(int block = 0; block < residualKronNBlocks; ++block){
        const arma::uvec & rows =
          residualKronRows[static_cast<std::size_t>(block)];
        const arma::uvec & columns =
          residualActiveColumns[static_cast<std::size_t>(block)];

        if(columns.n_elem == 0){
          continue;
        }

        // blockW (including any diagonal-H weighting) is constant across
        // REML iterations and was precomputed once before the main loop.
        const arma::mat & blockW =
          residualKronBlockWCache[static_cast<std::size_t>(block)];

        arma::mat blockRiW = RkronInv * blockW;

        if(useH){
          arma::vec hDiagonal(rows.n_elem);
          for(arma::uword localRow = 0; localRow < rows.n_elem; ++localRow){
            hDiagonal(localRow) = Hs(rows(localRow), rows(localRow));
          }
          blockRiW.each_col() %= hDiagonal;
        }

        const arma::mat blockCross = blockW.t() * blockRiW;

        for(arma::uword localCol = 0; localCol < columns.n_elem; ++localCol){
          const arma::uword globalCol = columns(localCol);
          for(arma::uword localRow = 0; localRow < columns.n_elem; ++localRow){
            const double crossValue = blockCross(localRow, localCol);
            if(crossValue != 0.0){
              Ctriplets.emplace_back(
                static_cast<int>(columns(localRow)),
                static_cast<int>(globalCol),
                crossValue
              );
            }
          }
        }

        for(arma::uword localRow = 0; localRow < rows.n_elem; ++localRow){
          for(arma::uword localCol = 0; localCol < columns.n_elem; ++localCol){
            const double precisionDesignValue = blockRiW(localRow, localCol);
            if(precisionDesignValue != 0.0){
              RiWspTriplets.emplace_back(
                static_cast<int>(rows(localRow)),
                static_cast<int>(columns(localCol)),
                precisionDesignValue
              );
            }
          }
        }
      }

      C.setFromTriplets(Ctriplets.begin(), Ctriplets.end());
      C.makeCompressed();
      RiWsp.setFromTriplets(RiWspTriplets.begin(), RiWspTriplets.end());
      RiWsp.makeCompressed();

      {
        const Eigen::Map<const Eigen::VectorXd> RiyEig(
          Riy.memptr(), static_cast<Eigen::Index>(Riy.n_elem)
        );
        const Eigen::VectorXd rhsMMEEig = W.transpose() * RiyEig;
        rhsMME.set_size(rhsMMEEig.size());
        std::copy(rhsMMEEig.data(), rhsMMEEig.data() + rhsMMEEig.size(), rhsMME.memptr());
      }
      RiWisSparse = true;
    }else{
      Eigen::MatrixXd Wdense = Eigen::MatrixXd(W);
      arma::mat Wd(Wdense.data(), Wdense.rows(), Wdense.cols());
      RiWdense = applyRiDense(Wd);
      arma::mat Cdense = Wd.t() * RiWdense;
      // Build C directly from the dense product in a single pass, instead
      // of round-tripping through an intermediate arma::sp_mat (which would
      // itself scan Cdense once before a second scan converted it to Eigen
      // triplets).
      {
        const arma::uword nEff = Cdense.n_rows;
        std::vector<Eigen::Triplet<double>> Ctriplets;
        Ctriplets.reserve(Cdense.n_elem);
        for(arma::uword col = 0; col < Cdense.n_cols; ++col){
          for(arma::uword row = 0; row < nEff; ++row){
            const double v = Cdense(row, col);
            if(v != 0.0){
              Ctriplets.emplace_back(
                static_cast<int>(row),
                static_cast<int>(col),
                v
              );
            }
          }
        }
        C.resize(static_cast<Eigen::Index>(nEff), static_cast<Eigen::Index>(nEff));
        C.setFromTriplets(Ctriplets.begin(), Ctriplets.end());
        C.makeCompressed();
      }
      rhsMME = arma::vec(Wd.t() * Riy);
    }

    // ------------------------------------------------------------
    // Construct inverses of random-effect covariance matrices
    // and add G^-1 blocks to M
    //
    // IMPORTANT:
    // lambda must remain in THIS scope because it is subsequently
    // used to construct Wu and the score equations.
    // ------------------------------------------------------------

    arma::field<arma::sp_mat> lambda(nReAl);
    std::vector<Eigen::Triplet<double>> Gtriplets;

    if(nZs > 0){

        for(int i = 0; i < nRe; ++i){

            arma::mat thetaSym =
                arma::symmatu(
                    theta(i)
                );

            arma::mat lambdaDense;

            bool lambdaOK =
                eigenSpdInverse(
                    thetaSym,
                    lambdaDense
                );

            if(!lambdaOK){

                // The covariance update machinery normally guarantees
                // positive definiteness.  Repair only when the fast
                // SPD inverse actually fails.
                arma::mat bend =
                    nearPDcpp(
                        thetaSym,
                        100,
                        1e-06,
                        1e-07
                    );

                lambdaOK =
                    eigenSpdInverse(
                        bend,
                        lambdaDense
                    );

                if(!lambdaOK){
                    Rcpp::stop(
                        "Unable to invert a random-effect covariance matrix "
                        "even after the nearPD fallback."
                    );
                }
            }

            lambda(i) =
                arma::sp_mat(
                    lambdaDense
                );

            // Add G_i^{-1} = lambda_i \kron A_i directly by blocks.
            // This avoids materialising the potentially very large
            // Kronecker product GI.
            const std::size_t iCache =
                static_cast<std::size_t>(i);

            if(
                lambdaDense.n_rows
                !=
                partitionStartCache[iCache].size()
                ||
                lambdaDense.n_cols
                !=
                partitionStartCache[iCache].size()
            ){
                Rcpp::stop(
                    "Random-effect covariance dimensions are inconsistent "
                    "with the cached MME partitions."
                );
            }

            for(arma::uword iRow = 0; iRow < lambdaDense.n_rows; ++iRow){
                const arma::uword rowStart =
                    partitionStartCache[iCache][static_cast<std::size_t>(iRow)];
                for(arma::uword iCol = 0; iCol < lambdaDense.n_cols; ++iCol){
                    const double coefficient = lambdaDense(iRow,iCol);
                    if(coefficient == 0.0){ continue; }
                    const arma::uword colStart =
                        partitionStartCache[iCache][static_cast<std::size_t>(iCol)];
                    for(arma::sp_mat::const_iterator relEntry = Ai(i).begin();
                        relEntry != Ai(i).end(); ++relEntry){
                      Gtriplets.emplace_back(
                        static_cast<int>(rowStart + relEntry.row()),
                        static_cast<int>(colStart + relEntry.col()),
                        coefficient * (*relEntry)
                      );
                    }
                }
            }
        }

        if(!Gtriplets.empty()){
          Eigen::SparseMatrix<double> Gcontribution(nEffects, nEffects);
          Gcontribution.setFromTriplets(Gtriplets.begin(), Gtriplets.end());
          Gcontribution.makeCompressed();
          C = C + Gcontribution;
          C.makeCompressed();
        }
    }


    // ------------------------------------------------------------
    // CHANGE #2: sparse symmetric factorisation of C
    //
    // We use Eigen::SimplicialLDLT because C is the symmetric
    // mixed-model coefficient matrix.  Unlike Armadillo's
    // spsolve_factoriser interface, this factorisation exposes the
    // diagonal D, so log|C| is available without making C dense.
    // ------------------------------------------------------------

    if(
        C.rows() > static_cast<Eigen::Index>(std::numeric_limits<int>::max())
        ||
        C.cols() > static_cast<Eigen::Index>(std::numeric_limits<int>::max())
    ){
      Rcpp::stop(
        "C is too large for the 32-bit Eigen sparse index type used in ai_mme_sp()."
      );
    }

    // Stage 2 Eigen migration: C is already native Eigen sparse, so there is
    // no per-iteration Armadillo<->Eigen copy left before factorisation.
    C.makeCompressed();

    double logDetC = 0.0;
    double logDetD = 0.0; // only computed/used when reml==false (Nu==0 => 0)
    // Diagnostic only: minimum LDLT pivot in direct mode; minimum diagonal
    // entry of C in PCG mode. This must not be used by the PCG algorithm.
    double minD = std::numeric_limits<double>::quiet_NaN();
    bool reuseCselectedTopology = false;

    // D = Z'R^{-1}Z + G^{-1}, the bottom-right (nX..nEffects-1) block of C,
    // re-indexed to 0..Nu-1. This block never involves X, so it is exactly
    // the same matrix regardless of the fixed-effects design - a single
    // pass over C's stored entries suffices. Shared by both the ldlt and
    // cholmod reml=FALSE (ML) branches below.
    auto buildDFromC = [&]() -> EigenSpMat {
      std::vector<Eigen::Triplet<double>> Dtriplets;
      for(int col = 0; col < static_cast<int>(C.outerSize()); ++col){
        if(col < nX){ continue; }
        for(EigenSpMat::InnerIterator it(C, col); it; ++it){
          if(it.row() < nX){ continue; }
          Dtriplets.emplace_back(
            static_cast<int>(it.row() - nX),
            col - nX,
            it.value()
          );
        }
      }
      EigenSpMat Dmat(Nu, Nu);
      Dmat.setFromTriplets(Dtriplets.begin(), Dtriplets.end());
      Dmat.makeCompressed();
      return Dmat;
    };

    if(solverName == "ldlt"){
      const bool sameCPattern =
        CsymbolicReady && eigenSparsePatternMatches(C, CouterPattern, CinnerPattern);

      if(!sameCPattern){
        Cfactor.analyzePattern(C);
        if(Cfactor.info() != Eigen::Success){
          Rcpp::stop("Sparse symbolic analysis of the MME coefficient matrix C failed.");
        }
        cacheEigenSparsePattern(C, CouterPattern, CinnerPattern);
        CsymbolicReady = true;
      }

      Cfactor.factorize(C);
      if(Cfactor.info() != Eigen::Success && sameCPattern){
        Cfactor.analyzePattern(C);
        if(Cfactor.info() == Eigen::Success){ Cfactor.factorize(C); }
      }
      if(Cfactor.info() != Eigen::Success){
        Rcpp::stop("Sparse LDLT factorisation of the MME coefficient matrix C failed.");
      }
      CnumericReady = true;

      const Eigen::VectorXd Dldlt = Cfactor.vectorD();
      if(Dldlt.size() > 0){ minD = Dldlt.minCoeff(); }
      for(Eigen::Index j = 0; j < Dldlt.size(); ++j){
        const double dj = Dldlt(j);
        if(!std::isfinite(dj) || dj <= 0.0){
          Rcpp::stop("Sparse LDLT produced a non-positive/non-finite pivot in C.");
        }
        logDetC += std::log(dj);
      }
      reuseCselectedTopology = sameCPattern && CselectedTopologyReady;

      if(!reml && Nu > 0){
        EigenSpMat Dmat = buildDFromC();

        Dfactor.analyzePattern(Dmat);
        if(Dfactor.info() != Eigen::Success){
          Rcpp::stop("Sparse symbolic analysis of the random-effects-only matrix D failed (reml=FALSE).");
        }
        Dfactor.factorize(Dmat);
        if(Dfactor.info() != Eigen::Success){
          Rcpp::stop("Sparse LDLT factorisation of the random-effects-only matrix D failed (reml=FALSE).");
        }

        const Eigen::VectorXd DdiagLDLT = Dfactor.vectorD();
        for(Eigen::Index j = 0; j < DdiagLDLT.size(); ++j){
          const double dj = DdiagLDLT(j);
          if(!std::isfinite(dj) || dj <= 0.0){
            Rcpp::stop("Sparse LDLT produced a non-positive/non-finite pivot in D (reml=FALSE).");
          }
          logDetD += std::log(dj);
        }
      }
    }else if(solverName == "cholmod"){
      // Supernodal (BLAS-3) direct factorisation via R's Matrix package.
      // No Takahashi selected inverse is available for a supernodal factor
      // in this phase; score/AI-matrix traces instead fall back to batched
      // direct solves (see sparseTraceInverseTimes's useGenericCSolve path).
      cholmodState.ensureStarted();
      cholmod_sparse Cview = eigenToCholmodSparseView(C);

      const bool sameCPattern =
        CholmodSymbolicReady && eigenSparsePatternMatches(C, CouterPattern, CinnerPattern);

      if(!sameCPattern || cholmodState.factor == nullptr){
        if(cholmodState.factor != nullptr){
          M_cholmod_free_factor(&cholmodState.factor, &cholmodState.common);
        }
        cholmodState.factor = M_cholmod_analyze(&Cview, &cholmodState.common);
        if(cholmodState.factor == nullptr || cholmodState.common.status != CHOLMOD_OK){
          Rcpp::stop("CHOLMOD symbolic analysis of the MME coefficient matrix C failed.");
        }
        cacheEigenSparsePattern(C, CouterPattern, CinnerPattern);
        CholmodSymbolicReady = true;
      }

      const int factorizeOk =
        M_cholmod_factorize(&Cview, cholmodState.factor, &cholmodState.common);
      if(!factorizeOk || cholmodState.common.status != CHOLMOD_OK){
        Rcpp::stop("CHOLMOD supernodal factorisation of the MME coefficient matrix C failed.");
      }
      CnumericReady = true;

      logDetC = M_cholmod_factor_ldetA(cholmodState.factor);
      if(!std::isfinite(logDetC)){
        Rcpp::stop("CHOLMOD produced a non-finite log-determinant for C.");
      }

      if(!reml && Nu > 0){
        EigenSpMat Dmat = buildDFromC();

        cholmodDState.ensureStarted();
        cholmod_sparse Dview = eigenToCholmodSparseView(Dmat);

        const bool sameDPattern =
          CholmodDSymbolicReady && eigenSparsePatternMatches(Dmat, DouterPattern, DinnerPattern);

        if(!sameDPattern || cholmodDState.factor == nullptr){
          if(cholmodDState.factor != nullptr){
            M_cholmod_free_factor(&cholmodDState.factor, &cholmodDState.common);
          }
          cholmodDState.factor = M_cholmod_analyze(&Dview, &cholmodDState.common);
          if(cholmodDState.factor == nullptr || cholmodDState.common.status != CHOLMOD_OK){
            Rcpp::stop("CHOLMOD symbolic analysis of the random-effects-only matrix D failed (reml=FALSE).");
          }
          cacheEigenSparsePattern(Dmat, DouterPattern, DinnerPattern);
          CholmodDSymbolicReady = true;
        }

        const int factorizeDOk =
          M_cholmod_factorize(&Dview, cholmodDState.factor, &cholmodDState.common);
        if(!factorizeDOk || cholmodDState.common.status != CHOLMOD_OK){
          Rcpp::stop("CHOLMOD supernodal factorisation of the random-effects-only matrix D failed (reml=FALSE).");
        }

        logDetD = M_cholmod_factor_ldetA(cholmodDState.factor);
        if(!std::isfinite(logDetD)){
          Rcpp::stop("CHOLMOD produced a non-finite log-determinant for D (reml=FALSE).");
        }
      }
    }else{
      // Genuine factorisation-free MME path: no analyzePattern(), factorize(),
      // vectorD(), matrixL(), or Takahashi call is made for C.
      CnumericReady = false;
      preparePCG(C);
      // Keep the historical verbose diagnostic column defined without
      // introducing a factorisation in PCG mode.
      if(C.rows() > 0){
        minD = std::numeric_limits<double>::infinity();
        for(Eigen::Index jj = 0; jj < C.rows(); ++jj){
          minD = std::min(minD, C.coeff(jj,jj));
        }
      }
      logDetC = pcgApproxLogDet(C);
    }

    // ------------------------------------------------------------
    // Solve C * bu = W' R^-1 y
    // ------------------------------------------------------------

    Eigen::Map<const Eigen::VectorXd> rhsBu(
      rhsMME.memptr(),
      static_cast<Eigen::Index>(rhsMME.n_elem)
    );

    Eigen::VectorXd buEig =
      solveCVector(
        rhsBu,
        "BLUE/BLUP calculation"
      );

    bu.set_size(nEffects);

    for(arma::uword j = 0; j < bu.n_elem; ++j){
      bu(j) = buEig(static_cast<Eigen::Index>(j));
    }

    // ------------------------------------------------------------
    // y' P y from the Schur complement
    //
    // M = [ C  r ]
    //     [ r' d ]
    //
    // y'Py = d - r' C^{-1} r
    // ------------------------------------------------------------

    const double yPy =
      yRiy
      -
      arma::dot(
        rhsMME,
        bu
      );

    if(!std::isfinite(yPy)){
      Rcpp::stop(
        "Non-finite y'Py obtained from the sparse Schur-complement calculation."
      );
    }

    // ###########################
    // # 1.1) calculate the log-likelihood
    // ###########################

    double llikp = 0.0;

    if(nZs > 0){
      for(int i = 0; i < nRe; ++i){

        double val;
        double sign;

        bool ok1 =
          arma::log_det(
            val,
            sign,
            theta(i)
          );

        if(ok1 == false){
          Rcpp::Rcout
            << "log determinant failed "
            << arma::endl;
        }

        llikp =
          llikp
          +
          (
            nUsTotal(i)
            *
            val
            *
            sign
          )
          +
          (
            logDetA(i)
            *
            theta(i).n_rows
          );
      }
    }

    llik(iIter) =
      (-0.5)
      *
      (
        llikp
        +
        (reml ? logDetC : logDetD)
        +
        logDetR
        +
        yPy
      );

    // ==========================================================
    // GLOBAL EXACT-LIKELIHOOD STEP ACCEPTANCE / BACKTRACKING
    // ==========================================================
    //
    // The current likelihood belongs to the CURRENT covPar/theta point.
    // If that point is a proposal from the previous accepted iteration,
    // compare it with the previous accepted likelihood before computing any
    // new score/information quantities.
    //
    // A likelihood decrease is NOT convergence.  It triggers global
    // step-halving in the complete free covariance-parameter vector.
    // ==========================================================

    if(!haveAcceptedLikelihood){

      acceptedLikelihood =
        llik(iIter);

      haveAcceptedLikelihood =
        true;

    }else if(pendingLineSearch){

      const double trialLikelihood =
        llik(iIter);

      const bool likelihoodAcceptable =
        std::isfinite(trialLikelihood)
        &&
        (
          trialLikelihood
          +
          likelihoodAcceptTolerance
          >=
          acceptedLikelihood
        );

      if(!likelihoodAcceptable){

        if(lineSearchHalvings < maxLineSearchHalvings){

          ++lineSearchHalvings;
          lineSearchAlpha *= 0.5;

          arma::vec trialParameters =
            lineSearchBase
            +
            lineSearchAlpha
            *
            (
              lineSearchTarget
              -
              lineSearchBase
            );

          // Respect model-defined fixed parameters exactly.
          for(int iStruct = 0; iStruct < nRRe; ++iStruct){

            const arma::uword g0 =
              static_cast<arma::uword>(
                nVcStart(iStruct) - 1
              );

            const arma::uword g1 =
              static_cast<arma::uword>(
                nVcEnd(iStruct) - 1
              );

            arma::uvec local =
              arma::regspace<arma::uvec>(
                g0,
                static_cast<arma::uword>(1),
                g1
              );

            covPar(iStruct) =
              trialParameters(local);

            theta(iStruct) =
              evaluateStructure(
                iStruct,
                covPar(iStruct)
              );
          }

          if(verbose){
            Rcpp::Rcout
              << "  REML likelihood decreased; global step halved to alpha="
              << lineSearchAlpha
              << " (trial "
              << lineSearchHalvings
              << "/"
              << maxLineSearchHalvings
              << ")"
              << arma::endl;
          }

          // Recompute this SAME optimizer iteration at the smaller trial.
          --iIter;
          continue;

        }else{

          // No improving step was found.  Restore the previous accepted point
          // and evaluate it once more so all final C/b/u objects correspond to
          // an accepted covariance parameter vector.
          for(int iStruct = 0; iStruct < nRRe; ++iStruct){

            const arma::uword g0 =
              static_cast<arma::uword>(
                nVcStart(iStruct) - 1
              );

            const arma::uword g1 =
              static_cast<arma::uword>(
                nVcEnd(iStruct) - 1
              );

            arma::uvec local =
              arma::regspace<arma::uvec>(
                g0,
                static_cast<arma::uword>(1),
                g1
              );

            covPar(iStruct) =
              lineSearchBase(local);

            theta(iStruct) =
              evaluateStructure(
                iStruct,
                covPar(iStruct)
              );
          }

          pendingLineSearch = false;
          lineSearchStalled = true;

          if(verbose){
            Rcpp::Rcout
              << "  No likelihood-improving step found after "
              << maxLineSearchHalvings
              << " halvings; restoring the previous accepted covariance parameters."
              << arma::endl;
          }

          --iIter;
          continue;
        }

      }else{

        acceptedLikelihood =
          trialLikelihood;

        pendingLineSearch =
          false;

        if(verbose && lineSearchHalvings > 0){
          Rcpp::Rcout
            << "  Global likelihood line search accepted alpha="
            << lineSearchAlpha
            << arma::endl;
        }
      }
    }

    // Trace/AI preparation is unnecessary for likelihood trials rejected above.
    if(solverName == "ldlt"){
      buildSelectedInverseSubset(
        Cfactor,
        "C",
        CselectedTopology,
        reuseCselectedTopology
      );
      CselectedTopologyReady = true;

      if(!reml && Nu > 0){
        // ML score/trace terms need D's own selected inverse (never the
        // REML-adjusted [C^{-1}]_uu block); rebuilt fresh every iteration
        // (no topology-reuse fast path yet - correctness first).
        buildSelectedInverseSubset(
          Dfactor,
          "D",
          DselectedTopology,
          false
        );
        DselectedTopologyReady = true;
      }
    }else if(solverName == "pcg"){
      preparePCGTraceProbes(C);
    }

    b = bu(bInd); // move BLUEs to a different vector
    if(nZs > 0){ // move BLUPs to a different vector
      u = bu(uInd);
    }
    // ============================================================
    // CHANGE 8: analytic sparse-factor differentiation for the AI matrix
    // ============================================================
    //
    // The previous implementation constructed one observation-space
    // working variate per variance component,
    //
    //        w_i = V_i P y,
    //
    // and then evaluated
    //
    //        AI_ij = w_i' P w_j.
    //
    // That formulation is exact, but for many RANDOM covariance
    // parameters it creates an nR x nVC working-variate matrix and sends
    // all of those columns through R^{-1}.
    //
    // Here we forward-differentiate the already-factorised mixed-model
    // equations instead.  If
    //
    //        C b = r,
    //
    // then for variance parameter theta_i
    //
    //        C b_i = r_i - C_i b.
    //
    // The SAME sparse LDLT factorisation of C is reused for all sensitivity
    // right-hand sides.  For random-effect covariance parameters r_i = 0
    // and C_i is confined to the corresponding G^{-1} block, so no
    // observation-space random working variate is needed.
    //
    // The AI matrix is recovered from one half of the exact second
    // derivative of y' P y.  Because the marginal covariance V is linear
    // in the variance parameters,
    //
    //        0.5 d^2(y'Py)/(d theta_i d theta_j)
    //          = y' P V_i P V_j P y
    //          = AI_ij.
    //
    // For residual covariance parameters we retain the compact residual
    // working variates only.  This is intentional: the number of residual
    // covariance parameters is normally small, and this avoids explicitly
    // constructing second derivatives of R^{-1}.  Cross random/residual
    // terms are obtained from the differentiated coefficient system.
    //
    // This is an exact analytic forward differentiation of the sparse
    // factorised MME solve; no finite differences or numerical derivatives
    // are used.
    // ============================================================

    arma::field<arma::sp_mat> uSinv(nReAl);

    // Number/index at which residual variance components begin.
    const arma::uword residualVcStart =
      static_cast<arma::uword>(
        nVcStart(nRRe - 1) - 1
      );

    const arma::uword nRandomVc =
      residualVcStart;


    // Forward-sensitivity right-hand sides:
    //
    //   random component i:    -(dC_i) b
    //   residual component i:  -W'Ri S_i Ri e
    //
    // Solving C * db = sensitivityRHS gives db/dtheta.
    arma::mat sensitivityRHS(
      nEffects,
      nVcTotal,
      arma::fill::zeros
    );

    // Cache the small covariance-space objects needed for the random/random
    // curvature correction:
    //
    // d Lambda_i = -Lambda B_i Lambda
    //
    // d2 Lambda_ij =
    //     Lambda B_j Lambda B_i Lambda
    //   + Lambda B_i Lambda B_j Lambda.
    std::vector<arma::mat> randomLambda(
      static_cast<std::size_t>(nRe)
    );

    std::vector<arma::mat> randomQuadraticBase(
      static_cast<std::size_t>(nRe)
    );

    std::vector< std::vector<arma::mat> > randomParameterBasis(
      static_cast<std::size_t>(nRe)
    );

    std::vector< std::vector<arma::mat> > randomDLambda(
      static_cast<std::size_t>(nRe)
    );

    if(nZs > 0){

      for(int iR = 0; iR < nRe; ++iR){

        const std::size_t iCache =
          static_cast<std::size_t>(iR);

        arma::mat partitionsP =
          partitions(iR);

        const arma::uword nLevels =
          static_cast<arma::uword>(
            partitionsP(0,1)
            -
            partitionsP(0,0)
            +
            1
          );

        const arma::uword nTraits =
          static_cast<arma::uword>(
            partitionsP.n_rows
          );

        arma::mat U(
          nLevels,
          nTraits,
          arma::fill::zeros
        );

        for(arma::uword iCol = 0; iCol < nTraits; ++iCol){

          const arma::uvec & idx =
            partitionIndexCache[
              iCache
            ][
              static_cast<std::size_t>(iCol)
            ];

          if(idx.n_elem != nLevels){
            Rcpp::stop(
              "Random-effect partition length mismatch in the "
              "factor-differentiation AI engine."
            );
          }

          U.col(iCol) =
            bu(idx);
        }

        arma::mat lambdaDense =
          arma::mat(
            lambda(iR)
          );

        randomLambda[iCache] =
          lambdaDense;

        // Keep u * Lambda for the existing score calculation below.
        arma::mat uSinvDense =
          U
          *
          lambdaDense;

        uSinv(iR) =
          arma::sp_mat(
            uSinvDense
          );

        // A * U is reused for every local covariance derivative.
        arma::mat AU =
          arma::mat(
            Ai(iR)
            *
            U
          );

        randomQuadraticBase[iCache] =
          U.t()
          *
          AU;

        const arma::uword localNvc =
          static_cast<arma::uword>(
            nVc(iR)
          );

        randomParameterBasis[iCache].resize(
          static_cast<std::size_t>(localNvc)
        );

        randomDLambda[iCache].resize(
          static_cast<std::size_t>(localNvc)
        );

        for(arma::uword localK = 0; localK < localNvc; ++localK){

          // Generic first derivative dSigma/dphi_k.  For legacy US/DIAG
          // models this is the familiar constant covariance-cell basis;
          // for AR1 it is the analytic derivative of sigma2*rho^|i-j|.
          const arma::mat & Bk =
            cachedCovarianceD1(
              iR,
              localK
            );

          randomParameterBasis[iCache][
            static_cast<std::size_t>(localK)
          ] =
            Bk;

          arma::mat dLambda =
            -lambdaDense
            *
            Bk
            *
            lambdaDense;

          dLambda =
            0.5
            *
            (
              dLambda
              +
              dLambda.t()
            );

          randomDLambda[iCache][
            static_cast<std::size_t>(localK)
          ] =
            dLambda;

          // (dC_k) b for this random covariance structure is
          //
          //      (dLambda_k kron A) b
          //
          // and can be evaluated blockwise as
          //
          //      (A U) dLambda_k'.
          arma::mat dCbuLocal =
            AU
            *
            dLambda.t();

          const arma::uword globalK =
            static_cast<arma::uword>(
              nVcStart(iR) - 1
            )
            +
            localK;

          for(arma::uword trait = 0; trait < nTraits; ++trait){

            const arma::uvec & idx =
              partitionIndexCache[
                iCache
              ][
                static_cast<std::size_t>(trait)
              ];

            for(arma::uword rr = 0; rr < idx.n_elem; ++rr){
              sensitivityRHS(
                idx(rr),
                globalK
              ) =
                -dCbuLocal(
                  rr,
                  trait
                );
            }
          }
        }
      }
    }

    // ------------------------------------------------------------
    // Residual sensitivities.
    // ------------------------------------------------------------
    arma::vec e;
    {
      const Eigen::Map<const Eigen::VectorXd> buEig(
        bu.memptr(), static_cast<Eigen::Index>(bu.n_elem)
      );
      const Eigen::VectorXd WbuEig = W * buEig;
      e = yDense - arma::vec(WbuEig.data(), WbuEig.size());
    }

    arma::vec Rie =
      applyRiVec(
        e
      );

    // Only residual-parameter working variates are retained.  AR1 has
    // two columns (sigma2, rho) no matter how many lag bases define R.
    arma::mat residualWorking(
      nR,
      nResidualPar,
      arma::fill::zeros
    );

    for(arma::uword iP = 0; iP < nResidualPar; ++iP){
      residualWorking.col(iP) =
        arma::vec(
          residualDerivativeBasis(iP)
          *
          Rie
        );
    }

    arma::mat RiResidualWorking =
      applyRiDense(
        residualWorking
      );

    arma::mat residualCrossRHS;
    {
      const Eigen::Map<const Eigen::MatrixXd> RiResidualWorkingEig(
        RiResidualWorking.memptr(),
        static_cast<Eigen::Index>(RiResidualWorking.n_rows),
        static_cast<Eigen::Index>(RiResidualWorking.n_cols)
      );
      const Eigen::MatrixXd residualCrossRHSEig = W.transpose() * RiResidualWorkingEig;
      residualCrossRHS.set_size(residualCrossRHSEig.rows(), residualCrossRHSEig.cols());
      std::copy(residualCrossRHSEig.data(), residualCrossRHSEig.data() + residualCrossRHSEig.size(), residualCrossRHS.memptr());
    }

    for(arma::uword iP = 0; iP < nResidualPar; ++iP){

      const arma::uword globalK =
        residualVcStart
        +
        iP;

      sensitivityRHS.col(globalK) =
        -residualCrossRHS.col(iP);
    }

    // ------------------------------------------------------------
    // Differentiate the sparse MME solve.
    //
    // C is already numerically factorised.  Solve all derivative RHSs
    // simultaneously using that same sparse LDLT factor.
    // ------------------------------------------------------------
    Eigen::Map<
      const Eigen::Matrix<
        double,
        Eigen::Dynamic,
        Eigen::Dynamic,
        Eigen::ColMajor
      >
    > sensitivityRHSEig(
      sensitivityRHS.memptr(),
      static_cast<Eigen::Index>(
        sensitivityRHS.n_rows
      ),
      static_cast<Eigen::Index>(
        sensitivityRHS.n_cols
      )
    );

    Eigen::MatrixXd dBuEig =
      solveCMatrix(
        sensitivityRHSEig,
        "factor-differentiation AI sensitivity equations"
      );

    // Zero-copy Armadillo view of Eigen's column-major result.  dBuEig owns
    // the memory and remains alive for the complete AI assembly below.
    arma::mat dBu(
      dBuEig.data(),
      static_cast<arma::uword>(dBuEig.rows()),
      static_cast<arma::uword>(dBuEig.cols()),
      false,
      true
    );

    // ------------------------------------------------------------
    // Assemble the exact Average Information matrix.
    // ------------------------------------------------------------
    avInf.zeros(
      nVcTotal,
      nVcTotal
    );

    // Random rows:
    //
    // AI_ij =
    //      b' C_i b_j
    //    + 0.5 b' C_ij b
    //
    // The first term for every j is obtained in one dense BLAS product,
    // because sensitivityRHS_i = -C_i b.
    if(nRandomVc > 0){

      arma::mat randomFirstTerm =
        (
          -sensitivityRHS.cols(
            static_cast<arma::uword>(0),
            nRandomVc - 1
          )
        ).t()
        *
        dBu;

      avInf.submat(
        static_cast<arma::uword>(0),
        static_cast<arma::uword>(0),
        nRandomVc - 1,
        static_cast<arma::uword>(nVcTotal - 1)
      ) =
        randomFirstTerm;

      // Add 0.5 b' C_ij b for pairs of covariance parameters belonging
      // to the SAME random covariance structure.  Cross-structure second
      // derivatives are exactly zero.
    #ifdef _OPENMP
      #pragma omp parallel for if(nRe > 1) schedule(static)
    #endif
      for(int iR = 0; iR < nRe; ++iR){

        const std::size_t iCache =
          static_cast<std::size_t>(iR);

        const arma::mat & lambdaDense =
          randomLambda[iCache];

        const arma::mat & quadBase =
          randomQuadraticBase[iCache];

        const arma::uword localNvc =
          static_cast<arma::uword>(
            nVc(iR)
          );

        const arma::uword globalStart =
          static_cast<arma::uword>(
            nVcStart(iR) - 1
          );

        if(randomStructurallyDiagonal[iCache]){
          // Closed-form fast path: for a structurally diagonal covariance
          // product, lambdaDense and every dSigma/dphi_k (Bi/Bj) are
          // diagonal matrices, so d2Lambda collapses to
          //   d2Lambda_kk = 2 * Bi_kk * Bj_kk * lambda_kk^3
          // and secondCorrection = sum_k Bi_kk*Bj_kk*lambda_kk^3*quadBase_kk.
          // This replaces localNvc^2 dense q x q matrix-chain products
          // (O(q^3) each) with localNvc^2 O(q) dot products.
          const arma::vec lambdaCubeDiag =
            arma::pow(lambdaDense.diag(), 3);
          const arma::vec quadDiag =
            quadBase.diag();

          std::vector<arma::vec> basisDiag(localNvc);
          for(arma::uword localK = 0; localK < localNvc; ++localK){
            basisDiag[localK] =
              randomParameterBasis[iCache][
                static_cast<std::size_t>(localK)
              ].diag();
          }

#ifdef _OPENMP
          #pragma omp parallel for if(nRe <= 1 && localNvc > 1) schedule(static)
#endif
          for(arma::uword localI = 0; localI < localNvc; ++localI){
            for(arma::uword localJ = 0; localJ < localNvc; ++localJ){

              const double secondCorrection =
                arma::accu(
                  basisDiag[localI]
                  %
                  basisDiag[localJ]
                  %
                  lambdaCubeDiag
                  %
                  quadDiag
                );

              avInf(
                globalStart + localI,
                globalStart + localJ
              ) +=
                secondCorrection;
            }
          }

          continue;
        }

      #ifdef _OPENMP
        #pragma omp parallel for if(nRe <= 1 && localNvc > 1) schedule(static)
      #endif
        for(arma::uword localI = 0; localI < localNvc; ++localI){

          const arma::mat & Bi =
            randomParameterBasis[iCache][
              static_cast<std::size_t>(localI)
            ];

          // lambdaDense * Bi * lambdaDense depends only on localI, so it
          // is hoisted out of the localJ loop: this halves the dense
          // matrix-chain work of the O(q^3) AI second-derivative term.
          const arma::mat Mi =
            lambdaDense
            *
            Bi
            *
            lambdaDense;

          for(arma::uword localJ = 0; localJ < localNvc; ++localJ){

            const arma::mat & Bj =
              randomParameterBasis[iCache][
                static_cast<std::size_t>(localJ)
              ];

            // Average-information metric in generic coordinates:
            // J' AI(Sigma) J.  We intentionally use first covariance
            // derivatives only.  The nonlinear d2Sigma term belongs to the
            // exact observed Hessian, not to this positive AI/Gauss-Newton
            // metric.
            //
            // d2Lambda = lambdaDense*Bj*Mi + (lambdaDense*Bj*Mi)'
            // is algebraically identical to the original four-term chain
            // (lambdaDense*Bj*lambdaDense*Bi*lambdaDense +
            //  lambdaDense*Bi*lambdaDense*Bj*lambdaDense) since Mi is
            // symmetric and Mi*Bj*lambdaDense = (lambdaDense*Bj*Mi)'.
            const arma::mat crossTerm =
              lambdaDense
              *
              Bj
              *
              Mi;

            arma::mat d2Lambda =
              crossTerm
              +
              crossTerm.t();

            d2Lambda =
              0.5
              *
              (
                d2Lambda
                +
                d2Lambda.t()
              );

            const double secondCorrection =
              0.5
              *
              arma::accu(
                d2Lambda
                %
                quadBase
              );

            avInf(
              globalStart + localI,
              globalStart + localJ
            ) +=
              secondCorrection;
          }
        }
      }

      // Numerical roundoff can make the random/random block very slightly
      // asymmetric even though the analytical expression is symmetric.
      arma::mat randomBlock =
        avInf.submat(
          static_cast<arma::uword>(0),
          static_cast<arma::uword>(0),
          nRandomVc - 1,
          nRandomVc - 1
        );

      randomBlock =
        0.5
        *
        (
          randomBlock
          +
          randomBlock.t()
        );

      avInf.submat(
        static_cast<arma::uword>(0),
        static_cast<arma::uword>(0),
        nRandomVc - 1,
        nRandomVc - 1
      ) =
        randomBlock;
    }

    // Residual/residual block:
    //
    // w_i' R^{-1} w_j - q_i' C^{-1} q_j
    //
    // and because db_j = -C^{-1} q_j,
    //
    // = w_i' R^{-1} w_j + q_i' db_j.
    arma::mat residualGram =
      residualWorking.t()
      *
      RiResidualWorking;

    arma::mat residualSensitivity =
      dBu.cols(
        residualVcStart,
        static_cast<arma::uword>(
          nVcTotal - 1
        )
      );

    arma::mat residualAI =
      residualGram
      +
      residualCrossRHS.t()
      *
      residualSensitivity;

    residualAI =
      0.5
      *
      (
        residualAI
        +
        residualAI.t()
      );

    avInf.submat(
      residualVcStart,
      residualVcStart,
      static_cast<arma::uword>(
        nVcTotal - 1
      ),
      static_cast<arma::uword>(
        nVcTotal - 1
      )
    ) =
      residualAI;

    // Random/residual cross block was already calculated in the random
    // rows via b' C_i b_j.  Copy it explicitly to the residual/random
    // block instead of halving it through a blanket symmetrisation.
    if(nRandomVc > 0 && nResidualPar > 0){

      avInf.submat(
        residualVcStart,
        static_cast<arma::uword>(0),
        static_cast<arma::uword>(
          nVcTotal - 1
        ),
        nRandomVc - 1
      ) =
        avInf.submat(
          static_cast<arma::uword>(0),
          residualVcStart,
          nRandomVc - 1,
          static_cast<arma::uword>(
            nVcTotal - 1
          )
        ).t();
    }

    if(!avInf.is_finite()){
      Rcpp::stop(
        "Non-finite Average Information matrix produced by the "
        "factor-differentiation engine."
      );
    }

    // ##########################
    // # 5) get 1st derivatives (dL/ds2i) from MME-version
    // # PAPER FORMULA (Lee and Van der Werf, 2006)
    // # dL/ds2u = -0.5 [(Nu/s2u) - (tr(AiCuu)/s4u) -  (e/s2e)'(Zu/s2u)]
    // # dL/ds2e = -0.5 [((Nr-Nb)/s2e) - [(Nu - (tr(AiCuu)/s2u))*(1/s2e)] - ... - (e/s2e)'(e/s2e)]
    // #
    // # PAPER FORMULA (Jensen and Madsen, 1997)
    // # dL/ds2u = (q.i * lambda) - (lambda * (T + S) * lambda)  Eq. 18
    // # dL/ds2e = tr(Rij*Ri) - tr(Ci*W'*Ri*Rij*Ri*W) - (e'*Ri*Rij*Ri*e)
    // ###########################
    
    // Score traces use the sparse inverse subset; no full C^{-1}
    // is formed during iterations.
    arma::field<arma::mat> emInfList(nRRe);
    arma::vec dLu;

    // ============================================================
    // CHANGE #3A: selected inverse blocks for random-effect traces
    // ============================================================
    if(nZs > 0){

      for(int iR = 0; iR < nRe; ++iR){

        arma::mat thetaCprov = thetaC[iR];

        arma::mat traces(
          lambda(iR).n_rows,
          lambda(iR).n_cols,
          arma::fill::zeros
        );

        arma::mat partitionsP = partitions(iR);

      #ifdef _OPENMP
        #pragma omp parallel for if(lambda(iR).n_cols > 1) schedule(static)
      #endif
        for(int iCol = 0; iCol < static_cast<int>(lambda(iR).n_cols); ++iCol){

          bool columnBlockNeeded =
            covType[static_cast<std::size_t>(iR)] != "legacy";

          if(!columnBlockNeeded){
            for(int iRow = 0; iRow < static_cast<int>(lambda(iR).n_rows); ++iRow){
              if(thetaCprov(iRow, iCol) > 0){
                columnBlockNeeded = true;
                break;
              }
            }
          }
          if(!columnBlockNeeded){ continue; }

          const arma::uword colStart =
            static_cast<arma::uword>(partitionsP(iCol, 0) - 1);
          const arma::uword colEnd =
            static_cast<arma::uword>(partitionsP(iCol, 1) - 1);
          const arma::uword blockWidth = colEnd - colStart + 1;

          bool fallbackBlockAvailable = false;
          Eigen::MatrixXd fallbackBlockSolution;
          Eigen::MatrixXd pcgAiZ;

          if(solverName == "pcg"){
            const Eigen::SparseMatrix<double> & AiEigen =
              AiEigenCache[static_cast<std::size_t>(iR)];
            if(AiEigen.cols() != static_cast<Eigen::Index>(blockWidth)){
              Rcpp::stop("Random-effect inverse block dimensions are inconsistent with Ai.");
            }
            pcgAiZ.noalias() =
              AiEigen
              *
              pcgTraceZ.middleRows(
                static_cast<Eigen::Index>(colStart),
                static_cast<Eigen::Index>(blockWidth)
              );
          }

          for(int iRow = 0; iRow < static_cast<int>(lambda(iR).n_rows); ++iRow){

            if(
                covType[static_cast<std::size_t>(iR)] == "legacy"
                &&
                thetaCprov(iRow, iCol) <= 0
            ){
              continue;
            }

            // For a structurally diagonal covariance product, lambda and
            // every dSigma/dphi_k are diagonal matrices, so the score's
            // lambda*traces*lambda term only ever reads traces' diagonal
            // entries. Off-diagonal (iRow != iCol) blocks are exactly
            // unused and skipping them avoids the expensive Ai/selected-
            // inverse work below for q^2-q of the q^2 (row,col) pairs.
            if(
                covType[static_cast<std::size_t>(iR)] != "legacy"
                &&
                randomStructurallyDiagonal[static_cast<std::size_t>(iR)]
                &&
                iRow != iCol
            ){
              continue;
            }

            const arma::uword rowStart =
              static_cast<arma::uword>(partitionsP(iRow, 0) - 1);
            const arma::uword rowEnd =
              static_cast<arma::uword>(partitionsP(iRow, 1) - 1);
            const arma::uword blockHeight = rowEnd - rowStart + 1;

            double trAiCuu = 0.0;

            if(solverName == "pcg"){
              if(pcgAiZ.rows() != static_cast<Eigen::Index>(blockHeight)){
                Rcpp::stop("Random-effect inverse block dimensions are inconsistent with Ai.");
              }
              trAiCuu =
                (
                  pcgAiZ.array()
                  *
                  pcgTraceX.middleRows(
                    static_cast<Eigen::Index>(rowStart),
                    static_cast<Eigen::Index>(blockHeight)
                  ).array()
                ).sum()
                /
                static_cast<double>(pcgTraceProbes);
            }else{
              bool subsetComplete = true;
              for(arma::sp_mat::const_iterator ait = Ai(iR).begin();
                  ait != Ai(iR).end(); ++ait){
                const arma::uword ar = ait.row();
                const arma::uword ac = ait.col();
                if(ar >= blockHeight || ac >= blockWidth){
                  Rcpp::stop("Random-effect inverse block dimensions are inconsistent with Ai.");
                }
                double cij = 0.0;
                const bool haveEntry = reml
                  ? getSelectedInverseOriginal(
                      CselectedTopology,
                      static_cast<int>(colStart + ac),
                      static_cast<int>(rowStart + ar),
                      cij)
                  : getSelectedInverseOriginal(
                      DselectedTopology,
                      static_cast<int>(colStart + ac - nX),
                      static_cast<int>(rowStart + ar - nX),
                      cij);
                if(!haveEntry){
                  subsetComplete = false;
                  break;
                }
                trAiCuu += (*ait) * cij;
              }

              if(!subsetComplete){
                if(!fallbackBlockAvailable){
#ifdef _OPENMP
                  #pragma omp critical(sommer_factor_solve)
#endif
                  {
                    if(reml){
                      Eigen::MatrixXd selectedRHS = Eigen::MatrixXd::Zero(
                        static_cast<Eigen::Index>(nEffects),
                        static_cast<Eigen::Index>(blockWidth));
                      for(arma::uword j = 0; j < blockWidth; ++j){
                        selectedRHS(static_cast<Eigen::Index>(colStart + j),
                                    static_cast<Eigen::Index>(j)) = 1.0;
                      }
                      fallbackBlockSolution = solveCMatrix(
                        selectedRHS,
                        "selected-block random-effect score-trace fallback");
                    }else{
                      Eigen::MatrixXd selectedRHS = Eigen::MatrixXd::Zero(
                        static_cast<Eigen::Index>(Nu),
                        static_cast<Eigen::Index>(blockWidth));
                      for(arma::uword j = 0; j < blockWidth; ++j){
                        selectedRHS(static_cast<Eigen::Index>(colStart + j - nX),
                                    static_cast<Eigen::Index>(j)) = 1.0;
                      }
                      if(solverName == "cholmod"){
                        fallbackBlockSolution = solveDMatrixCholmod(
                          selectedRHS,
                          "selected-block random-effect score-trace fallback");
                      }else{
                        fallbackBlockSolution = Dfactor.solve(selectedRHS);
                        if(Dfactor.info() != Eigen::Success){
                          Rcpp::stop("Sparse LDLT trace fallback solve failed in selected-block random-effect score-trace fallback (reml=FALSE).");
                        }
                      }
                    }
                  }
                  fallbackBlockAvailable = true;
                }
                arma::mat inverseBlock(blockHeight, blockWidth);
                const arma::uword rowOffset = reml ? 0 : nX;
                for(arma::uword rr = 0; rr < blockHeight; ++rr){
                  for(arma::uword cc = 0; cc < blockWidth; ++cc){
                    inverseBlock(rr,cc) = fallbackBlockSolution(
                      static_cast<Eigen::Index>(rowStart + rr - rowOffset),
                      static_cast<Eigen::Index>(cc));
                  }
                }
                // trace(Ai*inverseBlock) via a full dense product is O(q^3)
                // (q=blockHeight), dominant when Ai is a dense relationship
                // matrix. Ai is symmetric, so trace(Ai*B) == accu(Ai % B);
                // accumulate directly over Ai's stored entries instead,
                // which is O(nnz(Ai)) (== O(q^2) even when Ai is dense, but
                // avoids the extra BLAS-3 pass entirely).
                double tsum = 0.0;
                for(arma::sp_mat::const_iterator ait2 = Ai(iR).begin();
                    ait2 != Ai(iR).end(); ++ait2){
                  tsum += (*ait2) * inverseBlock(ait2.row(), ait2.col());
                }
                trAiCuu = tsum;
              }
            }

            traces(iRow, iCol) = trAiCuu;
          }
        }

        traces =
          arma::symmatu(traces);

        arma::mat dLuProv =
          (
            arma::as_scalar(nUsTotal(iR))
            *
            lambda(iR)
          )
          -
          (
            uSinv(iR).t()
            *
            Ai(iR)
            *
            uSinv(iR)
          )
          -
          (
            lambda(iR)
            *
            traces
            *
            lambda(iR)
          );

        if(covType[static_cast<std::size_t>(iR)] == "legacy"){

          emInfList(iR) =
            buildEmInformationDiagonal(
              covPar(iR),
              arma::as_scalar(nUsTotal(iR))
            );

          dLu =
            arma::join_cols(
              dLu,
              mat_to_vecCpp2(
                arma::mat(dLuProv),
                thetaCprov
              )
            );

        }else{

          arma::vec structuredScore(
            covPar(iR).n_elem,
            arma::fill::zeros
          );

          const arma::mat scoreMatrix =
            arma::mat(
              dLuProv
            );

          for(arma::uword k = 0; k < covPar(iR).n_elem; ++k){
            const arma::mat & dSigma =
              cachedCovarianceD1(
                iR,
                k
              );

            structuredScore(k) =
              arma::accu(
                scoreMatrix
                %
                dSigma
              );
          }

          dLu =
            arma::join_cols(
              dLu,
              structuredScore
            );

          // A true EM update in nonlinear AR1 coordinates is not the same
          // diagonal closed form used for variance-cell parameters.  Use the
          // positive diagonal of the exact AI block as the stabilizing
          // information contribution.
          const arma::uword g0 =
            static_cast<arma::uword>(
              nVcStart(iR) - 1
            );

          const arma::uword g1 =
            static_cast<arma::uword>(
              nVcEnd(iR) - 1
            );

          arma::vec aiDiag =
            arma::diagvec(
              avInf.submat(
                g0,
                g0,
                g1,
                g1
              )
            );

          for(arma::uword k = 0; k < aiDiag.n_elem; ++k){
            if(!std::isfinite(aiDiag(k)) || aiDiag(k) <= tolParInv){
              aiDiag(k) = 1.0;
            }
          }

          emInfList(iR) =
            arma::diagmat(
              aiDiag
            );
        }
      }
    }

    // ============================================================
    // Residual score traces in generic working-parameter coordinates
    // ============================================================
    arma::vec dLe(
      nResidualPar,
      arma::fill::zeros
    );

    // RiWsp is constant across the iP loop below; converting it to
    // Armadillo once here avoids repeating a full triplet-based
    // Eigen->Armadillo conversion for every residual parameter.
    const arma::sp_mat RiWspArmaShared =
      RiWisSparse ? eigenSparseToArmaGlobal(RiWsp) : arma::sp_mat();

    Eigen::MatrixXd pcgRiWZ;
    Eigen::MatrixXd pcgRiWX;
    if(solverName == "pcg" && Rdiag && RiWisSparse){
      pcgRiWZ.noalias() = RiWsp * pcgTraceZ;
      pcgRiWX.noalias() = RiWsp * pcgTraceX;
    }

    for(arma::uword iP = 0; iP < nResidualPar; ++iP){

      const arma::sp_mat & Sprov =
        residualDerivativeBasis(iP);

      arma::sp_mat traceBasisDeriv;

      if(useH){
        traceBasisDeriv =
          Hs.t()
          *
          Sprov
          *
          Hs;
      }else{
        traceBasisDeriv =
          Sprov;
      }

      // First trace: tr((dR/dphi_i) * effective R^{-1}).
      double traceSRi = 0.0;

      if(Rdiag){

        for(arma::sp_mat::const_iterator it = traceBasisDeriv.begin();
            it != traceBasisDeriv.end();
            ++it){

          if(it.row() == it.col()){
            traceSRi +=
              (*it)
              *
              RdiagInv(
                static_cast<arma::uword>(
                  it.row()
                )
              );
          }
        }

      }else if(Rkron){

        // Independent per-block extraction/accumulation: parallelize over
        // blocks, and use trace(A*B) = accu(A % B) (valid since RkronInv is
        // symmetric) instead of forming the full q x q product.
        #ifdef _OPENMP
          #pragma omp parallel for if(residualKronNBlocks > 1) schedule(static) reduction(+:traceSRi)
        #endif
        for(int b = 0; b < residualKronNBlocks; ++b){

          const arma::uvec & rows =
            residualKronRows[
              static_cast<std::size_t>(b)
            ];

          arma::mat local(
            residualKronBlockSize,
            residualKronBlockSize,
            arma::fill::zeros
          );

          for(arma::uword aa = 0; aa < rows.n_elem; ++aa){
            for(arma::uword bb = 0; bb < rows.n_elem; ++bb){
              local(aa,bb) =
                traceBasisDeriv(
                  rows(aa),
                  rows(bb)
                );
            }
          }

          traceSRi +=
            arma::accu(
              local
              %
              RkronInv
            );
        }

      }else{

        bool usedRTraceFallback = false;

        traceSRi =
          sparseTraceInverseTimes(
            traceBasisDeriv,
            Rfactor,
            Rselected,
            "residual trace tr((dR/dphi) R^{-1})",
            usedRTraceFallback,
            solverName == "cholmod" ? std::function<Eigen::MatrixXd(const Eigen::Ref<const Eigen::MatrixXd> &, const std::string &)>(solveRMatrixCholmod) : nullptr
          );

        if(verbose && usedRTraceFallback && iIter == 0){
          Rcpp::Rcout
            << "Residual R trace required exact sparse-solve fallback for parameter "
            << iP + 1
            << arma::endl;
        }
      }

      bool usedCTraceFallback = false;
      double traceCorrection = 0.0;

      if(solverName == "pcg" && Rdiag && RiWisSparse){
        for(arma::sp_mat::const_iterator it = Sprov.begin();
            it != Sprov.end(); ++it){
          if(it.row() == it.col()){
            traceCorrection +=
              (*it)
              *
              (
                pcgRiWZ.row(static_cast<Eigen::Index>(it.row())).array()
                *
                pcgRiWX.row(static_cast<Eigen::Index>(it.row())).array()
              ).sum();
          }
        }
        traceCorrection /= static_cast<double>(pcgTraceProbes);
      }else{
        arma::sp_mat Btrace;

        if(RiWisSparse){
          Btrace =
            RiWspArmaShared.t()
            *
            Sprov
            *
            RiWspArmaShared;
        }else{
          arma::mat BtraceDense =
            RiWdense.t()
            *
            arma::mat(
              Sprov
              *
              RiWdense
            );
          Btrace = arma::sp_mat(BtraceDense);
        }

        traceCorrection =
          solverName == "pcg"
          ? pcgTraceCInverseTimesSparse(Btrace)
          : (reml
              ? sparseTraceInverseTimes(
                  Btrace,
                  Cfactor,
                  CselectedTopology,
                  "residual trace tr(C^{-1} W'Ri(dR/dphi)RiW)",
                  usedCTraceFallback,
                  solverName == "cholmod" ? std::function<Eigen::MatrixXd(const Eigen::Ref<const Eigen::MatrixXd> &, const std::string &)>(solveCMatrix) : nullptr
                )
              : (Nu > 0
                  ? sparseTraceInverseTimes(
                      arma::sp_mat(Btrace.submat(
                        static_cast<arma::uword>(nX),
                        static_cast<arma::uword>(nX),
                        static_cast<arma::uword>(nEffects - 1),
                        static_cast<arma::uword>(nEffects - 1)
                      )),
                      Dfactor,
                      DselectedTopology,
                      "residual trace tr(D^{-1} Z'Ri(dR/dphi)RiZ)",
                      usedCTraceFallback,
                      solverName == "cholmod" ? std::function<Eigen::MatrixXd(const Eigen::Ref<const Eigen::MatrixXd> &, const std::string &)>(solveDMatrixCholmod) : nullptr
                    )
                  : 0.0)
            );
      }

      if(verbose && solverName == "ldlt" && usedCTraceFallback && iIter == 0){
        Rcpp::Rcout
          << "Residual C trace required exact sparse-solve fallback for parameter "
          << iP + 1
          << arma::endl;
      }

      const double residualQuadratic =
        arma::dot(
          Rie,
          arma::vec(
            Sprov
            *
            Rie
          )
        );

      dLe(iP) =
        traceSRi
        -
        traceCorrection
        -
        residualQuadratic;
    }

    {
      const arma::uword g0 =
        static_cast<arma::uword>(
          nVcStart(residualStruct) - 1
        );

      const arma::uword g1 =
        static_cast<arma::uword>(
          nVcEnd(residualStruct) - 1
        );

      arma::vec aiDiag =
        arma::diagvec(
          avInf.submat(
            g0,
            g0,
            g1,
            g1
          )
        );

      for(arma::uword k = 0; k < aiDiag.n_elem; ++k){
        if(!std::isfinite(aiDiag(k)) || aiDiag(k) <= tolParInv){
          aiDiag(k) = 1.0;
        }
      }

      emInfList(residualStruct) =
        arma::diagmat(
          aiDiag
        );
    }

    dLu = join_cols( dLu, dLe );// join the random and residual first derivatives in a single vector
    
    for(int i = 0; i < nRRe; ++i){
      emInf.submat(nVcStart(i)-1, nVcStart(i)-1, nVcEnd(i)-1, nVcEnd(i)-1 ) = emInfList(i);
    }
    
    // ###########################
    // # 6) update the variance paramters using the Newton method
    // # PAPER FORMULA (Lee and Van der Werf, 2006)
    // # theta.n+1 = theta.n + (AInfi * dL/ds2)
    // ###########################
    
    arma::vec thetaUnlisted, thetaCUnlisted;
    arma::vec parLowerUnlisted, parUpperUnlisted, parScaleUnlisted;

    for(int i = 0; i < nRRe; ++i){

      thetaUnlisted =
        arma::join_cols(
          thetaUnlisted,
          covPar(i)
        );

      thetaCUnlisted =
        arma::join_cols(
          thetaCUnlisted,
          covConstraint(i)
        );

      parLowerUnlisted =
        arma::join_cols(
          parLowerUnlisted,
          covLower(i)
        );

      parUpperUnlisted =
        arma::join_cols(
          parUpperUnlisted,
          covUpper(i)
        );

      parScaleUnlisted =
        arma::join_cols(
          parScaleUnlisted,
          covScale(i)
        );
    }
    // Rcpp::Rcout << "thetaCUnlisted" << thetaCUnlisted << arma::endl;
    // create the 'weight' EM information matrix (TO BE USED LATER WITHIN THE OPTIMIZATION)
    arma::vec v2(nVcTotal, arma::fill::ones) ;
    arma::mat weightEmInfMat = diagmat(v2) * arma::as_scalar(weightEmInf(iIter));
    arma::mat weightAiInfMat = diagmat(v2) * (1 - arma::as_scalar(weightEmInf(iIter)));
    // Joint information matrix and update
    //                  AVERAGE INFORMATION                         +       EXPECTATION MAXIMIZATION
    InfMat = (weightAiInfMat * avInf) + (weightEmInfMat * emInf);

    arma::vec weightedScore =
      arma::as_scalar(weightInf(iIter))
      *
      dLu;

    bool informationSolveOK =
      solveInformationSystem(
        delta,
        InfMat,
        weightedScore,
        "main variance-component update"
      );

    if(!informationSolveOK){
      Rcpp::stop(
        "Unable to solve the variance-component information system."
      );
    }

    arma::vec expectedNewTheta =
      thetaUnlisted
      -
      delta;
    
    // #######################
    // # 7) APPLY CONSTRAINTS to VC
    // # suggestions from Madsen and Jensen (1997) and Gilmour (2019)
    // #######################
    //
    // Robust boundary / positive-definiteness handling:
    //
    //  * Use one numerical floor consistently for positive variance
    //    components and for covariance-matrix PD checks.
    //  * Touching a boundary does NOT immediately remove a parameter
    //    from the information-system update.
    //  * A positive parameter is frozen at its current boundary value
    //    only after 3 boundary hits.
    //  * PD repairs are local to the failing covariance structure.
    //  * If every proposed repair fails, reject only that structure's
    //    update and retain its previous covariance matrix rather than
    //    aborting the whole REML fit.
    // #######################

    const double vcFloor =
      std::max(
        1.0e-8,
        tolParInv
      );

    // A covariance matrix only needs to be numerically inside the open
    // positive-definite cone.  This threshold is deliberately below the
    // enforced scalar variance floor.
    const double pdCheckFloor =
      std::max(
        1.0e-12,
        0.1 * vcFloor
      );

    // Helper to impose descriptor-defined scalar search-space bounds.
    // Legacy structures reproduce the previous behaviour; structured
    // models such as AR1 can supply their own bounds (e.g. |rho| < 1).
    auto applyVarianceBounds =
      [&](arma::vec & candidate,
          const arma::uvec & indices,
          const bool markBoundary){

        for(arma::uword k = 0; k < indices.n_elem; ++k){

          const arma::uword g =
            indices(k);

          if(thetaCUnlisted(g) == 3){
            continue;
          }

          double lower =
            parLowerUnlisted(g);

          double upper =
            parUpperUnlisted(g);

          if(thetaCUnlisted(g) == 1){
            lower =
              std::max(
                lower,
                vcFloor
              );
          }

          if(
              std::isfinite(lower)
              &&
              candidate(g) < lower
          ){
            candidate(g) =
              lower;

            if(markBoundary){
              toBoundary(iIter,g) =
                1;
            }
          }

          if(
              std::isfinite(upper)
              &&
              candidate(g) > upper
          ){
            candidate(g) =
              upper;

            if(markBoundary){
              toBoundary(iIter,g) =
                1;
            }
          }
        }
      };

    // A) Apply ordinary scalar bounds to the main proposal.
    arma::uvec allParameterIndices =
      arma::regspace<arma::uvec>(
        static_cast<arma::uword>(0),
        static_cast<arma::uword>(1),
        static_cast<arma::uword>(nVcTotal - 1)
      );

    applyVarianceBounds(
      expectedNewTheta,
      allParameterIndices,
      true
    );

    // Recompute cumulative boundary-hit counts once per iteration.
    // The previous implementation accidentally treated a single boundary
    // hit as constrained even though the comments intended a 3-hit rule.
    for(int j = 0; j < toBoundary.n_cols; ++j){
      sumToBoundary(j) =
        arma::accu(
          toBoundary.col(j)
        );
    }

    // Parameters explicitly fixed by the model specification.
    arma::uvec modelFixed =
      arma::find(
        thetaCUnlisted == 3
      );

    // Positive variance parameters are frozen only after 3 boundary hits.
    // Keep this state separate from thetaC so dynamically boundary-frozen
    // parameters are not accidentally interpreted through thetaF.
    arma::uvec boundaryForced =
      arma::find(
        (thetaCUnlisted == 1)
        %
        (sumToBoundary >= 3)
      );

    arma::uvec constrained;

    if(modelFixed.n_elem > 0 && boundaryForced.n_elem > 0){
      constrained =
        arma::unique(
          arma::join_cols(
            modelFixed,
            boundaryForced
          )
        );
    }else if(modelFixed.n_elem > 0){
      constrained =
        modelFixed;
    }else{
      constrained =
        boundaryForced;
    }

    // Complement of constrained parameters.
    arma::uvec constrainedMask(
      nVcTotal,
      arma::fill::zeros
    );

    if(constrained.n_elem > 0){
      constrainedMask(constrained).ones();
    }

    arma::uvec unconstrained =
      arma::find(
        constrainedMask == 0
      );

    // Model-defined fixed working parameters simply retain their current
    // accepted value.  thetaF/addScaleParam are no longer part of the
    // Henderson covariance interface.
    if(modelFixed.n_elem > 0){
      expectedNewTheta(modelFixed) =
        thetaUnlisted(modelFixed);
    }

    // Boundary-forced parameters stay exactly at their current accepted
    // value.  They are no longer allowed to move after the third hit.
    if(boundaryForced.n_elem > 0){
      expectedNewTheta(boundaryForced) =
        thetaUnlisted(boundaryForced);
    }

    // B) If there are constrained parameters, re-solve only the genuinely
    // free variance-component block.  Merely touching a boundary once or
    // twice no longer removes a parameter from this solve.
    arma::mat InfMat_uu;
    arma::vec dLu_uu, delta_uu;

    if(constrained.n_elem > 0){

      delta(constrained).zeros();

      if(unconstrained.n_elem > 0){

        InfMat_uu =
          InfMat(
            unconstrained,
            unconstrained
          );

        dLu_uu =
          dLu(unconstrained);

        arma::vec weightedScore_uu =
          arma::as_scalar(
            weightInf(iIter)
          )
          *
          dLu_uu;

        bool constrainedSolveOK =
          solveInformationSystem(
            delta_uu,
            InfMat_uu,
            weightedScore_uu,
            "free variance-component block"
          );

        if(!constrainedSolveOK){
          Rcpp::stop(
            "Unable to solve the free variance-component information block."
          );
        }

        delta(unconstrained) =
          delta_uu;

        expectedNewTheta(unconstrained) =
          thetaUnlisted(unconstrained)
          -
          delta_uu;

        applyVarianceBounds(
          expectedNewTheta,
          unconstrained,
          true
        );
      }
    }

    // C) Delta-change bookkeeping is intentionally delayed until AFTER
    //    local PD repair and global working-coordinate trust scaling, so the
    //    recorded delta is the actual proposal sent to the likelihood line
    //    search rather than the raw unconstrained AI proposal.

    // D) Positive-definiteness fallback PER covariance structure.
    //    A failing structure is repaired locally; valid covariance
    //    structures keep their original AI/EM update unchanged.
    for(int iStruct = 0; iStruct < nRRe; ++iStruct){

      arma::uvec structIdx =
        arma::regspace<arma::uvec>(
          static_cast<arma::uword>(
            nVcStart(iStruct) - 1
          ),
          static_cast<arma::uword>(1),
          static_cast<arma::uword>(
            nVcEnd(iStruct) - 1
          )
        );

      auto structureIsPD =
        [&](const arma::vec & candidate) -> bool {

          const arma::vec localPar =
            candidate(structIdx);

          // Respect descriptor bounds before even constructing Sigma.
          for(arma::uword k = 0; k < localPar.n_elem; ++k){

            const double lo =
              covLower(iStruct)(k);

            const double hi =
              covUpper(iStruct)(k);

            if(
                std::isfinite(lo)
                &&
                localPar(k) < lo
            ){
              return false;
            }

            if(
                std::isfinite(hi)
                &&
                localPar(k) > hi
            ){
              return false;
            }
          }

          arma::mat m;

          if(!tryEvaluateStructure(iStruct, localPar, m)){
            return false;
          }

          m =
            arma::symmatu(m);

          arma::vec ev;

          const bool ok =
            arma::eig_sym(
              ev,
              m
            );

          if(
              !ok
              ||
              ev.n_elem == 0
              ||
              !ev.is_finite()
          ){
            return false;
          }

          return
            ev.min()
            >
            pdCheckFloor;
        };

      if(structureIsPD(expectedNewTheta)){
        continue;
      }

      if(verbose){
        Rcpp::Rcout
          << "Covariance structure "
          << iStruct + 1
          << " is not PD; applying local fallback."
          << arma::endl;
      }

      // Identify genuinely free parameters belonging only to this
      // covariance structure.
      std::vector<arma::uword> localFreeStd;

      for(arma::uword a = 0; a < structIdx.n_elem; ++a){

        const arma::uword g =
          structIdx(a);

        if(constrainedMask(g) == 0){
          localFreeStd.push_back(g);
        }
      }

      arma::uvec localFree(
        localFreeStd.size()
      );

      for(std::size_t k = 0; k < localFreeStd.size(); ++k){
        localFree(k) =
          localFreeStd[k];
      }

      bool repaired =
        false;

      arma::vec bestCandidate =
        expectedNewTheta;

      // ------------------------------------------------------------
      // One-dimensional descriptor structures still use log(sigma2), so the
      // parameter itself is not a variance directly. Reverting unconditionally
      // to the previous accepted value produces zero forward progress every
      // iteration for a boundary-approaching variance component (true value
      // at/near zero): the identical oversized step gets proposed and
      // rejected again next iteration, freezing the parameter permanently.
      // Instead, backtrack along the proposed step, trying progressively
      // smaller fractions of it until the structure is PD again - the same
      // step-halving idea used below for multi-parameter structures. Falls
      // through to the shared final fallback (full revert) further below
      // only if no fraction of the step is repairable (e.g. a genuine
      // non-finite proposal).
      // ------------------------------------------------------------
      if(theta(iStruct).n_rows == 1 && covPar(iStruct).n_elem == 1){

        const arma::uword g0 = structIdx(0);
        const double previousValue = thetaUnlisted(g0);
        const double proposedValue = expectedNewTheta(g0);

        if(std::isfinite(proposedValue)){

          arma::vec candidate = expectedNewTheta;
          double stepFraction = 1.0;

          for(int halvingAttempt = 0; halvingAttempt < 10 && !repaired; ++halvingAttempt){
            stepFraction *= 0.5;
            candidate(g0) = previousValue + stepFraction * (proposedValue - previousValue);
            repaired = structureIsPD(candidate);
          }

          if(repaired){
            bestCandidate = candidate;
          }
        }
      }

      if(
          !repaired
          &&
          localFree.n_elem > 0
          &&
          theta(iStruct).n_rows > 1
      ){

        // ----------------------------------------------------------
        // First local fallback: 50/50 AI-EM information.
        // ----------------------------------------------------------
        arma::mat localInfo =
          0.5
          *
          avInf(
            localFree,
            localFree
          )
          +
          0.5
          *
          emInf(
            localFree,
            localFree
          );

        arma::vec localScore =
          arma::as_scalar(
            weightInf(iIter)
          )
          *
          dLu(localFree);

        arma::vec localDelta;

        bool localOK =
          solveInformationSystem(
            localDelta,
            localInfo,
            localScore,
            "per-structure 50/50 AI-EM PD fallback"
          );

        if(localOK){

          bestCandidate =
            expectedNewTheta;

          bestCandidate(localFree) =
            thetaUnlisted(localFree)
            -
            localDelta;

          applyVarianceBounds(
            bestCandidate,
            localFree,
            true
          );

          repaired =
            structureIsPD(
              bestCandidate
            );
        }

        // ----------------------------------------------------------
        // Second local fallback: pure EM block.
        // ----------------------------------------------------------
        if(!repaired){

          localInfo =
            emInf(
              localFree,
              localFree
            );

          localOK =
            solveInformationSystem(
              localDelta,
              localInfo,
              localScore,
              "per-structure pure-EM PD fallback"
            );

          if(localOK){

            bestCandidate =
              expectedNewTheta;

            bestCandidate(localFree) =
              thetaUnlisted(localFree)
              -
              localDelta;

            applyVarianceBounds(
              bestCandidate,
              localFree,
              true
            );

            repaired =
              structureIsPD(
                bestCandidate
              );
          }
        }

        // ----------------------------------------------------------
        // Third local fallback: step-halving from the current accepted
        // covariance structure toward the best local proposal.
        //
        // The current point is the safest anchor.  Thirty halvings are
        // cheap because covariance structures are tiny and provide a
        // much more robust approach to the open PD boundary.
        // ----------------------------------------------------------
        if(!repaired){

          const arma::vec target =
            bestCandidate;

          arma::vec base =
            expectedNewTheta;

          base(structIdx) =
            thetaUnlisted(structIdx);

          // Preserve any model-defined fixed values in the candidate.
          if(modelFixed.n_elem > 0){
            for(arma::uword kk = 0; kk < modelFixed.n_elem; ++kk){
              const arma::uword g = modelFixed(kk);

              if(
                  g >= structIdx.min()
                  &&
                  g <= structIdx.max()
              ){
                base(g) =
                  expectedNewTheta(g);
              }
            }
          }

          for(
              int half = 0;
              half < 30 && !repaired;
              ++half
          ){

            const double alpha =
              std::pow(
                0.5,
                static_cast<double>(half)
              );

            arma::vec trial =
              base;

            trial(localFree) =
              base(localFree)
              +
              alpha
              *
              (
                target(localFree)
                -
                base(localFree)
              );

            applyVarianceBounds(
              trial,
              localFree,
              false
            );

            if(structureIsPD(trial)){
              bestCandidate =
                trial;

              repaired =
                true;
            }
          }
        }

        // ----------------------------------------------------------
        // Fourth local fallback: strict near-PD projection.
        //
        // The replacement nearPDcpp() guarantees strict PD for the
        // projected matrix.  Only genuinely free entries are copied
        // back so model-defined fixed entries remain untouched.
        // ----------------------------------------------------------
        if(
            !repaired
            &&
            covType[static_cast<std::size_t>(iStruct)] == "legacy"
        ){



          arma::mat bad =
            vec_to_matCpp(
              bestCandidate(structIdx),
              thetaC[iStruct]
            );

          arma::mat repairedMat =
            nearPDcpp(
              arma::symmatu(bad),
              100,
              1e-06,
              1e-07
            );

          arma::vec repairedVec =
            mat_to_vecCpp2(
              repairedMat,
              thetaC[iStruct]
            );

          for(arma::uword k = 0; k < localFree.n_elem; ++k){

            const arma::uword g =
              localFree(k);

            const arma::uword localPos =
              g
              -
              static_cast<arma::uword>(
                nVcStart(iStruct) - 1
              );

            bestCandidate(g) =
              repairedVec(localPos);
          }

          applyVarianceBounds(
            bestCandidate,
            localFree,
            true
          );

          repaired =
            structureIsPD(
              bestCandidate
            );
        }
      }

      // ------------------------------------------------------------
      // Final safe fallback:
      //
      // Reject ONLY this covariance structure's proposed step and
      // retain its previous accepted value.  This allows the other
      // covariance structures to continue updating instead of
      // terminating an otherwise well-behaved REML fit.
      // ------------------------------------------------------------
      if(!repaired){

        arma::vec rejectedCandidate =
          expectedNewTheta;

        rejectedCandidate(structIdx) =
          thetaUnlisted(structIdx);

        if(structureIsPD(rejectedCandidate)){

          bestCandidate =
            rejectedCandidate;

          repaired =
            true;

          if(verbose){
            Rcpp::Rcout
              << "Covariance structure "
              << iStruct + 1
              << " update rejected; retaining previous PD value."
              << arma::endl;
          }
        }
      }

      // If the previous accepted covariance structure itself is not PD,
      // we have an internal/model-specification inconsistency rather than
      // merely a bad proposed update.
      if(!repaired){

        Rcpp::stop(
          "Covariance structure "
          +
          std::to_string(iStruct + 1)
          +
          " is not positive definite even at its previous accepted value."
        );
      }

      expectedNewTheta(structIdx) =
        bestCandidate(structIdx);

      delta(structIdx) =
        thetaUnlisted(structIdx)
        -
        expectedNewTheta(structIdx);
    }

    // ==========================================================
    // GLOBAL WORKING-COORDINATE TRUST SCALING
    // ==========================================================
    //
    // expectedNewTheta is already model-fixed / boundary / PD repaired.
    // Limit the largest transformed-coordinate move while preserving the
    // complete joint AI direction.  The exact likelihood line search above
    // remains the final acceptance criterion.
    // ==========================================================

    arma::vec proposalStep =
      expectedNewTheta
      -
      thetaUnlisted;

    double trustScale = 1.0;

    for(arma::uword kk = 0; kk < unconstrained.n_elem; ++kk){

      const arma::uword g =
        unconstrained(kk);

      const double a =
        std::abs(
          proposalStep(g)
        );

      const double cap =
        workingStepCap(g);

      if(
          std::isfinite(a)
          &&
          std::isfinite(cap)
          &&
          cap > 0.0
          &&
          a > cap
      ){
        trustScale =
          std::min(
            trustScale,
            cap / a
          );
      }
    }

    if(trustScale < 1.0){

      expectedNewTheta =
        thetaUnlisted
        +
        trustScale
        *
        (
          expectedNewTheta
          -
          thetaUnlisted
        );

      // Fixed / dynamically constrained parameters must remain exact.
      if(constrained.n_elem > 0){
        expectedNewTheta(constrained) =
          thetaUnlisted(constrained);
      }

      delta =
        thetaUnlisted
        -
        expectedNewTheta;

      if(verbose){
        Rcpp::Rcout
          << "  Working-coordinate trust scaling applied: alpha="
          << trustScale
          << arma::endl;
      }
    }

    // C) Quantify changes in the ACTUAL safeguarded proposal.
    if(iIter == 0){

      delta_minus1 =
        delta;

    }else{

      percDelta.col(iIter) =
        delta
        /
        delta_minus1;

      delta_minus1 =
        delta;
    }

    // #######################
    // # 8) Save CURRENT ACCEPTED parameter vector in the monitor
    // #######################
    //
    // llik(iIter), C, b, u, thetaUnlisted and this monitor column now all
    // describe the same accepted parameter point.  The new proposal is not
    // committed until after the stopping decision below.
    monitor.col(iIter) =
      thetaUnlisted;

    // #######################
    // # 9) Stopping criteria
    // #######################

    time_t now = time(0);
    tm *ltm = localtime(&now);
    seconds = difftime(now,before);
    time_t before = time(0);
    localtime(&before);

    normMonitor(0,iIter) =
      arma::norm(
        delta(unconstrained),
        1
      );

    normMonitor(1,iIter) =
      arma::norm(
        dLu(unconstrained),
        1
      );

    arma::vec stopDelta =
      delta(unconstrained);

    arma::vec stopTheta =
      thetaUnlisted(unconstrained);

    normMonitor(2,iIter) =
      arma::norm(
        stopDelta,
        2
      )
      /
      (
        arma::norm(
          stopTheta,
          2
        )
        +
        tolParInv
      );

    arma::uvec restrained =
      arma::find(
        toBoundary.row(iIter) > 0
      );

    if(verbose == true){

      if(iIter == 0){
        Rcpp::Rcout
          << "iteration   "
          << " LogLik   "
          << "  wall    "
          << "cpu(sec)   "
          << "restrained   "
          << "EM weight";

        if(solverName == "ldlt"){
          Rcpp::Rcout
            << "      pivot";
        }

        Rcpp::Rcout
          << arma::endl;
      }

      Rcpp::Rcout
        << "    "
        << iIter+1
        << "      "
        << llik(iIter)
        << "   "
        << ltm->tm_hour
        << ":"
        << ltm->tm_min
        << ":"
        << ltm->tm_sec
        << "      "
        << seconds
        << "           "
        << restrained.n_elem
        << "      "
        << arma::as_scalar(weightEmInf(iIter));

      if(solverName == "ldlt"){
        Rcpp::Rcout
          << "      "
          << minD;
      }

      Rcpp::Rcout
        << arma::endl;
    }

    dLuOut =
      dLu;

    // A failed global line search means the optimizer could not identify an
    // improving step from this accepted point.  Return the accepted solution
    // rather than cycling indefinitely.  This is a stall, not a likelihood
    // decrease falsely labelled as convergence.
    if(lineSearchStalled){

      convergence = false;

      monitor = monitor.cols(0,iIter);
      normMonitor = normMonitor.cols(0,iIter);
      percDelta = percDelta.cols(0,iIter);
      llik = llik.cols(0,iIter);

      if(verbose){
        Rcpp::Rcout
          << "Optimization stopped at the last accepted REML point because the global line search stalled."
          << arma::endl;
      }

      break;
    }

    bool shouldStop = false;

    if(iIter > 0){

      const double delta_llik =
        llik(iIter)
        -
        llik(iIter-1);

      const bool likelihoodConverged =
        std::abs(
          delta_llik
        )
        <
        tolParConvLL;

      const bool parameterConverged =
        normMonitor(2,iIter)
        <
        tolParConvNorm;

      if(
          likelihoodConverged
          ||
          parameterConverged
      ){
        convergence = true;
        shouldStop = true;
      }
    }

    // Never return an unevaluated proposal simply because nIters was reached.
    if(iIter == nIters - 1){
      shouldStop = true;
    }

    if(shouldStop){

      monitor = monitor.cols(0,iIter);
      normMonitor = normMonitor.cols(0,iIter);
      percDelta = percDelta.cols(0,iIter);
      llik = llik.cols(0,iIter);

      break;
    }

    // ==========================================================
    // Queue the new joint proposal for exact-likelihood checking.
    // ==========================================================

    lineSearchBase =
      thetaUnlisted;

    lineSearchTarget =
      expectedNewTheta;

    lineSearchAlpha =
      1.0;

    lineSearchHalvings =
      0;

    pendingLineSearch =
      true;

    for(int i = 0; i < nRRe; ++i){

      arma::uvec toFill =
        arma::regspace<arma::uvec>(
          nVcStart(i)-1,
          1,
          nVcEnd(i)-1
        );

      const arma::vec proposedParameters =
        expectedNewTheta(toFill);

      arma::mat proposedTheta;

      if(tryEvaluateStructure(i, proposedParameters, proposedTheta)){
        covPar(i) =
          proposedParameters;

        theta(i) =
          proposedTheta;
      }else{
        covPar(i) =
          lineSearchBase(toFill);

        lineSearchTarget(toFill) =
          lineSearchBase(toFill);

        if(verbose){
          Rcpp::Rcout
            << "Covariance structure "
            << i + 1
            << " became invalid after joint trust scaling; retaining its previous accepted value."
            << arma::endl;
        }
      }

      if(covType[static_cast<std::size_t>(i)] == "legacy"){

        arma::mat thetaCProvNew =
          vec_to_matCpp(
            thetaCUnlisted(toFill),
            thetaC[i]
          );

        thetaC(i) =
          thetaCProvNew;
      }
    }
  }// end of iterative optimization
  
  
  double AIC = (-2 * llik((llik.n_cols-1))) + (2 * nX);
  double BIC = (-2 * llik((llik.n_cols-1))) + (log(nR) * nX);
  // move constraints to vector form binding the columns
  arma::vec thetaCUnlistedFinal;
  for (int i = 0; i < nRRe; ++i) {
    thetaCUnlistedFinal =
      arma::join_cols(
        thetaCUnlistedFinal,
        covConstraint(i)
      );
  }
  // ============================================================
  // FINAL-OUTPUT COMPUTATIONS ONLY
  // ============================================================

  // Final inverse information matrix for theta_se.
  // It is formed once AFTER optimisation and is never used for updates.
  {
    arma::mat infoIdentity =
      arma::eye<arma::mat>(
        nVcTotal,
        nVcTotal
      );

    bool finalInfoInverseOK =
      arma::solve(
        InfMatInv,
        InfMat,
        infoIdentity,
        arma::solve_opts::likely_sympd
      );

    if(!finalInfoInverseOK){
      finalInfoInverseOK =
        arma::solve(
          InfMatInv,
          InfMat,
          infoIdentity
        );
    }

    if(!finalInfoInverseOK){
      Rcpp::stop(
        "Unable to obtain the final inverse information matrix for theta_se."
      );
    }
  }

  // ------------------------------------------------------------
  // Optional final C-inverse work controlled by computeCi:
  //
  //   0 = no extra C-related work after REML convergence
  //   1 = final Takahashi sparse inverse subset
  //   2 = final full C inverse
  //
  // In mode 0, Ci and uPevList are returned empty and C is not
  // factorised again after the iterative optimisation.
  // ------------------------------------------------------------

  // Pre-existing bug fix: arma::sp_mat::reset() clears BOTH content and
  // dimensions (to 0x0), so it must be followed by set_size() before the
  // computeCi>0 fill loops below index into Ci.
  Ci.reset();
  Ci.set_size(nEffects, nEffects);
  SelectedInverseSubset finalCselected;

  if(computeCi > 0){

    if(solverName == "ldlt" && !CnumericReady){
      Rcpp::stop("No final sparse LDLT factorisation of C is available.");
    }
    if(solverName == "pcg" && !CpcgReady){
      Rcpp::stop("No final PCG coefficient operator is available.");
    }
    if(solverName == "cholmod" && !CnumericReady){
      Rcpp::stop("No final CHOLMOD factorisation of C is available.");
    }

    // Reuse the numerical factorisation from the final REML iteration.
    // C has not changed since that factorisation.

    if(computeCi == 1){

      // CHOLMOD's supernodal factor has no Takahashi selected-inverse
      // equivalent implemented here. The REML iterations already ran with
      // CHOLMOD for speed; refactorise the converged C ONCE with Eigen's
      // LDLT purely for this final extraction, reusing the same
      // buildSelectedInverseSubset() used by solver="ldlt" unmodified.
      if(solverName == "cholmod"){
        Cfactor.analyzePattern(C);
        if(Cfactor.info() != Eigen::Success){
          Rcpp::stop("Sparse LDLT symbolic analysis failed while preparing the final Takahashi extraction.");
        }
        Cfactor.factorize(C);
        if(Cfactor.info() != Eigen::Success){
          Rcpp::stop("Sparse LDLT factorisation failed while preparing the final Takahashi extraction.");
        }
      }

      // Takahashi sparse inverse subset in the filled LDLT pattern.
      buildSelectedInverseSubset(
        Cfactor,
        "final C",
        finalCselected,
        false
      );

      // Reconstruct the selected inverse subset as a sparse Ci in the
      // ORIGINAL MME ordering.
      std::vector<int> permutedToOriginal(nEffects, -1);

      for(int original = 0; original < nEffects; ++original){
        const int permuted =
          finalCselected.originalToPermuted[original];

        if(permuted < 0 || permuted >= nEffects){
          Rcpp::stop(
            "Invalid permutation while reconstructing the final Takahashi inverse subset."
          );
        }

        permutedToOriginal[permuted] =
          original;
      }

      for(int pcol = 0; pcol < nEffects; ++pcol){

        const int originalCol =
          permutedToOriginal[pcol];

        if(originalCol < 0){
          Rcpp::stop(
            "Unable to invert the final LDLT permutation while reconstructing Ci."
          );
        }

        const std::vector<int> & rr =
          finalCselected.rows[pcol];

        const std::vector<double> & vv =
          finalCselected.values[pcol];

        for(std::size_t k = 0; k < rr.size(); ++k){

          const int originalRow =
            permutedToOriginal[rr[k]];

          if(originalRow < 0){
            Rcpp::stop(
              "Unable to map a Takahashi inverse entry back to original MME ordering."
            );
          }

          const double value =
            vv[k];

          Ci(
            static_cast<arma::uword>(originalRow),
            static_cast<arma::uword>(originalCol)
          ) = value;

          if(originalRow != originalCol){
            Ci(
              static_cast<arma::uword>(originalCol),
              static_cast<arma::uword>(originalRow)
            ) = value;
          }
        }
      }

    }else if(computeCi == 2){

      // Complete inverse from the final sparse LDLT factorisation.
      Eigen::MatrixXd finalIdentity =
        Eigen::MatrixXd::Identity(
          static_cast<Eigen::Index>(nEffects),
          static_cast<Eigen::Index>(nEffects)
        );

      Eigen::MatrixXd finalCiEig =
        solveCMatrix(finalIdentity, "final full inverse of C");

      arma::mat finalCiDense(
        static_cast<arma::uword>(nEffects),
        static_cast<arma::uword>(nEffects)
      );

      std::copy(
        finalCiEig.data(),
        finalCiEig.data() + finalCiEig.size(),
        finalCiDense.memptr()
      );

      finalCiDense =
        0.5
        *
        (
          finalCiDense
          +
          finalCiDense.t()
        );

      Ci =
        arma::sp_mat(
          finalCiDense
        );
    }
  }

  // bring back to original scale the variance components
  for (int i = 0; i < nRRe; ++i) {
    theta(i) = theta(i)*vary;
  }
  bu = bu*stdy;
  b = b*stdy;
  // Rcpp::Rcout << intercept << arma::endl;
  if(intercept==true){ // if true mu is required in the first position only
    b(0,0) = b(0,0) + muy;
    bu(0,0) = bu(0,0) + muy;
  }else{ // if false mu is required all over b
    bu.submat(0, 0, (b.n_rows-1), 0) = bu.submat(0, 0, (b.n_rows-1), 0) + muy;
    b = b + muy;
  }
  if(computeCi > 0 && Ci.n_elem > 0){ Ci = Ci*vary; }
  // Transform covariance-parameter uncertainty and monitor output from
  // unconstrained working coordinates to reported natural coordinates.
  auto reportParameters =
    [&](const int iStruct,
        const arma::vec & work,
        arma::vec * jacobianOut = nullptr) -> arma::vec {

      arma::vec out = work;
      arma::vec jac(work.n_elem, arma::fill::ones);

      // Product-level scale is the only covariance convention intentionally
      // owned by the solver: internal tau = log(sigma2 / var(y)).
      out(0) = std::exp(work(0)) * vary;
      jac(0) = out(0);

      Rcpp::List cs =
        covDescriptor[static_cast<std::size_t>(iStruct)];
      Rcpp::List factors = cs["factors"];

      for(int fidx = 0; fidx < factors.size(); ++fidx){
        Rcpp::List f = Rcpp::as<Rcpp::List>(factors[fidx]);
        const int start1 = Rcpp::as<int>(f["par_start"]);
        const int end1 = Rcpp::as<int>(f["par_end"]);
        if(end1 < start1){
          continue;
        }

        if(!f.containsElementNamed("report")){
          Rcpp::stop("CovarianceFactor is missing report specification.");
        }
        Rcpp::List report = Rcpp::as<Rcpp::List>(f["report"]);
        const std::string backend =
          Rcpp::as<std::string>(report["backend"]);

        const int nFactorPar = end1 - start1 + 1;

        if(backend == "builtin"){
          Rcpp::CharacterVector tr = report["transform"];
          Rcpp::NumericVector lo = report["lower"];
          Rcpp::NumericVector hi = report["upper"];
          if(tr.size() != nFactorPar || lo.size() != nFactorPar ||
             hi.size() != nFactorPar){
            Rcpp::stop("CovarianceFactor report specification has incompatible length.");
          }

          for(int local = 0; local < nFactorPar; ++local){
            const arma::uword k =
              static_cast<arma::uword>(start1 - 1 + local);
            const std::string code = Rcpp::as<std::string>(tr[local]);
            const double eta = work(k);

            if(code == "identity"){
              out(k) = eta;
              jac(k) = 1.0;

            }else if(code == "exp"){
              out(k) = std::exp(eta);
              jac(k) = out(k);

            }else if(code == "tanh"){
              out(k) = std::tanh(eta);
              jac(k) = 1.0 - out(k)*out(k);

            }else if(code == "bounded_logit"){
              const double lower = lo[local];
              const double upper = hi[local];
              if(!std::isfinite(lower) || !std::isfinite(upper) ||
                 !(lower < upper)){
                Rcpp::stop("Invalid bounded-logit reporting interval.");
              }
              const double logistic =
                eta >= 0.0
                ? 1.0 / (1.0 + std::exp(-eta))
                : std::exp(eta) / (1.0 + std::exp(eta));
              out(k) = lower + (upper-lower)*logistic;
              jac(k) = (upper-lower)*logistic*(1.0-logistic);

            }else{
              Rcpp::stop("Unknown CovarianceFactor reporting transform: " + code);
            }
          }

        }else{
          Rcpp::stop("Unknown CovarianceFactor reporting backend: " + backend);
        }
      }

      if(jacobianOut != nullptr){
        (*jacobianOut) = jac;
      }
      return out;
    };

  arma::field<arma::vec> covParOut(nRRe);

  arma::vec finalJacobian(
    nVcTotal,
    arma::fill::ones
  );

  arma::mat monitorOut =
    monitor;

  arma::uword offset =
    0;

  for(int i = 0; i < nRRe; ++i){

    arma::vec localJac;

    covParOut(i) =
      reportParameters(
        i,
        covPar(i),
        &localJac
      );

    const arma::uword localN =
      covPar(i).n_elem;

    finalJacobian.subvec(
      offset,
      offset + localN - 1
    ) =
      localJac;

    for(arma::uword iter = 0; iter < monitor.n_cols; ++iter){

      if(
          !std::isfinite(
            monitor(
              offset,
              iter
            )
          )
      ){
        continue;
      }

      arma::vec localWork =
        monitor.submat(
          offset,
          iter,
          offset + localN - 1,
          iter
        );

      monitorOut.submat(
        offset,
        static_cast<arma::uword>(iter),
        offset + localN - 1,
        static_cast<arma::uword>(iter)
      ) =
        reportParameters(
          i,
          localWork
        );
    }

    offset +=
      localN;
  }

  if(
      InfMatInv.n_rows == finalJacobian.n_elem
      &&
      InfMatInv.n_cols == finalJacobian.n_elem
  ){

    InfMatInv =
      arma::diagmat(
        finalJacobian
      )
      *
      InfMatInv
      *
      arma::diagmat(
        finalJacobian
      );
  }

  // dLuOut=dLuOut/vary;
  // move the effects from vector to a field with matrices
  arma::field<arma::mat> uList(nRe), uPevList(nRe); // store indices of the random effects

  for (int i = 0; i < nRe; ++i) {

    arma::mat partitionsP = partitions(i);

    arma::mat uMat(
      arma::as_scalar(partitionsP(0,1) - partitionsP(0,0) + 1),
      partitionsP.n_rows
    );

    arma::mat eMat;

    if(computeCi > 0){
      eMat.set_size(
        arma::as_scalar(partitionsP(0,1) - partitionsP(0,0) + 1),
        partitionsP.n_rows
      );
    }

    for (int j = 0; j < partitionsP.n_rows; ++j) {

      const arma::uword firstIndex =
        static_cast<arma::uword>(partitionsP(j,0)-1);

      const arma::uword lastIndex =
        static_cast<arma::uword>(partitionsP(j,1)-1);

      uMat.col(j) =
        bu.submat(
          firstIndex,
          0,
          lastIndex,
          0
        );

      if(computeCi == 1){

        // Diagonal entries are always available in the Takahashi subset.
        arma::vec pevDiag(lastIndex - firstIndex + 1);

        for(arma::uword k = firstIndex; k <= lastIndex; ++k){

          double cii = 0.0;

          bool found =
            getSelectedInverseOriginal(
              finalCselected,
              static_cast<int>(k),
              static_cast<int>(k),
              cii
            );

          if(!found){
            Rcpp::stop(
              "Final Takahashi inverse subset is missing a diagonal entry required for uPevList."
            );
          }

          // finalCselected is still on the scaled-response scale.
          pevDiag(k - firstIndex) =
            cii
            *
            vary;
        }

        eMat.col(j) =
          pevDiag;

      }else if(computeCi == 2){

        // Ci has already been rescaled to the original response scale.
        eMat.col(j) =
          arma::diagvec(
            Ci.submat(
              firstIndex,
              firstIndex,
              lastIndex,
              lastIndex
            )
          );
      }
    }

    uList(i) = uMat;

    if(computeCi > 0){
      uPevList(i) = eMat;
    }else{
      // Mode 0: list shape is retained, but no PEV/inverse work is done.
      uPevList(i) = arma::mat();
    }
  }

  // return results in a list form
  // W/C are native Eigen sparse matrices internally; bridge to arma::sp_mat
  // only here, for the R-facing Matrix output.
  const arma::sp_mat Wout = eigenSparseToArmaGlobal(W);
  const arma::sp_mat Cout = eigenSparseToArmaGlobal(C);
  return Rcpp::List::create(
    Rcpp::Named("llik") = llik,
    // Rcpp::Named("M") = M,
    Rcpp::Named("W") = Wout,
    Rcpp::Named("C") = Cout,
    Rcpp::Named("Cscale") = vary,
    Rcpp::Named("b") = b,
    Rcpp::Named("u") = u,
    Rcpp::Named("bu") = bu,
    Rcpp::Named("Ci") = Ci,
    // Only mode 2 (full explicit inverse) supports arbitrary D %*% Ci %*% t(D)
    // queries; mode 1 fills only the Takahashi selected-inverse pattern, same
    // as post_mme_Cinverse_cpp()'s CiComputed = (mode == 2) for consistency.
    Rcpp::Named("CiComputed") = (computeCi == 2),
    Rcpp::Named("CiMode") = computeCi,
    Rcpp::Named("solver") = solverName,
    Rcpp::Named("pcgTol") = pcgTol,
    Rcpp::Named("pcgMaxIters") = pcgMaxIters,
    Rcpp::Named("theta") = theta,
    Rcpp::Named("covPar") = covParOut,
    Rcpp::Named("covType") = Rcpp::wrap(covType),
    Rcpp::Named("theta_se") = InfMatInv,
    Rcpp::Named("InfMat") = InfMat, //InfMat,
    Rcpp::Named("monitor") = monitorOut,
    // Rcpp::Named("constraints") = thetaCUnlistedFinal,
    Rcpp::Named("uList") = uList,
    Rcpp::Named("uPevList") = uPevList,
    Rcpp::Named("AIC") = AIC,
    Rcpp::Named("BIC") = BIC,
    Rcpp::Named("convergence") = convergence,
    Rcpp::Named("partitions") = partitions,
    Rcpp::Named("percDelta") = percDelta,
    Rcpp::Named("normMonitor") = normMonitor,
    Rcpp::Named("toBoundary") = toBoundary,
    Rcpp::Named("dLu") = dLuOut,
    Rcpp::Named("Cchol") = arma::mat()
    // Rcpp::Named("PMWu_mat") = PMWu_mat, 
    // Rcpp::Named("MWuchol") = MWuchol
  );
  
}

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

static arma::mat directEvalNativeFactor(const Rcpp::List & f,
                                        const arma::vec & localPar,
                                        const std::string & op){

  const arma::uword q =
    static_cast<arma::uword>(Rcpp::as<int>(f["dim"]));

  if(op == "identity"){
    return arma::eye<arma::mat>(q,q);
  }

  if(op == "diag"){
    if(localPar.n_elem + 1 != q){
      Rcpp::stop("Malformed diagonal covariance factor.");
    }
    arma::vec d(q, arma::fill::ones);
    for(arma::uword k = 0; k < localPar.n_elem; ++k){
      d(k+1) = std::exp(localPar(k));
    }
    return arma::diagmat(d);
  }

  if(op == "ar1"){
    if(localPar.n_elem != 1){
      Rcpp::stop("Malformed AR1 covariance factor.");
    }
    const double rho = std::tanh(localPar(0));
    arma::mat K(q,q,arma::fill::zeros);
    for(arma::uword i = 0; i < q; ++i){
      for(arma::uword j = 0; j < q; ++j){
        const arma::uword d = (i > j ? i-j : j-i);
        K(i,j) = std::pow(rho, static_cast<double>(d));
      }
    }
    return K;
  }

  if(op == "us"){
    Rcpp::IntegerVector rr = f["us_row"];
    Rcpp::IntegerVector cc = f["us_col"];
    Rcpp::LogicalVector dd = f["us_diag"];

    if(localPar.n_elem != static_cast<arma::uword>(rr.size()) ||
       rr.size() != cc.size() || rr.size() != dd.size()){
      Rcpp::stop("Malformed unstructured covariance factor.");
    }

    arma::mat L(q,q,arma::fill::zeros);
    L(0,0) = 1.0;

    for(arma::uword k = 0; k < localPar.n_elem; ++k){
      const arma::uword i =
        static_cast<arma::uword>(rr[static_cast<int>(k)] - 1);
      const arma::uword j =
        static_cast<arma::uword>(cc[static_cast<int>(k)] - 1);
      if(dd[static_cast<int>(k)]){
        L(i,j) = std::exp(localPar(k));
      }else{
        L(i,j) = localPar(k);
      }
    }

    return L * L.t();
  }

  if(op == "cor_uniform"){
    if(q < 2 || localPar.n_elem != 1){
      Rcpp::stop("Malformed compound-symmetry/uniform-correlation factor.");
    }
    const double lo = -1.0 / static_cast<double>(q - 1);
    const double eta = localPar(0);
    const double s = eta >= 0.0 ? 1.0/(1.0+std::exp(-eta)) : std::exp(eta)/(1.0+std::exp(eta));
    const double rho = lo + (1.0-lo)*s;
    arma::mat K(q, q, arma::fill::value(rho));
    K.diag().ones();
    return K;
  }

  if(op == "corh"){
    if(q < 2 || localPar.n_elem != q){
      Rcpp::stop("Malformed heterogeneous uniform-correlation factor.");
    }
    const double lo = -1.0 / static_cast<double>(q - 1);
    const double eta = localPar(0);
    const double s = eta >= 0.0 ? 1.0/(1.0+std::exp(-eta)) : std::exp(eta)/(1.0+std::exp(eta));
    const double rho = lo + (1.0-lo)*s;
    arma::mat C(q, q, arma::fill::value(rho));
    C.diag().ones();
    arma::vec variances(q, arma::fill::ones);
    for(arma::uword k = 1; k < q; ++k){
      variances(k) = std::exp(localPar(k));
    }
    arma::vec sd = arma::sqrt(variances);
    return arma::diagmat(sd) * C * arma::diagmat(sd);
  }

  if(op == "arp"){
    const int order = Rcpp::as<int>(f["order"]);
    if(order < 1 || localPar.n_elem != static_cast<arma::uword>(order) ||
       q <= static_cast<arma::uword>(order)){
      Rcpp::stop("Malformed AR(p) covariance factor.");
    }
    arma::vec pacf(static_cast<arma::uword>(order), arma::fill::zeros);
    for(int j = 0; j < order; ++j){
      pacf(static_cast<arma::uword>(j)) = std::tanh(localPar(static_cast<arma::uword>(j)));
    }
    arma::vec phi;
    for(int m = 1; m <= order; ++m){
      arma::vec next(static_cast<arma::uword>(m), arma::fill::zeros);
      next(static_cast<arma::uword>(m-1)) = pacf(static_cast<arma::uword>(m-1));
      if(m > 1){
        for(int j = 0; j < m-1; ++j){
          next(static_cast<arma::uword>(j)) =
            phi(static_cast<arma::uword>(j)) -
            pacf(static_cast<arma::uword>(m-1)) * phi(static_cast<arma::uword>(m-2-j));
        }
      }
      phi = next;
    }
    arma::mat A(static_cast<arma::uword>(order), static_cast<arma::uword>(order), arma::fill::zeros);
    arma::vec b(static_cast<arma::uword>(order), arma::fill::zeros);
    for(int kk = 1; kk <= order; ++kk){
      A(static_cast<arma::uword>(kk-1), static_cast<arma::uword>(kk-1)) += 1.0;
      for(int jj = 1; jj <= order; ++jj){
        const int d = std::abs(kk - jj);
        const double pj = phi(static_cast<arma::uword>(jj-1));
        if(d == 0){
          b(static_cast<arma::uword>(kk-1)) += pj;
        }else{
          A(static_cast<arma::uword>(kk-1), static_cast<arma::uword>(d-1)) -= pj;
        }
      }
    }
    arma::vec rInitial;
    bool ok = arma::solve(rInitial, A, b);
    if(!ok || !rInitial.is_finite()){
      Rcpp::stop("Unable to solve Yule-Walker equations for AR(p) factor.");
    }
    arma::vec rho(q, arma::fill::zeros);
    rho(0) = 1.0;
    for(int h = 1; h <= order; ++h){
      rho(static_cast<arma::uword>(h)) = rInitial(static_cast<arma::uword>(h-1));
    }
    for(arma::uword h = static_cast<arma::uword>(order+1); h < q; ++h){
      double value = 0.0;
      for(int jj = 1; jj <= order; ++jj){
        value += phi(static_cast<arma::uword>(jj-1)) * rho(h - static_cast<arma::uword>(jj));
      }
      rho(h) = value;
    }
    arma::mat K(q, q, arma::fill::zeros);
    for(arma::uword i = 0; i < q; ++i){
      for(arma::uword j = 0; j < q; ++j){
        const arma::uword d = i > j ? i-j : j-i;
        K(i,j) = rho(d);
      }
    }
    return K;
  }

  if(op == "ma"){
    const int order = Rcpp::as<int>(f["order"]);
    if(order < 1 || localPar.n_elem != static_cast<arma::uword>(order) ||
       q <= static_cast<arma::uword>(order)){
      Rcpp::stop("Malformed MA(q) covariance factor.");
    }
    arma::vec coef(static_cast<arma::uword>(order + 1), arma::fill::zeros);
    coef(0) = 1.0;
    for(int j = 1; j <= order; ++j){
      coef(static_cast<arma::uword>(j)) = localPar(static_cast<arma::uword>(j-1));
    }
    arma::vec rho(q, arma::fill::zeros);
    double gamma0 = 0.0;
    for(int j = 0; j <= order; ++j){
      const double c = coef(static_cast<arma::uword>(j));
      gamma0 += c*c;
    }
    if(!std::isfinite(gamma0) || gamma0 <= 0.0){
      Rcpp::stop("Invalid MA covariance normalization.");
    }
    rho(0) = 1.0;
    for(int h = 1; h <= order; ++h){
      double gamma = 0.0;
      for(int j = 0; j <= order-h; ++j){
        gamma += coef(static_cast<arma::uword>(j)) * coef(static_cast<arma::uword>(j+h));
      }
      rho(static_cast<arma::uword>(h)) = gamma / gamma0;
    }
    arma::mat K(q, q, arma::fill::zeros);
    for(arma::uword i = 0; i < q; ++i){
      for(arma::uword j = 0; j < q; ++j){
        const arma::uword d = i > j ? i-j : j-i;
        K(i,j) = d <= static_cast<arma::uword>(order) ? rho(d) : 0.0;
      }
    }
    return K;
  }

  if(op == "corg"){
    Rcpp::IntegerVector rr = f["corg_row"];
    Rcpp::IntegerVector cc = f["corg_col"];
    if(localPar.n_elem != static_cast<arma::uword>(rr.size()) || rr.size() != cc.size()){
      Rcpp::stop("Malformed general-correlation factor.");
    }
    arma::mat A(q, q, arma::fill::eye);
    for(arma::uword k = 0; k < localPar.n_elem; ++k){
      const arma::uword i = static_cast<arma::uword>(rr[static_cast<int>(k)] - 1);
      const arma::uword j = static_cast<arma::uword>(cc[static_cast<int>(k)] - 1);
      A(i,j) = localPar(k);
    }
    arma::mat S = A * A.t();
    arma::vec sd = arma::sqrt(S.diag());
    arma::mat denom = sd * sd.t();
    arma::mat K = S % arma::pow(denom, -1.0);
    K.diag().ones();
    return 0.5 * (K + K.t());
  }

  if(op == "fa"){
    const int order = Rcpp::as<int>(f["order"]);
    const int nload = Rcpp::as<int>(f["fa_nload"]);
    Rcpp::IntegerVector rr = f["fa_row"];
    Rcpp::IntegerVector cc = f["fa_col"];
    Rcpp::LogicalVector dd = f["fa_diag"];
    if(order < 1 || nload < 1 || rr.size() != nload || cc.size() != nload ||
       dd.size() != nload ||
       localPar.n_elem != static_cast<arma::uword>(nload + static_cast<int>(q) - 1)){
      Rcpp::stop("Malformed factor-analytic covariance factor.");
    }
    arma::mat L(q, static_cast<arma::uword>(order), arma::fill::zeros);
    double referenceLogLoading = 0.0;
    for(int a = 0; a < nload; ++a){
      if(rr[a] == 1 && cc[a] == 1 && dd[a]){
        referenceLogLoading = localPar(static_cast<arma::uword>(a));
        break;
      }
    }
    const double twiceReferenceLogLoading = 2.0 * referenceLogLoading;
    const double logScale =
      twiceReferenceLogLoading > 0.0
      ? twiceReferenceLogLoading + std::log1p(std::exp(-twiceReferenceLogLoading))
      : std::log1p(std::exp(twiceReferenceLogLoading));
    const double inverseReferenceSd = std::exp(-0.5 * logScale);
    for(int a = 0; a < nload; ++a){
      const arma::uword i = static_cast<arma::uword>(rr[a] - 1);
      const arma::uword j = static_cast<arma::uword>(cc[a] - 1);
      const double value =
        dd[a]
        ? std::exp(localPar(static_cast<arma::uword>(a)) - 0.5 * logScale)
        : localPar(static_cast<arma::uword>(a)) * inverseReferenceSd;
      L(i,j) = value;
    }
    arma::vec psi(q, arma::fill::zeros);
    psi(0) = std::exp(-logScale);
    for(arma::uword i = 1; i < q; ++i){
      psi(i) = std::exp(localPar(static_cast<arma::uword>(nload) + i - 1) - logScale);
    }
    return L * L.t() + arma::diagmat(psi);
  }

  if(op == "ante"){
    const int ncoef = Rcpp::as<int>(f["ante_ncoef"]);
    Rcpp::IntegerVector rr = f["ante_row"];
    Rcpp::IntegerVector cc = f["ante_col"];
    if(ncoef < 1 || rr.size() != ncoef || cc.size() != ncoef ||
       localPar.n_elem != static_cast<arma::uword>(ncoef + static_cast<int>(q) - 1)){
      Rcpp::stop("Malformed antedependence covariance factor.");
    }
    arma::mat T(q, q, arma::fill::eye);
    for(int a = 0; a < ncoef; ++a){
      const arma::uword i = static_cast<arma::uword>(rr[a] - 1);
      const arma::uword j = static_cast<arma::uword>(cc[a] - 1);
      T(i,j) = -localPar(static_cast<arma::uword>(a));
    }
    arma::vec innovation(q, arma::fill::ones);
    for(arma::uword i = 1; i < q; ++i){
      innovation(i) = std::exp(localPar(static_cast<arma::uword>(ncoef) + i - 1));
    }
    arma::mat Ti = arma::inv(arma::trimatl(T));
    arma::mat M = Ti * arma::diagmat(innovation) * Ti.t();
    const double scale = M(0,0);
    return M / scale;
  }

  Rcpp::stop("Unsupported native covariance evaluator opcode: " + op);
  return arma::mat();
}

static arma::mat directEvalFactor(const Rcpp::List & f, const arma::vec & localPar){

  if(!f.containsElementNamed("evaluator")){
    Rcpp::stop("CovarianceFactor is missing evaluator specification.");
  }
  Rcpp::List spec = Rcpp::as<Rcpp::List>(f["evaluator"]);
  const std::string backend = Rcpp::as<std::string>(spec["backend"]);
  const arma::uword q = static_cast<arma::uword>(Rcpp::as<int>(f["dim"]));

  arma::mat K;

  if(backend == "native"){
    if(!spec.containsElementNamed("op")){
      Rcpp::stop("Native CovarianceFactor evaluator is missing op.");
    }
    K = directEvalNativeFactor(f, localPar, Rcpp::as<std::string>(spec["op"]));
  }else if(backend == "fixed"){
    if(!spec.containsElementNamed("matrix")){
      Rcpp::stop("Fixed CovarianceFactor evaluator is missing matrix.");
    }
    K = Rcpp::as<arma::mat>(spec["matrix"]);
  }else if(backend == "R"){
    if(!spec.containsElementNamed("fun") || Rf_isNull(spec["fun"])){
      Rcpp::stop("R CovarianceFactor evaluator is missing fun(par).");
    }
    Rcpp::Function fun = spec["fun"];
    SEXP ans = fun(Rcpp::wrap(localPar));
    K = Rcpp::as<arma::mat>(ans);
  }else{
    Rcpp::stop("Unknown CovarianceFactor evaluator backend: " + backend);
  }

  if(K.n_rows != q || K.n_cols != q || !K.is_finite()){
    Rcpp::stop("CovarianceFactor evaluator returned an invalid matrix.");
  }
  return 0.5 * (K + K.t());
}

// Central finite-difference derivative. Covariance-shaping factors are
// always small (levels/traits/lags, never records or coefficients), so
// re-evaluating the factor twice per free parameter is negligible.
static arma::mat directFactorD1(const Rcpp::List & f, const arma::vec & localPar,
                                const arma::uword k, const double relStep = 1.0e-6){

  if(k >= localPar.n_elem){
    Rcpp::stop("Invalid CovarianceFactor derivative parameter index.");
  }
  if(!f.containsElementNamed("derivative")){
    Rcpp::stop("CovarianceFactor is missing derivative specification.");
  }
  Rcpp::List spec = Rcpp::as<Rcpp::List>(f["derivative"]);
  const std::string backend = Rcpp::as<std::string>(spec["backend"]);

  if(backend == "none"){
    Rcpp::stop("Derivative requested for a CovarianceFactor with no parameters.");
  }
  if(backend == "R"){
    if(!spec.containsElementNamed("fun") || Rf_isNull(spec["fun"])){
      Rcpp::stop("R derivative specification is missing fun(par,k).");
    }
    Rcpp::Function dfun = spec["fun"];
    SEXP ans = dfun(Rcpp::wrap(localPar), static_cast<int>(k + 1));
    arma::mat D = Rcpp::as<arma::mat>(ans);
    return 0.5 * (D + D.t());
  }

  const double h = relStep * (1.0 + std::abs(localPar(k)));
  arma::vec plus = localPar;
  arma::vec minus = localPar;
  plus(k) += h;
  minus(k) -= h;
  return (directEvalFactor(f, plus) - directEvalFactor(f, minus)) / (2.0*h);
}

static arma::mat directEvaluateDescriptor(const Rcpp::List & cs, const arma::vec & par){

  if(par.n_elem < 1){
    Rcpp::stop("Covariance descriptor must contain log_sigma2.");
  }
  Rcpp::List factors = cs["factors"];
  arma::mat K(1,1,arma::fill::ones);

  for(int fidx = 0; fidx < factors.size(); ++fidx){
    Rcpp::List f = Rcpp::as<Rcpp::List>(factors[fidx]);
    const int start1 = f.containsElementNamed("par_start") ? Rcpp::as<int>(f["par_start"]) : 1;
    const int end1 = f.containsElementNamed("par_end") ? Rcpp::as<int>(f["par_end"]) : 0;
    arma::vec localPar;
    if(end1 >= start1){
      localPar = par.subvec(static_cast<arma::uword>(start1-1), static_cast<arma::uword>(end1-1));
    }
    K = arma::kron(K, directEvalFactor(f, localPar));
  }
  return std::exp(par(0)) * K;
}

static arma::mat directDescriptorD1(const Rcpp::List & cs, const arma::vec & par,
                                    const arma::uword k){

  if(k == 0){
    return directEvaluateDescriptor(cs, par);
  }
  Rcpp::List factors = cs["factors"];
  arma::mat K(1,1,arma::fill::ones);
  bool found = false;

  for(int fidx = 0; fidx < factors.size(); ++fidx){
    Rcpp::List f = Rcpp::as<Rcpp::List>(factors[fidx]);
    const int start1 = f.containsElementNamed("par_start") ? Rcpp::as<int>(f["par_start"]) : 1;
    const int end1 = f.containsElementNamed("par_end") ? Rcpp::as<int>(f["par_end"]) : 0;
    arma::vec localPar;
    if(end1 >= start1){
      localPar = par.subvec(static_cast<arma::uword>(start1-1), static_cast<arma::uword>(end1-1));
    }
    arma::mat piece;
    const int k1 = static_cast<int>(k) + 1;
    if(end1 >= start1 && k1 >= start1 && k1 <= end1){
      piece = directFactorD1(f, localPar, static_cast<arma::uword>(k1-start1));
      found = true;
    }else{
      piece = directEvalFactor(f, localPar);
    }
    K = arma::kron(K, piece);
  }
  if(!found){
    Rcpp::stop("Covariance derivative parameter does not belong to any factor.");
  }
  return std::exp(par(0)) * K;
}

// Transforms unconstrained working parameters (log_sigma2, atanh(rho), ...)
// to reported natural-scale values, with the elementwise Jacobian used for
// the delta-method transform of the parameter covariance matrix.
static arma::vec directReportParameters(const Rcpp::List & cs, const arma::vec & work,
                                        const double vary, arma::vec & jacobianOut){

  arma::vec out = work;
  arma::vec jac(work.n_elem, arma::fill::ones);

  out(0) = std::exp(work(0)) * vary;
  jac(0) = out(0);

  Rcpp::List factors = cs["factors"];
  for(int fidx = 0; fidx < factors.size(); ++fidx){
    Rcpp::List f = Rcpp::as<Rcpp::List>(factors[fidx]);
    const int start1 = Rcpp::as<int>(f["par_start"]);
    const int end1 = Rcpp::as<int>(f["par_end"]);
    if(end1 < start1){ continue; }
    if(!f.containsElementNamed("report")){
      Rcpp::stop("CovarianceFactor is missing report specification.");
    }
    Rcpp::List report = Rcpp::as<Rcpp::List>(f["report"]);
    const std::string backend = Rcpp::as<std::string>(report["backend"]);
    const int nFactorPar = end1 - start1 + 1;

    if(backend != "builtin"){
      Rcpp::stop("Unknown CovarianceFactor reporting backend: " + backend);
    }
    Rcpp::CharacterVector tr = report["transform"];
    Rcpp::NumericVector lo = report["lower"];
    Rcpp::NumericVector hi = report["upper"];
    if(tr.size() != nFactorPar || lo.size() != nFactorPar || hi.size() != nFactorPar){
      Rcpp::stop("CovarianceFactor report specification has incompatible length.");
    }
    for(int local = 0; local < nFactorPar; ++local){
      const arma::uword k = static_cast<arma::uword>(start1 - 1 + local);
      const std::string code = Rcpp::as<std::string>(tr[local]);
      const double eta = work(k);
      if(code == "identity"){
        out(k) = eta; jac(k) = 1.0;
      }else if(code == "exp"){
        out(k) = std::exp(eta); jac(k) = out(k);
      }else if(code == "tanh"){
        out(k) = std::tanh(eta); jac(k) = 1.0 - out(k)*out(k);
      }else if(code == "bounded_logit"){
        const double lower = lo[local], upper = hi[local];
        if(!std::isfinite(lower) || !std::isfinite(upper) || !(lower < upper)){
          Rcpp::stop("Invalid bounded-logit reporting interval.");
        }
        const double logistic = eta >= 0.0 ? 1.0/(1.0+std::exp(-eta)) : std::exp(eta)/(1.0+std::exp(eta));
        out(k) = lower + (upper-lower)*logistic;
        jac(k) = (upper-lower)*logistic*(1.0-logistic);
      }else{
        Rcpp::stop("Unknown CovarianceFactor reporting transform: " + code);
      }
    }
  }
  jacobianOut = jac;
  return out;
}

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
                              const bool & reml = true
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

  double vary = arma::var(arma::vec(arma::mat(y0).col(0)));
  if(!std::isfinite(vary) || vary <= 0.0){
    Rcpp::stop("Response variance must be positive and finite.");
  }
  const double muy = arma::mean(arma::vec(arma::mat(y0).col(0)));
  const double stdy = std::sqrt(vary);

  arma::vec y = (arma::vec(arma::mat(y0).col(0)) - muy) / stdy;

  bool intercept = false;
  if(X.n_cols > 0 && arma::accu(X.col(0)) == X.n_rows){
    intercept = true;
  }
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
    covPar(i) = Rcpp::as<arma::vec>(cs["par"]);
    Rcpp::LogicalVector freeR = cs["free"];
    if(freeR.size() != static_cast<int>(covPar(i).n_elem)){
      Rcpp::stop("covStruct$free and covStruct$par have inconsistent lengths.");
    }
    covFree(i).set_size(covPar(i).n_elem);
    for(arma::uword k = 0; k < covPar(i).n_elem; ++k){
      covFree(i)(k) = freeR[static_cast<int>(k)] ? 1.0 : 0.0;
    }
    covLower(i) = arma::vec(covPar(i).n_elem, arma::fill::value(-std::numeric_limits<double>::infinity()));
    covUpper(i) = arma::vec(covPar(i).n_elem, arma::fill::value(std::numeric_limits<double>::infinity()));

    covPar(i)(0) -= std::log(vary);
    covLower(i)(0) = std::log(std::max(1.0e-8, tolParInv)) - std::log(vary);

    theta(i) = directEvaluateDescriptor(cs, covPar(i));
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
  const int residualDim = static_cast<int>(theta(residualStruct).n_rows);
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

  // ------------------------------------------------------------
  // Main covariance assembly: V(par) and its derivatives.
  // ------------------------------------------------------------
  auto buildV = [&](const arma::field<arma::vec> & par) -> arma::mat {
    arma::mat V = buildR(directEvaluateDescriptor(covDescriptor[residualStruct], par(residualStruct)));
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
      return buildR(directDescriptorD1(covDescriptor[residualStruct], par(residualStruct), k));
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
    arma::mat Th = directEvaluateDescriptor(covDescriptor[iStruct], par);
    arma::vec eigval;
    bool ok = arma::eig_sym(eigval, arma::symmatu(Th));
    return ok && eigval.n_elem > 0 && eigval.min() > tolParInv;
  };

  time_t before = time(0);
  localtime(&before);

  // Constant across iterations: only free (non-fixed) parameters are
  // updated, and their bounds never change.
  arma::vec freeUnlisted(nVcTotal), lowerUnlisted(nVcTotal), upperUnlisted(nVcTotal);
  for(int i = 0; i < nRRe; ++i){
    freeUnlisted.subvec(nVcStart(i)-1, nVcEnd(i)-1) = covFree(i);
    lowerUnlisted.subvec(nVcStart(i)-1, nVcEnd(i)-1) = covLower(i);
    upperUnlisted.subvec(nVcStart(i)-1, nVcEnd(i)-1) = covUpper(i);
  }
  arma::uvec freeIdx = arma::find(freeUnlisted > 0.5);

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

  for(int iIter = 0; iIter < nIters; ++iIter){
    lastIter = iIter;

    arma::mat V = buildV(covPar);
    V = arma::symmatu(V);

    arma::mat Vi;
    bool okV = arma::inv_sympd(Vi, V);
    arma::mat D = arma::eye<arma::mat>(nR, nR);
    for(int jit = 0; !okV && jit < 4; ++jit){
      V = V + D * (tolParInv * std::pow(10.0, jit));
      okV = arma::inv_sympd(Vi, V);
    }
    if(!okV){
      Rcpp::stop("V is numerically singular. Try a larger tolParInv or check the supplied covariance structures.");
    }

    arma::mat VX = Vi * Xd;
    arma::mat Q = Xd.t() * VX;
    arma::mat Qinv;
    bool okQ = arma::inv_sympd(Qinv, Q);
    arma::mat Dx = arma::eye<arma::mat>(nX, nX);
    for(int jit = 0; !okQ && jit < 4; ++jit){
      Q = Q + Dx * (tolParInv * std::pow(10.0, jit));
      okQ = arma::inv_sympd(Qinv, Q);
    }
    if(!okQ){
      Rcpp::stop("X'V^{-1}X is numerically singular. Try a larger tolParInv or check the fixed-effects design for collinearity.");
    }

    arma::mat P = Vi - VX * Qinv * VX.t();
    P = arma::symmatu(P);
    Py = P * y;

    double valV, signV, valQ, signQ;
    arma::log_det(valV, signV, V);
    const double logDetV = valV;
    double logDetQ = 0.0;
    if(reml){
      arma::log_det(valQ, signQ, Q);
      logDetQ = valQ;
    }
    const double yPy = arma::as_scalar(y.t() * Py);
    const double n = static_cast<double>(nR);
    const double constTerm = reml ? (n - static_cast<double>(rankX)) * std::log(2.0*arma::datum::pi)
                                  : n * std::log(2.0*arma::datum::pi);
    llik(iIter) = -0.5 * (logDetV + logDetQ + yPy + constTerm) + llikScaleCorrection;

    if(iIter > 0){
      const double deltaLL = llik(iIter) - llik(iIter-1);
      if(deltaLL < tolParConvLL){
        convergence = true;
        monitor.col(iIter) = monitor.col(iIter-1);
        break;
      }
    }

    // Score and average information.
    arma::vec score(nVcTotal, arma::fill::zeros);
    arma::mat AI(nVcTotal, nVcTotal, arma::fill::zeros);
    std::vector<arma::vec> vList(nVcTotal), PvList(nVcTotal);

    int g = 0;
    for(int i = 0; i < nRRe; ++i){
      for(arma::uword k = 0; k < covPar(i).n_elem; ++k){
        arma::mat Dk = buildD(covPar, i, k);
        arma::vec vk = Dk * Py;
        const double traceTerm = arma::accu(P % Dk);
        score(g) = 0.5 * (arma::dot(Py, vk) - traceTerm);
        vList[g] = vk;
        PvList[g] = P * vk;
        g++;
      }
    }
    for(int a = 0; a < nVcTotal; ++a){
      for(int b = a; b < nVcTotal; ++b){
        AI(a,b) = 0.5 * arma::dot(vList[a], PvList[b]);
        AI(b,a) = AI(a,b);
      }
    }

    arma::vec emDiag = arma::diagvec(AI);
    for(int k = 0; k < nVcTotal; ++k){
      if(!std::isfinite(emDiag(k)) || emDiag(k) <= tolParInv){ emDiag(k) = 1.0; }
    }
    arma::mat InfMat = (1.0 - weightEmInf(iIter)) * AI + weightEmInf(iIter) * arma::diagmat(emDiag);

    arma::vec thetaUnlisted(nVcTotal);
    for(int i = 0; i < nRRe; ++i){
      thetaUnlisted.subvec(nVcStart(i)-1, nVcEnd(i)-1) = covPar(i);
    }

    arma::vec delta(nVcTotal, arma::fill::zeros);
    arma::vec weightedScore = weightInf(iIter) * score;
    if(freeIdx.n_elem > 0){
      arma::mat InfFree = InfMat.submat(freeIdx, freeIdx);
      arma::vec scoreFree = weightedScore(freeIdx);
      arma::vec deltaFree;
      bool solvedOK = arma::solve(deltaFree, InfFree, scoreFree, arma::solve_opts::likely_sympd);
      if(!solvedOK){ solvedOK = arma::solve(deltaFree, arma::pinv(InfFree), scoreFree); }
      if(!solvedOK){
        Rcpp::stop("Unable to solve the variance-component information system.");
      }
      delta(freeIdx) = deltaFree;
    }

    arma::vec candidate = thetaUnlisted + delta;
    for(int k = 0; k < nVcTotal; ++k){
      if(freeUnlisted(k) < 0.5){ candidate(k) = thetaUnlisted(k); continue; }
      if(candidate(k) < lowerUnlisted(k)){ candidate(k) = lowerUnlisted(k); }
      if(candidate(k) > upperUnlisted(k)){ candidate(k) = upperUnlisted(k); }
    }

    for(int halving = 0; halving < 10; ++halving){
      bool allPD = true;
      for(int i = 0; i < nRRe && allPD; ++i){
        arma::vec parCandidate = candidate.subvec(nVcStart(i)-1, nVcEnd(i)-1);
        if(!structureIsPD(i, parCandidate)){ allPD = false; }
      }
      if(allPD){ break; }
      candidate = thetaUnlisted + 0.5 * (candidate - thetaUnlisted);
      if(halving == 9){ candidate = thetaUnlisted; }
    }

    for(int i = 0; i < nRRe; ++i){
      covPar(i) = candidate.subvec(nVcStart(i)-1, nVcEnd(i)-1);
      theta(i) = directEvaluateDescriptor(covDescriptor[i], covPar(i));
    }
    monitor.col(iIter) = candidate;

    if(verbose){
      time_t now = time(0);
      tm *ltm = localtime(&now);
      double seconds = difftime(now, before);
      before = time(0);
      if(iIter == 0){ Rcpp::Rcout << "iteration   LogLik   wall     cpu(sec)" << arma::endl; }
      Rcpp::Rcout << "    " << iIter+1 << "      " << llik(iIter) << "   "
                  << ltm->tm_hour << ":" << ltm->tm_min << ":" << ltm->tm_sec
                  << "      " << seconds << arma::endl;
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

  // Clean (non-EM-blended) average information at the converged point, used
  // only for theta_se: the per-iteration InfMat blends in a diagonal EM
  // stabilizer for optimizer robustness, which must not leak into reported
  // standard errors.
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
    if(freeIdx.n_elem > 0){
      InfMatInvWorking.submat(freeIdx, freeIdx) = arma::pinv(AIFinal.submat(freeIdx, freeIdx));
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
        const int startA = static_cast<int>(partitions(i)(a,0)) - 1;
        const int endA = static_cast<int>(partitions(i)(a,1)) - 1;
        Zfull.cols(startA, endA) = Zdense[i][a];
        for(arma::uword b = 0; b < q; ++b){
          const int startB = static_cast<int>(partitions(i)(b,0)) - 1;
          const int endB = static_cast<int>(partitions(i)(b,1)) - 1;
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
        const int startA = static_cast<int>(partitions(i)(a,0)) - 1;
        const int endA = static_cast<int>(partitions(i)(a,1)) - 1;
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
  if(intercept){ beta(0) += muy; }
  VarBeta = VarBeta * vary;
  u = u * stdy;
  for(int i = 0; i < nRe; ++i){ uList(i) = uList(i) * stdy; }

  arma::field<arma::vec> covParOut(nRRe);
  arma::vec finalJacobian(nVcTotal, arma::fill::ones);
  arma::mat monitorOut = monitor;
  for(int i = 0; i < nRRe; ++i){
    arma::vec localJac;
    covParOut(i) = directReportParameters(covDescriptor[i], covPar(i), vary, localJac);
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
    Rcpp::Named("theta_se") = theta_se,
    Rcpp::Named("monitor") = monitorOut.cols(0, lastIter),
    Rcpp::Named("uList") = uList,
    Rcpp::Named("uPevList") = uPevList,
    Rcpp::Named("AIC") = AIC,
    Rcpp::Named("BIC") = BIC,
    Rcpp::Named("convergence") = convergence,
    Rcpp::Named("partitions") = partitions,
    Rcpp::Named("CiMode") = 0,
    Rcpp::Named("Ci") = arma::sp_mat()
  );
}


