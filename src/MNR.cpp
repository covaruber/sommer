// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#define ARMA_DONT_PRINT_ERRORS
#include "RcppArmadillo.h"
#include "stdlib.h" 

#include <progress.hpp>

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
arma::mat amat(const arma::mat & Xo, const bool & endelman, double minMAF) {
  
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
  
  if(endelman == true){ //  Endelman
    
    arma::rowvec ms012 = mean( X+1, 0 ); // means of columns
    arma::rowvec freq = ms012/2;
    double v = 2 * mean(freq % (1 - freq));
    
    arma::mat one(n, 1, arma::fill::ones);
    arma::mat freqmat = one * freq;
    arma::mat W = (X + 1) - (2 * freqmat);
    // 
    arma::mat K = W * W.t();
    A = K/v/p;
    
  }else{ // regular vanRaden
    
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
  int p = X.n_cols;
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
arma::rowvec scorecalc(const arma::mat & Mimv, 
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
  // make the test if the inversion went well
  arma::rowvec score(nt, arma::fill::zeros);
  if(Winv.n_rows > 0 && MAF > minMAF){ 
    arma::mat XZMimvVy = XZMimv.t() * (Vinv*Ymv); // XZM Vi y
    arma::colvec b = Winv * XZMimvVy; // (XZV-ZX)- XZV-y
    arma::colvec e = Ymv - (XZMimv * b); // Y - XB
    arma::mat mVar = (e.t() * (Vinv * e))/v2; // eV-e/(n-p) = variance
    double mVarAsDouble = mVar(0,0);
    arma::mat bVar = Winv * mVarAsDouble; // eVe * B
    // extract the right fixed effect
    arma::mat bn(b.n_rows,1); // empty vector 
    for (int i = 0; i < b.n_rows; ++i) {
      bn(i) = i; // fill it with their own position
    }
    arma::uvec ps = find(bn > (Xmv.n_cols-1)); // index for good markers
    arma::colvec bMarker = b(ps); // beta for marker
    arma::mat bMarkerVar = bVar(ps,ps); // var beta for marker
    arma::vec fStat = arma::pow(bMarker,2)/diagvec(bMarkerVar);//F statistic
    arma::vec x = v2/(v2 + v1 * fStat); // probabilty
    for (int j = 0; j < nt; ++j) {
      score(j) = x(j);
    }
    // score = x;//-1 * (log10(x)); // log10(pbeta(x, v2/2, v1/2)); 
  }
  
  return score; // return score, fStat, bMarker, R2
}

// [[Rcpp::depends(RcppProgress)]]
// [[Rcpp::export]]
arma::mat gwasForLoop(const arma::mat & M, // marker matrix 
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
  arma::mat scores(n_marker,nt, arma::fill::zeros); 
  
  // start for loop for each marker
  Progress p(n_marker, display_progress);
  for (int i = 0; i < n_marker; ++i) {
    if (Progress::check_abort() ){
      return 0;
    }
    p.increment(); 
    arma::mat Mi = M.col(i); // extract marker i
    arma::mat Mimv = arma::kron(Mi,Dnt); // kronecker for multivariate
    arma::rowvec prov = scorecalc(Mimv,Ymv, Zmv, Xmv, Vinv, nt, minMAF);
    // fill the vector in case of multiple traits
    for (int j = 0; j < nt; ++j) {
      scores(i,j) = prov(j);
    }
  }
  
  return scores;
}

// [[Rcpp::export]]
Rcpp::List MNR(const arma::mat & Y, const Rcpp::List & X,
               const Rcpp::List & Gx,
               const Rcpp::List & Z, const Rcpp::List & K,
               const Rcpp::List & R, const Rcpp::List & Ge,
               const Rcpp::List & GeI, const arma::vec & ws,
               int iters, double tolpar, double tolparinv, 
               const bool & ai, const bool & pev, 
               const bool & verbose,const bool & retscaled) {
  
  time_t before = time(0);
  localtime(&before);
  
  int n_fixed = X.size(); // define nre=number of fixed effects
  int n_random = Z.size(); // define nre=number of random effects
  int n_rcov = R.size(); // define nre=number of residual effects
  int n_re = n_random + n_rcov; // define nre=number of total random effects z+r
  int n_traits = Y.n_cols; // define n_traits=number of traits
  int no = Y.n_rows; // define n_traits=number of traits
  // ****************************************************
  // define ZKZ' and R
  // ****************************************************
  // calculate and concatenate ZKZ' and R
  arma::cube ZKZtR(no,no,n_re);
  
  for (int i = 0; i < n_re; ++i) {
    int irw = i - n_random;
    if(i < n_random && n_random > 0){
      arma::sp_mat zp = Rcpp::as<arma::sp_mat>(Z[i]);
      bool dcheck = isIdentity_mat(Rcpp::as<arma::mat>(K[i]));
      if(dcheck == true){ // if K[i] is diagonal
        if(zp.n_rows == zp.n_cols){//is a square matrix
          bool dcheck2 = isIdentity_spmat(zp);
          if(dcheck2 == true){ // if Z[i] is diagonal
            ZKZtR.slice(i) = Rcpp::as<arma::mat>(K[i]); 
          }else{ZKZtR.slice(i) = zp * zp.t(); }
        }else{
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
      double dcheck3 = accu(ws - 1);
      if(dcheck3 == 0){ // if W is diagonal no need to multiply
        ZKZtR.slice(i) = Rcpp::as<arma::sp_mat>(R[irw]);
      }else{ // if W (weights) is not diagonal then multiply Wis R Wis 
        // arma::vec ws = W.diag();// 1 / sqrt(diagvec(W));
        arma::vec ws2 = 1/sqrt(ws);
        arma::mat Wis = diagmat(ws2); // W inverse squared
        ZKZtR.slice(i) = Wis * Rcpp::as<arma::sp_mat>(R[irw]) * Wis;
      }
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
  int no_vc = 0; // to add and find out how many VC exist in total
  for (int i = 0; i < n_re; ++i) { // for each random effect fill the cube
    sigma.slice(i) = Rcpp::as<arma::mat>(Ge[i]); // take Ge for a random effect (initial VC values) and save them in a slice
    arma::vec oo = mat_to_vecCpp(sigma.slice(i),GeI[i]) ; // extract upper triangular from that slice in a vector form, pass the constraints as 2nd argument
    sigma_ut[i] = oo; // oo is sigma2 in vector form and stored in the list sigma_ut
    constraintsL[i] = mat_to_vecCpp(GeI[i],GeI[i]) ; // who are diagonal and non-diagonal VCs, pass constraints in list form
    no_vc = no_vc + oo.n_elem; // keep adding the #of VC
  }
  // sigma_ut_un will have all VC for all random effects in a single vector
  arma::vec sigma_ut_un; // vector to unlist the LIST of VC for all random effects
  arma::vec constraints; // vector to unlist constraints
  for(int i=0; i < n_re ; i++){ // for each random effect unlist
    sigma_ut_un = join_cols(sigma_ut_un,sigma_ut[i]); // column bind vectors so we end up with a very long vector with all VC
    constraints = join_cols(constraints,constraintsL[i]); // column bind vectors so we end up with a very long vector with all constraints
  }
  arma::vec sigmaF_ut_un = sigma_ut_un; // make a copy for fixed-value vc's when we use constraints
  arma::vec coef_ut_un = sigma_ut_un; // make a 2nd copy of the same vector for stabilization
  
  arma::vec taper(iters);  taper.fill(0.9); // create the weight vector for VC
  
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
  // Rcpp::List var_comp_ret(iters);
  // Rcpp::List PdViList(kk); // list to store the multivariate derivatives * P or PVi=P*dZKZ'/ds
  
  arma::vec v(nom, arma::fill::ones); // generate enough ones for an identity matrix of dimensions nt x nt
  arma::mat Vi(nom,nom); // V or phenotypic variance matrix
  arma::mat P(nom,nom); // to fill the projection matrix
  arma::sp_mat D = arma::speye<arma::sp_mat>(nom,nom);
  arma::vec seqrankX = seqCpp(0,rankX-1); // will be used to keep only the eigen values for indices 1 to rankX
  arma::vec seqkk = seqCpp(0,kk-1);
  arma::vec popo = arma::vec(rankX, arma::fill::zeros);
  for(int i=0; i < rankX; i++){popo(i) = 1;}
  arma::mat A(kk,kk,arma::fill::zeros); // to store second derivatives
  arma::mat A_svd;
  arma::vec eigval2; // will be used for the decomposition of P, within the algorithm
  arma::mat eigvec2; // will be used for the decomposition of P
  arma::mat sigma_store(sigma_ut_un.n_elem,iters); // to store variance comp through the different iterations
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
  int cycle, cycle2;
  double ldet, llik, llik0, delta_llik, checkP, seconds; // to store likelihoods and determinants
  // ###############
  // LOOP for cycles
  // ###############
  for(cycle=0; cycle < iters; cycle++){ // for each cycle
    
    if(cycle == 0){taper(cycle) = 0.5;} // taper 0.5 to relax change in variance component in the 1st iteration
    if(cycle == 1){taper(cycle) = 0.7;} // taper 0.5 to relax change in variance component in the 2nd iteration
    
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
            Rcpp::Rcout << "System is singular (V). Stopping the job. Try a bigger number of tolparinv." << arma::endl;
            return 0;
          }
        }
      }
    }
    // if last iteration let's make Xm in the opposite direction
    if(last_iteration == true){
      // Xm = arma::kron(X,Gx);
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
    tXVXVX = arma::solve(tXVX, VX.t()); // X'VXVX
    arma::solve(tXVXVX,tXVX,VX.t());
    if(tXVXVX.n_rows == 0){ // if fails try to invert with diag(1e-6)
      arma::solve(tXVXVX,tXVX + (D*(tolparinv)),VX.t());
      if(tXVXVX.n_rows == 0){// if fails try to invert with diag(1e-5)
        arma::solve(tXVXVX,tXVX + (D*(tolparinv*10)),VX.t());
        if(tXVXVX.n_rows == 0){ // if fails try to invert with diag(1e-4)
          arma::solve(tXVXVX,tXVX + (D*(tolparinv*100)),VX.t());
          if(tXVXVX.n_rows == 0){ // finally stop
            Rcpp::Rcout << "System is singular (tXVXVX). Aborting the job. Try a bigger number of tolparinv." << arma::endl;
            return 0;
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
      
      // calculate first derivatives (dL/ds)
      arma::vec ww(kk); // vector to store the product Y'PViPY - tr(PVi) = dL/ds
      arma::cube PdViList(nom,nom,kk); // list to store the multivariate derivatives * P or PVi=P*dZKZ'/ds
      for(int i=0; i < kk; i++){
        int re = re_mapper(i);
        arma::mat zkzp = ZKZtR.slice(re); // it repeats the same ZKZtR if is a vc for the same random effect
        arma::mat PdVi = P * kron(deriv_dummy.slice(i),zkzp); // multivariate dVi = dZKZ'/ds
        if(ai && cycle > 2){
          ww[i] = - (0.5 * arma::as_scalar(trace(PdVi))) + (0.5 * arma::as_scalar((Ysm.t() * PdVi * P * Ysm)));
        }else{
          ww[i] = arma::as_scalar(Ysm.t() * PdVi * P * Ysm) - accu(diagvec(PdVi));
        }
        PdViList.slice(i) = PdVi;
      }
      // theta(k) * dL/ds  ..... are scalar values
      ww = ww % var_components; // to be used later for updating the variance components
      
      // calculate second derivatives (AverageInformation)
      // Fisher's Information tr(PVi * PVi) .... A*=Vi=dV/ds .... [Vi Vj'] si sj ; TT is the list of derivatives for all random effects - trait combos
      for (int i = 0; i < kk; i++){
        for (int j = 0; j < kk; j++){
          if (i > j){}else{//only upper triangular
            if(ai && cycle > 2){ // if average informatio
              A(i,j) = 0.5 * arma::as_scalar(Ysm.t() * PdViList.slice(i) * P * PdViList.slice(j) * P * (P * Ysm)); // j is .t() ?
            }else{ // if newton raphson
              A(i,j) = accu(PdViList.slice(i) % PdViList.slice(j).t()) * arma::as_scalar(var_components(i)) * arma::as_scalar(var_components(j));
            }
          }
        }
      }
      A = arma::symmatu(A); // copy lower in upper triangular
      A_svd = arma::pinv(A, 1.490116e-08); // Inverse of Fishers
      
      if(A_svd.n_rows == 0){ // if fails 
        Rcpp::Rcout << "System is singular (A_svd). Aborting the job. Try a bigger number of tolparinv." << arma::endl;
        return 0;
      }
      
      // F- * sigma(k) * dL/ds
      arma::vec new_ww(kk);
      new_ww = A_svd * ww; //update variance components where: ww = theta(k) * dL/ds
      
      // ^^^^^^^^^^^^^^^^^^
      // ^^^^^^^^^^^^^^^^^^
      // parameter restrain
      // GeI values
      // 0 not estimated
      // 1 positive
      // 2 flexible
      // 3 fixed
      arma::vec coef_ut_unC = coef_ut_un + (taper(cycle) * new_ww);
      arma::uvec restrain = find(constraints == 1 && coef_ut_unC < 0);
      arma::vec cc = coef_ut_unC(restrain); // extract the ones that suppose to be positive
      // arma::vec cc2 = cc(find(cc < 0)); // identify var comp < 0 (1's)
      if(cc.n_elem > 0){ // we have to restrain
        // rest0 = '(';  rest1=cc.n_elem; rest2 = 'restrained)';
        arma::uvec no_restrain = find((constraints == 1 && coef_ut_unC > 0) || (constraints > 1)); // indices of columns that are OK to use
        arma::mat Ac = A.submat(no_restrain,no_restrain); // subset of A
        arma::mat Ac_svd = arma::inv(Ac); // Inverse of Fishers (subset of A)
        if(Ac_svd.n_rows == 0){ // if fails 
          Rcpp::Rcout << "System is singular (Ac_svd). Stopping the job. Try a bigger number of tolparinv." << arma::endl;
          return 0;
        }
        arma::vec wwc = ww(no_restrain); // subset of ww
        arma::vec newc_wwc = Ac_svd * wwc; //update variance components
        new_ww(no_restrain) = newc_wwc;
        new_ww(restrain) = new_ww(restrain)*0;
      }//else just keep going
      // end of parameter restrain
      // ^^^^^^^^^^^^^^^^^^
      // ^^^^^^^^^^^^^^^^^^
      
      // update
      // sigma + f[s*F-*dL/ds] ..... = coef + taper[x]
      coef_ut_un = coef_ut_un + (taper(cycle) * new_ww); 
      if(cc.n_elem > 0){
        coef_ut_un(restrain) = coef_ut_un(restrain)*0; // the ones that still go below zero and shouldn't let's fix them
      }
      // weight the projection matrix to provide stability
      sigmatwo(arma::find(pos == 0)) =  coef_ut_un(arma::find(pos == 0)); // index of pos
      sigmatwo(arma::find(pos == 1)) = exp(coef_ut_un(arma::find(pos == 1)));
      // the fixed paramters are forced to be the original value
      sigmatwo(find(constraints == 3)) = sigmaF_ut_un(find(constraints == 3));
      // the equalized paramters are forced to be equal
      // arma::uvec ko = find(constraints == 4);
      // if(ko.n_elem > 0){
      //   
      //   double meansigman = mean(sigmatwo(find(constraints == 4)));
      //   sigmatwo(find(constraints == 4)) = (sigmatwo(find(constraints == 4))*0) + meansigman;
      //   Rcpp::Rcout << meansigman << arma::endl;
      //   Rcpp::Rcout << sigmatwo(find(constraints == 4)) << arma::endl;
      //   // sigmatwo(find(constraints == 4))[1] = sigmatwo(find(constraints == 4))[1];
      // }
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
        arma::mat  FI = A/2;
        arma::vec myone(pos.n_elem,arma::fill::ones);
        arma::vec sp = ((sigmatwo - myone) % pos) + myone;
        arma::mat FI_c = FI / (sp * sp.t());
        sigma_cov = pinv(FI_c);
        if(sigma_cov.n_rows == 0){ // if fails 
          Rcpp::Rcout << "System is singular (sigma_cov). Aborting the job." << arma::endl;
          return 0;
        }
      }
      
    }else{// if we are in the last iteration now we calculate u, PEV, B, XB
      
      arma::inv(tXVXi,tXVX);
      if(tXVXi.n_rows == 0){ // if fails try to invert with diag(1e-6)
        arma::inv(tXVXi,tXVX+(D*(tolparinv)));
        if(tXVXi.n_rows == 0){// if fails try to invert with diag(1e-5)
          arma::inv(tXVXi,tXVX+(D*(tolparinv*10)));
          if(tXVXi.n_rows == 0){
            Rcpp::Rcout << "System is singular (tXVXi). Aborting the job. Try a bigger number of tolparinv." << arma::endl;
            return 0;
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
          arma::mat VarK = arma::kron(Rcpp::as<arma::mat>(K[i]),sigma.slice(i));
          arma::mat ZKfv = VarK * arma::kron(Zprov.t(),dD);
          U(i) = ZKfv * Vie; // Z' G Vi (Y - Xb)
          if(pev==true){
            VarU(i) = ZKfv * (P * ZKfv.t()); // var(u) = Z' G [Vi - (VX*tXVXVX)] G Z'
            PevU(i) = VarK - Rcpp::as<arma::mat>(VarU(i)); // PEV = G - var(u)
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
  arma::mat FISH = (sigma_cov % (dd*dd.t())) / (ee*ee.t());
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
    Rcpp::Named("monitor") = monitor2
  );
}
