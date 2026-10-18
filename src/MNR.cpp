#include "mme_sp2_support.h"
#include "mme_progress.h"
#include <Matrix/stubs.c>

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
  
  if(nishio == true){ // Nishio and Satoh (2014); Vitezica et al. (2013) orthogonal coding
    
    arma::rowvec p = mean( X+1, 0 )/2; // frequency of the allele counted by x=1
    arma::rowvec q = 1-p;
    arma::rowvec p2q = 2*(p%q);
    // w = -2q^2 (x=1), 2pq (x=0), -2p^2 (x=-1), i.e. w = -x^2 + (p-q)x + 2pq
    arma::mat W = X.each_row() % (p-q) - arma::square(X);
    W.each_row() += p2q;
    double varHW = arma::accu(arma::square(p2q));
    D = (W * W.t())/varHW;
    
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

arma::mat predict_mmes_vcov_cpp(Rcpp::List model, const arma::sp_mat & Dmat);

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

  if(model.containsElementNamed("coefficientMap") &&
     !Rf_isNull(model["coefficientMap"])){
    const arma::sp_mat mapping = Rcpp::as<arma::sp_mat>(model["coefficientMap"]);
    const arma::uword originalEffects = mapping.n_rows;
    arma::sp_mat Ci(originalEffects, originalEffects);
    if(mode == 2){
      const arma::sp_mat identity = arma::speye<arma::sp_mat>(originalEffects, originalEffects);
      Ci = arma::sp_mat(predict_mmes_vcov_cpp(model, identity));
    }else{
      for(arma::uword first = 0; first < originalEffects; first += 32){
        const arma::uword count = std::min<arma::uword>(32, originalEffects - first);
        arma::sp_mat selection(count, originalEffects);
        for(arma::uword row = 0; row < count; ++row){ selection(row, first + row) = 1.0; }
        const arma::mat block = predict_mmes_vcov_cpp(model, selection);
        for(arma::uword row = 0; row < count; ++row){ Ci(first + row, first + row) = block(row, row); }
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
  arma::sp_mat design = Dmat;
  if(model.containsElementNamed("coefficientMap") &&
     !Rf_isNull(model["coefficientMap"])){
    const arma::sp_mat mapping = Rcpp::as<arma::sp_mat>(model["coefficientMap"]);
    if(design.n_cols != mapping.n_rows){
      Rcpp::stop("D must have one column per original mixed-model coefficient.");
    }
    design = design * mapping;
  }

  if(C.n_rows != C.n_cols){
    Rcpp::stop("Stored C must be square.");
  }

  const int nEffects = static_cast<int>(C.n_rows);

  if(static_cast<int>(design.n_cols) != nEffects){
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

  const int k = static_cast<int>(design.n_rows);

  Eigen::MatrixXd rhs = Eigen::MatrixXd::Zero(nEffects, k);
  for(arma::sp_mat::const_iterator it = design.begin(); it != design.end(); ++it){
    rhs(static_cast<Eigen::Index>(it.col()), static_cast<Eigen::Index>(it.row())) = (*it);
  }

  Eigen::MatrixXd X = Cfactor.solve(rhs);

  if(Cfactor.info() != Eigen::Success){
    Rcpp::stop("Sparse LDLT solve failed while computing prediction variances.");
  }

  const arma::mat Xarma(X.data(), static_cast<arma::uword>(nEffects), static_cast<arma::uword>(k));

  arma::mat vcov = (design * Xarma) * Cscale;
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
              Inf(i,j) = 0.5 * arma::as_scalar(Ysm.t() * PdViList.slice(i) * PdViList.slice(j) * Py);
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

  if(op == "fa" || op == "rr"){
    return sommer::loadingFactorCovariance(f, localPar);
  }

  if(op == "ante"){
    const int ncoef = Rcpp::as<int>(f["ante_ncoef"]);
    Rcpp::IntegerVector rr = f["ante_row"];
    Rcpp::IntegerVector cc = f["ante_col"];
    if(ncoef < 1 || rr.size() != ncoef || cc.size() != ncoef ||
       localPar.n_elem != static_cast<arma::uword>(ncoef + static_cast<int>(q) - 1)){
      Rcpp::stop("Malformed antedependence covariance factor.");
    }
    return anteInnovationCovariance(q, rr, cc, localPar, sparseFactorKernelsEnabled(f));
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

  if(sparseFactorKernelsEnabled(f) && f.containsElementNamed("model") &&
     (Rcpp::as<std::string>(f["model"]) == "sar" || Rcpp::as<std::string>(f["model"]) == "car")){
    K = evaluateSparseSpatialFactor(f, localPar, false).covariance;
  }else if(backend == "native"){
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
  if(sparseFactorKernelsEnabled(f) && f.containsElementNamed("model") &&
     (Rcpp::as<std::string>(f["model"]) == "sar" || Rcpp::as<std::string>(f["model"]) == "car")){
    return evaluateSparseSpatialFactor(f, localPar, true).covarianceD1;
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

// Residual structures with correlated or section-specific factors are assembled
// directly in record space instead of over the full covariance-product space.
static bool directResidualCompact(const Rcpp::List & cs){
  Rcpp::List factors = cs["factors"];
  for(int fidx = 0; fidx < factors.size(); ++fidx){
    Rcpp::List f = Rcpp::as<Rcpp::List>(factors[fidx]);
    if(f.containsElementNamed("section_levels") ||
       !f.containsElementNamed("structurally_diagonal") ||
       !Rcpp::as<bool>(f["structurally_diagonal"])){
      return true;
    }
  }
  return false;
}

static int directSectionListed(const Rcpp::List & f, const int t){
  if(f.containsElementNamed("section_index")){
    Rcpp::IntegerVector index = f["section_index"];
    return index[t];
  }
  return t;
}

static arma::vec directSectionSlice(const Rcpp::List & f, const arma::vec & localPar, const int listed){
  const arma::uword npar = static_cast<arma::uword>(Rcpp::as<int>(f["section_npar"]));
  if(npar == 0 || listed < 0){
    return arma::vec();
  }
  const arma::uword first = static_cast<arma::uword>(listed) * npar;
  return localPar.subvec(first, first + npar - 1);
}

// Kronecker PD: sigma2 * prod_f min eig(F_f), taking the worst section per factor.
static bool directDescriptorPD(const Rcpp::List & cs, const arma::vec & par, const double floor){
  if(par.n_elem < 1 || !par.is_finite()){
    return false;
  }
  Rcpp::List factors = cs["factors"];
  double logMinimum = par(0);
  for(int fidx = 0; fidx < factors.size(); ++fidx){
    Rcpp::List f = Rcpp::as<Rcpp::List>(factors[fidx]);
    const int start1 = f.containsElementNamed("par_start") ? Rcpp::as<int>(f["par_start"]) : 1;
    const int end1 = f.containsElementNamed("par_end") ? Rcpp::as<int>(f["par_end"]) : 0;
    arma::vec localPar;
    if(end1 >= start1){
      localPar = par.subvec(static_cast<arma::uword>(start1-1), static_cast<arma::uword>(end1-1));
    }
    const int nSections = f.containsElementNamed("section_levels") ? Rf_length(f["section_levels"]) : 0;
    double factorMinimum = 1.0;
    bool any = false;
    for(int t = 0; t < std::max(nSections, 1); ++t){
      arma::mat K;
      try{
        if(nSections > 0){
          const int listed = directSectionListed(f, t);
          if(listed < 0){
            continue;
          }
          K = directEvalFactor(f, directSectionSlice(f, localPar, listed));
        }else{
          K = directEvalFactor(f, localPar);
        }
      }catch(...){
        return false;
      }
      arma::vec eigenvalues;
      if(!arma::eig_sym(eigenvalues, K) || !eigenvalues.is_finite() || eigenvalues.min() <= 0.0){
        return false;
      }
      factorMinimum = any ? std::min(factorMinimum, eigenvalues.min()) : eigenvalues.min();
      any = true;
    }
    logMinimum += std::log(factorMinimum);
  }
  return std::isfinite(logMinimum) && logMinimum > std::log(floor);
}

static arma::mat directEvaluateDescriptor(const Rcpp::List & cs, const arma::vec & par){

  if(par.n_elem < 1){
    Rcpp::stop("Covariance descriptor must contain log_sigma2.");
  }
  Rcpp::List factors = cs["factors"];
  arma::mat K(1,1,arma::fill::ones);

  for(int fidx = 0; fidx < factors.size(); ++fidx){
    Rcpp::List f = Rcpp::as<Rcpp::List>(factors[fidx]);
    if(f.containsElementNamed("section_levels")){
      Rcpp::stop("Internal error: dsumm() structures must be assembled in record space.");
    }
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

// [[Rcpp::export]]
arma::mat evaluate_covstruct_cpp(const Rcpp::List & covStruct, const arma::vec & par){
  return directEvaluateDescriptor(covStruct, par);
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
    covLower(i) = arma::vec(covPar(i).n_elem, arma::fill::value(-std::numeric_limits<double>::infinity()));
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
  const bool adaptive = covStructI.hasAttribute("adaptiveRecovery") &&
    Rcpp::as<bool>(covStructI.attr("adaptiveRecovery"));
  const std::string blendMethod = covStructI.hasAttribute("blend")
    ? Rcpp::as<std::string>(covStructI.attr("blend")) : "em";
  if(blendMethod != "em" && blendMethod != "pxem") {
    Rcpp::stop("Unknown direct-engine EM blending method.");
  }
  const bool parameterExpanded = blendMethod == "pxem";
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
    arma::mat feasible = scoreTangent;
    for(int activeSet = 0; weight < 1.0 && activeSet <= nVcTotal && feasible.n_cols > 0; ++activeSet){
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
      if(adaptive && (!solved || !direction.is_finite() ||
          arma::norm(normalized * direction - rhs, 2) >
            1.0e-6 * (1.0 + arma::norm(rhs, 2)) ||
          !sommer::reliableAIInformation(
          information, 1.0 / scaling))) {
        ++regularizedAttempts;
        const sommer::RegularizedAIDirection retry =
          sommer::regularizedAIDirection(AI, score, feasible, stepCaps);
        if(retry.reliable) {
          delta = retry.direction;
          regularized = true;
          regularizedRidges.push_back(retry.ridge);
          regularizedIterations.push_back(iIter + 1);
          break;
        }
        weight = 1.0;
        recoveryReason = "information-solve";
        break;
      }
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
          if(!fullEM && delta(k) != 0.0) {
            alpha = std::min(alpha, stepCaps(k) / std::abs(delta(k)));
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
