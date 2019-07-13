// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#define ARMA_DONT_PRINT_ERRORS
#include "RcppArmadillo.h"

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
          // if(i == j){}else{out2[counter]=false;}
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
  arma::vec Ym = vectorise(Y); // multivariate Y
  int nom = Ym.n_rows;
  
  arma::mat Xm;
  for (int i = 0; i < n_fixed; ++i) {
    if(i==0){
      Xm = kron(Rcpp::as<arma::mat>(Gx[i]), Rcpp::as<arma::mat>(X[i]));
    }else{
      Xm = arma::join_horiz( Xm , kron(Rcpp::as<arma::mat>(Gx[i]), Rcpp::as<arma::mat>(X[i])) );
    }
  }
  // Xm = arma::kron(Gx,X);
  
  // Xm = Xm.t();
  // arma::mat Xm = kron(Gx,X); // multivariate X
  arma::mat Ys = scaleCpp(Y); // scale Y using the scaleCpp function made
  arma::vec Ysm = vectorise(Ys); // multivariate Y
  // ****************************************************
  // initial VC
  // ****************************************************
  arma::mat base_var = cov(Y); // original variance
  arma::mat sc_var = cov(Ys); // scaled variance
  int rankX = Xm.n_rows - rank(Xm); // rank arma::uword
  // sigma
  arma::cube sigma(n_traits,n_traits,n_re);
  arma::cube sigma_scaled(n_traits,n_traits,n_re);
  
  arma::field<arma::vec> sigma_ut(n_re);
  arma::field<arma::vec> constraintsL(n_re);
  int no_vc = 0;
  for (int i = 0; i < n_re; ++i) { // scale the provided Ge
    // arma::mat provj = as<arma::mat>(Ge[i])/base_var;
    sigma.slice(i) = Rcpp::as<arma::mat>(Ge[i]); //provj % sc_var ; // scaled variance
    // sigmaddd[i] = provj % sc_var ; // scaled variance
    arma::vec oo = mat_to_vecCpp(sigma.slice(i),GeI[i]) ; // extract upper triangular in a vector form
    sigma_ut[i] = oo;
    constraintsL[i] = mat_to_vecCpp(GeI[i],GeI[i]) ; // who are diagonal and non-diagonal VCs
    no_vc = no_vc + oo.n_elem;
  }
  arma::vec sigma_ut_un;
  arma::vec constraints; 
  for(int i=0; i < n_re ; i++){
    sigma_ut_un = join_cols(sigma_ut_un,sigma_ut[i]);
    constraints = join_cols(constraints,constraintsL[i]);
  }
  arma::vec sigmaF_ut_un = sigma_ut_un; // make a copy for fixed-value vc's
  arma::vec coef_ut_un = sigma_ut_un; // make a 2nd copy
  
  arma::vec taper(iters);  taper.fill(0.9); 
  
  int  kk = sigma_ut_un.n_elem; // number of VCs
  arma::vec llstore(iters); // container for LL
  // std::vector<bool> pos(sigma_ut_un.n_elem, true); // create position vector
  arma::vec pos(sigma_ut_un.n_elem, arma::fill::zeros); // create position vector
  
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
  arma::mat P(nom,nom);
  arma::sp_mat D = arma::speye<arma::sp_mat>(nom,nom);
  arma::vec seqrankX = seqCpp(0,rankX-1); // will be used to keep only the eigen values for 1:rankX
  arma::vec seqkk = seqCpp(0,kk-1);
  arma::vec popo = arma::vec(rankX, arma::fill::zeros);
  for(int i=0; i < rankX; i++){popo(i) = 1;}
  arma::mat A(kk,kk,arma::fill::zeros);
  arma::mat A_svd;
  arma::vec eigval2; // will be used for the decomposition of P, within the algorithm
  arma::mat eigvec2; // will be used for the decomposition of P
  arma::mat sigma_store(sigma_ut_un.n_elem,iters);
  arma::mat llik_store(1,iters);
  
  arma::mat beta, fitted, residuals;
  // arma::field<arma::mat> VarG(n_random);
  Rcpp::List VarU(n_random);
  Rcpp::List PevU(n_random);
  Rcpp::List U(n_random);
  // arma::sp_mat dD = arma::speye<arma::sp_mat>(n_traits,n_traits);
  
  arma::vec vdD(n_traits,arma::fill::ones);
  arma::mat dD = arma::diagmat(vdD);
  arma::mat sigma_cov;
  arma::mat tXVXi; // var-cov fixed effects
  // arma::mat dD = arma::speye<arma::mat>(n_traits,n_traits);
  
  bool convergence = false;
  bool last_iteration = false;
  int cycle, cycle2;
  double ldet, llik, llik0, delta_llik, checkP, seconds; // to store likelihoods and determinants
  
  for(cycle=0; cycle < iters; cycle++){ // for each cycle
    
    if(cycle == 0){taper(cycle) = 0.5;}
    if(cycle == 1){taper(cycle) = 0.7;}
    
    for (int i = 0; i < n_re; ++i) { 
      sigma_ut[i] = mat_to_vecCpp(sigma.slice(i),Rcpp::as<arma::mat>(GeI[i])) ; // extract upper triangular in a vector form
    }
    arma::vec sigmatwo;
    for(int i=0; i < n_re ; i++){
      sigmatwo = join_cols(sigmatwo,sigma_ut[i]);
    }
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
    if(Vi.n_rows == 0){ // if fails try to invert with diag(1e-6)
      V = V + (D*tolparinv);
      arma::inv_sympd(Vi,V);
      if(Vi.n_rows == 0){// if fails try to invert with diag(1e-5)
        V = V + (D*(tolparinv*10));
        arma::inv_sympd(Vi,V);
        if(Vi.n_rows == 0){ // if fails try to invert with diag(1e-4)
          V = V + (D*(tolparinv*100));
          arma::inv_sympd(Vi,V);
          if(Vi.n_rows == 0){ // finally, if fails try to invert with diag(1e-3)
            V = V + (D*(tolparinv*1000));
            arma::inv_sympd(Vi,V);
            if(Vi.n_rows == 0){ // finally stop
              Rcpp::Rcout << "Sistem is singular. Stopping the job. Try a smaller number of tolparinv." << arma::endl;
            }
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
    arma::mat VX = Vi * Xm;
    arma::mat tXVX = Xm.t() * VX;
    
    arma::mat tXVXVX;
    tXVXVX = arma::solve(tXVX, VX.t());
    arma::solve(tXVXVX,tXVX,VX.t());
    if(tXVXVX.n_rows == 0){ // if fails try to invert with diag(1e-6)
      arma::solve(tXVXVX,tXVX + (D*(tolparinv)),VX.t());
      if(tXVXVX.n_rows == 0){// if fails try to invert with diag(1e-5)
        arma::solve(tXVXVX,tXVX + (D*(tolparinv*10)),VX.t());
        if(tXVXVX.n_rows == 0){ // if fails try to invert with diag(1e-4)
          arma::solve(tXVXVX,tXVX + (D*(tolparinv*100)),VX.t());
          if(tXVXVX.n_rows == 0){ // finally stop
            Rcpp::Rcout << "Sistem is singular. Stopping the job. Try a smaller number of tolparinv." << arma::endl;
          }
        }
      }
    }
    // try{
    //   tXVXVX = arma::solve(tXVX, VX.t());//, arma::solve_opts::fast); 
    // }catch(std::exception &ex) {
    //   // tXVXVX <- try(solve((t(X)%*%VX + (tolparinv * diag(dim(t(X)%*%VX)[2]))),t(VX)), silent = TRUE)
    //   tXVXVX = arma::solve(tXVX + (D*tolparinv), VX.t());
    // }
    // projection matrix
    P = Vi - (VX*tXVXVX); // set V to NULL
    
    if(last_iteration == false){
      
      arma::vec rss = Ysm.t() * (P * Ysm); // scalar RSS
      
      double rankXorss = arma::as_scalar(rankX/rss);
      double rssorankX = arma::as_scalar(rss/rankX);
      
      sigmatwo = sigmatwo * rssorankX;
      
      // weight the projection matrix to provide stability
      coef_ut_un(arma::find(pos == 0)) =  sigmatwo(arma::find(pos == 0)); // index of pos
      coef_ut_un(arma::find(pos == 1)) = log(sigmatwo(arma::find(pos == 1)));
      
      // calculate the log-likelihood
      P = P * rankXorss; 
      rss = rankX; 
      arma::eig_sym(eigval2, eigvec2, P);
      eigval2 = sort(eigval2,"descend"); // sort eigen vectors
      eigval2 = eigval2(arma::find(popo == 1));//(find(seqrankX < rankX)); // only take the values from 1 to
      checkP = eigval2.min();
      if(checkP < 0){ // if any eigen value is < 0 recalculate P
        P = P + (D * (tolpar - eigval2.min())) ;
        eigval2 = eigval2 + tolpar - eigval2.min();
      }
      ldet = accu(log(eigval2));
      llik = ldet/2 - (arma::as_scalar(rss)/2);
      
      if(cycle == 0){llik0 = llik;}
      delta_llik = llik - llik0;
      llik0 = llik;
      
      // use the stabilization
      arma::vec var_components(kk, arma::fill::ones);
      double check00 = accu(pos);
      if(check00 > 0){
        arma::uvec ind = find(pos == 1);
        var_components(ind) = sigmatwo(ind);
      }
      
      // calculate first derivatives
      arma::vec ww(kk); // vector to store the product Y'PViPY - tr(PVi)
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
      ww = ww % var_components;
      // calculate second derivatives
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
      
      // F- * sigma(k) * dL/ds
      arma::vec new_ww(kk);
      new_ww = A_svd * ww; //update variance components
      
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
        arma::vec wwc = ww(no_restrain); // subset of ww
        arma::vec newc_wwc = Ac_svd * wwc; //update variance components
        new_ww(no_restrain) = newc_wwc;
        new_ww(restrain) = new_ww(restrain)*0;
      }//else just keep going
      // end of parameter restrain
      // ^^^^^^^^^^^^^^^^^^
      // ^^^^^^^^^^^^^^^^^^
      
      coef_ut_un = coef_ut_un + (taper(cycle) * new_ww); // sigma + f[s*F-*dL/ds] ..... = coef + taper[x]
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
      }
      
    }else{// if we are in the last iteration now we calculate u, PEV, B, XB
      
      arma::inv(tXVXi,tXVX);
      if(tXVXi.n_rows == 0){ // if fails try to invert with diag(1e-6)
        arma::inv(tXVXi,tXVX+(D*(tolparinv)));
        if(tXVXi.n_rows == 0){// if fails try to invert with diag(1e-5)
          arma::inv(tXVXi,tXVX+(D*(tolparinv*10)));
          if(tXVXi.n_rows == 0){
            Rcpp::Rcout << "Sistem is singular. Stopping the job. Try a smaller number of tolparinv." << arma::endl;
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
