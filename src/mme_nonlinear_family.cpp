#include <RcppArmadillo.h>
#include <RcppEigen.h>
#include "mme_nonlinear_family.h"

// [[Rcpp::depends(RcppArmadillo, RcppEigen)]]

// [[Rcpp::export]]
double sommer_nonlinear_sparse_admission(const arma::sp_mat & design,
    const arma::sp_mat & precision,const arma::sp_mat & curvature,
    const double memoryMB=1000){
  return sommer::nonlinearSparseAdmission(design,precision,curvature,memoryMB);
}

// [[Rcpp::export]]
arma::sp_mat sommer_nonlinear_coordinate_map(const arma::sp_mat & random,
    const Rcpp::List & descriptors,const Rcpp::List & relationships,
    const double memoryMB=1000){
  sommer::NonlinearCovariance covariance(random,descriptors,relationships,memoryMB);
  std::vector<double> free;
  for(int term=0;term<descriptors.size();++term){
    const Rcpp::List descriptor=descriptors[term];
    const arma::vec parameters=Rcpp::as<arma::vec>(descriptor["par"]);
    const Rcpp::LogicalVector flags=descriptor["free"];
    for(int j=0;j<flags.size();++j) if(flags[j]) free.push_back(parameters(j));
  }
  std::string reason;
  if(!covariance.update(arma::vec(free),reason))
    Rcpp::stop("Invalid post-fit covariance: "+reason);
  const arma::uword q=covariance.latentStart(random.n_cols);
  if(24.0*random.n_cols*q>memoryMB*1048576)
    Rcpp::stop("Nonlinear coordinate-map workspace exceeds execution$denseMemoryMB.");
  arma::sp_mat map(random.n_cols,q);
  for(arma::uword j=0;j<q;++j){
    Rcpp::checkUserInterrupt();
    arma::vec unit(q,arma::fill::zeros); unit(j)=1;
    const arma::vec column=covariance.restore(unit);
    for(arma::uword i=0;i<column.n_elem;++i)
      if(column(i)!=0) map(i,j)=column(i);
  }
  return map;
}

// [[Rcpp::export]]
SEXP sommer_nonlinear_family_create(const arma::sp_mat & fixed,
    const arma::sp_mat & random,const arma::vec & response,
    const arma::vec & offset,const arma::vec & input,const std::string & curve,
    const int coordinate,const int curveStart,const Rcpp::List & descriptors,
    const Rcpp::List & relationships,const Rcpp::List & specifications,
    const arma::uvec & index,const arma::vec & trials,
    const double memoryMB=1000,const int maxIters=100,
    const double scoreTol=1e-8,const double stepTol=1e-8,
    Rcpp::Nullable<Rcpp::NumericMatrix> coordinateMap=R_NilValue){
  if(coordinate<0 || curveStart<0 ||
     (coordinate==0 && static_cast<arma::uword>(curveStart)!=random.n_cols))
    Rcpp::stop("Invalid nonlinear family coordinate map.");
  if(!std::isfinite(memoryMB) || memoryMB<=0 ||
     8.0*(150.0*response.n_elem+60.0*random.n_cols+8.0*random.n_nonzero)>
       memoryMB*1048576)
    Rcpp::stop("Nonlinear family workspace exceeds execution$denseMemoryMB.");
  return Rcpp::XPtr<sommer::NonlinearFamilyObjective>(
    new sommer::NonlinearFamilyObjective(fixed,random,response,offset,input,curve,
      coordinate,curveStart,descriptors,relationships,specifications,index,trials,
      memoryMB,maxIters,scoreTol,stepTol,coordinateMap.isNull() ? arma::mat() :
        Rcpp::as<arma::mat>(coordinateMap.get())),true);
}

// [[Rcpp::export]]
arma::vec sommer_nonlinear_curve_joint_eval(const std::string & curve,
    const arma::vec & input,const arma::vec & parameters,
    const arma::mat & deviations){
  if(deviations.n_rows!=input.n_elem || deviations.n_cols!=parameters.n_elem ||
     parameters.n_elem!=sommer::nonlinearParameterCount(curve) ||
     !deviations.is_finite())
    Rcpp::stop("Invalid joint curve prediction dimensions or deviations.");
  arma::vec result(input.n_elem);
  for(arma::uword i=0;i<input.n_elem;++i){
    if(i%1024==0) Rcpp::checkUserInterrupt();
    sommer::CurveDerivatives value;
    if(!sommer::nonlinearCurveDerivatives(curve,input(i),
        parameters+deviations.row(i).t(),value))
      Rcpp::stop("Non-finite joint curve prediction.");
    result(i)=value.mean;
  }
  return result;
}

// [[Rcpp::export]]
Rcpp::List sommer_nonlinear_family_eval(SEXP context,const arma::vec & parameters,
    const bool gradient=true){
  Rcpp::XPtr<sommer::NonlinearFamilyObjective> objective(context);
  return objective->evaluate(parameters,gradient);
}

// [[Rcpp::export]]
Rcpp::List sommer_nonlinear_family_conditional(SEXP context,const arma::vec & parameters,
    const arma::vec & mode){
  Rcpp::XPtr<sommer::NonlinearFamilyObjective> objective(context);
  const auto value=objective->conditional(parameters,mode);
  return Rcpp::List::create(Rcpp::_["valid"]=value.valid,Rcpp::_["reason"]=value.reason,
    Rcpp::_["objective"]=value.valid ? value.objective : R_PosInf,
    Rcpp::_["score"]=value.score,Rcpp::_["hessian"]=value.hessian,
    Rcpp::_["predictor"]=value.mean);
}
