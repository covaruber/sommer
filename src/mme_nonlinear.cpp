#include <RcppArmadillo.h>
#include "mme_nonlinear_gaussian.h"
#include "mme_nonlinear_laplace.h"

// [[Rcpp::export]]
Rcpp::List sommer_nonlinear_curve_eval(const std::string & curve,
    const arma::vec & input, const arma::vec & parameters, const bool gradient=true){
  sommer::NonlinearCurveEvaluation result;
  if(!sommer::nonlinearCurve(curve,input,parameters,result,gradient))
    Rcpp::stop("Nonlinear curve evaluation produced non-finite values or derivatives.");
  return Rcpp::List::create(Rcpp::_["mean"]=result.mean,
    Rcpp::_["jacobian"]=result.jacobian);
}

// [[Rcpp::export]]
SEXP sommer_nonlinear_gaussian_create(const arma::sp_mat & fixed,
    const arma::sp_mat & random, const arma::vec & response,
    const arma::vec & offset, const arma::vec & input, const std::string & curve,
    const Rcpp::List & descriptors, const Rcpp::List & relationships,
    const double memoryMB=1000){
  return Rcpp::XPtr<sommer::NonlinearGaussianObjective>(
    new sommer::NonlinearGaussianObjective(fixed,random,response,offset,input,
      curve,descriptors,relationships,memoryMB), true);
}

// [[Rcpp::export]]
Rcpp::List sommer_nonlinear_gaussian_eval(SEXP context,
    const arma::vec & parameters, const bool gradient=true){
  Rcpp::XPtr<sommer::NonlinearGaussianObjective> objective(context);
  return objective->evaluate(parameters,gradient);
}

// [[Rcpp::export]]
Rcpp::List sommer_nonlinear_curve_derivatives(const std::string & curve,
    const double input, const arma::vec & parameters){
  sommer::CurveDerivatives value;
  if(!sommer::nonlinearCurveDerivatives(curve,input,parameters,value))
    Rcpp::stop("Invalid nonlinear derivative inputs or non-finite result.");
  const arma::uword k=parameters.n_elem;
  arma::vec first(k); arma::mat second(k,k); arma::cube third(k,k,k);
  for(arma::uword i=0;i<k;++i){
    first(i)=value.d[i];
    for(arma::uword j=0;j<k;++j){
      second(i,j)=value.d2[i][j];
      for(arma::uword l=0;l<k;++l) third(i,j,l)=value.d3[i][j][l];
    }
  }
  return Rcpp::List::create(Rcpp::_["mean"]=value.mean,Rcpp::_["first"]=first,
    Rcpp::_["second"]=second,Rcpp::_["third"]=third);
}

// [[Rcpp::export]]
SEXP sommer_nonlinear_laplace_create(const arma::sp_mat & fixed,
    const arma::sp_mat & random, const arma::vec & response,
    const arma::vec & offset, const arma::vec & input, const std::string & curve,
    const int coordinate, const int curveStart,
    const Rcpp::List & descriptors, const Rcpp::List & relationships,
    const double memoryMB=1000, const int maxIters=100,
    const double scoreTol=1e-8, const double stepTol=1e-8){
  if(coordinate<1 || curveStart<0) Rcpp::stop("Invalid random curve coordinate map.");
  return Rcpp::XPtr<sommer::NonlinearLaplaceObjective>(
    new sommer::NonlinearLaplaceObjective(fixed,random,response,offset,input,
      curve,coordinate-1,curveStart,descriptors,relationships,memoryMB,
      maxIters,scoreTol,stepTol),true);
}

// [[Rcpp::export]]
Rcpp::List sommer_nonlinear_laplace_eval(SEXP context,
    const arma::vec & parameters, const bool gradient=true){
  Rcpp::XPtr<sommer::NonlinearLaplaceObjective> objective(context);
  return objective->evaluate(parameters,gradient);
}

// [[Rcpp::export]]
Rcpp::List sommer_nonlinear_conditional_eval(SEXP context,
    const arma::vec & parameters, const arma::vec & mode){
  Rcpp::XPtr<sommer::NonlinearLaplaceObjective> objective(context);
  const auto value=objective->conditional(parameters,mode);
  return Rcpp::List::create(Rcpp::_["valid"]=value.valid,
    Rcpp::_["reason"]=value.reason,Rcpp::_["objective"]=value.valid ? value.objective : R_PosInf,
    Rcpp::_["score"]=value.score,Rcpp::_["hessian"]=value.hessian,
    Rcpp::_["mean"]=value.mean);
}

// [[Rcpp::export]]
arma::vec sommer_nonlinear_curve_random_eval(const std::string & curve,
    const arma::vec & input, const arma::vec & parameters, const int coordinate,
    const arma::vec & deviations){
  if(coordinate<1 || static_cast<arma::uword>(coordinate)>parameters.n_elem ||
     parameters.n_elem!=sommer::nonlinearParameterCount(curve) ||
     deviations.n_elem!=input.n_elem || !deviations.is_finite())
    Rcpp::stop("Invalid random curve prediction dimensions or deviations.");
  arma::vec result(input.n_elem);
  for(arma::uword i=0;i<input.n_elem;++i){
    if(i%1024==0) Rcpp::checkUserInterrupt();
    arma::vec theta=parameters;
    theta(coordinate-1)+=deviations(i);
    sommer::CurveDerivatives value;
    if(!sommer::nonlinearCurveDerivatives(curve,input(i),theta,value))
      Rcpp::stop("Non-finite random curve prediction.");
    result(i)=value.mean;
  }
  return result;
}
