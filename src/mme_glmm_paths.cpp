#include <Rcpp.h>
#include "mme_glmm_paths.h"

// [[Rcpp::export(rng = false)]]
Rcpp::List sommer_glmm_path_plan(const std::string & estimator,
    const bool maximumLikelihood, const bool ordinal, const bool nominal,
    const bool mixture, const bool independentResiduals){
  if(estimator != "pql" && estimator != "laplace")
    Rcpp::stop("Unknown GLMM estimator.");
  const auto plan = sommer::planGLMMPath({
    estimator == "pql" ? sommer::GLMMEstimator::PQL : sommer::GLMMEstimator::Laplace,
    maximumLikelihood, ordinal, nominal, mixture, independentResiduals});
  return Rcpp::List::create(
    Rcpp::_["requested"]=sommer::glmmPathName(plan.requested),
    Rcpp::_["eligible"]=plan.eligible,
    Rcpp::_["selected"]=sommer::glmmPathName(plan.selected),
    Rcpp::_["reason"]=sommer::glmmPathReasonCode(plan.reason),
    Rcpp::_["message"]=sommer::glmmPathReasonMessage(plan.reason),
    Rcpp::_["likelihood"]=nominal ? "category-block" : "scalar",
    Rcpp::_["composition"]=nominal ? (ordinal ? "ordinal-nominal" : "nominal") :
      (ordinal ? "ordinal" : "scalar"));
}
