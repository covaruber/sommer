#ifndef SOMMER_MME_FACTORS_H
#define SOMMER_MME_FACTORS_H

#include <RcppArmadillo.h>
#include "mme_paths.h"

namespace sommer {

inline PrecisionFamily precisionFamily(const Rcpp::List & factor) {
  if(!factor.containsElementNamed("representation")) { return PrecisionFamily::Dense; }
  const Rcpp::List representation = factor["representation"];
  const std::string kind = representation.containsElementNamed("kind")
    ? Rcpp::as<std::string>(representation["kind"]) : "dense";
  const std::string form = representation.containsElementNamed("form")
    ? Rcpp::as<std::string>(representation["form"]) : "";
  if(kind == "diagonal") { return PrecisionFamily::Diagonal; }
  if(kind == "lowrank_diagonal") { return PrecisionFamily::LowRankDiagonal; }
  if(kind == "compound_symmetry") { return PrecisionFamily::CompoundSymmetry; }
  if(form == "ar1") { return PrecisionFamily::AR1; }
  if(form == "stationary_ar") { return PrecisionFamily::StationaryAR; }
  if(form == "modified_cholesky") { return PrecisionFamily::ModifiedCholesky; }
  if(form == "normalized_normal_equations") { return PrecisionFamily::NormalEquations; }
  if(form == "normalized_affine_precision") { return PrecisionFamily::Affine; }
  return PrecisionFamily::Dense;
}

inline Rcpp::List factorCapability(const Rcpp::List & factor, const char * name) {
  if(!factor.containsElementNamed("representation")) {
    Rcpp::stop("Covariance factor is missing its representation.");
  }
  const Rcpp::List representation = factor["representation"];
  if(!representation.containsElementNamed(name) || Rf_isNull(representation[name])) {
    Rcpp::stop(std::string("Covariance factor is missing its ") + name + " capability.");
  }
  return Rcpp::as<Rcpp::List>(representation[name]);
}

struct FactorPrecisionState {
  arma::mat precision;
  std::vector<arma::mat> precisionD1;
  std::vector<double> logDetD1;
  std::vector<arma::mat> factorCovariance;
  std::vector<arma::mat> factorPrecision;
  std::vector<arma::sp_mat> sparseFactorPrecision;
  std::vector<Rcpp::List> factorDescriptor;
  std::vector<double> factorLogDet;
  std::vector<arma::vec> factorPar;
  std::vector<int> factorStart1;
  std::vector<int> factorEnd1;
  std::vector<arma::mat> localCovarianceD1;
  std::vector<arma::mat> localPrecisionD1;
  std::vector<arma::sp_mat> sparseLocalPrecisionD1;
  std::vector<int> derivativeFactor;
  arma::mat lowRankLoadings;
  arma::vec lowRankSpecific;
  arma::vec lowRankPrecisionDiagonal;
  int lowRankFactorIndex = -1;
  double lowRankCoreLogDet = 0.0;
  double inverseScale = 1.0;
  double logDet = 0.0;
  arma::uword totalDim = 0;
  bool factorWise = false;
};

}

#endif