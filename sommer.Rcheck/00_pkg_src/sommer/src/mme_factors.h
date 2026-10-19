#ifndef SOMMER_MME_FACTORS_H
#define SOMMER_MME_FACTORS_H

#include <RcppArmadillo.h>
#include "mme_paths.h"

namespace sommer {

inline bool zeroSpecificFactor(const Rcpp::List & factor) {
  if(factor.containsElementNamed("zero_specific") &&
     Rcpp::as<bool>(factor["zero_specific"])) { return true; }
  if(!factor.containsElementNamed("representation")) { return false; }
  const Rcpp::List representation = factor["representation"];
  if(!representation.containsElementNamed("loadings") ||
     Rf_isNull(representation["loadings"])) { return false; }
  const Rcpp::List loadings = representation["loadings"];
  return loadings.containsElementNamed("zero_specific") &&
    Rcpp::as<bool>(loadings["zero_specific"]);
}

inline bool heterogeneousSpecificFactor(const Rcpp::List & factor) {
  if(factor.containsElementNamed("representation")) {
    const Rcpp::List representation = factor["representation"];
    if(representation.containsElementNamed("loadings") &&
       !Rf_isNull(representation["loadings"])) {
      const Rcpp::List loadings = representation["loadings"];
      if(loadings.containsElementNamed("heterogeneous_specific")) {
        return Rcpp::as<bool>(loadings["heterogeneous_specific"]);
      }
    }
  }
  return !zeroSpecificFactor(factor);
}

inline arma::mat normalizedFactorLoadings(const Rcpp::List & factor,
                                         const arma::vec & par) {
  const bool rr = zeroSpecificFactor(factor);
  const bool rrMetadata = factor.containsElementNamed("rr_nload");
  const int q = Rcpp::as<int>(factor["dim"]);
  const int rank = Rcpp::as<int>(factor["order"]);
  const int count = Rcpp::as<int>(factor[rrMetadata ? "rr_nload" : "fa_nload"]);
  const Rcpp::IntegerVector rows = factor[rrMetadata ? "rr_row" : "fa_row"];
  const Rcpp::IntegerVector columns = factor[rrMetadata ? "rr_col" : "fa_col"];
  const Rcpp::LogicalVector diagonal = factor[rrMetadata ? "rr_diag" : "fa_diag"];
  const bool heterogeneous = !rr && heterogeneousSpecificFactor(factor);
  if(q < 1 || rank < 1 || rank > q || count < 1 ||
     rows.size() != count || columns.size() != count || diagonal.size() != count ||
     par.n_elem != static_cast<arma::uword>(count + (heterogeneous ? q - 1 : 0))) {
    Rcpp::stop("Malformed factor loading metadata.");
  }
  double reference = 0.0;
  for(int a = 0; a < count; ++a) {
    if(rows[a] == 1 && columns[a] == 1) { reference = par(a); }
  }
  const double twiceReference = 2.0 * reference;
  const double logScale = rr ? twiceReference :
    std::max(0.0, twiceReference) + std::log1p(std::exp(-std::abs(twiceReference)));
  arma::mat L(q, rank, arma::fill::zeros);
  for(int a = 0; a < count; ++a) {
    if(rows[a] < 1 || rows[a] > q || columns[a] < 1 || columns[a] > rank) {
      Rcpp::stop("Invalid factor loading coordinate.");
    }
    L(rows[a] - 1, columns[a] - 1) = diagonal[a]
      ? std::exp(par(a) - 0.5 * logScale)
      : par(a) * std::exp(-0.5 * logScale);
  }
  return L;
}

inline arma::mat loadingFactorCovariance(const Rcpp::List & factor,
                                        const arma::vec & par) {
  const arma::mat L = normalizedFactorLoadings(factor, par);
  arma::mat covariance = L * L.t();
  if(zeroSpecificFactor(factor)) { return covariance; }
  const bool rrMetadata = factor.containsElementNamed("rr_nload");
  const int count = Rcpp::as<int>(factor[rrMetadata ? "rr_nload" : "fa_nload"]);
  const Rcpp::IntegerVector rows = factor[rrMetadata ? "rr_row" : "fa_row"];
  const Rcpp::IntegerVector columns = factor[rrMetadata ? "rr_col" : "fa_col"];
  double reference = 0.0;
  for(int a = 0; a < count; ++a) {
    if(rows[a] == 1 && columns[a] == 1) { reference = par(a); }
  }
  const double twiceReference = 2.0 * reference;
  const double logScale = std::max(0.0, twiceReference) +
    std::log1p(std::exp(-std::abs(twiceReference)));
  covariance.diag() += std::exp(-logScale);
  if(heterogeneousSpecificFactor(factor)) {
    for(arma::uword row = 1; row < covariance.n_rows; ++row) {
      covariance(row, row) += std::exp(par(count + row - 1) - logScale) -
        std::exp(-logScale);
    }
  }
  return covariance;
}

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
