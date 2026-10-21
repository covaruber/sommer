// [[Rcpp::depends(RcppArmadillo)]]
#include "../../src/mme_ai_recovery.h"

// [[Rcpp::export]]
Rcpp::List regularized_ai_probe(const arma::mat & information,
    const arma::vec & score, const arma::mat & tangent, const arma::vec & caps) {
  const auto result = sommer::regularizedAIDirection(information, score, tangent, caps);
  return Rcpp::List::create(Rcpp::Named("reliable") = result.reliable,
    Rcpp::Named("direction") = result.direction, Rcpp::Named("ridge") = result.ridge,
    Rcpp::Named("ordinaryReliable") = sommer::reliableAIInformation(
      tangent.t() * information * tangent, sommer::aiTrustScale(tangent, caps)));
}

// [[Rcpp::export]]
Rcpp::List reduced_ai_recovery_probe(const arma::mat & information,
    const arma::vec & score, const arma::mat & basis,
    const arma::mat & feasible, const arma::vec & caps) {
  const arma::mat reduced = basis.t() * information * basis;
  const arma::mat projected = sommer::aiInformationInTangent(reduced, basis, feasible);
  const auto full = sommer::regularizedAIDirection(information, score, feasible, caps);
  const auto packed = sommer::regularizedReducedAIDirection(
    projected, feasible.t() * score, feasible, caps);
  return Rcpp::List::create(Rcpp::Named("projected") = projected,
    Rcpp::Named("full") = full.direction, Rcpp::Named("packed") = packed.direction,
    Rcpp::Named("fullReliable") = full.reliable, Rcpp::Named("packedReliable") = packed.reliable);
}

// [[Rcpp::export]]
Rcpp::List feasible_ai_probe(const arma::mat & information, const arma::vec & score,
    const arma::mat & tangent, const arma::vec & point, const arma::vec & lower,
    const arma::vec & caps, bool reduced = false) {
  const arma::mat packed = reduced ? arma::mat(tangent.t() * information * tangent) : information;
  const auto result = sommer::feasibleAIDirection(packed, score,
    reduced ? tangent : arma::mat(), tangent, point, lower, caps);
  const arma::vec candidate = point + result.direction;
  return Rcpp::List::create(Rcpp::Named("reliable") = result.reliable,
    Rcpp::Named("direction") = result.direction,
    Rcpp::Named("cappedDirection") = result.cappedDirection,
    Rcpp::Named("ridge") = result.ridge,
    Rcpp::Named("attempts") = result.regularizedAttempts,
    Rcpp::Named("nextTangent") = sommer::varianceFloorTangent(tangent, candidate, lower,
      arma::vec(tangent * (tangent.t() * score))));
}
