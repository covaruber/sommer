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
