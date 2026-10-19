#include <Rcpp.h>
#include "../../src/mme_control_r.h"

// [[Rcpp::export]]
std::string control_mode_probe(const Rcpp::List & descriptors,
                              const std::string & fallback = "auto") {
  sommer::control::Mode mode = sommer::control::Mode::Auto;
  if(!sommer::control::parseMode(fallback, mode)) { Rcpp::stop("Invalid test default."); }
  const Rcpp::Function getOption = Rcpp::Environment::base_env()["getOption"];
  return sommer::control::modeName(sommer::control::Transport(descriptors).experimentalMode(
    "testMode", "sommer.mme.testMode", mode, getOption));
}

// [[Rcpp::export]]
Rcpp::List control_proposal_probe(const Rcpp::List & descriptors) {
  const sommer::control::Transport transport(descriptors);
  const sommer::control::Method method = transport.method();
  return Rcpp::List::create(
    Rcpp::Named("hybrid") = method.curvature.method == sommer::control::Curvature::HybridAI,
    Rcpp::Named("aitken") = method.proposal.extrapolation == sommer::control::Extrapolation::Aitken,
    Rcpp::Named("accelerated") = method.proposal.accelerationRequested,
    Rcpp::Named("blend") = sommer::control::blendName(transport.blend("Unknown EM blending method.")),
    Rcpp::Named("adaptive") = transport.logical("adaptiveRecovery"));
}

// [[Rcpp::export]]
std::string control_memory_probe(const Rcpp::List & descriptors) {
  const Rcpp::Function getOption = Rcpp::Environment::base_env()["getOption"];
  return sommer::control::memoryName(sommer::control::Transport(descriptors).memory(getOption));
}
