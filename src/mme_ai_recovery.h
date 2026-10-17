#ifndef SOMMER_MME_AI_RECOVERY_H
#define SOMMER_MME_AI_RECOVERY_H

#include <RcppArmadillo.h>
#include <limits>
#include <utility>

namespace sommer {

struct RegularizedAIDirection {
  arma::vec direction;
  double ridge = 0.0;
  bool reliable = false;
};

inline arma::vec aiTrustScale(const arma::mat & tangent, const arma::vec & caps) {
  const arma::mat metric = tangent.each_col() / caps;
  return 1.0 / arma::sqrt(arma::sum(arma::square(metric), 0).t());
}

inline bool reliableAIInformation(const arma::mat & information, const arma::vec & scale) {
  if(information.n_rows == 0 || !scale.is_finite()) { return false; }
  const arma::mat units = arma::diagmat(scale);
  const double reciprocalCondition = arma::rcond(units * information * units);
  // A small backward residual alone can hide loss of most forward precision.
  return std::isfinite(reciprocalCondition) &&
    reciprocalCondition >= std::sqrt(std::numeric_limits<double>::epsilon());
}

// Solve in trust-cap units so a nearly flat log-variance cannot shrink every
// other coordinate through the global trust multiplier.
inline RegularizedAIDirection regularizedAIDirection(
    const arma::mat & information, const arma::vec & score,
    const arma::mat & tangent, const arma::vec & caps) {
  RegularizedAIDirection result;
  if(tangent.n_cols == 0) { return result; }
  const arma::vec scale = aiTrustScale(tangent, caps);
  if(!scale.is_finite()) { return result; }
  const arma::mat scaledTangent = tangent * arma::diagmat(scale);
  const arma::mat curvature = arma::symmatu(
    scaledTangent.t() * information * scaledTangent);
  const arma::vec projected = scaledTangent.t() * score;
  const double magnitude = arma::abs(curvature.diag()).max();
  if(!std::isfinite(magnitude) || magnitude <= 0.0) { return result; }
  double ridge = std::max(magnitude * 1.0e-6, std::numeric_limits<double>::epsilon());
  for(int attempt = 0; attempt < 12; ++attempt, ridge *= 10.0) {
    arma::mat system = curvature;
    system.diag() += ridge;
    arma::mat upper;
    arma::vec step;
    if(!arma::chol(upper, system) ||
       !arma::solve(step, system, projected,
         arma::solve_opts::likely_sympd + arma::solve_opts::no_approx) ||
       !step.is_finite() ||
       arma::norm(system * step - projected, 2) >
         1.0e-6 * (1.0 + arma::norm(projected, 2))) { continue; }
    arma::vec direction = scaledTangent * step;
    if(!direction.is_finite() || arma::dot(direction, score) <= 0.0) { continue; }
    if(arma::max(arma::abs(direction) / caps) > 1.0) { continue; }
    result.direction = std::move(direction);
    result.ridge = ridge;
    result.reliable = true;
    return result;
  }
  return result;
}

} // namespace sommer

#endif
