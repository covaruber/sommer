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

inline arma::mat varianceFloorTangent(
    const arma::mat & tangent, const arma::vec & point,
    const arma::vec & lower, const arma::vec & direction) {
  std::vector<arma::uword> active;
  for(arma::uword j = 0; j < point.n_elem; ++j) {
    if(std::isfinite(lower(j)) &&
       point(j) <= lower(j) + 1.0e-10 * (1.0 + std::abs(lower(j))) &&
       direction(j) < 0.0 && arma::any(tangent.row(j) != 0.0)) {
      active.push_back(j);
    }
  }
  if(active.empty()) { return tangent; }
  return tangent * arma::null(tangent.rows(arma::uvec(active)));
}

inline arma::mat aiInformationInTangent(
    const arma::mat & information, const arma::mat & basis,
    const arma::mat & tangent) {
  if(basis.n_rows == 0) { return tangent.t() * information * tangent; }
  if(tangent.n_cols == 0) { return arma::mat(0, 0); }
  arma::mat coordinates;
  if(!arma::solve(coordinates, basis, tangent, arma::solve_opts::no_approx) ||
     !coordinates.is_finite() ||
     arma::norm(basis * coordinates - tangent, "fro") >
       1.0e-8 * (1.0 + arma::norm(tangent, "fro"))) {
    Rcpp::stop("AI feasible directions are outside the reduced information tangent.");
  }
  return coordinates.t() * information * coordinates;
}

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
inline RegularizedAIDirection regularizedAIScaledDirection(
    const arma::mat & curvature, const arma::vec & projected,
    const arma::mat & scaledTangent, const arma::vec & caps,
    const arma::vec & fullScore = arma::vec()) {
  RegularizedAIDirection result;
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
    const double descent = fullScore.n_elem > 0
      ? arma::dot(direction, fullScore) : arma::dot(step, projected);
    if(!direction.is_finite() || descent <= 0.0) { continue; }
    if(arma::max(arma::abs(direction) / caps) > 1.0) { continue; }
    result.direction = std::move(direction);
    result.ridge = ridge;
    result.reliable = true;
    return result;
  }
  return result;
}

inline RegularizedAIDirection regularizedReducedAIDirection(
    const arma::mat & information, const arma::vec & projectedScore,
    const arma::mat & tangent, const arma::vec & caps) {
  if(tangent.n_cols == 0) { return RegularizedAIDirection(); }
  const arma::vec scale = aiTrustScale(tangent, caps);
  if(!scale.is_finite()) { return RegularizedAIDirection(); }
  return regularizedAIScaledDirection(
    arma::symmatu(arma::diagmat(scale) * information * arma::diagmat(scale)),
    scale % projectedScore, tangent * arma::diagmat(scale), caps);
}

inline RegularizedAIDirection regularizedAIDirection(
    const arma::mat & information, const arma::vec & score,
    const arma::mat & tangent, const arma::vec & caps) {
  if(tangent.n_cols == 0) { return RegularizedAIDirection(); }
  const arma::vec scale = aiTrustScale(tangent, caps);
  if(!scale.is_finite()) { return RegularizedAIDirection(); }
  const arma::mat scaledTangent = tangent * arma::diagmat(scale);
  return regularizedAIScaledDirection(
    arma::symmatu(scaledTangent.t() * information * scaledTangent),
    scaledTangent.t() * score, scaledTangent, caps, score);
}

struct FeasibleAIDirection {
  arma::vec direction;
  arma::vec cappedDirection;
  bool reliable = false;
  int regularizedAttempts = 0;
  double ridge = 0.0;
};

inline FeasibleAIDirection feasibleAIDirection(
    const arma::mat & information, const arma::vec & score,
    const arma::mat & informationTangent, const arma::mat & tangent,
    const arma::vec & point, const arma::vec & lower, const arma::vec & caps) {
  FeasibleAIDirection result;
  result.direction.zeros(point.n_elem);
  result.cappedDirection.zeros(point.n_elem);
  const arma::mat base = varianceFloorTangent(tangent, point, lower,
    arma::vec(tangent * (tangent.t() * score)));
  arma::mat feasible = base;
  arma::vec boundary(point.n_elem, arma::fill::zeros);
  std::vector<arma::uword> active;
  for(arma::uword iteration = 0; iteration <= point.n_elem; ++iteration) {
    const arma::mat curvature = aiInformationInTangent(information, informationTangent, feasible);
    arma::vec projected = feasible.t() * score;
    if(!active.empty() && feasible.n_cols > 0) {
      arma::mat joint = arma::join_rows(feasible, boundary);
      const arma::mat cross = aiInformationInTangent(information, informationTangent, joint);
      projected -= cross.submat(0, feasible.n_cols, feasible.n_cols - 1, feasible.n_cols);
    }
    arma::vec coordinates;
    bool solved = feasible.n_cols == 0;
    bool reliable = solved;
    if(feasible.n_cols > 0) {
      const arma::vec units = arma::sqrt(arma::clamp(arma::abs(curvature.diag()),
        std::numeric_limits<double>::min(), std::numeric_limits<double>::max()));
      const arma::mat normalized = curvature / (units * units.t());
      const arma::vec rhs = projected / units;
      arma::mat upper;
      solved = arma::chol(upper, normalized) &&
        arma::solve(coordinates, normalized, rhs,
        arma::solve_opts::likely_sympd + arma::solve_opts::no_approx) &&
        coordinates.is_finite() &&
        arma::norm(normalized * coordinates - rhs, 2) <= 1.0e-6 * (1.0 + arma::norm(rhs, 2));
      if(solved) { coordinates /= units; }
      reliable = solved && reliableAIInformation(curvature, aiTrustScale(feasible, caps)) &&
        (arma::dot(coordinates, projected) > 0.0 || arma::norm(projected, 2) == 0.0);
    } else {
      coordinates.set_size(0);
    }
    arma::vec remaining = solved ? arma::vec(feasible * coordinates)
                                : arma::vec(point.n_elem, arma::fill::zeros);
    arma::vec direction = boundary + remaining;
    // A finite Newton direction can identify a boundary even when tiny
    // log-variance curvature makes the full information solve unreliable.
    arma::uword crossing = point.n_elem;
    double fraction = 1.0;
    if(solved) {
      for(arma::uword j = 0; j < point.n_elem; ++j) {
        if(std::isfinite(lower(j)) && direction(j) < 0.0 &&
           point(j) + direction(j) < lower(j) - 1.0e-10 &&
           arma::any(feasible.row(j) != 0.0)) {
          const arma::vec axis = feasible.row(j).t();
          const double length = arma::dot(axis, axis);
          const double marginalScore = arma::dot(axis, projected) / length;
          const double marginalCurvature = arma::as_scalar(axis.t() * curvature * axis) /
            (length * length);
          // Do not pin a well-supported variance just because a correlated
          // joint Newton step overshoots far from the optimum.
          if(marginalScore >= 0.0 || marginalCurvature <= 0.0 ||
             -marginalScore / marginalCurvature < point(j) - lower(j)) { continue; }
          const double candidate = (point(j) - lower(j)) / -direction(j);
          if(candidate < fraction) { fraction = candidate; crossing = j; }
        }
      }
    }
    if(crossing < point.n_elem) {
      active.push_back(crossing);
      const arma::uvec indices(active);
      const arma::mat constraints = base.rows(indices);
      const arma::vec targets = lower(indices) - point(indices);
      const arma::vec offset = arma::pinv(constraints) * targets;
      if(!offset.is_finite() ||
         arma::norm(constraints * offset - targets, 2) > 1.0e-8 * (1.0 + arma::norm(targets, 2))) {
        Rcpp::stop("Inconsistent variance-floor constraints in the AI tangent.");
      }
      boundary = base * offset;
      feasible = base * arma::null(constraints);
      continue;
    }
    double trust = 1.0;
    for(arma::uword j = 0; j < point.n_elem; ++j) {
      if(remaining(j) != 0.0) { trust = std::min(trust, caps(j) / std::abs(remaining(j))); }
    }
    if(!reliable || trust < 0.1) {
      ++result.regularizedAttempts;
      const RegularizedAIDirection retry =
        regularizedReducedAIDirection(curvature, projected, feasible, caps);
      if(!retry.reliable) { return result; }
      remaining = retry.direction;
      direction = boundary + remaining;
      result.ridge = retry.ridge;
    }
    if(!direction.is_finite() || arma::dot(direction, score) < 0.0) { return result; }
    result.direction = direction;
    result.cappedDirection = remaining;
    result.reliable = true;
    return result;
  }
  return result;
}

} // namespace sommer

#endif
