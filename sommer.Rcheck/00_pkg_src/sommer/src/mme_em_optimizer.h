#ifndef SOMMER_MME_EM_OPTIMIZER_H
#define SOMMER_MME_EM_OPTIMIZER_H

#include <RcppArmadillo.h>
#include <functional>

namespace sommer {

struct EMMinimum {
  arma::vec point;
  double value = 0.0;
  bool converged = false;
  int evaluations = 0;
};

// Every accepted numerical M-step decreases -2Q, including early termination.
inline EMMinimum minimizeEM(
    const arma::vec & start,
    const std::function<double(const arma::vec &)> & objective) {
  EMMinimum out;
  out.point = start;
  auto evaluate = [&](const arma::vec & point) {
    ++out.evaluations;
    Rcpp::checkUserInterrupt();
    return objective(point);
  };
  out.value = evaluate(start);
  if(!std::isfinite(out.value)) {
    Rcpp::stop("The EM objective is non-finite at the accepted parameter vector.");
  }
  if(start.is_empty()) { out.converged = true; return out; }
  auto gradient = [&](const arma::vec & point, double value) {
    arma::vec g(point.n_elem);
    for(arma::uword j = 0; j < point.n_elem; ++j) {
      double h = 1.0e-5 * (1.0 + std::abs(point(j)));
      double plus = std::numeric_limits<double>::infinity();
      double minus = plus;
      for(int refinement = 0; refinement < 10; ++refinement) {
        arma::vec trial = point;
        trial(j) += h;
        plus = evaluate(trial);
        trial(j) = point(j) - h;
        minus = evaluate(trial);
        if(std::isfinite(plus) || std::isfinite(minus)) { break; }
        if(refinement < 9) { h *= 0.25; }
      }
      if(std::isfinite(plus) && std::isfinite(minus)) {
        g(j) = (plus - minus) / (2.0 * h);
      } else if(std::isfinite(plus)) {
        g(j) = (plus - value) / h;
      } else if(std::isfinite(minus)) {
        g(j) = (value - minus) / h;
      } else {
        Rcpp::stop("The EM objective has no finite numerical derivative in a free direction.");
      }
    }
    return g;
  };
  arma::mat inverse = arma::eye<arma::mat>(start.n_elem, start.n_elem);
  arma::vec g = gradient(out.point, out.value);
  for(int iteration = 0; iteration < 60; ++iteration) {
    if(arma::abs(g).max() <= 1.0e-7 * (1.0 + std::abs(out.value))) {
      out.converged = true;
      break;
    }
    arma::vec direction = -inverse * g;
    if(arma::dot(direction, g) >= 0.0) {
      inverse.eye();
      direction = -g;
    }
    const double largest = arma::abs(direction).max();
    if(largest > 1.0) { direction /= largest; }
    double alpha = 1.0;
    arma::vec candidate;
    double next = out.value;
    bool accepted = false;
    for(int search = 0; search < 30; ++search) {
      candidate = out.point + alpha * direction;
      next = evaluate(candidate);
      if(std::isfinite(next) &&
         next <= out.value + 1.0e-4 * alpha * arma::dot(g, direction)) {
        accepted = true;
        break;
      }
      alpha *= 0.5;
    }
    if(!accepted) { break; }
    const arma::vec nextGradient = gradient(candidate, next);
    const arma::vec step = candidate - out.point;
    const arma::vec change = nextGradient - g;
    const double curvature = arma::dot(step, change);
    if(curvature > 1.0e-12 * arma::norm(step) * arma::norm(change)) {
      const double rho = 1.0 / curvature;
      const arma::mat identity = arma::eye<arma::mat>(start.n_elem, start.n_elem);
      const arma::mat transfer = identity - rho * step * change.t();
      inverse = transfer * inverse * transfer.t() + rho * step * step.t();
    } else {
      inverse.eye();
    }
    out.point = std::move(candidate);
    out.value = next;
    g = nextGradient;
  }
  return out;
}

} // namespace sommer

#endif
