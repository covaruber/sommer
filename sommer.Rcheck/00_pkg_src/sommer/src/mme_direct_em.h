#ifndef SOMMER_MME_DIRECT_EM_H
#define SOMMER_MME_DIRECT_EM_H

#include "mme_em_optimizer.h"
#include "mme_recovery.h"
#include "mme_ai_recovery.h"

namespace sommer {

struct DirectEMTerm {
  arma::mat map;
  arma::mat covariance;
};

// Keep zero-specific loadings in the design, not in a singular prior density.
inline DirectEMTerm directEMTerm(const Rcpp::List & descriptor,
                                 const arma::vec & parameters) {
  DirectEMTerm out;
  out.map.ones(1, 1);
  out.covariance.ones(1, 1);
  const Rcpp::List factors = descriptor["factors"];
  bool reducedRank = false;
  for(int axis = 0; axis < factors.size(); ++axis) {
    reducedRank = reducedRank || zeroSpecificFactor(Rcpp::as<Rcpp::List>(factors[axis]));
  }
  if(!reducedRank) {
    out.covariance = directEvaluateDescriptor(descriptor, parameters);
    out.map.eye(out.covariance.n_rows, out.covariance.n_rows);
    return out;
  }
  for(int axis = 0; axis < factors.size(); ++axis) {
    const Rcpp::List factor = factors[axis];
    const int dimension = Rcpp::as<int>(factor["dim"]);
    const int first = Rcpp::as<int>(factor["par_start"]) - 1;
    const int last = Rcpp::as<int>(factor["par_end"]) - 1;
    const arma::vec local = last >= first
      ? arma::vec(parameters.subvec(first, last)) : arma::vec();
    if(zeroSpecificFactor(factor)) {
      const arma::mat loading = normalizedFactorLoadings(factor, local);
      out.map = arma::kron(out.map, loading);
      out.covariance = arma::kron(out.covariance,
        arma::eye<arma::mat>(loading.n_cols, loading.n_cols));
    } else {
      out.map = arma::kron(out.map, arma::eye<arma::mat>(dimension, dimension));
      out.covariance = arma::kron(out.covariance, directEvalFactor(factor, local));
    }
  }
  out.covariance *= std::exp(parameters(0));
  return out;
}

inline arma::mat directEMKernel(const arma::field<arma::mat> & kernels,
                                const arma::mat & coefficients) {
  arma::mat out(kernels(0, 0).n_rows, kernels(0, 0).n_cols, arma::fill::zeros);
  for(arma::uword a = 0; a < coefficients.n_rows; ++a) {
    for(arma::uword b = 0; b < coefficients.n_cols; ++b) {
      if(coefficients(a, b) != 0.0) { out += coefficients(a, b) * kernels(a, b); }
    }
  }
  return out;
}

struct DirectEMProposal {
  arma::vec point;
  arma::vec ordinary;
  arma::vec expansion;
  double qGain = 0.0;
  double pxGain = 0.0;
  bool closed = false;
  bool converged = false;
  bool pxConverged = false;
  int evaluations = 0;
  int ordinaryEvaluations = 0;
  int pxEvaluations = 0;
};

inline DirectEMProposal directEMProposal(
    const std::vector<Rcpp::List> & descriptors,
    const arma::field<arma::vec> & parameters,
    const arma::vec & freeFlags, const arma::vec & lower, const arma::vec & upper,
    const arma::mat & tangent, const bool constrained,
    const std::vector<arma::field<arma::mat>> & kernels,
    const std::vector<arma::uword> & widths,
    const arma::mat & X, const arma::vec & y,
    const arma::mat & Vi, const arma::mat & P, const bool reml,
    const std::function<arma::mat(const arma::vec &)> & residual,
    const bool expanded) {
  DirectEMProposal out;
  const int terms = parameters.n_elem, randomTerms = terms - 1;
  const arma::uword n = y.n_elem;
  std::vector<arma::uword> offsets(terms + 1, 0);
  for(int term = 0; term < terms; ++term) {
    offsets[term + 1] = offsets[term] + parameters(term).n_elem;
  }
  arma::vec point(offsets.back());
  for(int term = 0; term < terms; ++term) {
    point.subvec(offsets[term], offsets[term + 1] - 1) = parameters(term);
  }
  const arma::vec r = P * y;
  const arma::mat A = reml ? P : Vi;
  const arma::mat VX = Vi * X;
  arma::mat betaVariance;
  if(!arma::inv_sympd(betaVariance, X.t() * VX)) {
    Rcpp::stop("Direct EM fixed-effect information is numerically singular.");
  }
  const arma::vec beta = betaVariance * (VX.t() * y);
  std::vector<DirectEMTerm> current(randomTerms);
  std::vector<arma::mat> moments(randomTerms);
  arma::mat G(n, n, arma::fill::zeros);
  bool movingMaps = false;
  for(int term = 0; term < randomTerms; ++term) {
    current[term] = directEMTerm(descriptors[term], parameters(term));
    const DirectEMTerm & state = current[term];
    const arma::uword q = state.covariance.n_rows;
    movingMaps = movingMaps || state.map.n_rows != state.map.n_cols;
    const Rcpp::List factors = descriptors[term]["factors"];
    for(int axis = 0; axis < factors.size(); ++axis) {
      movingMaps = movingMaps || zeroSpecificFactor(Rcpp::as<Rcpp::List>(factors[axis]));
    }
    arma::mat contraction(q, q, arma::fill::zeros);
    for(arma::uword a = 0; a < q; ++a) {
      for(arma::uword b = a; b < q; ++b) {
        const arma::mat B = directEMKernel(kernels[term],
          arma::mat(state.map.col(a) * state.map.col(b).t()));
        contraction(a, b) = arma::as_scalar(r.t() * B * r) - arma::accu(A % B);
        contraction(b, a) = contraction(a, b);
      }
    }
    // Contract posterior coefficient uncertainty without forming coefficient-sized matrices.
    moments[term] = arma::symmatu(widths[term] * state.covariance +
      state.covariance * contraction * state.covariance);
    G += directEMKernel(kernels[term],
      arma::mat(state.map * state.covariance * state.map.t()));
  }
  const arma::mat R = residual(parameters(randomTerms));
  const arma::mat residualCovariance = arma::symmatu(R - R * A * R);
  const arma::vec randomMean = G * r;
  const arma::vec error = y - X * beta - randomMean;
  const arma::mat baseResidualMoment = arma::symmatu(error * error.t() +
    residualCovariance);
  arma::vec cachedResidualParameters;
  double cachedResidualObjective = 0.0;
  auto valid = [&](const arma::vec & candidate) {
    return candidate.is_finite() && arma::all(candidate >= lower) &&
      arma::all(candidate <= upper);
  };
  auto objective = [&](const arma::vec & candidate, const arma::vec & expansion) {
    ++out.evaluations;
    Rcpp::checkUserInterrupt();
    const double invalid = std::numeric_limits<double>::infinity();
    if(!valid(candidate) || !expansion.is_finite()) { return invalid; }
    const bool px = !expansion.is_empty();
    const arma::vec amplitudes = arma::exp(expansion);
    if(px && (!amplitudes.is_finite() || arma::any(amplitudes <= 0.0))) { return invalid; }
    const bool changedDesign = movingMaps || px;
    arma::mat T, candidateG;
    if(changedDesign) {
      T.zeros(n, n);
      candidateG.zeros(n, n);
    }
    double value = 0.0;
    for(int term = 0; term < randomTerms; ++term) {
      const arma::vec local = candidate.subvec(offsets[term], offsets[term + 1] - 1);
      DirectEMTerm state = directEMTerm(descriptors[term], local);
      if(px) {
        state.map *= amplitudes(term);
        state.covariance /= amplitudes(term) * amplitudes(term);
      }
      arma::mat chol;
      if(!state.covariance.is_finite() || !arma::chol(chol, state.covariance)) { return invalid; }
      const arma::mat solved = arma::solve(arma::trimatu(chol),
        arma::solve(arma::trimatl(chol.t()), moments[term]));
      value += 2.0 * widths[term] * arma::accu(arma::log(chol.diag())) + arma::trace(solved);
      if(changedDesign) {
        const DirectEMTerm & old = current[term];
        // Freeze the E-step prior; only candidate loadings/expansion change the design.
        T += directEMKernel(kernels[term],
          arma::mat(state.map * old.covariance * old.map.t()));
        candidateG += directEMKernel(kernels[term],
          arma::mat(state.map * old.covariance * state.map.t()));
      }
    }
    const arma::vec local = candidate.subvec(offsets[randomTerms], offsets[terms] - 1);
    if(!changedDesign && cachedResidualParameters.n_elem == local.n_elem &&
       arma::all(cachedResidualParameters == local)) {
      return value + cachedResidualObjective;
    }
    const arma::mat candidateR = residual(local);
    arma::mat chol;
    if(!candidateR.is_finite() || !arma::chol(chol, candidateR)) { return invalid; }
    arma::mat moment = baseResidualMoment;
    if(changedDesign || !reml) {
      const arma::vec mean = changedDesign ? arma::vec(T * r) : randomMean;
      arma::mat covariance = residualCovariance;
      if(changedDesign) {
        covariance = candidateG - T * Vi * T.t();
        if(reml) {
          // Integrating flat-prior fixed effects adds their joint posterior uncertainty.
          const arma::mat fixed = X - T * VX;
          covariance += fixed * betaVariance * fixed.t();
        }
      }
      arma::vec fixedMean = beta;
      if(!reml) {
        const arma::mat RX = arma::solve(arma::trimatu(chol),
          arma::solve(arma::trimatl(chol.t()), X));
        if(!arma::solve(fixedMean, X.t() * RX, RX.t() * (y - mean),
                       arma::solve_opts::likely_sympd + arma::solve_opts::no_approx)) {
          return invalid;
        }
      }
      const arma::vec shifted = y - X * fixedMean - mean;
      moment = arma::symmatu(covariance + shifted * shifted.t());
    }
    const arma::mat solved = arma::solve(arma::trimatu(chol),
      arma::solve(arma::trimatl(chol.t()), moment));
    const double residualObjective =
      2.0 * arma::accu(arma::log(chol.diag())) + arma::trace(solved);
    if(!changedDesign) {
      cachedResidualParameters = local;
      cachedResidualObjective = residualObjective;
    }
    return value + residualObjective;
  };
  const auto ordinaryObjective = [&](const arma::vec & candidate) {
    return objective(candidate, arma::vec());
  };
  const double base = ordinaryObjective(point);
  if(!std::isfinite(base)) {
    Rcpp::stop("Direct EM Q is non-finite at the accepted parameters.");
  }
  arma::vec candidate = point;
  bool closed = !constrained && !movingMaps;
  for(int term = 0; term < terms && closed; ++term) {
    const arma::uword first = offsets[term], last = offsets[term + 1] - 1;
    const arma::uvec free = arma::find(freeFlags.subvec(first, last) > 0.5);
    if(free.is_empty()) { continue; }
    if(free.n_elem == 1 && free(0) == 0) {
      arma::mat covariance, moment;
      double dimension;
      if(term < randomTerms) {
        covariance = current[term].covariance;
        moment = moments[term];
        dimension = widths[term] * covariance.n_rows;
      } else {
        covariance = R;
        moment = baseResidualMoment;
        if(!reml) {
          arma::mat RX;
          arma::vec fixedMean;
          if(!arma::solve(RX, R, X, arma::solve_opts::likely_sympd) ||
             !arma::solve(fixedMean, X.t() * RX, RX.t() * (y - randomMean),
                          arma::solve_opts::likely_sympd)) {
            Rcpp::stop("Direct EM ML fixed-effect M-step failed.");
          }
          const arma::vec shifted = y - X * fixedMean - randomMean;
          moment = arma::symmatu(residualCovariance + shifted * shifted.t());
        }
        dimension = n;
      }
      arma::mat solved;
      if(!arma::solve(solved, covariance, moment, arma::solve_opts::likely_sympd)) {
        Rcpp::stop("Direct EM scalar M-step failed.");
      }
      const double ratio = arma::trace(solved) / dimension;
      if(!std::isfinite(ratio) || ratio <= 0.0) { closed = false; break; }
      candidate(first) = std::max(lower(first), point(first) + std::log(ratio));
    } else if(term < randomTerms && free.n_elem == last - first + 1) {
      const Rcpp::List factors = descriptors[term]["factors"];
      if(factors.size() != 1) { closed = false; break; }
      const Rcpp::List factor = factors[0];
      const std::string model = factor.containsElementNamed("model")
        ? Rcpp::as<std::string>(factor["model"]) : "";
      if(model != "diag" && model != "us") { closed = false; break; }
      arma::mat covariance = moments[term] / widths[term];
      const double scale = covariance(0, 0);
      if(scale <= 0.0) { closed = false; break; }
      candidate(first) = std::log(scale);
      if(model == "diag") {
        for(arma::uword j = 1; j < covariance.n_rows; ++j) {
          candidate(first + j) = std::log(covariance(j, j) / scale);
        }
      } else {
        arma::mat chol;
        if(!arma::chol(chol, covariance / scale)) { closed = false; break; }
        const arma::mat loading = chol.t();
        const Rcpp::IntegerVector rows = factor["us_row"], columns = factor["us_col"];
        const Rcpp::LogicalVector diagonal = factor["us_diag"];
        for(int j = 0; j < rows.size(); ++j) {
          const double entry = loading(rows[j] - 1, columns[j] - 1);
          candidate(first + 1 + j) = diagonal[j] ? std::log(entry) : entry;
        }
      }
    } else {
      closed = false;
    }
  }
  if(closed && valid(candidate) && ordinaryObjective(candidate) <= base) {
    out.point = candidate;
    out.closed = true;
    out.converged = true;
  } else {
    const auto reducedObjective = [&](const arma::vec & reduced) {
      return ordinaryObjective(point + tangent * reduced);
    };
    const EMMinimum minimum = minimizeEM(
      arma::vec(tangent.n_cols, arma::fill::zeros), reducedObjective);
    out.point = point + tangent * minimum.point;
    out.converged = minimum.converged;
  }
  out.ordinary = out.point;
  const double ordinary = ordinaryObjective(out.point);
  out.ordinaryEvaluations = out.evaluations;
  out.expansion.zeros(randomTerms);
  if(expanded && randomTerms > 0) {
    const auto pxObjective = [&](const arma::vec & reduced) {
      return objective(out.ordinary + tangent * reduced.head(tangent.n_cols),
                       reduced.tail(randomTerms));
    };
    const EMMinimum minimum = minimizeEM(
      arma::vec(tangent.n_cols + randomTerms, arma::fill::zeros), pxObjective);
    out.point = out.ordinary + tangent * minimum.point.head(tangent.n_cols);
    out.expansion = minimum.point.tail(randomTerms);
    out.pxGain = 0.5 * (ordinary - minimum.value);
    out.qGain = 0.5 * (base - minimum.value);
    out.pxConverged = minimum.converged;
    out.pxEvaluations = out.evaluations - out.ordinaryEvaluations;
  } else {
    out.qGain = 0.5 * (base - ordinary);
  }
  if(!std::isfinite(out.qGain) || out.qGain < -1.0e-8 * (1.0 + std::abs(base))) {
    Rcpp::stop("Direct EM M-step decreased its expected complete-data log-likelihood.");
  }
  return out;
}

} // namespace sommer

#endif
