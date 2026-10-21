#ifndef SOMMER_MME_SCALE_H
#define SOMMER_MME_SCALE_H

#include "mme_ai_recovery.h"

namespace sommer {

using CommonScaleFit = std::function<Rcpp::List(const Rcpp::List &, int)>;

inline Rcpp::List fitCommonScale(
    const Rcpp::List & descriptors, int iterations,
    double likelihoodTolerance, double parameterTolerance, double inverseTolerance,
    const arma::vec & steps, int inverseMode, bool verbose, bool reml,
    const CommonScaleFit & conditionalFit) {
  if(hasZeroSpecific(descriptors)){
    Rcpp::stop("Common-scale profiling does not support latent reduced-rank covariance.");
  }
  const int terms = descriptors.size();
  std::vector<arma::uword> offsets(terms + 1, 0);
  for(int term = 0; term < terms; ++term){
    const Rcpp::List descriptor = descriptors[term];
    offsets[term + 1] = offsets[term] + Rcpp::as<arma::vec>(descriptor["par"]).n_elem;
  }
  const arma::uword count = offsets.back(), reference = offsets[terms - 1];
  arma::vec parameters(count), common(count, arma::fill::zeros), freeFlags(count);
  for(int term = 0; term < terms; ++term){
    const Rcpp::List descriptor = descriptors[term];
    parameters.subvec(offsets[term], offsets[term + 1] - 1) =
      Rcpp::as<arma::vec>(descriptor["par"]);
    const Rcpp::LogicalVector free = descriptor["free"];
    for(arma::uword k = offsets[term]; k < offsets[term + 1]; ++k){
      freeFlags(k) = free[k - offsets[term]] == TRUE;
    }
    if(freeFlags(offsets[term]) == 0.0){
      Rcpp::stop("Common-scale profiling requires every product-level variance scale to be free.");
    }
    common(offsets[term]) = 1.0;
  }
  const arma::uvec free = arma::find(freeFlags > 0.5);
  const VarianceParameterMap map = readVarianceParameterMap(descriptors, count);
  arma::mat tangent(count, free.n_elem, arma::fill::zeros);
  for(arma::uword k = 0; k < free.n_elem; ++k){ tangent(free(k), k) = 1.0; }
  if(map.active){ tangent = map.T; }
  arma::vec scaleCoordinates;
  if(!arma::solve(scaleCoordinates, tangent, common, arma::solve_opts::no_approx) ||
     !scaleCoordinates.is_finite() ||
     arma::norm(tangent * scaleCoordinates - common, 2) > 1.0e-8){
    Rcpp::stop("Common-scale profiling requires constraints that preserve a common variance scale.");
  }
  const arma::uvec anchors = arma::find(arma::abs(tangent.row(reference).t()) > 0.0);
  if(anchors.n_elem == 0){ Rcpp::stop("Common-scale profiling has no free residual reference scale."); }
  const arma::uword anchor = anchors(0);
  arma::mat shape(count, tangent.n_cols - 1);
  for(arma::uword k = 0, column = 0; k < tangent.n_cols; ++k){
    if(k == anchor){ continue; }
    shape.col(column++) = tangent.col(k) -
      tangent.col(anchor) * (tangent(reference, k) / tangent(reference, anchor));
  }
  const arma::vec varianceLower = covarianceVarianceLowerBounds(descriptors, inverseTolerance);
  const arma::vec caps = covarianceWorkingStepCaps(descriptors);
  const bool retainInformation = !descriptors.hasAttribute("covarianceInformation") ||
    Rcpp::as<std::string>(descriptors.attr("covarianceInformation")) == "full";
  parameters -= common * parameters(reference);
  int evaluations = 0, halvings = 0, regularizedSteps = 0, boundaryRejections = 0;
  double columns = 0.0, coefficientEvaluations = 0.0, profileColumns = 0.0;
  double scoreNorm = std::numeric_limits<double>::infinity();
  const double scoreTolerance = std::sqrt(std::max(likelihoodTolerance, parameterTolerance));
  arma::vec history, scaleHistory, phaseSeconds;
  arma::mat monitor, norms;
  Rcpp::List fit;
  Rcpp::RObject ldltCache = R_NilValue, cholmodCache = R_NilValue;
  const auto observe = [&](Rcpp::List value){
    const Rcpp::List diagnostics = value["engineDiagnostics"];
    coefficientEvaluations += Rcpp::as<double>(diagnostics["coefficientEvaluations"]);
    columns += Rcpp::as<double>(diagnostics["nativeAISensitivityColumns"]);
    if(diagnostics.containsElementNamed("profile") && !Rf_isNull(diagnostics["profile"])){
      const Rcpp::List profile = diagnostics["profile"];
      const arma::vec seconds = Rcpp::as<arma::vec>(profile["seconds"]);
      if(phaseSeconds.is_empty()){ phaseSeconds.zeros(seconds.n_elem); }
      phaseSeconds += seconds;
    }
    if(value.containsElementNamed(".ldltCache")){ ldltCache = value[".ldltCache"]; }
    if(value.containsElementNamed(".cholmodCache")){ cholmodCache = value[".cholmodCache"]; }
    return value;
  };
  const auto covarianceAt = [&](const arma::vec & point){
    Rcpp::List covariance = Rcpp::clone(descriptors);
    covariance.attr("commonScale") = "off";
    covariance.attr("directionBlend") = false;
    covariance.attr("emInformationOnly") = true;
    covariance.attr("retainLDLTCache") = true;
    covariance.attr("retainCholmodCache") = true;
    if(!ldltCache.isNULL()){ covariance.attr("ldltCache") = ldltCache; }
    if(!cholmodCache.isNULL()){ covariance.attr("cholmodCache") = cholmodCache; }
    for(int term = 0; term < terms; ++term){
      Rcpp::List descriptor = covariance[term];
      descriptor["par"] = point.subvec(offsets[term], offsets[term + 1] - 1);
      covariance[term] = descriptor;
    }
    return covariance;
  };
  const auto evaluate = [&](const arma::vec & point){
    ++evaluations;
    Rcpp::List covariance = covarianceAt(point);
    covariance.attr("aiTangent") = shape;
    covariance.attr("profileCommonScale") = true;
    covariance.attr("covarianceInformation") = "none";
    Rcpp::List value = observe(conditionalFit(covariance, 0));
    const Rcpp::List diagnostics = value["engineDiagnostics"];
    profileColumns += Rcpp::as<double>(diagnostics["nativeAISensitivityColumns"]);
    return value;
  };
  const auto reportedPoint = [&](const arma::vec & point, const Rcpp::List & value) -> arma::vec {
    const double scale = Rcpp::as<double>(value[".commonScale"]);
    return point + common * std::log(scale);
  };
  const auto insideScaleBounds = [&](const arma::vec & point){
    if(!point.is_finite()){ return false; }
    for(arma::uword k = 0; k < count; ++k){
      if(point(k) <= varianceLower(k)){ return false; }
    }
    return true;
  };
  fit = evaluate(parameters);
  if(!insideScaleBounds(reportedPoint(parameters, fit))){
    Rcpp::stop("Common-scale profiling reaches a variance floor at the initial shape; refit with commonScale='off'.");
  }
  MMEProgress progress(terms);
  bool converged = false, stalled = false;
  double previousLikelihood = -std::numeric_limits<double>::infinity();
  for(int iteration = 0; iteration < iterations; ++iteration){
    Rcpp::checkUserInterrupt();
    const arma::vec score = -0.5 * Rcpp::as<arma::vec>(fit["dLu"]);
    const arma::mat information = 0.5 * Rcpp::as<arma::mat>(fit["InfMat"]);
    const arma::vec projected = shape.t() * score;
    scoreNorm = 0.0;
    for(arma::uword k = 0; k < projected.n_elem; ++k){
      scoreNorm = std::max(scoreNorm, std::abs(projected(k)) /
        std::sqrt(std::max(inverseTolerance, std::abs(information(k, k)))));
    }
    const arma::vec point = reportedPoint(parameters, fit);
    const arma::vec likelihoods = Rcpp::as<arma::vec>(fit["llik"]);
    const double likelihood = likelihoods(0);
    history.insert_rows(history.n_elem, arma::vec({likelihood}));
    scaleHistory.insert_rows(scaleHistory.n_elem, arma::vec({std::exp(point(reference))}));
    const arma::mat raw = Rcpp::as<arma::mat>(fit["monitor"]);
    arma::vec reported = raw.col(0);
    const double scale = Rcpp::as<double>(fit[".commonScale"]);
    for(int term = 0; term < terms; ++term){ reported(offsets[term]) *= scale; }
    monitor.insert_cols(monitor.n_cols, reported);
    arma::vec reducedStep(shape.n_cols, arma::fill::zeros);
    bool regularized = false;
    bool stepFailed = false;
    if(shape.n_cols > 0){
      const bool reliable = reliableAIInformation(information, aiTrustScale(shape, caps)) &&
        arma::solve(reducedStep, information, projected,
          arma::solve_opts::likely_sympd + arma::solve_opts::no_approx) &&
        reducedStep.is_finite() && arma::dot(reducedStep, projected) > 0.0;
      if(!reliable){
        const RegularizedAIDirection retry = regularizedReducedAIDirection(
          information, projected, shape, caps);
        if(!retry.reliable){
          stalled = true;
          stepFailed = true;
        }else{
          if(!arma::solve(reducedStep, shape, retry.direction, arma::solve_opts::no_approx)){
            Rcpp::stop("Unable to recover profiled AI shape coordinates.");
          }
          regularized = true;
          ++regularizedSteps;
        }
      }
    }
    const arma::vec direction = shape * reducedStep;
    const double relativeStep = arma::norm(direction, 2) / (1.0 + arma::norm(parameters, 2));
    norms.insert_cols(norms.n_cols, arma::vec({arma::norm(direction, 1), scoreNorm, relativeStep}));
    if(verbose){
      if(iteration == 0){
        Rcpp::Rcout << "MME path: exact common-scale profile AI"
                    << " [" << mmePathStatus(fit, false) << "]\n";
        progress.header();
      }
      progress.row(iteration + 1, likelihood, previousLikelihood,
        regularized ? "profile-AI-reg" : "profile-AI", "off", NA_REAL, NA_REAL);
    }
    if(scoreNorm <= scoreTolerance &&
       (shape.n_cols == 0 || relativeStep < parameterTolerance ||
        (iteration > 0 && std::abs(likelihood - previousLikelihood) < likelihoodTolerance))){
      converged = true;
      stalled = false;
      break;
    }
    if(stepFailed){ break; }
    if(iteration == iterations - 1){ break; }
    double alpha = steps(std::min<arma::uword>(iteration, steps.n_elem - 1));
    for(arma::uword k = 0; k < count; ++k){
      if(std::abs(direction(k)) > 0.0){ alpha = std::min(alpha, caps(k) / std::abs(direction(k))); }
    }
    bool accepted = false;
    for(int attempt = 0; attempt < 25; ++attempt){
      const arma::vec candidate = parameters + alpha * direction;
      Rcpp::List trial;
      try {
        trial = evaluate(candidate);
      } catch(const CovarianceFactorizationError &){
        ++halvings;
        alpha *= 0.5;
        continue;
      }
      if(!insideScaleBounds(reportedPoint(candidate, trial))){
        ++boundaryRejections;
        ++halvings;
        alpha *= 0.5;
        continue;
      }
      const arma::vec trialLikelihood = Rcpp::as<arma::vec>(trial["llik"]);
      if(trialLikelihood(0) >= likelihood - 1.0e-10){
        parameters = candidate;
        fit = trial;
        accepted = true;
        break;
      }
      ++halvings;
      alpha *= 0.5;
    }
    if(!accepted){ stalled = true; break; }
    previousLikelihood = likelihood;
  }
  const arma::vec finalPoint = reportedPoint(parameters, fit);
  const double finalScale = std::exp(finalPoint(reference));
  const double degrees = Rcpp::as<double>(fit[".commonScaleDF"]);
  Rcpp::List covariance = covarianceAt(finalPoint);
  if(!retainInformation && tangent.n_cols < count){ covariance.attr("aiTangent") = tangent; }
  Rcpp::List result = observe(conditionalFit(covariance, inverseMode));
  const arma::vec returnedLikelihood = Rcpp::as<arma::vec>(result["llik"]);
  if(std::abs(returnedLikelihood(0) - history(history.n_elem - 1)) >
     1.0e-8 * (1.0 + std::abs(returnedLikelihood(0)))){
    Rcpp::stop("Common-scale final reporting does not match the accepted profiled likelihood.");
  }
  const arma::vec finalScore = -0.5 * Rcpp::as<arma::vec>(result["dLu"]);
  const arma::mat finalInformation = 0.5 * Rcpp::as<arma::mat>(result["InfMat"]);
  const arma::mat basis = result.containsElementNamed(".aiTangent")
    ? Rcpp::as<arma::mat>(result[".aiTangent"]) : arma::mat();
  const arma::mat curvature = aiInformationInTangent(finalInformation, basis, tangent);
  const arma::vec finalProjected = tangent.t() * finalScore;
  double finalNorm = 0.0;
  for(arma::uword k = 0; k < finalProjected.n_elem; ++k){
    finalNorm = std::max(finalNorm, std::abs(finalProjected(k)) /
      std::sqrt(std::max(inverseTolerance, std::abs(curvature(k, k)))));
  }
  converged = converged && finalNorm <= scoreTolerance;
  if(stalled){
    Rcpp::warning("Common-scale profile AI stalled at the last accepted point. "
      "Convergence was not certified; refit with commonScale='off' for ordinary AI/EM recovery.");
  }
  result["llik"] = arma::mat(history.t());
  result["monitor"] = monitor;
  result["normMonitor"] = norms;
  result["convergence"] = converged && !stalled;
  result["percDelta"] = arma::mat(count, monitor.n_cols, arma::fill::zeros);
  result["toBoundary"] = arma::mat(monitor.n_cols, count, arma::fill::zeros);
  Rcpp::List diagnostics = result["engineDiagnostics"];
  diagnostics["optimizer"] = "common-scale-profile-AI";
  diagnostics["coefficientEvaluations"] = coefficientEvaluations;
  diagnostics["commonScale"] = Rcpp::List::create(
    Rcpp::Named("active") = true, Rcpp::Named("fallbackReason") = "",
    Rcpp::Named("likelihoodMode") = reml ? "REML" : "ML",
    Rcpp::Named("degreesOfFreedom") = degrees,
    Rcpp::Named("scale") = finalScale, Rcpp::Named("scaleHistory") = scaleHistory,
    Rcpp::Named("independentParameters") = static_cast<double>(tangent.n_cols),
    Rcpp::Named("shapeParameters") = static_cast<double>(shape.n_cols),
    Rcpp::Named("profileEvaluations") = evaluations,
    Rcpp::Named("sensitivityColumns") = columns,
    Rcpp::Named("profileSensitivityColumns") = profileColumns,
    Rcpp::Named("finalInformationEvaluations") = 1,
    Rcpp::Named("lineSearchHalvings") = halvings,
    Rcpp::Named("varianceFloorRejections") = boundaryRejections,
    Rcpp::Named("regularizedSteps") = regularizedSteps,
    Rcpp::Named("finalProjectedScoreNorm") = finalNorm,
    Rcpp::Named("scoreTolerance") = scoreTolerance,
    Rcpp::Named("stopReason") = converged ? "converged" : stalled ? "stalled" : "iteration-limit");
  if(!phaseSeconds.is_empty()){
    Rcpp::List profile = diagnostics["profile"];
    Rcpp::NumericVector seconds = profile["seconds"];
    std::copy(phaseSeconds.begin(), phaseSeconds.end(), seconds.begin());
    profile["seconds"] = seconds;
    profile["secondsByIteration"] = R_NilValue;
    profile["scope"] = "all common-scale trials and final reporting; counters describe final evaluation";
    diagnostics["profile"] = profile;
  }
  result["engineDiagnostics"] = diagnostics;
  if(!retainInformation){ result["InfMat"] = R_NilValue; result["theta_se"] = R_NilValue; }
  if(result.containsElementNamed(".aiTangent")){
    const Rcpp::CharacterVector names = result.names();
    for(int k = 0; k < names.size(); ++k){
      if(Rcpp::as<std::string>(names[k]) == ".aiTangent"){ result.erase(k); break; }
    }
  }
  return result;
}

}
#endif
