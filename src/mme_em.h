#ifndef SOMMER_MME_EM_H
#define SOMMER_MME_EM_H
#include "mme_recovery.h"
#include "mme_progress.h"
#include "mme_ai_recovery.h"
#include "mme_em_optimizer.h"
#include "mme_control_r.h"

namespace sommer {

struct CovarianceFactorizationError : std::runtime_error {
  using std::runtime_error::runtime_error;
};

using EMConditionalFit = std::function<Rcpp::List(
  const Rcpp::List &, const arma::vec &, const Rcpp::List &, int,
  const RRConditionalCallback *)>;

inline double projectedEMScoreNorm(
    const arma::mat & tangent, const arma::vec & point,
    const arma::vec & lower,
    const arma::vec & score, const arma::mat & information, double tolerance,
    const arma::mat & informationTangent = arma::mat()) {
  const arma::mat feasible = varianceFloorTangent(tangent, point, lower,
    arma::vec(tangent * (tangent.t() * score)));
  const arma::vec projected = feasible.t() * score;
  const arma::mat curvature = aiInformationInTangent(information, informationTangent, feasible);
  double norm = 0.0;
  for(arma::uword j = 0; j < projected.n_elem; ++j) {
    double floor = tolerance;
    for(arma::uword k = 0; k < point.n_elem; ++k) {
      if(feasible(k, j) != 0.0 && std::isfinite(lower(k)) &&
         point(k) <= lower(k) + 1.0e-10 * (1.0 + std::abs(lower(k)))) {
        floor = std::numeric_limits<double>::min();
        break;
      }
    }
    norm = std::max(norm, std::abs(projected(j)) /
      std::sqrt(std::max(floor, std::abs(curvature(j, j)))));
  }
  return norm;
}

inline Rcpp::List fitEMBlend(
    const arma::sp_mat & X, const Rcpp::List & Z, const arma::vec & Zind,
    const Rcpp::List & Ai, const arma::sp_mat & H, const bool useH,
    const arma::uvec & residualBlock, const arma::uvec & residualIndex,
    const Rcpp::List & descriptors, const int iterations,
    const double likelihoodTolerance, const double parameterTolerance,
    const double inverseTolerance, const arma::vec & weights,
    const arma::vec & steps, const int inverseMode, const bool verbose, const bool reml,
    const EMConditionalFit & conditionalFit) {
  const int terms = descriptors.size();
  const control::Transport controls(descriptors);
  control::Proposal proposal;
  proposal.blend = controls.blend("Unknown EM direction-blending method.");
  const std::string method = control::blendName(proposal.blend);
  const bool parameterExpanded = proposal.blend == control::Blend::PXEM;
  proposal.adaptiveRecovery = controls.logical("adaptiveRecovery");
  const bool adaptive = proposal.adaptiveRecovery;
  const double scoreTolerance = std::sqrt(std::max(likelihoodTolerance, parameterTolerance));
  AIRecoveryPolicy recoveryPolicy;
  int previousHalvings = 0, recoverySteps = 0, acceptedRecoveries = 0;
  bool previousStepWasAI = false;
  std::vector<std::string> recoveryReasons;
  std::vector<int> recoveryIterations;
  int regularizedAttempts = 0, regularizedAccepted = 0;
  std::vector<int> regularizedIterations;
  std::vector<double> regularizedRidges;
  std::string stopReason = "iteration-limit";
  double currentScoreNorm = std::numeric_limits<double>::infinity();
  const bool lowRank = hasZeroSpecific(descriptors);
  const std::string reducedMode = descriptors.hasAttribute("reducedAI")
    ? Rcpp::as<std::string>(descriptors.attr("reducedAI")) : "off";
  const bool retainInformation = !descriptors.hasAttribute("covarianceInformation") ||
    Rcpp::as<std::string>(descriptors.attr("covarianceInformation")) == "full";
  if(reducedMode == "force" && lowRank) {
    Rcpp::stop("reducedAI='force' does not support latent reduced-rank fitting.");
  }
  std::vector<arma::uword> offsets(terms + 1, 0);
  for(int term = 0; term < terms; ++term) {
    const Rcpp::List descriptor = descriptors[term];
    offsets[term + 1] = offsets[term] + Rcpp::as<arma::vec>(descriptor["par"]).n_elem;
  }
  const arma::uword count = offsets.back();
  const arma::vec stepCaps = covarianceWorkingStepCaps(descriptors);
  const double logFloor = std::log(std::max(1.0e-8, inverseTolerance));
  arma::vec varianceLower = covarianceVarianceLowerBounds(descriptors, inverseTolerance);
  arma::vec parameters(count), freeFlags(count, arma::fill::zeros);
  for(int term = 0; term < terms; ++term) {
    const Rcpp::List descriptor = descriptors[term];
    const arma::vec par = Rcpp::as<arma::vec>(descriptor["par"]);
    const Rcpp::LogicalVector free = descriptor["free"];
    parameters.subvec(offsets[term], offsets[term + 1] - 1) = par;
    for(arma::uword j = 0; j < par.n_elem; ++j) {
      freeFlags(offsets[term] + j) = free[j] == TRUE;
    }
  }
  const arma::uvec free = arma::find(freeFlags > 0.5);
  for(arma::uword j = 0; j < count; ++j) {
    varianceLower(j) = std::min(varianceLower(j), parameters(j));
  }
  const VarianceParameterMap map = readVarianceParameterMap(descriptors, count);
  arma::mat tangent(count, free.n_elem, arma::fill::zeros);
  for(arma::uword j = 0; j < free.n_elem; ++j) { tangent(free(j), j) = 1.0; }
  if(map.active) {
    tangent = map.T;
    for(arma::uword j = 0; j < count; ++j) {
      if(freeFlags(j) < 0.5 && arma::any(tangent.row(j) != 0.0)) {
        Rcpp::stop("The EM constraint map contains a free direction through a fixed parameter.");
      }
    }
    const arma::uvec used = arma::find(arma::sum(arma::abs(tangent), 0).t() > 0.0);
    tangent = tangent.cols(used);
  }
  const bool reducedAI = reducedMode != "off" && !lowRank && map.active &&
    tangent.n_cols < count;
  int reducedEvaluations = 0, fullInformationEvaluations = 0;
  double reducedSensitivityColumns = 0.0, skippedSensitivityColumns = 0.0;
  double reducedSensitivityEntries = 0.0, reducedAssemblyEntries = 0.0;
  std::vector<bool> variableTerm(terms, false);
  for(int term = 0; term < terms; ++term) {
    variableTerm[term] = arma::accu(arma::abs(
      tangent.rows(offsets[term], offsets[term + 1] - 1))) > 0.0;
  }
  const double memoryMB = descriptors.hasAttribute("denseMemoryMB")
    ? Rcpp::as<double>(descriptors.attr("denseMemoryMB")) : 1000.0;
  const double budget = memoryMB * 1.0e6 / sizeof(double);
  auto requireWorkspace = [&](double entries) {
    if(entries > budget) {
      Rcpp::stop("Exact EM moments exceed denseMemoryMB. Increase the engine memory budget; "
                 "no approximate EM moments or diagonal AI substitution are used.");
    }
  };
  CovarianceFactorEngine engine;
  int aiEvaluations = 0, emEvaluations = 0, qEvaluations = 0;
  int closedFormSteps = 0, generalizedSteps = 0, numericalConvergences = 0;
  int momentSolves = 0, likelihoodEvaluations = 0, lineHalvings = 0;
  int emFallbacks = 0;
  int pxEvaluations = 0, pxEMFallbacks = 0;
  int latentMomentPreparations = 0;
  arma::uword latentMomentBases = 0;
  arma::uword randomMomentPairs = 0, randomOffDiagonalPairsSkipped = 0;
  int scalarResidualTraceSteps = 0, scalarResidualContractions = 0;
  int invalidTrials = 0;
  std::map<std::string, int> invalidTrialReasons;
  arma::uword sensitivityColumns = 0;
  double coefficientEvaluations = 0.0;
  arma::vec profileSeconds, recordedProfile;
  arma::mat profileIterations;
  std::map<std::string, double> profileCounters;
  std::map<std::string, double> engineCounters, gpuCounters;
  Rcpp::RObject ldltCache = R_NilValue, cholmodCache = R_NilValue;
  struct Timing { double calls = 0.0, elapsed = 0.0, cpu = 0.0; };
  std::map<std::string, Timing> factorTimings;
  arma::vec likelihoodHistory, weightHistory, qGains, mstepConvergenceHistory;
  arma::vec pxGains, pxConvergenceHistory;
  arma::ivec pxEvaluationHistory;
  arma::mat expansionHistory;
  arma::ivec mstepEvaluationHistory;
  arma::mat monitor, norms;
  MMEProgress progress(terms);
  Rcpp::List result;
  bool converged = false, stalled = false;
  bool iterationNeedsAI = false;
  Rcpp::List pendingAI;
  bool havePendingAI = false, finalInformationReady = false;
  int nativeAIEvaluations = 0, latentAIEvaluations = 0, acceptedAIFitReuses = 0;
  double previousLikelihood = -std::numeric_limits<double>::infinity();
  auto observe = [&](Rcpp::List fit) {
    const Rcpp::List diagnostics = fit["engineDiagnostics"];
    if(fit.containsElementNamed(".aiTangent")) {
      ++reducedEvaluations;
      const arma::mat basis = Rcpp::as<arma::mat>(fit[".aiTangent"]);
      reducedSensitivityColumns += Rcpp::as<double>(diagnostics["nativeAISensitivityColumns"]);
      skippedSensitivityColumns += count - basis.n_cols;
      reducedSensitivityEntries = std::max(reducedSensitivityEntries,
        Rcpp::as<double>(diagnostics["nativeAISensitivityEntries"]));
      reducedAssemblyEntries = std::max(reducedAssemblyEntries,
        Rcpp::as<double>(diagnostics["nativeAIAssemblyEntries"]));
    } else if(Rcpp::as<double>(diagnostics["nativeAISensitivityColumns"]) > 0.0) {
      ++fullInformationEvaluations;
    }
    coefficientEvaluations += Rcpp::as<double>(diagnostics["coefficientEvaluations"]);
    for(const char * name : {"CsymbolicAnalyses", "DsymbolicAnalyses",
                            "DselectedTopologyBuilds", "DselectedTopologyReuses",
                            "cholmodDSelectedBuilds",
                            "gridFactorizations", "gridFactorReuses",
                            "factorEvaluationReuses", "sparseSpatialFactorizationReuses",
                            "emInverseTraceCalls", "emInverseTraceFallbackPanels",
                            "diagonalResidualTraceCalls", "diagonalResidualTraceSolveColumns",
                            "emDesignTraceCalls"}) {
      if(diagnostics.containsElementNamed(name)) {
        engineCounters[name] += Rcpp::as<double>(diagnostics[name]);
      }
    }
    if(fit.containsElementNamed(".ldltCache") && !Rf_isNull(fit[".ldltCache"])) {
      ldltCache = fit[".ldltCache"];
    }
    if(fit.containsElementNamed(".cholmodCache") && !Rf_isNull(fit[".cholmodCache"])) {
      cholmodCache = fit[".cholmodCache"];
    }
    if(fit.containsElementNamed("gpuDiagnostics")) {
      const Rcpp::List gpu = fit["gpuDiagnostics"];
      for(const char * name : {"choleskyCalls", "triangularSolveCalls", "productCalls"}) {
        gpuCounters[name] += Rcpp::as<double>(gpu[name]);
      }
    }
    if(diagnostics.containsElementNamed("profile") && !Rf_isNull(diagnostics["profile"])) {
      const Rcpp::List profile = diagnostics["profile"];
      if(profile.containsElementNamed("seconds")) {
        const arma::vec seconds = Rcpp::as<arma::vec>(profile["seconds"]);
        if(profileSeconds.is_empty()) { profileSeconds.zeros(seconds.n_elem); }
        profileSeconds += seconds;
        const Rcpp::List counters = profile["counters"];
        for(const char * name : {"numericFactorCalls", "numericFactorSeconds",
                                "sensitivityRHSColumns", "traceFallbackSolves", "traceBatchSolves",
                                "traceDerivativePreparationSeconds",
                                "traceRandomPrecisionDerivativeSeconds",
                                "traceFactorCovarianceDerivativeSeconds",
                                "traceFactorPrecisionDerivativeSeconds",
                                "traceRandomPrecisionOtherSeconds",
                                "traceResidualPrecisionDerivativeSeconds",
                                "traceGridDerivativeBasisSeconds",
                                "traceResidualDerivativeBasisSeconds",
                                "traceDerivativeOtherSeconds",
                                "traceSelectedInverseSeconds",
                                "scoreTraceFallbackSeconds", "aiAssemblyWallSeconds",
                                "aiRandomFirstTermSeconds", "aiRandomCurvatureSeconds",
                                "aiResidualBlockSeconds"}) {
          if(counters.containsElementNamed(name)) {
            profileCounters[name] += Rcpp::as<double>(counters[name]);
          }
        }
      }
    }
    const Rcpp::List factorization = fit["factorizationDiagnostics"];
    const Rcpp::List timings = factorization["timings"];
    const Rcpp::CharacterVector names = timings.names();
    for(int j = 0; j < timings.size(); ++j) {
      const Rcpp::List timing = timings[j];
      Timing & total = factorTimings[Rcpp::as<std::string>(names[j])];
      total.calls += Rcpp::as<double>(timing["calls"]);
      total.elapsed += Rcpp::as<double>(timing["elapsedSeconds"]);
      total.cpu += Rcpp::as<double>(timing["cpuSeconds"]);
    }
    return fit;
  };
  auto covarianceAt = [&](const arma::vec & point) {
    Rcpp::List covariance = Rcpp::clone(descriptors);
    covariance.attr("directionBlend") = false;
    covariance.attr("retainLDLTCache") = true;
    covariance.attr("retainCholmodCache") = true;
    if(!ldltCache.isNULL()) { covariance.attr("ldltCache") = ldltCache; }
    if(!cholmodCache.isNULL()) { covariance.attr("cholmodCache") = cholmodCache; }
    for(int term = 0; term < terms; ++term) {
      Rcpp::List descriptor = covariance[term];
      descriptor["par"] = point.subvec(offsets[term], offsets[term + 1] - 1);
      covariance[term] = descriptor;
    }
    return covariance;
  };
  auto validCandidate = [&](const arma::vec & point) {
    if(!point.is_finite()) { return false; }
    Rcpp::List candidateCovariance = descriptors;
    if(lowRank) {
      candidateCovariance = latentCovarianceModel(descriptors, Z, Zind, point, offsets).covariance;
    }
    for(arma::uword j = 0; j < count; ++j) {
      if(point(j) < varianceLower(j) - 1.0e-10) {
        return false;
      }
    }
    for(int term = 0; term < terms; ++term) {
      const double scale = std::exp(point(offsets[term]));
      if(!std::isfinite(scale) || scale <= 0.0) { return false; }
      const Rcpp::List descriptor = candidateCovariance[term];
      const Rcpp::List factors = descriptor["factors"];
      const arma::vec par = point.subvec(offsets[term], offsets[term + 1] - 1);
      for(int axis = 0; axis < factors.size(); ++axis) {
        const Rcpp::List factor = factors[axis];
        const int first = Rcpp::as<int>(factor["par_start"]) - 1;
        const int last = Rcpp::as<int>(factor["par_end"]) - 1;
        const arma::vec local = last >= first ? arma::vec(par.subvec(first, last)) : arma::vec();
        const int sections = engine.factorSections(factor);
        for(int section = 0; section < std::max(sections, 1); ++section) {
          const arma::mat value = sections ? engine.evalSectionFactor(factor, local, section)
                                          : engine.evalFactor(factor, local);
          arma::mat upper;
          if(!value.is_finite() || !arma::chol(upper, value)) { return false; }
        }
      }
    }
    return true;
  };
  auto evaluate = [&](const arma::vec & point, int ci,
                      const RRConditionalCallback * callback,
                      bool nativeInformation = false) {
    ++likelihoodEvaluations;
    Rcpp::List covariance = covarianceAt(point);
    Rcpp::List design = Z;
    arma::vec index = Zind;
    if(lowRank) {
      const LatentCovarianceModel latent = latentCovarianceModel(covariance, Z, Zind, point, offsets);
      covariance = latent.covariance;
      design = latent.design;
      index = latent.index;
    }
    const bool nativeAI = !lowRank && (nativeInformation || (callback && iterationNeedsAI));
    covariance.attr("rrConditional") = !nativeAI;
    covariance.attr("emInformationOnly") = nativeAI;
    if(nativeAI && reducedAI) {
      covariance.attr("aiTangent") = tangent;
      covariance.attr("covarianceInformation") = "none";
    }
    covariance.attr("emMomentsOnly") = !nativeAI && (!callback || !iterationNeedsAI);
    // Only native-information fits can be returned without a final refit.
    covariance.attr("omitCoefficientMatrix") = !nativeInformation;
    if(callback && !nativeAI) {
      for(int term = 0; term < terms; ++term) {
        Rcpp::List descriptor = covariance[term];
        const arma::vec par = Rcpp::as<arma::vec>(descriptor["par"]);
        descriptor["free"] = Rcpp::LogicalVector(par.n_elem, false);
        covariance[term] = descriptor;
      }
    }
    // Conditional likelihood trials need no score, traces, or AI.
    return observe(conditionalFit(design, index, covariance, ci, callback));
  };

  for(int iteration = 0; iteration < iterations; ++iteration) {
    Rcpp::checkUserInterrupt();
    double weight = adaptive ? 0.0 : weights(std::min<arma::uword>(iteration, weights.n_elem - 1));
    const double damping = steps(std::min<arma::uword>(iteration, steps.n_elem - 1));
    iterationNeedsAI = weight < 1.0 && (adaptive || iteration < iterations - 1) && tangent.n_cols > 0;
    arma::vec emProposal = parameters, aiDirection(count, arma::fill::zeros);
    arma::vec ordinaryProposal = parameters;
    arma::vec logExpansion(terms - 1, arma::fill::zeros);
    arma::vec score(count, arma::fill::zeros);
    double qGain = 0.0;
    double mstepConverged = NA_REAL;
    int mstepEvaluations = 0;
    double pxGain = 0.0, pxConverged = NA_REAL;
    int pxObjectiveEvaluations = 0;
    bool closedFormAccepted = false;
    bool recoveryMomentsReady = false;
    bool regularizedAI = false;
    arma::vec cappedAIDirection(count, arma::fill::zeros);
    std::string recoveryReason;
    const auto requestRecovery = [&](const std::string & reason) {
      weight = 1.0;
      recoveryReason = reason;
      recoveryPolicy.reset();
    };
    const RRConditionalCallback callback = [&](const RRConditionalOperators & operators) {
      if((iteration == iterations - 1 && !adaptive) || tangent.n_cols == 0) { return; }
      if(weight < 1.0) {
          RRAIResult ai;
          if(!lowRank) {
            ai.score = operators.score;
            ai.information = operators.information;
            ai.sensitivityColumns = operators.sensitivityColumns;
            ++nativeAIEvaluations;
          } else {
            ai = reducedRankInformation(operators, X, Z, Zind, Ai,
              residualBlock, residualIndex, H, useH, descriptors, parameters, offsets, reml);
            ++latentAIEvaluations;
          }
          ++aiEvaluations;
          sensitivityColumns += ai.sensitivityColumns;
          score = ai.score;
          if(adaptive) {
            if(!score.is_finite() || !ai.information.is_finite()) {
              Rcpp::stop("Non-finite AI score or information in adaptive recovery.");
            }
            currentScoreNorm = projectedEMScoreNorm(tangent, parameters, varianceLower,
              score, ai.information, inverseTolerance, operators.informationTangent);
            if(iteration == iterations - 1) { return; }
            const FeasibleAIDirection proposal = feasibleAIDirection(
              ai.information, score, operators.informationTangent, tangent,
              parameters, varianceLower, stepCaps);
            regularizedAttempts += proposal.regularizedAttempts;
            if(!proposal.reliable) {
              requestRecovery("information-solve");
              return;
            }
            aiDirection = proposal.direction;
            cappedAIDirection = proposal.cappedDirection;
            regularizedAI = proposal.ridge > 0.0;
            if(regularizedAI) {
              regularizedIterations.push_back(iteration + 1);
              regularizedRidges.push_back(proposal.ridge);
            }
            return;
          }
          arma::mat feasible = tangent;
          for(int activeSet = 0; activeSet <= terms && feasible.n_cols > 0; ++activeSet) {
            const arma::mat information = aiInformationInTangent(
              ai.information, operators.informationTangent, feasible);
            const arma::vec projected = feasible.t() * score;
            arma::vec direction;
            if(!arma::solve(direction, information, projected, arma::solve_opts::likely_sympd)) {
              direction = arma::pinv(information) * projected;
            }
            if(!direction.is_finite()) { Rcpp::stop("Non-finite AI direction in EM blending."); }
            aiDirection = feasible * direction;
            const arma::mat next = varianceFloorTangent(
              feasible, parameters, varianceLower, aiDirection);
            if(next.n_cols == feasible.n_cols) { break; }
            feasible = next;
            aiDirection.zeros();
          }
      }
      // Native AI consumption can omit moment operators; recovery obtains them
      // through the same conditional backend only after a trigger.
      if(adaptive && !operators.solve) { return; }
      if(weight == 0.0) { return; }
      recoveryMomentsReady = true;
      ++emEvaluations;
      const double logScale = std::log(operators.scale);
      LatentCovarianceModel current = latentCovarianceModel(descriptors, Z, Zind, parameters, offsets);
      std::vector<arma::mat> randomMoments(terms - 1);
      std::vector<arma::uword> coefficientOffsets(terms, X.n_cols);
      double storedRandomEntries = 0.0;
      arma::uword maximumSolveWidth = 0;
      const Rcpp::List residualDescriptor = descriptors[terms - 1];
      const bool scalarResidual = Rcpp::as<int>(residualDescriptor["dim"]) == 1 &&
        offsets[terms] - offsets[terms - 1] == 1 &&
        (variableTerm.back() || lowRank || parameterExpanded) &&
        static_cast<bool>(operators.posteriorRandomTraces);
      double inversePriorTrace = 0.0;
      for(int term = 0; term < terms - 1; ++term) {
        const arma::uword q = current.maps[term].n_cols;
        const arma::uword m = current.widths[term];
        coefficientOffsets[term + 1] = coefficientOffsets[term] + q * m;
        if(!variableTerm[term] && !parameterExpanded && !scalarResidual) { continue; }
        storedRandomEntries += static_cast<double>(q) * q;
        maximumSolveWidth = std::max(maximumSolveWidth, std::min<arma::uword>(32, m));
        requireWorkspace(5.0 * storedRandomEntries + 2.0 * m * q +
          4.0 * operators.C.n_rows * maximumSolveWidth);
        arma::mat means(m, q);
        for(arma::uword j = 0; j < q; ++j) {
          means.col(j) = operators.mean.subvec(coefficientOffsets[term] + j * m,
                                               coefficientOffsets[term] + (j + 1) * m - 1);
        }
        const arma::sp_mat precision = Rcpp::as<arma::sp_mat>(Ai[term]);
        const Rcpp::List descriptor = current.covariance[term];
        const Rcpp::List factors = descriptor["factors"];
        bool diagonal = true;
        for(int axis = 0; axis < factors.size(); ++axis) {
          const Rcpp::List factor = factors[axis];
          diagonal = diagonal && factor.containsElementNamed("structurally_diagonal") &&
            Rcpp::as<bool>(factor["structurally_diagonal"]);
        }
        randomMomentPairs += diagonal ? q : q * (q + 1) / 2;
        if(diagonal) { randomOffDiagonalPairsSkipped += q * (q - 1) / 2; }
        arma::mat moment(q, q, arma::fill::zeros);
        const arma::mat related = precision * means;
        if(diagonal) {
          moment.diag() = arma::sum(means % related, 0).t();
        } else {
          moment = means.t() * related;
        }
        arma::mat inverseMoment(q, q, arma::fill::zeros);
        if(operators.posteriorRandomTraces) {
          inverseMoment = operators.posteriorRandomTraces(term, diagonal);
        } else for(arma::uword j = 0; j < q; ++j) {
          for(arma::uword first = 0; first < m; first += 32) {
            const arma::uword panel = std::min<arma::uword>(32, m - first);
            arma::mat basis(operators.C.n_rows, panel, arma::fill::zeros);
            for(arma::uword k = 0; k < panel; ++k) {
              basis(coefficientOffsets[term] + j * m + first + k, k) = 1.0;
            }
            const arma::mat solved = operators.posteriorSolve(basis);
            ++momentSolves;
            for(arma::uword k = 0; k < panel; ++k) {
              for(arma::sp_mat::const_col_iterator entry = precision.begin_col(first + k);
                  entry != precision.end_col(first + k); ++entry) {
                for(arma::uword i = diagonal ? j : 0; i <= j; ++i) {
                  inverseMoment(i, j) += *entry * solved(coefficientOffsets[term] + i * m + entry.row(), k);
                }
              }
            }
          }
        }
        inverseMoment = arma::symmatu(inverseMoment);
        randomMoments[term] = moment + inverseMoment;
        if(scalarResidual) {
          arma::vec par = parameters.subvec(offsets[term], offsets[term + 1] - 1);
          par(0) -= logScale;
          const arma::mat covariance = engine.evaluateDescriptor(descriptor, par);
          inversePriorTrace += diagonal
            ? arma::accu(inverseMoment.diag() / covariance.diag())
            : arma::trace(arma::solve(covariance, inverseMoment));
        }
      }
      struct ResidualMoment {
        arma::uvec index;
        std::vector<arma::uvec> rows;
        arma::mat moment;
      };
      std::vector<ResidualMoment> residualMoments;
      std::map<std::vector<arma::uword>, std::size_t> groups;
      double storedEntries = 0.0;
      const bool residualNeeded = variableTerm.back() || lowRank || parameterExpanded;
      for(arma::uword block = 1; residualNeeded && block <= residualBlock.max(); ++block) {
        const arma::uvec rows = arma::find(residualBlock == block);
        if(rows.is_empty()) { continue; }
        const arma::uvec index = residualIndex.elem(rows) - 1;
        const std::vector<arma::uword> key(index.begin(), index.end());
        auto found = groups.find(key);
        if(found == groups.end()) {
          storedEntries += static_cast<double>(rows.n_elem) * rows.n_elem;
          requireWorkspace(5.0 * (storedRandomEntries + storedEntries) +
            4.0 * operators.C.n_rows * std::min<arma::uword>(32, rows.n_elem));
          ResidualMoment state;
          state.index = index;
          residualMoments.push_back(std::move(state));
          found = groups.emplace(key, residualMoments.size() - 1).first;
        }
        residualMoments[found->second].rows.push_back(rows);
      }
      for(const ResidualMoment & state : residualMoments) {
        const arma::uword width = state.index.n_elem == 1
          ? std::min<std::size_t>(32, state.rows.size())
          : std::min<arma::uword>(32, state.index.n_elem);
        maximumSolveWidth = std::max(maximumSolveWidth, width);
      }
      requireWorkspace(5.0 * (storedRandomEntries + storedEntries) +
        4.0 * operators.C.n_rows * maximumSolveWidth);
      arma::vec transformedY = operators.y;
      arma::sp_mat transformedDesign = operators.W;
      if(useH) {
        transformedY = operators.residualTransform * transformedY;
        transformedDesign = operators.residualTransform * transformedDesign;
      }
      auto selectRows = [](const arma::sp_mat & design, const arma::uvec & rows) {
        arma::sp_mat local(rows.n_elem, design.n_cols);
        for(arma::uword row = 0; row < rows.n_elem; ++row) {
          for(arma::sp_mat::const_row_iterator entry = design.begin_row(rows(row));
              entry != design.end_row(rows(row)); ++entry) {
            local(row, entry.col()) = *entry;
          }
        }
        return local;
      };
      auto residualMoment = [&](const ResidualMoment & state, const arma::sp_mat & design) {
        arma::mat moment(state.index.n_elem, state.index.n_elem, arma::fill::zeros);
        if(scalarResidual && !lowRank && state.rows.size() == transformedY.n_elem &&
           operators.posteriorTrace) {
          double entries = 0.0;
          for(arma::uword row = 0; row < design.n_rows; ++row) {
            double width = 0.0;
            for(arma::sp_mat::const_row_iterator entry = design.begin_row(row);
                entry != design.end_row(row); ++entry) { ++width; }
            entries += width * width;
          }
          const double available = budget - 5.0 * (storedRandomEntries + storedEntries);
          if(entries <= WorkspacePolicy::sparseEntryLimit() && 4.0 * entries <= available) {
            const arma::vec error = transformedY - design * operators.mean;
            const arma::sp_mat gram = design.t() * design;
            moment(0, 0) = arma::dot(error, error) + operators.posteriorTrace(gram);
            ++scalarResidualContractions;
            return moment;
          }
        }
        if(state.index.n_elem == 1) {
          for(std::size_t first = 0; first < state.rows.size(); first += 32) {
            const arma::uword panel = std::min<std::size_t>(32, state.rows.size() - first);
            arma::uvec rows(panel);
            for(arma::uword j = 0; j < panel; ++j) { rows(j) = state.rows[first + j](0); }
            const arma::sp_mat local = selectRows(design, rows);
            const arma::vec error = transformedY.elem(rows) - local * operators.mean;
            const arma::mat solved = operators.posteriorSolve(arma::mat(local.t()));
            ++momentSolves;
            moment(0, 0) += arma::dot(error, error) + arma::accu(arma::mat(local) % solved.t());
          }
          return moment;
        }
        for(const arma::uvec & rows : state.rows) {
          const arma::sp_mat localSparse = selectRows(design, rows);
          const arma::vec error = transformedY.elem(rows) - localSparse * operators.mean;
          moment += error * error.t();
          for(arma::uword first = 0; first < rows.n_elem; first += 32) {
            const arma::uword panel = std::min<arma::uword>(32, rows.n_elem - first);
            const arma::mat rhs = arma::mat(localSparse.rows(first, first + panel - 1).t());
            const arma::mat solved = operators.posteriorSolve(rhs);
            ++momentSolves;
            moment.cols(first, first + panel - 1) += localSparse * solved;
          }
        }
        return arma::mat(0.5 * (moment + moment.t()));
      };
      if(!lowRank) {
        for(ResidualMoment & state : residualMoments) {
          if(scalarResidual && state.index.n_elem == 1 &&
             state.rows.size() == transformedY.n_elem) {
            arma::vec par = parameters.subvec(offsets[terms - 1], count - 1);
            par(0) -= logScale;
            const arma::mat covariance = engine.evaluateDescriptor(residualDescriptor, par);
            const arma::vec error = transformedY - transformedDesign * operators.mean;
            state.moment.set_size(1, 1);
            // REML integrates fixed effects; ML conditions on their current GLS estimate.
            state.moment(0, 0) = arma::dot(error, error) +
              covariance(0, 0) * (operators.C.n_rows - (reml ? 0 : X.n_cols) - inversePriorTrace);
            ++scalarResidualTraceSteps;
          } else {
            state.moment = residualMoment(state, transformedDesign);
          }
        }
      }
      struct ExpansionMoment {
        std::vector<arma::mat> linear;
        std::vector<arma::mat> quadratic;
      };
      std::vector<ExpansionMoment> expansionMoments;
      struct LatentMomentBasis { int term; int row; int column; };
      std::vector<LatentMomentBasis> latentBasis;
      std::vector<arma::mat> latentMoments;
      double latentMomentEntries = 0.0;
      if(lowRank) {
        latentBasis.push_back({-1, -1, -1});
        for(int term = 0; term < terms - 1; ++term) {
          const Rcpp::List descriptor = descriptors[term];
          const Rcpp::List factors = descriptor["factors"];
          bool latent = false;
          arma::mat support(1, 1, arma::fill::ones);
          for(int axis = 0; axis < factors.size(); ++axis) {
            const Rcpp::List factor = factors[axis];
            const int dimension = Rcpp::as<int>(factor["dim"]);
            arma::mat axisSupport = arma::eye<arma::mat>(dimension, dimension);
            if(zeroSpecificFactor(factor)) {
              latent = true;
              axisSupport.zeros(dimension, Rcpp::as<int>(factor["order"]));
              const bool rrMetadata = factor.containsElementNamed("rr_nload");
              const Rcpp::IntegerVector rows = factor[rrMetadata ? "rr_row" : "fa_row"];
              const Rcpp::IntegerVector columns = factor[rrMetadata ? "rr_col" : "fa_col"];
              for(int j = 0; j < rows.size(); ++j) {
                axisSupport(rows[j] - 1, columns[j] - 1) = 1.0;
              }
            }
            support = arma::kron(support, axisSupport);
          }
          if(latent) {
            for(arma::uword column = 0; column < support.n_cols; ++column) {
              for(arma::uword row = 0; row < support.n_rows; ++row) {
                if(support(row, column) != 0.0) {
                  latentBasis.push_back({term, static_cast<int>(row), static_cast<int>(column)});
                }
              }
            }
          } else {
            latentBasis.push_back({term, -1, -1});
          }
        }
        const arma::uword bases = latentBasis.size();
        latentMomentEntries = static_cast<double>(bases) * bases * storedEntries;
        // The design is linear in normalized loading entries. Cache its exact
        // quadratic posterior moments once, not coefficient solves for every Q probe.
        requireWorkspace(5.0 * (storedRandomEntries + storedEntries) +
          latentMomentEntries +
          4.0 * bases * operators.C.n_rows * maximumSolveWidth +
          4.0 * bases * (operators.W.n_nonzero + operators.W.n_rows));
        std::vector<arma::sp_mat> designs;
        for(const LatentMomentBasis & basis : latentBasis) {
          arma::sp_mat design(operators.W.n_rows, operators.W.n_cols);
          if(basis.term < 0) {
            if(X.n_cols > 0) { design.cols(0, X.n_cols - 1) = X; }
          } else if(basis.row < 0) {
            const arma::uword first = coefficientOffsets[basis.term];
            const arma::uword last = coefficientOffsets[basis.term + 1] - 1;
            design.cols(first, last) = operators.W.cols(first, last);
          } else {
            const arma::uvec blocks = arma::find(Zind == basis.term + 1);
            const arma::uword first = coefficientOffsets[basis.term] +
              basis.column * current.widths[basis.term];
            design.cols(first, first + current.widths[basis.term] - 1) =
              Rcpp::as<arma::sp_mat>(Z[blocks(basis.row)]);
          }
          designs.push_back(useH
            ? arma::sp_mat(operators.residualTransform * design) : design);
        }
        for(const ResidualMoment & state : residualMoments) {
          const arma::uword dimension = state.index.n_elem;
          arma::mat polynomial(dimension * dimension, bases * bases, arma::fill::zeros);
          auto accumulate = [&](const arma::uvec & rows, bool scalarBatch) {
            std::vector<arma::sp_mat> local;
            arma::mat errors(rows.n_elem, bases);
            for(arma::uword j = 0; j < bases; ++j) {
              local.push_back(selectRows(designs[j], rows));
              errors.col(j) = -local[j] * operators.mean;
            }
            errors.col(0) += transformedY.elem(rows);
            for(arma::uword j = 0; j < bases; ++j) {
              for(arma::uword first = 0; first < rows.n_elem; first += 32) {
                const arma::uword panel = std::min<arma::uword>(32, rows.n_elem - first);
                const arma::mat solved = operators.posteriorSolve(
                  arma::mat(local[j].rows(first, first + panel - 1).t()));
                ++momentSolves;
                for(arma::uword k = 0; k < bases; ++k) {
                  arma::mat contraction = local[k] * solved;
                  contraction += errors.col(k) * errors.col(j).subvec(first, first + panel - 1).t();
                  if(scalarBatch) {
                    polynomial(0, k + j * bases) += arma::trace(contraction);
                  } else {
                    for(arma::uword column = 0; column < panel; ++column) {
                      polynomial.col(k + j * bases).subvec(
                        (first + column) * dimension, (first + column + 1) * dimension - 1) +=
                        contraction.col(column);
                    }
                  }
                }
              }
            }
          };
          if(dimension == 1) {
            for(std::size_t first = 0; first < state.rows.size(); first += 32) {
              const arma::uword panel = std::min<std::size_t>(32, state.rows.size() - first);
              arma::uvec rows(panel);
              for(arma::uword j = 0; j < panel; ++j) { rows(j) = state.rows[first + j](0); }
              accumulate(rows, true);
            }
          } else {
            for(const arma::uvec & rows : state.rows) { accumulate(rows, false); }
          }
          latentMoments.push_back(std::move(polynomial));
        }
        ++latentMomentPreparations;
        latentMomentBases = std::max(latentMomentBases, bases);
      }
      auto expandedObjective = [&](const arma::vec & point, const arma::vec & expansion) {
        ++qEvaluations;
        if(!point.is_finite() || !expansion.is_finite()) {
          return std::numeric_limits<double>::infinity();
        }
        const bool expanded = !expansion.is_empty();
        const arma::vec amplitudes = arma::exp(expansion);
        if(!amplitudes.is_finite() || arma::any(amplitudes <= 0.0)) {
          return std::numeric_limits<double>::infinity();
        }
        for(arma::uword j = 0; j < count; ++j) {
          if(point(j) < varianceLower(j) - 1.0e-10) {
            return std::numeric_limits<double>::infinity();
          }
        }
        LatentCovarianceModel candidate;
        candidate.covariance = descriptors;
        arma::vec latentProducts;
        if(lowRank) {
          candidate.covariance = current.covariance;
          std::vector<arma::mat> maps;
          for(int term = 0; term < terms - 1; ++term) {
            maps.push_back(reducedRankCoefficientMap(
              Rcpp::as<Rcpp::List>(descriptors[term]),
              point.subvec(offsets[term], offsets[term + 1] - 1)));
          }
          arma::vec coefficients(latentBasis.size(), arma::fill::ones);
          for(std::size_t j = 1; j < latentBasis.size(); ++j) {
            const LatentMomentBasis & basis = latentBasis[j];
            coefficients(j) = (expanded ? amplitudes(basis.term) : 1.0) *
              (basis.row < 0 ? 1.0 : maps[basis.term](basis.row, basis.column));
          }
          latentProducts = arma::vectorise(coefficients * coefficients.t());
        }
        double value = 0.0;
        for(int term = 0; term < terms - 1; ++term) {
          if(!variableTerm[term] && !parameterExpanded) { continue; }
          arma::vec par = point.subvec(offsets[term], offsets[term + 1] - 1);
          par(0) -= logScale;
          if(expanded) { par(0) -= 2.0 * expansion(term); }
          const arma::mat covariance = engine.evaluateDescriptor(
            Rcpp::as<Rcpp::List>(candidate.covariance[term]), par);
          arma::mat upper;
          if(!covariance.is_finite() || !arma::chol(upper, covariance)) {
            return std::numeric_limits<double>::infinity();
          }
          const arma::mat solved = arma::solve(arma::trimatu(upper),
            arma::solve(arma::trimatl(upper.t()), randomMoments[term]));
          value += 2.0 * current.widths[term] * arma::accu(arma::log(upper.diag())) +
            arma::trace(solved);
        }
        arma::vec par = point.subvec(offsets[terms - 1], count - 1);
        par(0) -= logScale;
        const Rcpp::List descriptor = descriptors[terms - 1];
        for(std::size_t group = 0; group < residualMoments.size(); ++group) {
          const ResidualMoment & state = residualMoments[group];
          const arma::mat covariance = engine.evaluateObservedDescriptor(
            descriptor, par, state.index, par.n_elem);
          arma::mat upper;
          if(!covariance.is_finite() || !arma::chol(upper, covariance)) {
            return std::numeric_limits<double>::infinity();
          }
          arma::mat moment = lowRank
            ? arma::reshape(latentMoments[group] * latentProducts, state.index.n_elem, state.index.n_elem)
            : state.moment;
          if(expanded && !lowRank) {
            const ExpansionMoment & polynomial = expansionMoments[group];
            const arma::vec delta = amplitudes - 1.0;
            std::size_t pair = 0;
            for(int j = 0; j < terms - 1; ++j) {
              moment += delta(j) * polynomial.linear[j];
              for(int k = j; k < terms - 1; ++k) {
                moment += delta(j) * delta(k) * polynomial.quadratic[pair++];
              }
            }
          }
          const arma::mat solved = arma::solve(arma::trimatu(upper),
            arma::solve(arma::trimatl(upper.t()), moment));
          value += 2.0 * state.rows.size() * arma::accu(arma::log(upper.diag())) +
            arma::trace(solved);
        }
        return value;
      };
      auto objective = [&](const arma::vec & point) {
        return expandedObjective(point, arma::vec());
      };
      const double baseQ = objective(parameters);
      bool closed = !map.active && !lowRank;
      arma::vec candidate = parameters;
      for(int term = 0; term < terms; ++term) {
        const arma::uword first = offsets[term], last = offsets[term + 1] - 1;
        const arma::uvec localFree = arma::find(freeFlags.subvec(first, last) > 0.5);
        if(localFree.is_empty()) { continue; }
        const Rcpp::List descriptor = descriptors[term];
        arma::vec par = parameters.subvec(first, last);
        par(0) -= logScale;
        if(localFree.n_elem == 1 && localFree(0) == 0 && !lowRank && !map.active) {
          double contraction = 0.0, dimension = 0.0;
          if(term < terms - 1) {
            const arma::mat covariance = engine.evaluateDescriptor(descriptor, par);
            contraction = arma::trace(arma::solve(covariance, randomMoments[term]));
            dimension = static_cast<double>(current.widths[term]) * covariance.n_rows;
          } else {
            for(const ResidualMoment & state : residualMoments) {
              const arma::mat covariance = engine.evaluateObservedDescriptor(
                descriptor, par, state.index, par.n_elem);
              contraction += arma::trace(arma::solve(covariance, state.moment));
              dimension += state.index.n_elem * state.rows.size();
            }
          }
          candidate(first) += std::log(contraction / dimension);
          candidate(first) = std::max(candidate(first), std::min(parameters(first), logFloor));
        } else if(localFree.n_elem == last - first + 1 && !lowRank && !map.active) {
          const Rcpp::List factors = descriptor["factors"];
          if(factors.size() != 1) { closed = false; continue; }
          const Rcpp::List factor = factors[0];
          const std::string model = factor.containsElementNamed("model")
            ? Rcpp::as<std::string>(factor["model"]) : "";
          if(model != "us" && model != "diag") { closed = false; continue; }
          arma::mat covariance;
          if(term < terms - 1) {
            covariance = randomMoments[term] / current.widths[term];
          } else {
            const arma::uword dimension = Rcpp::as<int>(descriptor["dim"]);
            requireWorkspace(5.0 * (storedRandomEntries + storedEntries +
              static_cast<double>(dimension) * dimension) +
              4.0 * operators.C.n_rows * maximumSolveWidth);
            covariance.zeros(dimension, dimension);
            arma::vec counts(dimension, arma::fill::zeros);
            bool complete = true;
            for(const ResidualMoment & state : residualMoments) {
              complete = complete && state.index.n_elem == dimension;
              covariance.submat(state.index, state.index) += state.moment;
              counts.elem(state.index) += state.rows.size();
            }
            if(counts.min() <= 0.0 || (model == "us" && !complete)) {
              closed = false;
              continue;
            }
            if(model == "diag") {
              covariance = arma::diagmat(covariance.diag() / counts);
            } else {
              covariance /= counts(0);
            }
          }
          const double scale = covariance(0, 0);
          candidate(first) = std::log(scale) + logScale;
          if(model == "diag") {
            for(arma::uword j = 1; j < covariance.n_rows; ++j) {
              candidate(first + j) = std::log(covariance(j, j) / scale);
            }
          } else {
            arma::mat upper;
            if(!arma::chol(upper, covariance / scale)) { closed = false; continue; }
            const arma::mat lower = upper.t();
            const Rcpp::IntegerVector rows = factor["us_row"], columns = factor["us_col"];
            const Rcpp::LogicalVector diagonal = factor["us_diag"];
            for(int j = 0; j < rows.size(); ++j) {
              const double entry = lower(rows[j] - 1, columns[j] - 1);
              candidate(first + 1 + j) = diagonal[j] ? std::log(entry) : entry;
            }
          }
        } else {
          closed = false;
        }
      }
      const double candidateQ = closed ? objective(candidate) : baseQ;
      if(closed && std::isfinite(candidateQ) && candidateQ <= baseQ) {
        emProposal = candidate;
        closedFormAccepted = true;
        ++closedFormSteps;
      } else {
        requireWorkspace(5.0 * (storedRandomEntries + storedEntries) + latentMomentEntries +
          4.0 * operators.C.n_rows * maximumSolveWidth +
          4.0 * tangent.n_cols * tangent.n_cols);
        const auto reducedObjective = [&](const arma::vec & reduced) {
          return objective(parameters + tangent * reduced);
        };
        const EMMinimum minimum = minimizeEM(arma::vec(tangent.n_cols, arma::fill::zeros),
                                             reducedObjective);
        emProposal = parameters + tangent * minimum.point;
        mstepConverged = minimum.converged;
        mstepEvaluations = minimum.evaluations;
        ++generalizedSteps;
        if(minimum.converged) { ++numericalConvergences; }
      }
      ordinaryProposal = emProposal;
      double proposalQ = objective(emProposal);
      if(parameterExpanded && terms > 1) {
        ++pxEvaluations;
        const arma::uword expansionCount = terms - 1;
        const double polynomialCount = expansionCount +
          static_cast<double>(expansionCount) * (expansionCount + 1) / 2.0;
        const double optimizerDimension = tangent.n_cols + expansionCount;
        requireWorkspace(5.0 * (storedRandomEntries + storedEntries) + latentMomentEntries +
          (lowRank ? 0.0 : (polynomialCount + expansionCount + 4.0) * storedEntries) +
          4.0 * operators.C.n_rows * maximumSolveWidth +
          4.0 * optimizerDimension * optimizerDimension);
        if(!lowRank) {
          expansionMoments.resize(residualMoments.size());
          auto shiftedDesign = [&](int first, int second, double shift) {
            arma::sp_mat design = transformedDesign;
            for(int term = 0; term < terms - 1; ++term) {
              if(term != first && term != second) { continue; }
              design.cols(coefficientOffsets[term], coefficientOffsets[term + 1] - 1) *= 1.0 + shift;
            }
            return design;
          };
          for(std::size_t group = 0; group < residualMoments.size(); ++group) {
            const ResidualMoment & state = residualMoments[group];
            ExpansionMoment & polynomial = expansionMoments[group];
            std::vector<arma::mat> plus(expansionCount);
            for(int j = 0; j < terms - 1; ++j) {
              plus[j] = residualMoment(state, shiftedDesign(j, j, 1.0));
              const arma::mat minus = residualMoment(state, shiftedDesign(j, j, -1.0));
              polynomial.linear.push_back(0.5 * (plus[j] - minus));
            }
            for(int j = 0; j < terms - 1; ++j) {
              const arma::mat diagonal = plus[j] - state.moment - polynomial.linear[j];
              polynomial.quadratic.push_back(diagonal);
              for(int k = j + 1; k < terms - 1; ++k) {
                polynomial.quadratic.push_back(
                  residualMoment(state, shiftedDesign(j, k, 1.0)) -
                  plus[j] - plus[k] + state.moment);
              }
            }
          }
        }
        // Optimize reduced, admissible covariance parameters and expansion jointly.
        // G* = G(phi)/a^2 makes fixed parameters and ties constraints on phi, not a.
        const auto pxObjective = [&](const arma::vec & reduced) {
          const arma::vec point = ordinaryProposal +
            tangent * reduced.head(tangent.n_cols);
          return expandedObjective(point, reduced.tail(expansionCount));
        };
        const EMMinimum minimum = minimizeEM(
          arma::vec(tangent.n_cols + expansionCount, arma::fill::zeros), pxObjective);
        emProposal = ordinaryProposal + tangent * minimum.point.head(tangent.n_cols);
        logExpansion = minimum.point.tail(expansionCount);
        pxGain = 0.5 * (proposalQ - minimum.value);
        proposalQ = minimum.value;
        pxConverged = minimum.converged;
        pxObjectiveEvaluations = minimum.evaluations;
        closedFormAccepted = false;
      }
      qGain = 0.5 * (baseQ - proposalQ);
      if(!std::isfinite(qGain) || qGain < -1.0e-8 * (1.0 + std::abs(baseQ))) {
        Rcpp::stop("The EM M-step decreased its expected complete-data log-likelihood.");
      }
    };
    const auto consumeNativeAI = [&]() {
      RRConditionalOperators operators;
      operators.score = -0.5 * Rcpp::as<arma::vec>(result["dLu"]);
      operators.information = 0.5 * Rcpp::as<arma::mat>(result["InfMat"]);
      if(result.containsElementNamed(".aiTangent")) {
        operators.informationTangent = Rcpp::as<arma::mat>(result[".aiTangent"]);
      }
      const Rcpp::List diagnostics = result["engineDiagnostics"];
      operators.sensitivityColumns = Rcpp::as<double>(diagnostics["nativeAISensitivityColumns"]);
      callback(operators);
    };
    if(havePendingAI) {
      result = pendingAI;
      pendingAI = Rcpp::List();
      havePendingAI = false;
      ++acceptedAIFitReuses;
      if(iterationNeedsAI) { consumeNativeAI(); }
      finalInformationReady = iteration == iterations - 1 || (adaptive && inverseMode == 0);
    } else if(!lowRank && weight == 0.0 &&
              (iterationNeedsAI || iteration == iterations - 1)) {
      finalInformationReady = iteration == iterations - 1 || (adaptive && inverseMode == 0);
      result = evaluate(parameters, finalInformationReady ? inverseMode : 0, nullptr, true);
      if(iterationNeedsAI) { consumeNativeAI(); }
    } else {
      result = evaluate(parameters, 0, &callback);
      finalInformationReady = false;
    }
    if(!profileSeconds.is_empty()) {
      if(recordedProfile.is_empty()) { recordedProfile.zeros(profileSeconds.n_elem); }
      profileIterations.insert_rows(profileIterations.n_rows,
        arma::rowvec((profileSeconds - recordedProfile).t()));
      recordedProfile = profileSeconds;
    }
    const arma::vec likelihoods = Rcpp::as<arma::vec>(result["llik"]);
    const double likelihood = likelihoods(likelihoods.n_elem - 1);
    if(!std::isfinite(likelihood)) { Rcpp::stop("Non-finite accepted likelihood in EM blending."); }
    if(iteration > 0 && likelihood < previousLikelihood - 1.0e-10) {
      Rcpp::stop("An EM proposal increased Q but decreased the observed likelihood; "
                 "the fit was not accepted.");
    }
    if(adaptive && iteration < iterations - 1 && recoveryReason.empty()) {
      const char * reason = recoveryPolicy.observe(previousStepWasAI, previousHalvings,
        likelihood - previousLikelihood, likelihoodTolerance, currentScoreNorm <= scoreTolerance);
      if(reason) { requestRecovery(reason); }
    }
    const auto prepareRecovery = [&]() {
      iterationNeedsAI = false;
      result = evaluate(parameters, 0, &callback);
      finalInformationReady = false;
      ++recoverySteps;
      recoveryReasons.push_back(recoveryReason);
      recoveryIterations.push_back(iteration + 1);
    };
    if(adaptive && weight == 1.0) {
      if(!recoveryMomentsReady) { prepareRecovery(); }
      else {
        ++recoverySteps;
        recoveryReasons.push_back(recoveryReason);
        recoveryIterations.push_back(iteration + 1);
      }
    }
    likelihoodHistory.insert_rows(likelihoodHistory.n_elem, arma::vec({likelihood}));
    weightHistory.insert_rows(weightHistory.n_elem, arma::vec({weight}));
    qGains.insert_rows(qGains.n_elem, arma::vec({qGain}));
    mstepConvergenceHistory.insert_rows(mstepConvergenceHistory.n_elem,
      arma::vec({mstepConverged}));
    mstepEvaluationHistory.insert_rows(mstepEvaluationHistory.n_elem,
      arma::ivec({mstepEvaluations}));
    if(parameterExpanded) {
      pxGains.insert_rows(pxGains.n_elem, arma::vec({pxGain}));
      pxConvergenceHistory.insert_rows(pxConvergenceHistory.n_elem, arma::vec({pxConverged}));
      pxEvaluationHistory.insert_rows(pxEvaluationHistory.n_elem, arma::ivec({pxObjectiveEvaluations}));
      expansionHistory.insert_cols(expansionHistory.n_cols, logExpansion);
    }
    const arma::mat reported = Rcpp::as<arma::mat>(result["monitor"]);
    monitor.insert_cols(monitor.n_cols, reported.col(0));
    arma::vec direction = weight * (emProposal - parameters) + (1.0 - weight) * aiDirection;
    const double relativeStep = arma::norm(direction, 2) / (1.0 + arma::norm(parameters, 2));
    norms.insert_cols(norms.n_cols, arma::vec({arma::norm(direction, 1), arma::norm(score, 1), relativeStep}));
    if(verbose) {
      if(iteration == 0) {
        const Rcpp::List diagnostics = result["engineDiagnostics"];
        Rcpp::Rcout << (adaptive
          ? "MME path: AI-first with on-demand EM/PX-GEM recovery"
          : parameterExpanded
          ? "MME path: scalar PX-GEM + AI direction blending"
          : "MME path: exact EM/GEM + AI direction blending")
          << " [" << mmePathStatus(result, lowRank) << "]\n";
        if(Rcpp::as<bool>(diagnostics["blockSchurActive"])) {
          const std::string coefficientEngine = Rcpp::as<bool>(diagnostics["factorSchurActive"])
            ? "Latent-factor Schur" : Rcpp::as<bool>(diagnostics["blockChainActive"])
            ? "Bordered block-chain" : Rcpp::as<bool>(diagnostics["blockBandedActive"])
            ? "Block-banded Schur" : "Dense block-Schur";
          Rcpp::Rcout << coefficientEngine << " engine active ("
            << Rcpp::as<int>(diagnostics["blockSchurGroups"]) << " random-effect blocks, "
            << Rcpp::as<int>(diagnostics["blockSchurBorder"]) << " border effects).\n";
        }
        progress.header();
      }
    }
    double progressTrust = NA_REAL, progressLine = NA_REAL;
    std::string progressUpdate = regularizedAI ? "AI-reg"
      : MMEProgress::update(weight, parameterExpanded);
    const double progressPrevious = previousLikelihood;
    const auto printProgress = [&]() {
      if(verbose) {
        progress.row(iteration + 1, likelihood, progressPrevious, progressUpdate,
          regularizedAI && weight == 0.0 ? "info-solve"
            : MMEProgress::policy(adaptive, parameterExpanded, weight, recoveryReason),
          progressTrust, progressLine);
      }
    };
    if(tangent.n_cols == 0) {
      converged = true; stopReason = "converged"; printProgress(); break;
    }
    if(iteration > 0 &&
       (adaptive
         ? currentScoreNorm <= scoreTolerance &&
           (relativeStep < parameterTolerance ||
            std::abs(likelihood - previousLikelihood) < likelihoodTolerance)
         : relativeStep < parameterTolerance &&
           std::abs(likelihood - previousLikelihood) < likelihoodTolerance)) {
      converged = true;
      stopReason = "converged";
      printProgress();
      break;
    }
    if(iteration == iterations - 1) { printProgress(); break; }
    previousLikelihood = likelihood;
    // An undamped increasing-Q EM step needs no duplicate likelihood-only fit.
    // Its next E-step evaluates and checks the observed likelihood above.
    if(!adaptive && weight == 1.0 && damping == 1.0 &&
       closedFormAccepted && validCandidate(emProposal)) {
      parameters = emProposal;
      progressTrust = 1.0;
      progressLine = 1.0;
      printProgress();
      continue;
    }
    bool accepted = false;
    const int halvingsBeforeSearch = lineHalvings;
    const int phases = adaptive && weight == 0.0 && currentScoreNorm > scoreTolerance ? 2 : 1;
    for(int phase = 0; phase < phases && !accepted; ++phase) {
      if(phase == 1) {
        requestRecovery(regularizedAI ? "information-solve" : "line-search");
        prepareRecovery();
        direction = emProposal - parameters;
        weightHistory(weightHistory.n_elem - 1) = weight;
        qGains(qGains.n_elem - 1) = qGain;
        mstepConvergenceHistory(mstepConvergenceHistory.n_elem - 1) = mstepConverged;
        mstepEvaluationHistory(mstepEvaluationHistory.n_elem - 1) = mstepEvaluations;
        norms(0, norms.n_cols - 1) = arma::norm(direction, 1);
        norms(2, norms.n_cols - 1) = arma::norm(direction, 2) / (1.0 + arma::norm(parameters, 2));
        if(parameterExpanded) {
          pxGains(pxGains.n_elem - 1) = pxGain;
          pxConvergenceHistory(pxConvergenceHistory.n_elem - 1) = pxConverged;
          pxEvaluationHistory(pxEvaluationHistory.n_elem - 1) = pxObjectiveEvaluations;
          expansionHistory.col(expansionHistory.n_cols - 1) = logExpansion;
        }
      }
      const bool mixed = weight > 0.0 && weight < 1.0;
      const int attempts = (mixed ? 2 : 1) + (parameterExpanded && weight > 0.0 && terms > 1 ? 1 : 0);
      for(int attempt = 0; attempt < attempts; ++attempt) {
      if(attempt == 1 && mixed) { direction = emProposal - parameters; ++emFallbacks; }
      else if(attempt > 0) { direction = ordinaryProposal - parameters; ++pxEMFallbacks; }
      progressUpdate = regularizedAI && weight == 0.0 ? "AI-reg"
        : MMEProgress::update(attempt > 0 ? 1.0 : weight,
            parameterExpanded && (attempt == 0 || (attempt == 1 && mixed)));
      double alpha = damping;
      if(weight < 1.0 && attempt == 0) {
        for(arma::uword j = 0; j < count; ++j) {
          const double step = adaptive && weight == 0.0 ? cappedAIDirection(j) : direction(j);
          if(std::abs(step) > 0.0) {
            alpha = std::min(alpha, stepCaps(j) / std::abs(step));
          }
        }
      }
      for(arma::uword index = 0; index < count; ++index) {
        if(std::isfinite(varianceLower(index)) && direction(index) < 0.0) {
          alpha = std::min(alpha,
            std::max(0.0, (parameters(index) - varianceLower(index)) / -direction(index)));
        }
      }
      progressTrust = damping > 0.0 ? alpha / damping : 1.0;
      const double initialAlpha = alpha;
      for(int search = 0; search < 25; ++search) {
        progressLine = initialAlpha > 0.0 ? alpha / initialAlpha : 0.0;
        const arma::vec candidate = parameters + alpha * direction;
        if(!validCandidate(candidate)) {
          alpha *= 0.5;
          ++lineHalvings;
          continue;
        }
        Rcpp::List trial;
        try {
          const bool nextNativeAI = !lowRank &&
            weights(std::min<arma::uword>(iteration + 1, weights.n_elem - 1)) == 0.0;
          trial = evaluate(candidate,
            nextNativeAI && iteration + 1 == iterations - 1 ? inverseMode : 0,
            nullptr, nextNativeAI);
        } catch(const CovarianceFactorizationError & error) {
          ++invalidTrials;
          ++invalidTrialReasons[error.what()];
          ++coefficientEvaluations;
          alpha *= 0.5;
          ++lineHalvings;
          continue;
        }
        const arma::vec trialLikelihoods = Rcpp::as<arma::vec>(trial["llik"]);
        const double trialLikelihood = trialLikelihoods(trialLikelihoods.n_elem - 1);
        const bool improves = adaptive && (weight > 0.0 || regularizedAI)
          ? trialLikelihood > likelihood + 1.0e-10
          : trialLikelihood >= likelihood - 1.0e-10;
        if(std::isfinite(trialLikelihood) && improves) {
          parameters = candidate;
          havePendingAI = !lowRank &&
            weights(std::min<arma::uword>(iteration + 1, weights.n_elem - 1)) == 0.0;
          if(havePendingAI) { pendingAI = trial; }
          accepted = true;
          break;
        }
        alpha *= 0.5;
        ++lineHalvings;
      }
      if(accepted) { break; }
      }
    }
    previousHalvings = lineHalvings - halvingsBeforeSearch;
    previousStepWasAI = weight == 0.0;
    if(accepted && regularizedAI && weight == 0.0) {
      ++regularizedAccepted;
    }
    if(adaptive && weight > 0.0) {
      previousHalvings = 0;
      if(accepted) { ++acceptedRecoveries; }
    }
    printProgress();
    if(!accepted) {
      stalled = true;
      stopReason = adaptive && weight > 0.0 ? "recovery-stalled" : "line-search-stalled";
      break;
    }
  }
  Rcpp::List finalCovariance = covarianceAt(parameters);
  finalCovariance.attr("emInformationOnly") = true;
  if(reducedAI && !retainInformation) {
    finalCovariance.attr("aiTangent") = tangent;
  }
  if(retainInformation && result.containsElementNamed(".aiTangent")) {
    finalInformationReady = false;
  }
  if(!finalInformationReady) {
    result = observe(conditionalFit(Z, Zind, finalCovariance, inverseMode, nullptr));
    ++likelihoodEvaluations;
  }
  if(profileIterations.n_rows > 0) {
    profileIterations.row(profileIterations.n_rows - 1) +=
      (profileSeconds - recordedProfile).t();
  }
  // The latent RR return uses +2 score; native conditional fits use -2 score.
  const arma::vec finalScore = (lowRank ? -1.0 : 1.0) *
    Rcpp::as<arma::vec>(result["dLu"]);
  const arma::mat finalInformation = Rcpp::as<arma::mat>(result["InfMat"]);
  const arma::mat finalInformationTangent = result.containsElementNamed(".aiTangent")
    ? Rcpp::as<arma::mat>(result[".aiTangent"]) : arma::mat();
  double projectedScoreNorm = 0.0;
  const arma::mat finalTangent = varianceFloorTangent(tangent, parameters, varianceLower,
    arma::vec(-tangent * (tangent.t() * finalScore)));
  if(finalTangent.n_cols > 0) {
    const arma::vec projectedScore = finalTangent.t() * finalScore;
    const arma::mat projectedInformation = aiInformationInTangent(
      finalInformation, finalInformationTangent, finalTangent);
    for(arma::uword j = 0; j < finalTangent.n_cols; ++j) {
      projectedScoreNorm = std::max(projectedScoreNorm, std::abs(projectedScore(j)) /
        std::sqrt(std::max(inverseTolerance, std::abs(projectedInformation(j, j)))));
    }
  }
  if(adaptive) {
    projectedScoreNorm = projectedEMScoreNorm(tangent, parameters, varianceLower,
      arma::vec(-0.5 * finalScore), arma::mat(0.5 * finalInformation),
      inverseTolerance, finalInformationTangent);
  }
  converged = converged && projectedScoreNorm <=
    std::sqrt(std::max(likelihoodTolerance, parameterTolerance));
  if(adaptive && !converged && stopReason == "converged") {
    stopReason = "score-not-certified";
  }
  result["monitor"] = monitor;
  result["llik"] = arma::mat(likelihoodHistory.t());
  result["normMonitor"] = norms;
  result["convergence"] = converged && !stalled;
  result["percDelta"] = arma::mat(count, monitor.n_cols, arma::fill::zeros);
  result["toBoundary"] = arma::mat(monitor.n_cols, count, arma::fill::zeros);
  Rcpp::List diagnostics = result["engineDiagnostics"];
  diagnostics["reducedAI"] = Rcpp::List::create(
    Rcpp::Named("active") = reducedAI,
    Rcpp::Named("fullParameters") = static_cast<double>(count),
    Rcpp::Named("independentParameters") = static_cast<double>(tangent.n_cols),
    Rcpp::Named("evaluations") = reducedEvaluations,
    Rcpp::Named("fullInformationEvaluations") = fullInformationEvaluations,
    Rcpp::Named("sensitivityColumns") = reducedSensitivityColumns,
    Rcpp::Named("sensitivityColumnsSkipped") = skippedSensitivityColumns,
    Rcpp::Named("maxSensitivityEntries") = reducedSensitivityEntries,
    Rcpp::Named("maxAssemblyEntries") = reducedAssemblyEntries,
    Rcpp::Named("fallbackReason") = reducedMode != "off" && !reducedAI
      ? (lowRank ? "latent-reduced-rank" : "no-independent-reduction") : "");
  diagnostics["optimizer"] = adaptive ? "AI-first-recovery"
    : parameterExpanded ? "PX-GEM-AI-direction" : "EM/GEM-AI-direction";
  diagnostics["likelihoodMode"] = reml ? "REML" : "ML";
  if(adaptive) {
    diagnostics["recovery"] = Rcpp::List::create(
      Rcpp::Named("method") = method,
      Rcpp::Named("likelihoodMode") = reml ? "REML" : "ML",
      Rcpp::Named("steps") = recoverySteps,
      Rcpp::Named("acceptedSteps") = acceptedRecoveries,
      Rcpp::Named("iterations") = recoveryIterations,
      Rcpp::Named("reasons") = recoveryReasons,
      Rcpp::Named("regularizedAttempts") = regularizedAttempts,
      Rcpp::Named("regularizedAcceptedSteps") = regularizedAccepted,
      Rcpp::Named("regularizedIterations") = regularizedIterations,
      Rcpp::Named("regularizedRidges") = regularizedRidges,
      Rcpp::Named("stopReason") = stopReason,
      Rcpp::Named("finalProjectedScoreNorm") = projectedScoreNorm,
      Rcpp::Named("scoreTolerance") = scoreTolerance);
    if(stalled) {
      Rcpp::warning("AI-first optimization stalled at the last accepted parameters (%s). "
                   "Consider different starting values; convergence was not certified.",
                   stopReason.c_str());
    }
  }
  diagnostics["coefficientEvaluations"] = coefficientEvaluations;
  for(const auto & counter : engineCounters) { diagnostics[counter.first] = counter.second; }
  diagnostics["lineSearchHalvings"] = lineHalvings;
  if(!profileSeconds.is_empty()) {
    Rcpp::List profile = diagnostics["profile"];
    Rcpp::NumericVector seconds = profile["seconds"];
    std::copy(profileSeconds.begin(), profileSeconds.end(), seconds.begin());
    profile["seconds"] = seconds;
    Rcpp::List counters = profile["counters"];
    for(const auto & entry : profileCounters) { counters[entry.first] = entry.second; }
    counters["iterations"] = static_cast<double>(monitor.n_cols);
    counters["sensitivityRHSColumns"] = profileCounters["sensitivityRHSColumns"] +
      (lowRank ? static_cast<double>(sensitivityColumns) : 0.0);
    profile["counters"] = counters;
    const Rcpp::CharacterVector phaseNames = seconds.names();
    for(int j = 0; j < phaseNames.size(); ++j) {
      const std::string phase = Rcpp::as<std::string>(phaseNames[j]);
      if(phase == "setup" || phase == "final") { profileIterations.col(j).zeros(); }
    }
    Rcpp::NumericMatrix byIteration = Rcpp::wrap(profileIterations);
    Rcpp::colnames(byIteration) = phaseNames;
    profile["secondsByIteration"] = byIteration;
    profile["scope"] = "successful conditional fits and final inference";
    diagnostics["profile"] = profile;
  }
  diagnostics["blend"] = Rcpp::List::create(
    Rcpp::Named("method") = method,
    Rcpp::Named("likelihoodMode") = reml ? "REML" : "ML",
    Rcpp::Named("weights") = weightHistory,
    Rcpp::Named("emEvaluations") = emEvaluations,
    Rcpp::Named("aiEvaluations") = aiEvaluations,
    Rcpp::Named("iterationSensitivityColumns") = static_cast<double>(sensitivityColumns),
    Rcpp::Named("finalAIRefreshes") = 1,
    Rcpp::Named("closedFormSteps") = closedFormSteps,
    Rcpp::Named("generalizedSteps") = generalizedSteps,
    Rcpp::Named("numericalMstepConvergences") = numericalConvergences,
    Rcpp::Named("numericalMstepConverged") = mstepConvergenceHistory,
    Rcpp::Named("numericalMstepEvaluations") = mstepEvaluationHistory,
    Rcpp::Named("qEvaluations") = qEvaluations,
    Rcpp::Named("qGains") = qGains,
    Rcpp::Named("momentPanelSolves") = momentSolves,
    Rcpp::Named("likelihoodEvaluations") = likelihoodEvaluations,
    Rcpp::Named("lineSearchHalvings") = lineHalvings,
    Rcpp::Named("emFallbacks") = emFallbacks,
    Rcpp::Named("lineSearchStalled") = stalled,
    Rcpp::Named("finalProjectedScoreNorm") = projectedScoreNorm,
    Rcpp::Named("finalInformationExact") = true);
  Rcpp::List blendDiagnostics = diagnostics["blend"];
  blendDiagnostics["latentMomentPreparations"] = latentMomentPreparations;
  blendDiagnostics["latentMomentBases"] = static_cast<double>(latentMomentBases);
  blendDiagnostics["randomMomentPairs"] = static_cast<double>(randomMomentPairs);
  blendDiagnostics["randomOffDiagonalPairsSkipped"] = static_cast<double>(randomOffDiagonalPairsSkipped);
  blendDiagnostics["scalarResidualTraceSteps"] = scalarResidualTraceSteps;
  blendDiagnostics["scalarResidualContractions"] = scalarResidualContractions;
  blendDiagnostics["nativeAIEvaluations"] = nativeAIEvaluations;
  blendDiagnostics["latentAIEvaluations"] = latentAIEvaluations;
  blendDiagnostics["acceptedAIFitReuses"] = acceptedAIFitReuses;
  if(parameterExpanded) {
    blendDiagnostics["pxEvaluations"] = pxEvaluations;
    blendDiagnostics["pxExpansionTerms"] = terms - 1;
    blendDiagnostics["pxLogExpansion"] = expansionHistory;
    blendDiagnostics["pxQGains"] = pxGains;
    blendDiagnostics["pxMstepConverged"] = pxConvergenceHistory;
    blendDiagnostics["pxMstepEvaluations"] = pxEvaluationHistory;
    blendDiagnostics["pxEMFallbacks"] = pxEMFallbacks;
  }
  Rcpp::IntegerVector invalidCounts(invalidTrialReasons.size());
  Rcpp::CharacterVector invalidReasons(invalidTrialReasons.size());
  int reasonIndex = 0;
  for(const auto & reason : invalidTrialReasons) {
    invalidReasons[reasonIndex] = reason.first;
    invalidCounts[reasonIndex++] = reason.second;
  }
  invalidCounts.names() = invalidReasons;
  blendDiagnostics["invalidTrialEvaluations"] = invalidTrials;
  blendDiagnostics["invalidTrialReasons"] = invalidCounts;
  diagnostics["blend"] = blendDiagnostics;
  if(invalidTrials > 0) {
    Rcpp::warning("EM/AI blending rejected %d numerically unfactorizable trial fits. "
                 "Reasons and counts are recorded in engineDiagnostics$blend$invalidTrialReasons.",
                 invalidTrials);
  }
  result["engineDiagnostics"] = diagnostics;
  Rcpp::List factorization = result["factorizationDiagnostics"];
  Rcpp::List timings = factorization["timings"];
  for(const auto & entry : factorTimings) {
    Rcpp::List timing = timings[entry.first];
    timing["calls"] = entry.second.calls;
    timing["elapsedSeconds"] = entry.second.elapsed;
    timing["cpuSeconds"] = entry.second.cpu;
    timing["averageCPUs"] = entry.second.elapsed > 0.0
      ? entry.second.cpu / entry.second.elapsed : NA_REAL;
    timings[entry.first] = timing;
  }
  factorization["timings"] = timings;
  result["factorizationDiagnostics"] = factorization;
  if(result.containsElementNamed("gpuDiagnostics")) {
    Rcpp::List gpu = result["gpuDiagnostics"];
    for(const auto & counter : gpuCounters) { gpu[counter.first] = counter.second; }
    result["gpuDiagnostics"] = gpu;
  }
  if(result.containsElementNamed(".aiTangent")) {
    const Rcpp::CharacterVector names = result.names();
    for(int j = 0; j < names.size(); ++j) {
      if(Rcpp::as<std::string>(names[j]) == ".aiTangent") { result.erase(j); break; }
    }
  }
  return result;
}

}

#endif
