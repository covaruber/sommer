#include <RcppArmadillo.h>
#include "mme_cholmod.h"
#include <Eigen/Sparse>
#include <cmath>
#include <limits>
#include <vector>
#include "mme_glmm_family.h"
#include "mme_glmm_precision.h"

namespace {

double stableLogLikelihood(const arma::vec & response,
                           const arma::vec & eta,
                           const sommer::GLMMFamilyModel & family,
                           arma::vec * mean,
                           arma::vec * weight,
                           arma::vec * score){
  const auto evaluation = family.evaluate(response, eta);
  *mean = evaluation.mean;
  *weight = evaluation.curvature;
  *score = evaluation.score;
  return evaluation.logDensity;
}

bool factorLogDeterminant(sommer::CholmodState & factor,
                          double & logDeterminant){
  return sommer::laplaceFactorLogDeterminant(factor, logDeterminant);
}

} // namespace

// Evaluate the conditional random-effect mode and the Laplace log likelihood
// for a supplied fixed-effect vector and random-effect precision matrix.
Rcpp::List laplaceModeEvaluation(const arma::sp_mat & fixedDesign,
                                    const arma::sp_mat & randomDesign,
                                    const arma::sp_mat & randomPrecision,
                                    const arma::vec & response,
                                    const arma::vec & offset,
                                    const arma::vec & fixedEffects,
                                    const arma::vec & modeStart,
                                    const sommer::GLMMFamilyModel & family,
                                    const double logDetPrecision,
                                    const int maxModeIterations,
                                    const double tolerance,
                                    sommer::CholmodState & factor,
                                    int * symbolicAnalyses,
                                    const bool verifyPrecision){
  const arma::uword observationCount = response.n_elem;
  const arma::uword randomCount = randomDesign.n_cols;
  if(fixedDesign.n_rows != observationCount ||
     randomDesign.n_rows != observationCount ||
     offset.n_elem != observationCount ||
     fixedEffects.n_elem != fixedDesign.n_cols ||
     modeStart.n_elem != randomCount ||
     randomPrecision.n_rows != randomCount ||
     randomPrecision.n_cols != randomCount){
    Rcpp::stop("Native Laplace inputs have incompatible dimensions.");
  }
  if(maxModeIterations < 1 || !std::isfinite(tolerance) || tolerance <= 0.0){
    Rcpp::stop("Native Laplace mode controls must be positive.");
  }
  if(!response.is_finite() || !offset.is_finite() ||
     !fixedDesign.is_finite() || !randomDesign.is_finite() ||
     !randomPrecision.is_finite() || !fixedEffects.is_finite() ||
     !modeStart.is_finite() || !std::isfinite(logDetPrecision)){
    Rcpp::stop("Native Laplace inputs must be finite.");
  }
  family.validateResponse(response);
  if(observationCount == 0 || randomCount == 0){
    Rcpp::stop("Native Laplace evaluation requires observations and random effects.");
  }
  if(!randomPrecision.is_symmetric(1e-12)){
    Rcpp::stop("Native Laplace random precision must be symmetric.");
  }

  int iterations = 0;
  auto invalidEvaluation = [&iterations](const std::string & reason){
    return Rcpp::List::create(Rcpp::_["valid"]=false,
      Rcpp::_["reason"]=reason, Rcpp::_["iterations"]=iterations);
  };
  if(verifyPrecision){
    sommer::CholmodState precisionFactor;
    const Eigen::SparseMatrix<double> eigenPrecision =
      sommer::laplaceEigenSparse(randomPrecision);
    if(!precisionFactor.analyzePattern(eigenPrecision) ||
       !precisionFactor.factorize(eigenPrecision)){
      return invalidEvaluation("non-positive-definite random precision");
    }
    double verifiedLogDetPrecision = NA_REAL;
    if(!factorLogDeterminant(precisionFactor, verifiedLogDetPrecision)){
      return invalidEvaluation("invalid random precision determinant");
    }
    if(std::abs(verifiedLogDetPrecision - logDetPrecision) >
       1e-8 * (1.0 + std::abs(verifiedLogDetPrecision))){
      Rcpp::stop("Native Laplace precision log determinant is inconsistent with its matrix.");
    }
  }

  arma::vec mode = modeStart;
  arma::vec mean(observationCount);
  arma::vec weight(observationCount);
  arma::vec score(observationCount);
  arma::vec eta = offset + fixedDesign * fixedEffects + randomDesign * mode;
  double conditionalLogLikelihood =
    stableLogLikelihood(response, eta, family, &mean, &weight, &score);
  if(!std::isfinite(conditionalLogLikelihood)){
    return invalidEvaluation("non-finite conditional likelihood");
  }
  double penalizedLogLikelihood =
    conditionalLogLikelihood - 0.5 * arma::as_scalar(mode.t() * randomPrecision * mode);
  bool converged = false;
  double hessianLogDeterminant = NA_REAL;
  double finalModeScore = NA_REAL;

  for(int iteration = 0; iteration < maxModeIterations; ++iteration){
    Rcpp::checkUserInterrupt();
    iterations = iteration + 1;
    arma::sp_mat weightedDesign = randomDesign;
    for(auto cell = weightedDesign.begin(); cell != weightedDesign.end(); ++cell){
      *cell *= weight(cell.row());
    }
    const arma::sp_mat hessian =
      randomPrecision + randomDesign.t() * weightedDesign;
    const arma::vec modeScore = randomDesign.t() * score - randomPrecision * mode;
    if(!modeScore.is_finite()){
      return invalidEvaluation("non-finite mode score");
    }
    finalModeScore = arma::norm(modeScore, "inf");
    if(finalModeScore <= tolerance * (1.0 + arma::norm(mode, "inf"))){
      converged = true;
      const Eigen::SparseMatrix<double> eigenHessian =
        sommer::laplaceEigenSparse(hessian);
      if(!factor.analyzePattern(eigenHessian, symbolicAnalyses) ||
         !factor.factorize(eigenHessian) ||
         !factorLogDeterminant(factor, hessianLogDeterminant)){
        return invalidEvaluation("non-positive-definite mode Hessian");
      }
      break;
    }
    const Eigen::SparseMatrix<double> eigenHessian = sommer::laplaceEigenSparse(hessian);
    if(!factor.analyzePattern(eigenHessian, symbolicAnalyses) || !factor.factorize(eigenHessian)){
      return invalidEvaluation("non-positive-definite mode Hessian");
    }
    Eigen::MatrixXd rightHandSide(static_cast<int>(randomCount), 1);
    for(arma::uword index = 0; index < randomCount; ++index){
      rightHandSide(static_cast<int>(index), 0) = modeScore(index);
    }
    const Eigen::MatrixXd solved = factor.solve(rightHandSide,
      "Native Laplace mode solve failed.");
    arma::vec step(randomCount);
    for(arma::uword index = 0; index < randomCount; ++index){
      step(index) = solved(static_cast<int>(index), 0);
    }
    if(!step.is_finite()){
      return invalidEvaluation("non-finite Newton step");
    }

    bool accepted = false;
    double stepScale = 1.0;
    arma::vec candidateMean(observationCount);
    arma::vec candidateWeight(observationCount);
    arma::vec candidateScore(observationCount);
    while(stepScale >= 1.0 / 1024.0){
      const arma::vec candidateMode = mode + stepScale * step;
      const arma::vec candidateEta =
        offset + fixedDesign * fixedEffects + randomDesign * candidateMode;
      const double candidateLikelihood = stableLogLikelihood(
        response, candidateEta, family, &candidateMean, &candidateWeight,
        &candidateScore);
      const double candidatePenalizedLikelihood = candidateLikelihood -
        0.5 * arma::as_scalar(candidateMode.t() * randomPrecision * candidateMode);
      const double likelihoodRoundoff = 64.0 * std::numeric_limits<double>::epsilon() *
        (1.0 + std::abs(penalizedLogLikelihood));
      bool improvesMode = false;
      if(std::isfinite(candidatePenalizedLikelihood) &&
         candidatePenalizedLikelihood >= penalizedLogLikelihood - likelihoodRoundoff){
        if(candidatePenalizedLikelihood > penalizedLogLikelihood + likelihoodRoundoff){
          improvesMode = true;
        }else{
          // At likelihood rounding precision, require a better mode score.
          const arma::vec candidateModeScore =
            randomDesign.t() * candidateScore - randomPrecision * candidateMode;
          const double candidateScoreNorm = arma::norm(candidateModeScore, "inf");
          improvesMode = candidateModeScore.is_finite() &&
            (candidateScoreNorm < finalModeScore ||
             candidateScoreNorm <= tolerance * (1.0 + arma::norm(candidateMode, "inf")));
        }
      }
      if(improvesMode){
        mode = candidateMode;
        eta = candidateEta;
        mean = candidateMean;
        weight = candidateWeight;
        score = candidateScore;
        conditionalLogLikelihood = candidateLikelihood;
        penalizedLogLikelihood = candidatePenalizedLikelihood;
        accepted = true;
        break;
      }
      stepScale *= 0.5;
    }
    if(!accepted){
      if(arma::norm(modeScore, "inf") <= 10.0 * tolerance *
         (1.0 + arma::norm(mode, "inf"))){
        converged = true;
        if(!factorLogDeterminant(factor, hessianLogDeterminant)){
          return invalidEvaluation("invalid mode Hessian determinant");
        }
        break;
      }
      return invalidEvaluation("mode line search failed");
    }
  }
  if(!converged){
    return invalidEvaluation("mode iteration limit reached");
  }
  const double logLikelihood = conditionalLogLikelihood -
    0.5 * arma::as_scalar(mode.t() * randomPrecision * mode) +
    0.5 * logDetPrecision - 0.5 * hessianLogDeterminant;
  if(!std::isfinite(logLikelihood)){
    return invalidEvaluation("non-finite Laplace likelihood");
  }
  const auto finalFamily = family.evaluate(response, eta);
  arma::sp_mat weightedRandomDesign = randomDesign;
  for(auto cell = weightedRandomDesign.begin();
      cell != weightedRandomDesign.end(); ++cell){
    *cell *= finalFamily.curvature(cell.row());
  }
  const arma::sp_mat finalModeHessian =
    randomPrecision + randomDesign.t() * weightedRandomDesign;
  return Rcpp::List::create(
    Rcpp::_["valid"]=true,
    Rcpp::_["reason"]="converged",
    Rcpp::_["logLik"]=logLikelihood,
    Rcpp::_["mode"]=mode,
    Rcpp::_["fitted"]=mean,
    Rcpp::_["linearPredictors"]=eta,
    Rcpp::_["modeHessian"]=finalModeHessian,
    Rcpp::_["observationCurvature"]=finalFamily.curvature,
    Rcpp::_["iterations"]=iterations,
    Rcpp::_["modeScore"]=finalModeScore,
    Rcpp::_["logDetHessian"]=hessianLogDeterminant);
}

// [[Rcpp::export]]
Rcpp::List sommer_laplace_glmm_eval(const arma::sp_mat & fixedDesign,
    const arma::sp_mat & randomDesign, const arma::sp_mat & randomPrecision,
    const arma::vec & response, const arma::vec & offset,
    const arma::vec & fixedEffects, const arma::vec & modeStart,
    const std::string & family, const double logDetPrecision,
    const int maxModeIterations = 100, const double tolerance = 1e-8){
  sommer::CholmodState modeFactor;
  const sommer::GLMMFamilyModel familyModel(family, response.n_elem);
  return laplaceModeEvaluation(fixedDesign, randomDesign, randomPrecision,
    response, offset, fixedEffects, modeStart, familyModel, logDetPrecision,
    maxModeIterations, tolerance, modeFactor, nullptr, true);
}

namespace {

class LaplaceObjective {
public:
  arma::sp_mat fixedDesign;
  arma::sp_mat randomDesign;
  arma::vec response;
  arma::vec offset;
  sommer::GLMMFamilyModel family;
  sommer::LaplacePrecision precision;
  sommer::CholmodState modeFactor;
  int symbolicAnalyses = 0;
  int evaluations = 0;
  arma::vec cachedParameters;
  Rcpp::List cachedResult;
  bool cachedMode = false;
  int cachedMaxModeIterations = 0;
  double cachedTolerance = 0.0;
  int modeEvaluations = 0;
  double memoryMB;

  LaplaceObjective(const arma::sp_mat & suppliedFixed,
      const arma::sp_mat & suppliedRandom, const arma::vec & suppliedResponse,
      const arma::vec & suppliedOffset, const std::string & suppliedFamily,
      const Rcpp::List & descriptors, const Rcpp::List & relationships,
      const double suppliedMemory, const Rcpp::List & familySpecs,
      const arma::uvec & familyIndex, const arma::vec & trials)
    : fixedDesign(suppliedFixed), randomDesign(suppliedRandom),
      response(suppliedResponse), offset(suppliedOffset),
      family(suppliedFamily, suppliedResponse.n_elem, familySpecs, familyIndex, trials),
      precision(descriptors, relationships), memoryMB(suppliedMemory) {
    if(!std::isfinite(memoryMB) || memoryMB <= 0){
      Rcpp::stop("Native Laplace memory budget must be positive and finite.");
    }
    family.validateResponse(response);
  }

  arma::mat solve(const arma::mat & rightHandSide){
    Eigen::Map<const Eigen::MatrixXd> mapped(rightHandSide.memptr(),
      rightHandSide.n_rows, rightHandSide.n_cols);
    const Eigen::MatrixXd result = modeFactor.solve(mapped,
      "Native Laplace derivative solve failed.");
    return arma::mat(result.data(), result.rows(), result.cols());
  }

  arma::vec gradient(const Rcpp::List & modeResult){
    const arma::uword randomCount = randomDesign.n_cols;
    const arma::uword observationCount = response.n_elem;
    // Bound dense solve panels rather than constructing H^-1 or Z H^-1 Z'.
    const double columnBytes = 8.0 * (3.0 * randomCount + observationCount);
    const arma::uword panelSize = static_cast<arma::uword>(std::min(32.0,
      std::floor(memoryMB * 1e6 / columnBytes)));
    if(panelSize == 0){
      Rcpp::stop("Native Laplace derivative solve exceeds the dense memory budget.");
    }
    const arma::vec mode = Rcpp::as<arma::vec>(modeResult["mode"]);
    const arma::vec predictor = Rcpp::as<arma::vec>(modeResult["linearPredictors"]);
    const auto familyResult = family.evaluate(response, predictor);
    arma::vec leverage(observationCount, arma::fill::zeros);
    for(arma::uword first = 0; first < observationCount; first += panelSize){
      Rcpp::checkUserInterrupt();
      const arma::uword last = std::min(observationCount - 1, first + panelSize - 1);
      const arma::mat designPanel(randomDesign.rows(first, last).t());
      const arma::mat solved = solve(designPanel);
      leverage.subvec(first, last) = arma::sum(designPanel % solved, 0).t();
    }
    const arma::vec correction = 0.5 * familyResult.curvatureD1 % leverage;
    const arma::vec correctionSolve = solve(arma::mat(randomDesign.t() * correction));
    const arma::vec fixedScore = fixedDesign.t() * (familyResult.score -
      correction + familyResult.curvature % (randomDesign * correctionSolve));
    const auto derivatives = precision.derivatives();
    arma::vec result(fixedScore.n_elem + derivatives.size() + family.parameterCount());
    if(fixedScore.n_elem) result.head(fixedScore.n_elem) = fixedScore;
    for(std::size_t parameterIndex = 0; parameterIndex < derivatives.size();
        ++parameterIndex){
      Rcpp::checkUserInterrupt();
      const auto & derivative = derivatives[parameterIndex];
      const arma::uword blockEnd = derivative.start + derivative.block.n_rows - 1;
      const arma::vec blockMode = mode.subvec(derivative.start, blockEnd);
      const arma::vec precisionMode = derivative.block * blockMode;
      double traceValue = 0.0;
      for(arma::uword first = 0; first < derivative.block.n_cols; first += panelSize){
        const arma::uword last = std::min(derivative.block.n_cols - 1,
          first + panelSize - 1);
        arma::mat rightHandSide(randomCount, last - first + 1, arma::fill::zeros);
        rightHandSide.rows(derivative.start, blockEnd) =
          arma::mat(derivative.block.cols(first, last));
        const arma::mat solved = solve(rightHandSide);
        for(arma::uword columnIndex = first; columnIndex <= last; ++columnIndex){
          traceValue += solved(derivative.start + columnIndex, columnIndex - first);
        }
      }
      result(fixedScore.n_elem + parameterIndex) =
        -0.5 * arma::dot(blockMode, precisionMode) +
        0.5 * derivative.logDet - 0.5 * traceValue +
        arma::dot(correctionSolve.subvec(derivative.start, blockEnd), precisionMode);
    }
    arma::uword outputPosition = fixedScore.n_elem + derivatives.size();
    const arma::vec modeCorrection = randomDesign * correctionSolve;
    for(std::size_t familyPosition = 0; familyPosition < family.specifications.size();
        ++familyPosition){
      if(!family.specifications[familyPosition].free) continue;
      const arma::uvec records = arma::find(family.index == familyPosition);
      const arma::vec rowScore = familyResult.parameterScore.elem(records) -
        0.5 * leverage.elem(records) % familyResult.curvatureParameter.elem(records) -
        modeCorrection.elem(records) % familyResult.scoreParameter.elem(records);
      const auto & specification = family.specifications[familyPosition];
      if(specification.zeroDesign.n_cols){
        const arma::vec zeroScore = specification.zeroDesign.rows(records).t() * rowScore;
        result.subvec(outputPosition, outputPosition + zeroScore.n_elem - 1) = zeroScore;
        outputPosition += zeroScore.n_elem;
      }else result(outputPosition++) = arma::accu(rowScore);
    }
    if(!result.is_finite()) Rcpp::stop("Native Laplace gradient is non-finite.");
    return result;
  }
};

} // namespace

// [[Rcpp::export]]
SEXP sommer_laplace_objective_create(const arma::sp_mat & fixedDesign,
    const arma::sp_mat & randomDesign, const arma::vec & response,
    const arma::vec & offset, const std::string & family,
    const Rcpp::List & descriptors, const Rcpp::List & relationships,
    const double memoryMB = 1000,
    Rcpp::Nullable<Rcpp::List> familySpecs = R_NilValue,
    Rcpp::Nullable<Rcpp::IntegerVector> familyIndex = R_NilValue,
    Rcpp::Nullable<Rcpp::NumericVector> trials = R_NilValue){
  arma::uvec nativeIndex;
  if(familyIndex.isNotNull()){
    const Rcpp::IntegerVector suppliedIndex(familyIndex);
    nativeIndex.set_size(suppliedIndex.size());
    for(int recordPosition = 0; recordPosition < suppliedIndex.size(); ++recordPosition){
      if(suppliedIndex[recordPosition] == NA_INTEGER || suppliedIndex[recordPosition] < 1)
        Rcpp::stop("Native family indices must be positive integers.");
      nativeIndex(recordPosition) = suppliedIndex[recordPosition] - 1;
    }
  }
  const Rcpp::List specifications = familySpecs.isNotNull() ?
    Rcpp::List(familySpecs) : Rcpp::List();
  const arma::vec nativeTrials = trials.isNotNull() ?
    Rcpp::as<arma::vec>(Rcpp::NumericVector(trials)) : arma::vec();
  return Rcpp::XPtr<LaplaceObjective>(new LaplaceObjective(fixedDesign,
    randomDesign, response, offset, family, descriptors, relationships, memoryMB,
    specifications, nativeIndex, nativeTrials), true);
}

// [[Rcpp::export]]
Rcpp::List sommer_laplace_objective_eval(SEXP context,
    const arma::vec & parameters, const int maxModeIterations = 100,
    const double tolerance = 1e-8, const bool gradient = true){
  Rcpp::XPtr<LaplaceObjective> objective(context);
  const arma::uword covarianceCount = objective->precision.parameterCount();
  const arma::uword nuisanceCount = objective->family.parameterCount();
  if(parameters.n_elem != objective->fixedDesign.n_cols + covarianceCount + nuisanceCount){
    Rcpp::stop("Native Laplace parameter vector has incompatible length.");
  }
  if(maxModeIterations < 1 || !std::isfinite(tolerance) || tolerance <= 0.0){
    Rcpp::stop("Native Laplace mode controls must be positive.");
  }
  std::string reason;
  if(!parameters.is_finite()){
    return Rcpp::List::create(Rcpp::_["valid"]=false,
      Rcpp::_["reason"]="non-finite parameter proposal");
  }
  ++objective->evaluations;
  const arma::uword fixedCount = objective->fixedDesign.n_cols;
  const bool reuseMode = objective->cachedMode &&
    maxModeIterations == objective->cachedMaxModeIterations &&
    tolerance == objective->cachedTolerance &&
    parameters.n_elem == objective->cachedParameters.n_elem &&
    arma::all(parameters == objective->cachedParameters);
  if(!objective->family.update(parameters.tail(nuisanceCount))){
    objective->cachedMode = false;
    return Rcpp::List::create(Rcpp::_["valid"]=false,
      Rcpp::_["reason"]="family shape/dispersion outside numerical range");
  }
  const arma::vec covarianceParameters = covarianceCount ?
    parameters.subvec(fixedCount, fixedCount + covarianceCount - 1) : arma::vec();
  if(!objective->precision.update(covarianceParameters, reason)){
    objective->cachedMode = false;
    return Rcpp::List::create(Rcpp::_["valid"]=false, Rcpp::_["reason"]=reason);
  }
  Rcpp::List result;
  if(reuseMode){
    result = Rcpp::clone(objective->cachedResult);
  }else{
    objective->cachedMode = false;
    ++objective->modeEvaluations;
    result = laplaceModeEvaluation(objective->fixedDesign,
      objective->randomDesign, objective->precision.matrix, objective->response,
      objective->offset, parameters.head(fixedCount),
      arma::vec(objective->randomDesign.n_cols, arma::fill::zeros), objective->family,
      objective->precision.logDet, maxModeIterations, tolerance,
      objective->modeFactor, &objective->symbolicAnalyses, false);
  }
  if(Rcpp::as<bool>(result["valid"]) && gradient &&
     !result.containsElementNamed("gradient")){
    result["gradient"] = objective->gradient(result);
  }
  objective->cachedMode = Rcpp::as<bool>(result["valid"]);
  objective->cachedParameters = parameters;
  objective->cachedMaxModeIterations = maxModeIterations;
  objective->cachedTolerance = tolerance;
  objective->cachedResult = Rcpp::clone(result);
  result["work"] = Rcpp::List::create(
    Rcpp::_["evaluations"]=objective->evaluations,
    Rcpp::_["modeEvaluations"]=objective->modeEvaluations,
    Rcpp::_["symbolicAnalyses"]=objective->symbolicAnalyses,
    Rcpp::_["relationshipFactorizations"]=objective->precision.relationshipFactorizations,
    Rcpp::_["descriptorUpdates"]=objective->precision.descriptorUpdates,
    Rcpp::_["descriptorReuses"]=objective->precision.descriptorReuses);
  return result;
}

// [[Rcpp::export]]
Rcpp::List sommer_glmm_family_eval(const arma::vec & response,
    const arma::vec & predictor, const Rcpp::List & familySpecs,
    const Rcpp::IntegerVector & familyIndex, const arma::vec & trials){
  if(response.n_elem != predictor.n_elem || familyIndex.size() != static_cast<int>(response.n_elem))
    Rcpp::stop("Native family evaluation dimensions must match.");
  arma::uvec nativeIndex(response.n_elem);
  for(int recordPosition = 0; recordPosition < familyIndex.size(); ++recordPosition){
    if(familyIndex[recordPosition] == NA_INTEGER || familyIndex[recordPosition] < 1)
      Rcpp::stop("Native family indices must be positive integers.");
    nativeIndex(recordPosition) = familyIndex[recordPosition] - 1;
  }
  sommer::GLMMFamilyModel family("", response.n_elem, familySpecs, nativeIndex, trials);
  family.validateResponse(response);
  const auto evaluation = family.evaluate(response, predictor);
  return Rcpp::List::create(Rcpp::_["logDensity"]=evaluation.logDensity,
    Rcpp::_["mean"]=evaluation.mean, Rcpp::_["score"]=evaluation.score,
    Rcpp::_["curvature"]=evaluation.curvature,
    Rcpp::_["curvatureD1"]=evaluation.curvatureD1,
    Rcpp::_["curvatureD2"]=evaluation.curvatureD2,
    Rcpp::_["parameterScore"]=evaluation.parameterScore,
    Rcpp::_["scoreParameter"]=evaluation.scoreParameter,
    Rcpp::_["curvatureParameter"]=evaluation.curvatureParameter);
}
