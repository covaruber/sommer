#ifndef SOMMER_MME_GLMM_FAMILY_H
#define SOMMER_MME_GLMM_FAMILY_H

#include <RcppArmadillo.h>
#include <algorithm>
#include <cmath>
#include <limits>

namespace sommer {

// Truncated Taylor coefficients in eta (order four) and log shape/dispersion
// (order one). Mixed coefficients supply score and curvature derivatives.
struct GLMMTaylor {
  double coefficient[5][2] = {};
  GLMMTaylor(double value = 0.0){ coefficient[0][0] = value; }
};

inline GLMMTaylor operator+(const GLMMTaylor & left, const GLMMTaylor & right){
  GLMMTaylor result;
  for(int etaOrder = 0; etaOrder <= 4; ++etaOrder)
    for(int shapeOrder = 0; shapeOrder <= 1; ++shapeOrder)
      result.coefficient[etaOrder][shapeOrder] =
        left.coefficient[etaOrder][shapeOrder] + right.coefficient[etaOrder][shapeOrder];
  return result;
}

inline GLMMTaylor operator-(const GLMMTaylor & left, const GLMMTaylor & right){
  GLMMTaylor result;
  for(int etaOrder = 0; etaOrder <= 4; ++etaOrder)
    for(int shapeOrder = 0; shapeOrder <= 1; ++shapeOrder)
      result.coefficient[etaOrder][shapeOrder] =
        left.coefficient[etaOrder][shapeOrder] - right.coefficient[etaOrder][shapeOrder];
  return result;
}

inline GLMMTaylor operator*(const GLMMTaylor & left, const GLMMTaylor & right){
  GLMMTaylor result;
  for(int etaOrder = 0; etaOrder <= 4; ++etaOrder)
    for(int shapeOrder = 0; shapeOrder <= 1; ++shapeOrder)
      for(int etaPart = 0; etaPart <= etaOrder; ++etaPart)
        for(int shapePart = 0; shapePart <= shapeOrder; ++shapePart)
          result.coefficient[etaOrder][shapeOrder] +=
            left.coefficient[etaPart][shapePart] *
            right.coefficient[etaOrder - etaPart][shapeOrder - shapePart];
  return result;
}

inline GLMMTaylor glmmExp(const GLMMTaylor & argument){
  const double baseValue = std::exp(argument.coefficient[0][0]);
  if(baseValue == 0.0) return GLMMTaylor();
  GLMMTaylor increment = argument - GLMMTaylor(argument.coefficient[0][0]);
  GLMMTaylor result(1.0), power(1.0);
  double inverseFactorial = 1.0;
  for(int order = 1; order <= 5; ++order){
    power = power * increment;
    inverseFactorial /= order;
    result = result + GLMMTaylor(inverseFactorial) * power;
  }
  return GLMMTaylor(baseValue) * result;
}

inline GLMMTaylor glmmLog(const GLMMTaylor & argument){
  const double baseValue = argument.coefficient[0][0];
  GLMMTaylor increment = GLMMTaylor(1.0 / baseValue) *
    (argument - GLMMTaylor(baseValue));
  GLMMTaylor result(std::log(baseValue)), power(1.0);
  for(int order = 1; order <= 5; ++order){
    power = power * increment;
    result = result + GLMMTaylor((order % 2 ? 1.0 : -1.0) / order) * power;
  }
  return result;
}

inline GLMMTaylor glmmLog1p(const GLMMTaylor & argument){
  GLMMTaylor result = glmmLog(GLMMTaylor(1.0) + argument);
  result.coefficient[0][0] = std::log1p(argument.coefficient[0][0]);
  return result;
}

inline GLMMTaylor glmmLogNormalCDF(double predictor){
  GLMMTaylor result(R::pnorm(predictor, 0.0, 1.0, true, true));
  double hazard[4] = {};
  hazard[0] = std::exp(R::dnorm(predictor, 0.0, 1.0, true) -
    result.coefficient[0][0]);
  for(int order = 0; order < 3; ++order){
    double squaredHazard = 0.0;
    for(int part = 0; part <= order; ++part)
      squaredHazard += hazard[part] * hazard[order - part];
    hazard[order + 1] = (-predictor * hazard[order] -
      (order > 0 ? hazard[order - 1] : 0.0) - squaredHazard) / (order + 1);
  }
  for(int order = 1; order <= 4; ++order)
    result.coefficient[order][0] = hazard[order - 1] / order;
  return result;
}

inline GLMMTaylor glmmShapeLGamma(const GLMMTaylor & argument){
  GLMMTaylor result(R::lgammafn(argument.coefficient[0][0]));
  result.coefficient[0][1] = R::digamma(argument.coefficient[0][0]) *
    argument.coefficient[0][1];
  return result;
}

inline GLMMTaylor glmmOrdinalLogCDF(double value, const std::string & link,
    double direction){
  if(link == "probit"){
    GLMMTaylor result = glmmLogNormalCDF(value);
    for(int order = 1; order <= 4; order += 2)
      result.coefficient[order][0] *= direction;
    return result;
  }
  GLMMTaylor argument(value);
  argument.coefficient[1][0] = direction;
  return (value >= 0.0 ? GLMMTaylor() : argument) -
    glmmLog1p(glmmExp(GLMMTaylor(value >= 0.0 ? -1.0 : 1.0) * argument));
}

inline GLMMTaylor glmmOrdinalLogDensity(double outcome, double predictor,
    const arma::vec & thresholds, const std::string & link){
  const arma::uword category = static_cast<arma::uword>(outcome) - 1;
  if(category == 0)
    return glmmOrdinalLogCDF(thresholds(0) - predictor, link, -1.0);
  if(category == thresholds.n_elem)
    return glmmOrdinalLogCDF(predictor - thresholds(category - 1), link, 1.0);
  const double lower = thresholds(category - 1) - predictor;
  const double upper = thresholds(category) - predictor;
  const double width = thresholds(category) - thresholds(category - 1);
  if(link == "logit")
    return glmmOrdinalLogCDF(upper, link, -1.0) +
      glmmOrdinalLogCDF(-lower, link, 1.0) +
      GLMMTaylor(std::log(-std::expm1(-width)));
  const double midpoint = lower + width / 2.0;
  if(width * (1.0 + std::abs(midpoint)) < 0.01){
    // Integrated normal-density series avoids cancellation of narrow-interval
    // derivatives. The omitted eighth-order term is below double precision.
    GLMMTaylor center(midpoint);
    center.coefficient[1][0] = -1.0;
    const GLMMTaylor square = center * center;
    const double width2 = width * width;
    const GLMMTaylor correction = GLMMTaylor(width2 / 24.0) *
      (square - GLMMTaylor(1.0)) + GLMMTaylor(width2 * width2 / 1920.0) *
      (square * square - GLMMTaylor(6.0) * square + GLMMTaylor(3.0)) +
      GLMMTaylor(width2 * width2 * width2 / 322560.0) *
      (square * square * square - GLMMTaylor(15.0) * square * square +
        GLMMTaylor(45.0) * square - GLMMTaylor(15.0));
    return GLMMTaylor(std::log(width) - 0.5 * std::log(2.0 * arma::datum::pi)) -
      GLMMTaylor(0.5) * square + glmmLog1p(correction);
  }
  // Use survival probabilities in the upper tail, avoiding CDFs rounded to one.
  const GLMMTaylor larger = lower >= 0.0 ?
    glmmOrdinalLogCDF(-lower, link, 1.0) :
    glmmOrdinalLogCDF(upper, link, -1.0);
  const GLMMTaylor smaller = lower >= 0.0 ?
    glmmOrdinalLogCDF(-upper, link, 1.0) :
    glmmOrdinalLogCDF(lower, link, -1.0);
  const GLMMTaylor logRatio = smaller - larger;
  GLMMTaylor difference = GLMMTaylor(1.0) - glmmExp(logRatio);
  difference.coefficient[0][0] = -std::expm1(logRatio.coefficient[0][0]);
  return larger + glmmLog(difference);
}

inline GLMMTaylor glmmOrdinalCutpointScore(double threshold, double predictor,
    const std::string & link, const GLMMTaylor & logDensity){
  const double value = threshold - predictor;
  GLMMTaylor logPDF;
  if(link == "logit"){
    logPDF = glmmOrdinalLogCDF(value, link, -1.0) +
      glmmOrdinalLogCDF(-value, link, 1.0);
  }else{
    GLMMTaylor argument(value);
    argument.coefficient[1][0] = -1.0;
    logPDF = GLMMTaylor(-0.5 * std::log(2.0 * arma::datum::pi)) -
      GLMMTaylor(0.5) * argument * argument;
  }
  return glmmExp(logPDF - logDensity);
}

struct GLMMFamilySpec {
  std::string family;
  std::string link;
  double logParameter = 0.0;
  bool free = false;
  arma::mat zeroDesign;
  arma::vec zeroCoefficients;
  arma::vec zeroOffset;
  arma::vec thresholds;
};

struct GLMMFamilyEvaluation {
  double logDensity = 0.0;
  arma::vec logDensityRows;
  arma::vec mean;
  arma::vec score;
  arma::vec curvature;
  arma::vec curvatureD1;
  arma::vec curvatureD2;
  arma::vec parameterScore;
  arma::vec scoreParameter;
  arma::vec curvatureParameter;
};

class GLMMFamilyModel {
public:
  std::vector<GLMMFamilySpec> specifications;
  arma::uvec index;
  arma::vec trials;

  GLMMFamilyModel(const std::string & family, arma::uword observations,
      const Rcpp::List & suppliedSpecs = Rcpp::List(),
      const arma::uvec & suppliedIndex = arma::uvec(),
      const arma::vec & suppliedTrials = arma::vec()){
    if(suppliedSpecs.size() == 0){
      specifications.push_back({family, family == "binomial" ? "logit" : "log",
        0.0, false});
    }else{
      for(int familyPosition = 0; familyPosition < suppliedSpecs.size(); ++familyPosition){
        const Rcpp::List specification = suppliedSpecs[familyPosition];
        specifications.push_back({Rcpp::as<std::string>(specification["family"]),
          Rcpp::as<std::string>(specification["link"]),
          Rcpp::as<double>(specification["logParameter"]),
          Rcpp::as<bool>(specification["free"])});
        if(specification.containsElementNamed("thresholds"))
          specifications.back().thresholds =
            Rcpp::as<arma::vec>(specification["thresholds"]);
        if(specification.containsElementNamed("zeroDesign")){
          auto & stored = specifications.back();
          stored.zeroDesign = Rcpp::as<arma::mat>(specification["zeroDesign"]);
          stored.zeroCoefficients = Rcpp::as<arma::vec>(specification["zeroCoefficients"]);
          stored.zeroOffset = Rcpp::as<arma::vec>(specification["zeroOffset"]);
          if(stored.family != "ziPoisson" || !stored.free ||
             stored.zeroDesign.n_rows != observations ||
             stored.zeroDesign.n_cols != stored.zeroCoefficients.n_elem ||
             stored.zeroDesign.n_cols == 0 || stored.zeroOffset.n_elem != observations ||
             !stored.zeroDesign.is_finite() || !stored.zeroCoefficients.is_finite() ||
             !stored.zeroOffset.is_finite())
            Rcpp::stop("Invalid native zero-component design or coefficients.");
        }
      }
    }
    for(const auto & specification : specifications){
      const bool binomial = specification.family == "binomial" &&
        (specification.link == "logit" || specification.link == "probit" ||
         specification.link == "cloglog");
      const bool logLinked = (specification.family == "poisson" ||
        specification.family == "ziPoisson" || specification.family == "nbinom2" ||
        specification.family == "Gamma") &&
        specification.link == "log";
      const bool gaussian = specification.family == "gaussian" && specification.link == "identity";
      const bool ordinal = specification.family == "ordinal" &&
        (specification.link == "logit" || specification.link == "probit");
      if(!binomial && !logLinked && !gaussian && !ordinal)
        Rcpp::stop("Unsupported native Laplace family/link combination.");
      if(ordinal && (specification.thresholds.n_elem == 0 ||
         !specification.thresholds.is_finite() ||
         arma::any(arma::diff(specification.thresholds) <= 0.0)))
        Rcpp::stop("Native ordinal evaluation requires finite, strictly increasing thresholds.");
      if(!std::isfinite(specification.logParameter) ||
         std::abs(specification.logParameter) >
           std::log(std::numeric_limits<double>::max()) / 4.0 ||
         ((specification.family == "binomial" || specification.family == "poisson") &&
          specification.free))
        Rcpp::stop("Invalid native family shape/dispersion specification.");
    }
    index = suppliedIndex.n_elem ? suppliedIndex : arma::uvec(observations, arma::fill::zeros);
    trials = suppliedTrials.n_elem ? suppliedTrials : arma::vec(observations, arma::fill::ones);
    if(index.n_elem != observations || trials.n_elem != observations ||
       !trials.is_finite() || arma::any(trials <= 0) ||
       arma::any(trials != arma::floor(trials)) ||
       (observations && index.max() >= specifications.size()))
      Rcpp::stop("Native family indices and trials have incompatible or invalid values.");
    for(std::size_t familyPosition = 0; familyPosition < specifications.size(); ++familyPosition){
      if(specifications[familyPosition].free &&
         !arma::any(index == familyPosition))
        Rcpp::stop("Estimated family parameters require observations of that family.");
    }
  }

  arma::uword parameterCount() const {
    arma::uword result = 0;
    for(const auto & specification : specifications)
      if(specification.free) result += specification.family == "ordinal" ?
        specification.thresholds.n_elem : (specification.zeroDesign.n_cols ?
          specification.zeroDesign.n_cols : 1);
    return result;
  }

  bool update(const arma::vec & parameters){
    if(parameters.n_elem != parameterCount())
      Rcpp::stop("Native family parameter vector has incompatible length.");
    arma::uword parameterPosition = 0;
    for(auto & specification : specifications){
      if(specification.free){
        if(specification.family == "ordinal"){
          arma::vec thresholds(specification.thresholds.n_elem);
          thresholds(0) = parameters(parameterPosition++);
          if(!std::isfinite(thresholds(0))) return false;
          for(arma::uword cutpoint = 1; cutpoint < thresholds.n_elem; ++cutpoint){
            const double gap = std::exp(parameters(parameterPosition++));
            thresholds(cutpoint) = thresholds(cutpoint - 1) + gap;
            if(!std::isfinite(gap) || !std::isfinite(thresholds(cutpoint)) ||
               thresholds(cutpoint) <= thresholds(cutpoint - 1)) return false;
          }
          specification.thresholds = thresholds;
        }else if(specification.zeroDesign.n_cols){
          const arma::uword count = specification.zeroDesign.n_cols;
          specification.zeroCoefficients = parameters.subvec(parameterPosition,
            parameterPosition + count - 1);
          parameterPosition += count;
          if(!specification.zeroCoefficients.is_finite()) return false;
        }else specification.logParameter = parameters(parameterPosition++);
      }
      if(!std::isfinite(specification.logParameter) ||
         std::abs(specification.logParameter) >
           std::log(std::numeric_limits<double>::max()) / 4.0) return false;
    }
    return true;
  }

  void validateResponse(const arma::vec & response) const {
    if(response.n_elem != index.n_elem || !response.is_finite())
      Rcpp::stop("Native Laplace responses must be finite with compatible dimensions.");
    for(arma::uword recordPosition = 0; recordPosition < response.n_elem; ++recordPosition){
      const auto & specification = specifications[index(recordPosition)];
      const double outcome = response(recordPosition);
      if(specification.family == "binomial" &&
         (outcome < 0 || outcome > trials(recordPosition) || outcome != std::floor(outcome)))
        Rcpp::stop("Native Laplace binomial responses must be integer successes between zero and trials (0/1 for binary data).");
      if((specification.family == "poisson" || specification.family == "ziPoisson" ||
          specification.family == "nbinom2") &&
         (outcome < 0 || outcome != std::floor(outcome)))
        Rcpp::stop("Native Laplace count responses must be non-negative integers.");
      if(specification.family == "Gamma" && outcome <= 0)
        Rcpp::stop("Native Laplace Gamma responses must be strictly positive.");
      if(specification.family == "ordinal" &&
         (outcome < 1 || outcome > specification.thresholds.n_elem + 1 ||
          outcome != std::floor(outcome)))
        Rcpp::stop("Native ordinal responses must be integer category codes between 1 and the number of thresholds plus one.");
      if(specification.family != "binomial" && trials(recordPosition) != 1)
        Rcpp::stop("Trials are only supported for binomial observations.");
    }
  }

  GLMMFamilyEvaluation evaluate(const arma::vec & response, const arma::vec & predictor,
                                const bool retainLogDensityRows = false) const {
    GLMMFamilyEvaluation result;
    if (retainLogDensityRows) result.logDensityRows.set_size(response.n_elem);
    result.mean.set_size(response.n_elem);
    result.score.set_size(response.n_elem);
    result.curvature.set_size(response.n_elem);
    result.curvatureD1.set_size(response.n_elem);
    result.curvatureD2.set_size(response.n_elem);
    result.parameterScore.zeros(response.n_elem);
    result.scoreParameter.zeros(response.n_elem);
    result.curvatureParameter.zeros(response.n_elem);
    for(arma::uword recordPosition = 0; recordPosition < response.n_elem; ++recordPosition){
      const double linearValue = predictor(recordPosition);
      const double outcome = response(recordPosition);
      const auto & specification = specifications[index(recordPosition)];
      if(!std::isfinite(linearValue)){
        result.logDensity = -INFINITY;
        return result;
      }
      GLMMTaylor eta(linearValue);
      eta.coefficient[1][0] = 1.0;
      const double parameterValue = specification.zeroDesign.n_cols ?
        arma::dot(specification.zeroDesign.row(recordPosition),
          specification.zeroCoefficients) + specification.zeroOffset(recordPosition) :
        specification.logParameter;
      if(!std::isfinite(parameterValue)){
        result.logDensity = -INFINITY;
        return result;
      }
      GLMMTaylor logParameter(parameterValue);
      logParameter.coefficient[0][1] = 1.0;
      GLMMTaylor density;
      if(specification.family == "binomial"){
        GLMMTaylor logSuccess, logFailure;
        if(specification.link == "logit"){
          const GLMMTaylor smallTail = glmmExp(GLMMTaylor(linearValue >= 0 ? -1.0 : 1.0) * eta);
          const GLMMTaylor logDenominator = glmmLog1p(smallTail);
          logSuccess = (linearValue >= 0 ? GLMMTaylor() : eta) - logDenominator;
          logFailure = (linearValue >= 0 ? GLMMTaylor() - eta : GLMMTaylor()) - logDenominator;
          result.mean(recordPosition) = linearValue >= 0 ?
            1.0 / (1.0 + std::exp(-linearValue)) :
            std::exp(linearValue) / (1.0 + std::exp(linearValue));
        }else if(specification.link == "probit"){
          logSuccess = glmmLogNormalCDF(linearValue);
          logFailure = glmmLogNormalCDF(-linearValue);
          for(int derivativeOrder = 1; derivativeOrder <= 4; derivativeOrder += 2)
            logFailure.coefficient[derivativeOrder][0] *= -1.0;
          result.mean(recordPosition) = R::pnorm(linearValue, 0.0, 1.0, true, false);
        }else{
          if(linearValue > std::log(std::numeric_limits<double>::max())){
            result.logDensity = -INFINITY;
            return result;
          }
          const GLMMTaylor intensity = glmmExp(eta);
          logFailure = GLMMTaylor() - intensity;
          if(intensity.coefficient[0][0] < std::log(2.0)){
            GLMMTaylor successProbability = GLMMTaylor(1.0) - glmmExp(logFailure);
            successProbability.coefficient[0][0] = -std::expm1(-intensity.coefficient[0][0]);
            logSuccess = glmmLog(successProbability);
          }else{
            logSuccess = glmmLog1p(GLMMTaylor() - glmmExp(logFailure));
          }
          result.mean(recordPosition) = -std::expm1(-intensity.coefficient[0][0]);
        }
        density = GLMMTaylor(R::lchoose(trials(recordPosition), outcome));
        if(outcome > 0) density = density + GLMMTaylor(outcome) * logSuccess;
        if(outcome < trials(recordPosition))
          density = density + GLMMTaylor(trials(recordPosition) - outcome) * logFailure;
      }else if(specification.family == "ordinal"){
        density = glmmOrdinalLogDensity(outcome, linearValue,
          specification.thresholds, specification.link);
        double mean = 1.0;
        for(const double threshold : specification.thresholds){
          const double value = linearValue - threshold;
          mean += specification.link == "probit" ?
            R::pnorm(value, 0.0, 1.0, true, false) :
            (value >= 0.0 ? 1.0 / (1.0 + std::exp(-value)) :
              std::exp(value) / (1.0 + std::exp(value)));
        }
        result.mean(recordPosition) = mean;
      }else if(specification.family == "gaussian"){
        const GLMMTaylor discrepancy = GLMMTaylor(outcome) - eta;
        density = GLMMTaylor(-0.5) * (GLMMTaylor(std::log(2.0 * arma::datum::pi)) +
          logParameter + discrepancy * discrepancy * glmmExp(GLMMTaylor() - logParameter));
        result.mean(recordPosition) = linearValue;
      }else{
        if(linearValue > std::log(std::numeric_limits<double>::max())){
          result.logDensity = -INFINITY;
          return result;
        }
        const GLMMTaylor intensity = glmmExp(eta);
        result.mean(recordPosition) = intensity.coefficient[0][0];
        if(specification.family == "poisson"){
          density = GLMMTaylor(outcome) * eta - intensity -
            GLMMTaylor(R::lgammafn(outcome + 1.0));
        }else if(specification.family == "ziPoisson"){
          GLMMTaylor logZeroProbability;
          GLMMTaylor logCountProbability;
          if(parameterValue >= 0.0){
            const GLMMTaylor tail =
              glmmLog1p(glmmExp(GLMMTaylor() - logParameter));
            logZeroProbability = GLMMTaylor() - tail;
            logCountProbability = GLMMTaylor() - logParameter - tail;
          }else{
            const GLMMTaylor tail = glmmLog1p(glmmExp(logParameter));
            logZeroProbability = logParameter - tail;
            logCountProbability = GLMMTaylor() - tail;
          }
          const double logCountMean =
            logCountProbability.coefficient[0][0] + linearValue;
          if(logCountMean > std::log(std::numeric_limits<double>::max())){
            result.logDensity = -INFINITY;
            return result;
          }
          result.mean(recordPosition) = std::exp(logCountMean);
          if(outcome == 0.0){
            const GLMMTaylor zeroFromInflation = logZeroProbability;
            const GLMMTaylor zeroFromCount = logCountProbability - intensity;
            const double shift = std::max(zeroFromInflation.coefficient[0][0],
              zeroFromCount.coefficient[0][0]);
            density = GLMMTaylor(shift) +
              glmmLog(glmmExp(zeroFromInflation - GLMMTaylor(shift)) +
                glmmExp(zeroFromCount - GLMMTaylor(shift)));
          }else{
            density = logCountProbability + GLMMTaylor(outcome) * eta -
              intensity - GLMMTaylor(R::lgammafn(outcome + 1.0));
          }
        }else if(specification.family == "nbinom2"){
          const GLMMTaylor shape = glmmExp(logParameter);
          // log(theta + mu) via log-sum-exp avoids overflow in the sum.
          const GLMMTaylor logTotal = linearValue >= specification.logParameter ?
            eta + glmmLog1p(glmmExp(logParameter - eta)) :
            logParameter + glmmLog1p(glmmExp(eta - logParameter));
          density = glmmShapeLGamma(GLMMTaylor(outcome) + shape) -
            glmmShapeLGamma(shape) - GLMMTaylor(R::lgammafn(outcome + 1.0)) +
            shape * (logParameter - logTotal) + GLMMTaylor(outcome) * (eta - logTotal);
          density.coefficient[0][0] = R::dnbinom_mu(outcome,
            shape.coefficient[0][0], intensity.coefficient[0][0], true);
        }else{
          const GLMMTaylor shape = glmmExp(GLMMTaylor() - logParameter);
          density = (shape - GLMMTaylor(1.0)) * GLMMTaylor(std::log(outcome)) -
            GLMMTaylor(outcome) * glmmExp(GLMMTaylor() - eta - logParameter) -
            glmmShapeLGamma(shape) - shape * (eta + logParameter);
          density.coefficient[0][0] = R::dgamma(outcome, shape.coefficient[0][0],
            intensity.coefficient[0][0] * std::exp(specification.logParameter), true);
        }
      }
      result.logDensity += density.coefficient[0][0];
      if (retainLogDensityRows) result.logDensityRows(recordPosition) = density.coefficient[0][0];
      result.score(recordPosition) = density.coefficient[1][0];
      result.curvature(recordPosition) = -2.0 * density.coefficient[2][0];
      result.curvatureD1(recordPosition) = -6.0 * density.coefficient[3][0];
      result.curvatureD2(recordPosition) = -24.0 * density.coefficient[4][0];
      result.parameterScore(recordPosition) = density.coefficient[0][1];
      result.scoreParameter(recordPosition) = density.coefficient[1][1];
      result.curvatureParameter(recordPosition) = -2.0 * density.coefficient[2][1];
    }
    return result;
  }
};

inline GLMMFamilyEvaluation evaluateGLMMFamily(const arma::vec & response,
    const arma::vec & predictor, const std::string & family){
  return GLMMFamilyModel(family, response.n_elem).evaluate(response, predictor);
}

} // namespace sommer

#endif
