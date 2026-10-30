#ifndef SOMMER_MME_GLMM_MULTINOMIAL_H
#define SOMMER_MME_GLMM_MULTINOMIAL_H

#include <RcppArmadillo.h>
#include <algorithm>
#include <cmath>

namespace sommer {

struct MultinomialEvaluation {
  double logDensity;
  arma::vec probability;
  arma::vec score;
  arma::mat curvature;
  arma::mat curvatureD1;
  arma::mat curvatureD2;
};

inline arma::vec multinomialCurvatureContraction(const arma::vec & probability,
    const arma::mat & leverage){
  const arma::uword dimension = probability.n_elem - 1;
  arma::vec result(dimension, arma::fill::zeros);
  const arma::vec nonreference = probability.head(dimension);
  // Contract the symmetric third cumulant without a category-cubed tensor.
  for(arma::uword category = 0; category <= dimension; ++category){
    arma::vec centered = -nonreference;
    if(category < dimension){
      double complement = 0.0;
      for(arma::uword other = 0; other <= dimension; ++other){
        if(other != category) complement += probability(other);
      }
      centered(category) = complement;
    }
    result += 0.5 * probability(category) *
      arma::dot(centered, leverage * centered) * centered;
  }
  return result;
}

inline MultinomialEvaluation multinomialEvaluation(const arma::vec & predictor,
    const arma::uword category, const arma::vec & direction,
    const int derivativeOrder = 2){
  const arma::uword dimension = predictor.n_elem;
  const arma::uword categories = dimension + 1;
  MultinomialEvaluation result;
  const double maximum = std::max(0.0, predictor.max());
  const arma::uword dominant = maximum > 0.0 ? predictor.index_max() : dimension;
  arma::vec centered(categories);
  centered.head(dimension) = predictor - maximum;
  centered(dimension) = -maximum;
  const arma::vec weights = arma::exp(centered);
  double remainingWeight = 0.0;
  for(arma::uword position = 0; position < categories; ++position){
    if(position != dominant) remainingWeight += weights(position);
  }
  const double logNormalizer = std::log1p(remainingWeight);
  result.probability = arma::exp(centered - logNormalizer);
  result.logDensity = centered(category) - logNormalizer;
  if(!std::isfinite(result.logDensity)){
    Rcpp::stop("Multinomial log probability exceeds the finite numerical range.");
  }
  const arma::vec probability = result.probability.head(dimension);
  result.score = -probability;
  if(category < dimension){
    double complement = 0.0;
    for(arma::uword position = 0; position < categories; ++position){
      if(position != category) complement += result.probability(position);
    }
    result.score(category) = complement;
  }
  if(derivativeOrder < 0) return result;
  result.curvature.zeros(dimension, dimension);
  if(derivativeOrder > 0) result.curvatureD1.zeros(dimension, dimension);
  if(derivativeOrder > 1) result.curvatureD2.zeros(dimension, dimension);
  const double meanDirection = arma::dot(probability, direction);
  arma::vec covarianceDirection(dimension, arma::fill::zeros);
  double varianceDirection = 0.0;
  for(arma::uword position = 0; position < categories; ++position){
    arma::vec centeredIndicator = -probability;
    if(position < dimension) centeredIndicator(position) += 1.0;
    const double centeredDirection =
      (position < dimension ? direction(position) : 0.0) - meanDirection;
    const arma::mat contribution = result.probability(position) *
      (centeredIndicator * centeredIndicator.t());
    result.curvature += contribution;
    if(derivativeOrder == 0) continue;
    result.curvatureD1 += centeredDirection * contribution;
    if(derivativeOrder > 1){
      result.curvatureD2 += centeredDirection * centeredDirection * contribution;
    }
    covarianceDirection += result.probability(position) * centeredDirection *
      centeredIndicator;
    varianceDirection += result.probability(position) *
      centeredDirection * centeredDirection;
  }
  // Third/fourth categorical cumulants contract Hessian derivatives with a
  // supplied direction, avoiding observation-by-category^4 tensor storage.
  if(derivativeOrder > 1){
    result.curvatureD2 -= varianceDirection * result.curvature +
      2.0 * (covarianceDirection * covarianceDirection.t());
  }
  if(!result.probability.is_finite() || !result.score.is_finite() ||
     !result.curvature.is_finite() || !result.curvatureD1.is_finite() ||
     !result.curvatureD2.is_finite()){
    Rcpp::stop("Multinomial block derivatives exceed the finite numerical range.");
  }
  return result;
}

} // namespace sommer

#endif
