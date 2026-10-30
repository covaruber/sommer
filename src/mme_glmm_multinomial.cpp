#include "mme_glmm_multinomial.h"

// [[Rcpp::export(rng = false)]]
Rcpp::List sommer_multinomial_family_eval(const arma::vec & response,
    const arma::mat & predictor, const arma::mat & direction,
    const double denseMemoryMB = 1000.0){
  const arma::uword observations = predictor.n_rows;
  const arma::uword dimension = predictor.n_cols;
  if(observations == 0 || dimension == 0 ||
     response.n_elem != observations ||
     direction.n_rows != observations || direction.n_cols != dimension){
    Rcpp::stop("Multinomial responses, predictors and directions require compatible nonempty dimensions.");
  }
  if(!response.is_finite() || !predictor.is_finite() || !direction.is_finite()){
    Rcpp::stop("Multinomial responses, predictors and directions must be finite.");
  }
  if(!std::isfinite(denseMemoryMB) || denseMemoryMB <= 0.0){
    Rcpp::stop("Multinomial denseMemoryMB must be positive and finite.");
  }
  const double d = static_cast<double>(dimension);
  const double n = static_cast<double>(observations);
  const double bytes = 8.0 * (2.0 * n * (3.0 * d * d + 2.0 * d + 2.0) +
    12.0 * d * d + 12.0 * d + 8.0);
  if(!std::isfinite(bytes) || bytes > denseMemoryMB * 1048576.0){
    Rcpp::stop("Multinomial block evaluation exceeds denseMemoryMB.");
  }
  const double categories = d + 1.0;
  for(arma::uword row = 0; row < observations; ++row){
    if(response(row) < 1.0 || response(row) > categories ||
       response(row) != std::floor(response(row))){
      Rcpp::stop("Multinomial response codes must be integers from 1 through the declared category count.");
    }
  }
  arma::vec logDensity(observations);
  arma::mat probability(observations, dimension + 1);
  arma::mat score(observations, dimension);
  arma::cube curvature(dimension, dimension, observations);
  arma::cube curvatureD1(dimension, dimension, observations);
  arma::cube curvatureD2(dimension, dimension, observations);
  for(arma::uword row = 0; row < observations; ++row){
    if(row % 128 == 0) Rcpp::checkUserInterrupt();
    const auto value = sommer::multinomialEvaluation(predictor.row(row).t(),
      static_cast<arma::uword>(response(row) - 1.0), direction.row(row).t());
    logDensity(row) = value.logDensity;
    probability.row(row) = value.probability.t();
    score.row(row) = value.score.t();
    curvature.slice(row) = value.curvature;
    curvatureD1.slice(row) = value.curvatureD1;
    curvatureD2.slice(row) = value.curvatureD2;
  }
  const double total = arma::accu(logDensity);
  if(!std::isfinite(total)){
    Rcpp::stop("Multinomial total log likelihood exceeds the finite numerical range.");
  }
  return Rcpp::List::create(Rcpp::_["logDensity"]=total,
    Rcpp::_["logDensityRows"]=logDensity, Rcpp::_["probability"]=probability,
    Rcpp::_["score"]=score, Rcpp::_["curvature"]=curvature,
    Rcpp::_["curvatureD1"]=curvatureD1, Rcpp::_["curvatureD2"]=curvatureD2);
}
