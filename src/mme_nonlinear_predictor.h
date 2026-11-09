#ifndef SOMMER_MME_NONLINEAR_PREDICTOR_H
#define SOMMER_MME_NONLINEAR_PREDICTOR_H

namespace sommer {

struct NonlinearCurveEvaluation {
  arma::vec mean;
  arma::mat jacobian;
};

inline arma::uword nonlinearParameterCount(const std::string & curve){
  if(curve == "logistic") return 3;
  if(curve == "exponential" || curve == "saturation") return 2;
  Rcpp::stop("Unsupported native nonlinear curve.");
  return 0;
}

// Derivatives are with respect to log-positive and unconstrained coordinates.
inline bool nonlinearCurve(const std::string & curve, const arma::vec & input,
    const arma::vec & parameters, NonlinearCurveEvaluation & result,
    const bool derivatives){
  const arma::uword count = nonlinearParameterCount(curve);
  if(parameters.n_elem != count || !input.is_finite())
    Rcpp::stop("Native nonlinear curve has incompatible or non-finite inputs.");
  if(curve != "logistic" && arma::any(input < 0))
    Rcpp::stop("Exponential and saturation curve inputs must be non-negative.");
  if(!parameters.is_finite()) return false;
  const double amplitude = std::exp(parameters(0));
  const double shape = std::exp(parameters(1));
  if(!std::isfinite(amplitude) || amplitude <= 0 ||
     !std::isfinite(shape) || shape <= 0) return false;
  result.mean.set_size(input.n_elem);
  if(derivatives) result.jacobian.set_size(input.n_elem, count);
  for(arma::uword row = 0; row < input.n_elem; ++row){
    if(row % 4096 == 0) Rcpp::checkUserInterrupt();
    const double x = input(row);
    double mean, shapeDerivative, locationDerivative = 0;
    if(curve == "exponential"){
      const double scaled = shape * x;
      if(!std::isfinite(scaled)) return false;
      mean = amplitude * std::exp(-scaled);
      shapeDerivative = -scaled * mean;
    }else if(curve == "logistic"){
      const double scaled = shape * (x - parameters(2));
      if(!std::isfinite(scaled)) return false;
      const double tail = std::exp(-std::abs(scaled));
      const double probability = scaled >= 0 ? 1/(1+tail) : tail/(1+tail);
      const double slope = amplitude * tail / ((1+tail)*(1+tail));
      mean = amplitude * probability;
      shapeDerivative = slope * scaled;
      locationDerivative = -slope * shape;
    }else{
      // Ratios avoid overflow in K+x for large finite input values.
      const double fraction = x > shape ? 1/(1+shape/x) :
        (x/shape)/(1+x/shape);
      mean = amplitude * fraction;
      shapeDerivative = -mean * (1-fraction);
    }
    result.mean(row) = mean;
    if(derivatives){
      result.jacobian(row, 0) = mean;
      result.jacobian(row, 1) = shapeDerivative;
      if(count == 3) result.jacobian(row, 2) = locationDerivative;
    }
  }
  return result.mean.is_finite() &&
    (!derivatives || result.jacobian.is_finite());
}

} // namespace sommer
#endif
