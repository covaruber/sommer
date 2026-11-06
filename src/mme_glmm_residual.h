#ifndef SOMMER_MME_GLMM_RESIDUAL_H
#define SOMMER_MME_GLMM_RESIDUAL_H

#include "mme_direct_residual.h"
#include <memory>

namespace sommer {

class GaussianObservationResidual {
  std::unique_ptr<DirectResidualCovariance> covariance;
  arma::vec parameters;
  arma::uvec free;
  arma::vec cached;
  bool ready = false;
public:
  arma::uvec rows;
  arma::mat precision;
  double logDet = 0.0;

  GaussianObservationResidual(const Rcpp::List & model, const double memoryMB){
    if(!model.size()) return;
    if(!std::isfinite(memoryMB) || memoryMB <= 0)
      Rcpp::stop("Gaussian residual memory budget must be positive and finite.");
    rows = Rcpp::as<arma::uvec>(model["rows"]) - 1;
    if(!rows.n_elem || 8.0 * rows.n_elem * rows.n_elem * 16 > memoryMB * 1048576.0)
      Rcpp::stop("Gaussian observation residual workspace exceeds denseMemoryMB.");
    const Rcpp::List descriptor = model["descriptor"];
    if(Rf_length(model["block"]) != static_cast<R_xlen_t>(rows.n_elem))
      Rcpp::stop("Gaussian residual mappings must match Gaussian rows.");
    parameters = Rcpp::as<arma::vec>(descriptor["par"]);
    free = arma::find(Rcpp::as<arma::uvec>(descriptor["free"]) != 0);
    covariance.reset(new DirectResidualCovariance(descriptor,
      Rcpp::as<arma::uvec>(model["block"]), Rcpp::as<arma::uvec>(model["index"]),
      Rcpp::as<arma::sp_mat>(model["W"]), Rcpp::as<bool>(model["weighted"])));
  }

  bool active() const { return static_cast<bool>(covariance); }
  arma::uword parameterCount() const { return free.n_elem; }

  bool update(const arma::vec & values){
    if(values.n_elem != free.n_elem)
      Rcpp::stop("Gaussian residual parameter vector has incompatible length.");
    if(!active()) return true;
    if(ready && arma::all(values == cached)) return true;
    ready = false;
    parameters.elem(free) = values;
    if(!parameters.is_finite()) return false;
    const arma::mat matrix = covariance->evaluate(parameters);
    arma::mat upper;
    if(!matrix.is_finite() || !arma::chol(upper, matrix)) return false;
    logDet = 2 * arma::accu(arma::log(upper.diag()));
    arma::mat identity(matrix.n_rows, matrix.n_rows, arma::fill::eye);
    arma::mat intermediate;
    if(!arma::solve(intermediate, arma::trimatl(upper.t()), identity,
         arma::solve_opts::fast) ||
       !arma::solve(precision, arma::trimatu(upper), intermediate,
         arma::solve_opts::fast)){
      return false;
    }
    precision = arma::symmatu(precision);
    ready = precision.is_finite();
    cached = values;
    return ready;
  }

  arma::mat derivative(arma::uword coordinate) const {
    return covariance->evaluate(parameters, free(coordinate));
  }

  bool initialize(){ return update(parameters.elem(free)); }
  arma::mat matrix() const { return covariance->evaluate(parameters); }
};

} // namespace sommer
#endif
