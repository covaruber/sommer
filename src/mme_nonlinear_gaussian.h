#ifndef SOMMER_MME_NONLINEAR_GAUSSIAN_H
#define SOMMER_MME_NONLINEAR_GAUSSIAN_H

#include "mme_glmm_system.h"
#include "mme_nonlinear_predictor.h"
#include "mme_nonlinear_covariance.h"
#include "mme_nonlinear_memory.h"

namespace sommer {

class NonlinearGaussianObjective {
  arma::sp_mat fixed, random;
  arma::vec response, offset, input;
  std::string curve;
  NonlinearCovariance covarianceModel;
  LaplacePrecision & precision;
  RandomLaplaceSystem system;
  double memoryMB;
  int evaluations = 0;

public:
  NonlinearGaussianObjective(const arma::sp_mat & x, const arma::sp_mat & z,
      const arma::vec & y, const arma::vec & o, const arma::vec & a,
      const std::string & function, const Rcpp::List & descriptors,
      const Rcpp::List & relationships, double memory)
    : fixed(x), random(z), response(y), offset(o), input(a), curve(function),
      covarianceModel(z,descriptors,relationships,memory),
      precision(covarianceModel.prior()), memoryMB(memory) {
    nonlinearParameterCount(curve);
    if(!std::isfinite(memoryMB) || memoryMB <= 0)
      Rcpp::stop("Nonlinear memory budget must be positive and finite.");
    if(!y.n_elem || x.n_rows != y.n_elem || z.n_rows != y.n_elem ||
       o.n_elem != y.n_elem || a.n_elem != y.n_elem ||
       !x.is_finite() || !z.is_finite() || !y.is_finite() ||
       !o.is_finite() || !a.is_finite())
      Rcpp::stop("Nonlinear Gaussian inputs must have compatible dimensions and finite values.");
    if(8.0*(40.0*y.n_elem+12.0*z.n_cols) > memoryMB*1048576.0)
      Rcpp::stop("Nonlinear Gaussian workspace exceeds execution$denseMemoryMB.");
  }

  Rcpp::List evaluate(const arma::vec & parameters, bool gradient){
    ++evaluations;
    const arma::uword p = fixed.n_cols, k = nonlinearParameterCount(curve);
    const arma::uword c = covarianceModel.parameterCount(), n = response.n_elem;
    if(parameters.n_elem != p+k+c+1)
      Rcpp::stop("Nonlinear Gaussian parameter vector has incompatible length.");
    auto invalid = [](const std::string & reason){
      return Rcpp::List::create(Rcpp::_["valid"]=false, Rcpp::_["reason"]=reason);
    };
    if(!parameters.is_finite()) return invalid("non-finite parameter proposal");
    arma::vec nonlinearMean;
    arma::mat nonlinearJacobian;
    NonlinearCurveEvaluation predictor;
    if(!nonlinearCurve(curve,input,parameters.subvec(p,p+k-1),predictor,gradient))
      return invalid("non-finite nonlinear mean or derivative");
    nonlinearMean=predictor.mean;
    if(gradient) nonlinearJacobian=predictor.jacobian;
    const double logVariance = parameters(parameters.n_elem-1);
    const double variance = std::exp(logVariance);
    const double weight = 1/variance;
    if(!std::isfinite(variance) || variance <= 0 || !std::isfinite(weight))
      return invalid("residual variance outside numerical range");
    arma::vec covariance(c);
    if(c) covariance = parameters.subvec(p+k,p+k+c-1);
    std::string reason;
    if(!covarianceModel.update(covariance, reason)) return invalid(reason);
    random=covarianceModel.design;
    const arma::uword q=random.n_cols;
    if(precision.matrix.n_rows != q)
      Rcpp::stop("Nonlinear random design and covariance layout are incompatible.");
    arma::vec mean = offset + nonlinearMean;
    if(p) mean += fixed * parameters.head(p);
    const arma::vec residual = response-mean;
    if(q==0){
      const arma::vec alpha=weight*residual;
      const double likelihood=-0.5*(n*(std::log(2*arma::datum::pi)+logVariance)+
        weight*arma::dot(residual,residual));
      if(!alpha.is_finite() || !std::isfinite(likelihood))
        return invalid("non-finite Gaussian likelihood or score");
      Rcpp::List result=Rcpp::List::create(Rcpp::_["valid"]=true,
        Rcpp::_["reason"]="evaluated",Rcpp::_["logLik"]=likelihood,
        Rcpp::_["mode"]=arma::vec(),Rcpp::_["fitted"]=mean,
        Rcpp::_["latentMode"]=arma::vec(),Rcpp::_["fixedMean"]=mean,
        Rcpp::_["modeScore"]=0.0);
      if(gradient){
        arma::vec score(parameters.n_elem,arma::fill::zeros);
        if(p) score.head(p)=fixed.t()*alpha;
        score.subvec(p,p+k-1)=nonlinearJacobian.t()*alpha;
        score(score.n_elem-1)=0.5*(-static_cast<double>(n)+
          weight*arma::dot(residual,residual));
        if(!score.is_finite()) return invalid("non-finite Gaussian score");
        result["gradient"]=score;
      }
      result["work"]=Rcpp::List::create(Rcpp::_["evaluations"]=evaluations,
        Rcpp::_["factorizations"]=0,Rcpp::_["contractionSolves"]=0,
        Rcpp::_["descriptorUpdates"]=0);
      return result;
    }
    arma::sp_mat curvature(n,n);
    curvature.diag().fill(weight);
    nonlinearSparseAdmission(random,precision.matrix,curvature,memoryMB);
    system.prepare(precision.matrix, random, precision.logDet);
    if(!system.factorize(curvature) || !system.determinant())
      return invalid("non-positive-definite Gaussian random system");
    const arma::vec mode = system.solve(arma::mat(weight*(random.t()*residual)));
    const arma::vec error = residual-random*mode;
    const arma::vec alpha = weight*error;
    const double quadratic = weight*arma::dot(error,error) +
      arma::dot(mode, precision.matrix*mode);
    const double likelihood = -0.5*(n*(std::log(2*arma::datum::pi)+logVariance) +
      system.logDetRatio + quadratic);
    if(!mode.is_finite() || !alpha.is_finite() || !std::isfinite(likelihood))
      return invalid("non-finite Gaussian likelihood or mode");
    Rcpp::List result = Rcpp::List::create(Rcpp::_["valid"]=true,
      Rcpp::_["reason"]="evaluated", Rcpp::_["logLik"]=likelihood,
      Rcpp::_["mode"]=covarianceModel.restore(mode), Rcpp::_["fitted"]=response-error,
      Rcpp::_["latentMode"]=mode,
      Rcpp::_["fixedMean"]=mean,
      Rcpp::_["modeScore"]=arma::norm(random.t()*alpha-precision.matrix*mode,"inf"));
    if(gradient){
      arma::vec score(parameters.n_elem, arma::fill::zeros);
      if(p) score.head(p) = fixed.t()*alpha;
      score.subvec(p,p+k-1) = nonlinearJacobian.t()*alpha;
      const double available = memoryMB*1048576.0 -
        8.0*(40.0*n+12.0*q);
      const double panelColumns = std::floor(available/(8.0*(3.0*q+n)));
      if(panelColumns < 1)
        Rcpp::stop("Nonlinear score solve workspace exceeds execution$denseMemoryMB.");
      const arma::uword panelSize = static_cast<arma::uword>(std::min(32.0,panelColumns));
      const auto derivatives = covarianceModel.derivatives();
      for(arma::uword coordinate = 0; coordinate < derivatives.size(); ++coordinate){
        Rcpp::checkUserInterrupt();
        const auto & derivative = derivatives[coordinate];
        const arma::uword end = derivative.start+derivative.block.n_rows-1;
        const arma::vec localMode = mode.subvec(derivative.start,end);
        double trace = 0;
        for(arma::uword first = 0; derivative.block.n_nonzero &&
            first < derivative.block.n_cols; first += panelSize){
          const arma::uword last = std::min(first+panelSize-1,derivative.block.n_cols-1);
          arma::mat rhs(q,last-first+1,arma::fill::zeros);
          rhs.rows(derivative.start,end) = arma::mat(derivative.block.cols(first,last));
          const arma::mat solved = system.solve(rhs);
          ++system.contractionSolves;
          for(arma::uword column = first; column <= last; ++column)
            trace += solved(derivative.start+column,column-first);
        }
        score(p+k+coordinate) = 0.5*(derivative.logDet-trace -
          arma::dot(localMode,derivative.block*localMode));
        if(covarianceModel.hasReducedRank()){
          const arma::sp_mat change=covarianceModel.designDerivative(coordinate);
          if(change.n_nonzero){
            double designTrace=0;
            for(arma::uword first=0;first<n;first+=panelSize){
              const arma::uword last=std::min(first+panelSize-1,n-1);
              const arma::mat panel(random.rows(first,last).t());
              const arma::mat solved=system.solve(panel);
              const arma::mat moving(change.rows(first,last).t());
              designTrace+=arma::accu(moving%solved);
            }
            score(p+k+coordinate)+=arma::dot(alpha,change*mode)-weight*designTrace;
          }
        }
      }
      double leverage = 0;
      for(arma::uword first = 0; first < n; first += panelSize){
        Rcpp::checkUserInterrupt();
        const arma::uword last = std::min(first+panelSize-1,n-1);
        const arma::mat panel(random.rows(first,last).t());
        const arma::mat solved = system.solve(panel);
        ++system.contractionSolves;
        leverage += arma::accu(panel % solved);
      }
      score(score.n_elem-1) = 0.5*(-static_cast<double>(n) +
        weight*(arma::dot(error,error)+leverage));
      if(!score.is_finite()) return invalid("non-finite Gaussian score");
      result["gradient"] = score;
    }
    result["work"] = Rcpp::List::create(Rcpp::_["evaluations"]=evaluations,
      Rcpp::_["factorizations"]=system.factorizations,
      Rcpp::_["contractionSolves"]=system.contractionSolves,
      Rcpp::_["descriptorUpdates"]=precision.descriptorUpdates);
    return result;
  }
};

} // namespace sommer
#endif
