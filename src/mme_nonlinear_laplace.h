#ifndef SOMMER_MME_NONLINEAR_LAPLACE_H
#define SOMMER_MME_NONLINEAR_LAPLACE_H

#include "mme_nonlinear_mode.h"
#include "mme_nonlinear_covariance.h"

namespace sommer {

class NonlinearLaplaceObjective {
  NonlinearModeModel model;
  NonlinearCovariance covarianceModel;
  LaplacePrecision & precision;
  NonlinearModeSystem system;
  double memoryMB, scoreTol, stepTol;
  int maxIters, evaluations=0;
  arma::uword q;
  arma::uword fullCurveStart;

  static const arma::sp_mat & admit(const arma::sp_mat & x,
      const arma::sp_mat & z, arma::uword n, double memory){
    if(!std::isfinite(memory) || memory<=0)
      Rcpp::stop("Invalid nonlinear mode memory budget.");
    if(8.0*(100.0*n+40.0*z.n_cols+8.0*z.n_nonzero+
        2.0*x.n_nonzero)>memory*1048576.0)
      Rcpp::stop("Nonlinear Laplace mode workspace exceeds execution$denseMemoryMB.");
    return x;
  }

  double traceDerivative(const LaplacePrecisionDerivative & derivative,
      arma::uword panelSize){
    if(!derivative.block.n_nonzero) return 0;
    double trace=0;
    const arma::uword end=derivative.start+derivative.block.n_rows-1;
    for(arma::uword first=0;first<derivative.block.n_cols;first+=panelSize){
      Rcpp::checkUserInterrupt();
      const arma::uword last=std::min(first+panelSize-1,derivative.block.n_cols-1);
      arma::mat rhs(q,last-first+1,arma::fill::zeros);
      rhs.rows(derivative.start,end)=arma::mat(derivative.block.cols(first,last));
      const arma::mat solved=system.solve(rhs);
      for(arma::uword column=first;column<=last;++column)
        trace+=solved(derivative.start+column,column-first);
    }
    return trace;
  }

public:
  NonlinearLaplaceObjective(const arma::sp_mat & x, const arma::sp_mat & z,
      const arma::vec & y, const arma::vec & o, const arma::vec & input,
      const std::string & curve, arma::uword coordinate, arma::uword start,
      const Rcpp::List & descriptors, const Rcpp::List & relationships,
      double memory, int iterations, double scoreTolerance, double stepTolerance)
    : model(admit(x,z,y.n_elem,memory),z,y,o,input,curve,coordinate,start,true,memory),
      covarianceModel(z,descriptors,relationships,memory),
      precision(covarianceModel.prior()),
      memoryMB(memory), scoreTol(scoreTolerance), stepTol(stepTolerance),
      maxIters(iterations), q(covarianceModel.latentStart(z.n_cols)),
      fullCurveStart(start) {
    if(!std::isfinite(memory) || memory<=0 || iterations<1 ||
       !std::isfinite(scoreTol) || scoreTol<=0 ||
       !std::isfinite(stepTol) || stepTol<=0)
      Rcpp::stop("Invalid nonlinear mode controls or memory budget.");
  }

  NonlinearModeValue conditional(const arma::vec & parameters, const arma::vec & mode){
    const arma::uword fixed=model.fixed.n_cols+nonlinearParameterCount(model.curve);
    const arma::uword c=covarianceModel.parameterCount();
    if(parameters.n_elem!=fixed+c+1 || mode.n_elem!=q)
      Rcpp::stop("Nonlinear Laplace parameter or mode length is incompatible.");
    arma::vec covariance(c);
    if(c) covariance=parameters.subvec(fixed,fixed+c-1);
    std::string reason;
    if(!covarianceModel.update(covariance,reason)){
      NonlinearModeValue result; result.reason=reason; return result;
    }
    if(covarianceModel.hasReducedRank())
      model.setDesign(covarianceModel.design,covarianceModel.latentStart(fullCurveStart));
    if(precision.matrix.n_rows!=q)
      Rcpp::stop("Nonlinear random covariance and mode dimensions differ.");
    return model.evaluate(parameters,mode,precision.matrix);
  }

  Rcpp::List evaluate(const arma::vec & parameters, bool gradient){
    ++evaluations;
    const arma::uword p=model.fixed.n_cols, k=nonlinearParameterCount(model.curve);
    const arma::uword c=covarianceModel.parameterCount(), n=model.response.n_elem;
    if(parameters.n_elem!=p+k+c+1)
      Rcpp::stop("Nonlinear Laplace parameter vector has incompatible length.");
    auto invalid=[](const std::string & reason){
      return Rcpp::List::create(Rcpp::_["valid"]=false,Rcpp::_["reason"]=reason);
    };
    if(!parameters.is_finite()) return invalid("non-finite outer parameters");
    const double weight=std::exp(-parameters(parameters.n_elem-1));
    if(!std::isfinite(weight) || weight<=0)
      return invalid("residual variance outside numerical range");
    arma::vec zero(q,arma::fill::zeros);
    const auto initial=conditional(parameters,zero);
    if(!initial.valid) return invalid(initial.reason);
    NonlinearModeValue best;
    int failed=0, competing=0;
    std::vector<std::string> reasons;
    for(int start=0;start<3;++start){
      arma::vec mode(q,arma::fill::zeros);
      if(start) mode.tail(q-model.curveStart).fill(start==1 ? 0.5 : -0.5);
      auto candidate=model.optimize(parameters,mode,precision.matrix,system,
        maxIters,scoreTol,stepTol);
      if(!candidate.valid){ ++failed; reasons.push_back(candidate.reason); continue; }
      if(best.valid && arma::norm(best.mode-candidate.mode,"inf")>1e-4) ++competing;
      if(!best.valid || candidate.objective<best.objective) best=std::move(candidate);
    }
    if(!best.valid){
      std::string reason="all deterministic mode starts failed";
      for(const auto & failure:reasons) reason+="; "+failure;
      return invalid(reason);
    }
    if(!system.factorize(best.hessian,best.logDet))
      return invalid("accepted mode Hessian factorization failed");
    const double likelihood=-best.objective-0.5*n*
      (std::log(2*arma::datum::pi)+parameters(parameters.n_elem-1))+
      0.5*precision.logDet-0.5*best.logDet;
    if(!std::isfinite(likelihood)) return invalid("non-finite Laplace likelihood");
    Rcpp::List result=Rcpp::List::create(Rcpp::_["valid"]=true,
      Rcpp::_["reason"]="evaluated",Rcpp::_["logLik"]=likelihood,
      Rcpp::_["mode"]=covarianceModel.restore(best.mode),Rcpp::_["fitted"]=best.mean,
      Rcpp::_["latentMode"]=best.mode,
      Rcpp::_["modeScore"]=best.scoreNorm,Rcpp::_["modeIterations"]=best.iterations,
      Rcpp::_["failedStarts"]=failed,Rcpp::_["competingModes"]=competing,
      Rcpp::_["failedStartReasons"]=reasons,Rcpp::_["logDetHessian"]=best.logDet);
    if(gradient){
      const double available=memoryMB*1048576.0-
        8.0*(100.0*n+40.0*q+8.0*model.incidence.n_nonzero+
          8.0*model.additive.n_nonzero);
      const double columns=std::floor(available/(8.0*(6.0*q+4.0*n)));
      if(columns<1) Rcpp::stop("Nonlinear Laplace score panels exceed execution$denseMemoryMB.");
      const arma::uword panelSize=static_cast<arma::uword>(std::min(32.0,columns));
      arma::vec jj(n,arma::fill::zeros), bb(n,arma::fill::zeros), bj(n,arma::fill::zeros);
      for(arma::uword first=0;first<n;first+=panelSize){
        Rcpp::checkUserInterrupt();
        const arma::uword last=std::min(first+panelSize-1,n-1);
        const arma::mat j(best.jacobian.rows(first,last).t());
        const arma::mat b(model.incidence.rows(first,last).t());
        const arma::mat solvedJ=system.solve(j), solvedB=system.solve(b);
        jj.subvec(first,last)=arma::sum(j%solvedJ,0).t();
        bb.subvec(first,last)=arma::sum(b%solvedB,0).t();
        bj.subvec(first,last)=arma::sum(b%solvedJ,0).t();
      }
      // Total dH includes movement of the conditional mode.
      auto modeTrace=[&](const arma::vec & movement){
        const arma::vec b=model.incidence*movement, j=best.jacobian*movement;
        return weight*arma::accu(2*best.second%b%bj+
          (best.second%j-best.residual%best.third%b)%bb);
      };
      arma::vec score(parameters.n_elem,arma::fill::zeros);
      for(arma::uword coordinate=0;coordinate<p+k;++coordinate){
        arma::vec derivative(n), mixed(n,arma::fill::zeros), mixedSecond(n,arma::fill::zeros);
        if(coordinate<p) derivative=arma::vec(model.fixed.col(coordinate));
        else{
          derivative=best.fixedDerivative.col(coordinate-p);
          mixed=best.mixedFirst.col(coordinate-p);
          mixedSecond=best.mixedSecond.col(coordinate-p);
        }
        const arma::vec partialScore=weight*(model.incidence.t()*
          (mixed%best.residual)-best.jacobian.t()*derivative);
        const arma::vec movement=system.solve(arma::mat(partialScore));
        const double trace=weight*arma::accu(2*mixed%bj+
          (derivative%best.second-best.residual%mixedSecond)%bb)+modeTrace(movement);
        score(coordinate)=weight*arma::dot(best.residual,derivative)-0.5*trace;
      }
      const auto derivatives=covarianceModel.derivatives();
      for(arma::uword coordinate=0;coordinate<derivatives.size();++coordinate){
        const auto & derivative=derivatives[coordinate];
        const arma::uword end=derivative.start+derivative.block.n_rows-1;
        const arma::vec local=best.mode.subvec(derivative.start,end);
        const arma::vec applied=derivative.block*local;
        arma::vec rhs(q,arma::fill::zeros);
        rhs.subvec(derivative.start,end)=-applied;
        double partialTrace=traceDerivative(derivative,panelSize);
        double envelope=0.5*(derivative.logDet-arma::dot(local,applied));
        if(covarianceModel.hasReducedRank()){
          const arma::sp_mat change=covarianceModel.designDerivative(coordinate);
          arma::sp_mat a=change;
          arma::sp_mat b(a.n_rows,a.n_cols);
          for(auto cell=change.begin();cell!=change.end();++cell){
            if(cell.col()<model.curveStart) continue;
            b(cell.row(),cell.col())=*cell;
            a(cell.row(),cell.col())=0;
          }
          if(a.n_nonzero || b.n_nonzero){
            const arma::vec deviation=b*best.mode;
            const arma::vec meanChange=a*best.mode+best.first%deviation;
            arma::sp_mat jacobianChange=a;
            for(auto cell=b.begin();cell!=b.end();++cell)
              jacobianChange(cell.row(),cell.col())+=(*cell)*best.first(cell.row());
            for(auto cell=model.incidence.begin();cell!=model.incidence.end();++cell)
              jacobianChange(cell.row(),cell.col())+=
                (*cell)*best.second(cell.row())*deviation(cell.row());
            rhs+=weight*(jacobianChange.t()*best.residual-best.jacobian.t()*meanChange);
            envelope+=weight*arma::dot(best.residual,meanChange);
            double contraction=0;
            for(arma::uword first=0;first<n;first+=panelSize){
              Rcpp::checkUserInterrupt();
              const arma::uword last=std::min(first+panelSize-1,n-1);
              const arma::mat j(best.jacobian.rows(first,last).t());
              const arma::mat moving(jacobianChange.rows(first,last).t());
              contraction+=2*arma::accu(moving%system.solve(j));
              const arma::mat base(model.incidence.rows(first,last).t());
              const arma::mat movingB(b.rows(first,last).t());
              const arma::vec cross=arma::sum(movingB%system.solve(base),0).t();
              contraction-=2*arma::dot(best.residual.subvec(first,last)%
                best.second.subvec(first,last),cross);
            }
            partialTrace+=weight*(contraction+arma::accu(
              (meanChange%best.second-best.residual%best.third%deviation)%bb));
          }
        }
        const arma::vec movement=system.solve(arma::mat(rhs));
        score(p+k+coordinate)=envelope-0.5*(partialTrace+modeTrace(movement));
      }
      const arma::vec movement=system.solve(arma::mat(-weight*
        (best.jacobian.t()*best.residual)));
      const double trace=-weight*arma::accu(jj-best.residual%best.second%bb)+
        modeTrace(movement);
      score(score.n_elem-1)=-0.5*n+0.5*weight*arma::dot(best.residual,best.residual)-
        0.5*trace;
      if(!score.is_finite()) return invalid("non-finite analytic Laplace score");
      result["gradient"]=score;
    }
    result["work"]=Rcpp::List::create(Rcpp::_["evaluations"]=evaluations,
      Rcpp::_["factorizations"]=system.factorizations,Rcpp::_["solves"]=system.solves,
      Rcpp::_["rejectedCurvatures"]=system.rejectedCurvatures,
      Rcpp::_["descriptorUpdates"]=precision.descriptorUpdates);
    return result;
  }
};

} // namespace sommer
#endif
