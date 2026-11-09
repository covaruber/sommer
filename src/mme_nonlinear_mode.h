#ifndef SOMMER_MME_NONLINEAR_MODE_H
#define SOMMER_MME_NONLINEAR_MODE_H

#include "mme_glmm_precision.h"
#include "mme_nonlinear_derivatives.h"
#include "mme_glmm_family.h"
#include "mme_nonlinear_memory.h"

namespace sommer {

struct NonlinearModeValue {
  bool valid = false;
  std::string reason = "not evaluated";
  double objective = 0, logDet = 0, scoreNorm = 0;
  int iterations = 0;
  arma::vec mode, mean, residual, score, first, second, third;
  arma::mat fixedDerivative, mixedFirst, mixedSecond;
  arma::sp_mat jacobian, hessian;
  GLMMFamilyEvaluation likelihood;
  std::vector<CurveDerivatives> derivatives;
};

class NonlinearModeSystem {
  CholmodState factor;
public:
  int factorizations = 0, solves = 0, rejectedCurvatures = 0;
  bool factorize(const arma::sp_mat & hessian, double & logDet){
    const auto eigen = laplaceEigenSparse(hessian);
    ++factorizations;
    factor.ensureStarted();
    // Indefinite trial curvature is expected during safeguarded Newton steps.
    factor.common.print=0;
    factor.common.error_handler=nullptr;
    const bool success=factor.analyzePattern(eigen) && factor.factorize(eigen) &&
      laplaceFactorLogDeterminant(factor, logDet);
    if(factor.common.status!=CHOLMOD_OK && factor.common.status!=CHOLMOD_NOT_POSDEF)
      Rcpp::stop("Nonlinear mode CHOLMOD failure (status %d).",factor.common.status);
    if(!success) ++rejectedCurvatures;
    return success;
  }
  arma::mat solve(const arma::mat & rhs){
    ++solves;
    Eigen::Map<const Eigen::MatrixXd> mapped(rhs.memptr(), rhs.n_rows, rhs.n_cols);
    const Eigen::MatrixXd value = factor.solve(mapped, "Nonlinear mode solve failed.");
    return arma::mat(value.data(), value.rows(), value.cols());
  }
};

class NonlinearModeModel {
public:
  arma::sp_mat fixed, additive, incidence;
  arma::vec response, offset, input;
  std::string curve;
  arma::uword coordinate, curveStart;
  double memoryMB;
  GLMMFamilyModel * family = nullptr;
  bool joint = false;
  std::vector<arma::sp_mat> bases, originalBases;
  arma::sp_mat originalAdditive;

  NonlinearModeModel(const arma::sp_mat & x, const arma::sp_mat & z,
      const arma::vec & y, const arma::vec & o, const arma::vec & t,
      const std::string & function, arma::uword selected, arma::uword start,
      bool curveRandom=true, double memory=1000,
      const arma::mat & coordinateMap=arma::mat())
    : fixed(x), response(y), offset(o), input(t), curve(function),
      coordinate(selected), curveStart(start), memoryMB(memory) {
    if(!y.n_elem || x.n_rows != y.n_elem || z.n_rows != y.n_elem ||
       o.n_elem != y.n_elem || t.n_elem != y.n_elem || !z.n_cols ||
       start > z.n_cols || (curveRandom && start==z.n_cols) ||
       selected >= nonlinearParameterCount(curve) ||
       !x.is_finite() || !z.is_finite() || !y.is_finite() ||
       !o.is_finite() || !t.is_finite())
      Rcpp::stop("Nonlinear random inputs have invalid dimensions or values.");
    additive = z;
    incidence = arma::sp_mat(z.n_rows,z.n_cols);
    if(coordinateMap.n_elem){
      const arma::uword k=nonlinearParameterCount(curve);
      if(coordinateMap.n_rows!=z.n_cols || coordinateMap.n_cols!=k ||
         !coordinateMap.is_finite() || arma::any(arma::vectorise(coordinateMap)<0) ||
         arma::any(arma::vectorise(coordinateMap)>1) ||
         arma::any(arma::sum(coordinateMap,1)>1))
        Rcpp::stop("Invalid joint curve coordinate map.");
      joint=true;
      bases.assign(k,arma::sp_mat(z.n_rows,z.n_cols));
      for(auto cell=z.begin();cell!=z.end();++cell){
        for(arma::uword j=0;j<k;++j) if(coordinateMap(cell.col(),j)!=0){
          if(coordinateMap(cell.col(),j)!=1 || *cell!=1)
            Rcpp::stop("Joint random curve loadings must be unit incidence.");
          bases[j](cell.row(),cell.col())=*cell;
          incidence(cell.row(),cell.col())=*cell;
          additive(cell.row(),cell.col())=0;
        }
      }
      for(arma::uword j=0;j<k;++j) if(arma::accu(coordinateMap.col(j))>0){
        arma::vec counts(z.n_rows,arma::fill::zeros);
        for(auto cell=bases[j].begin();cell!=bases[j].end();++cell)
          counts(cell.row())+=*cell;
        if(arma::any(counts!=1))
          Rcpp::stop("Joint curves require one incidence per selected parameter.");
      }
      originalBases=bases;
      originalAdditive=additive;
      return;
    }
    arma::uvec counts(y.n_elem, arma::fill::zeros);
    for(auto cell=z.begin(); cell!=z.end(); ++cell){
      if(cell.col() < start) continue;
      if(*cell != 1) Rcpp::stop("Random curve loadings must be unit incidence.");
      incidence(cell.row(),cell.col()) = 1;
      ++counts(cell.row());
      additive(cell.row(),cell.col()) = 0;
    }
    if(curveRandom && arma::any(counts != 1))
      Rcpp::stop("Random curves require exactly one group coordinate per observation.");
  }

  void setDesign(const arma::sp_mat & design, arma::uword start){
    curveStart=start;
    additive=design;
    incidence=arma::sp_mat(design.n_rows,design.n_cols);
    for(auto cell=design.begin();cell!=design.end();++cell){
      if(cell.col()<start) continue;
      incidence(cell.row(),cell.col())=*cell;
      additive(cell.row(),cell.col())=0;
    }
  }

  NonlinearModeValue evaluate(const arma::vec & parameters, const arma::vec & mode,
      const arma::sp_mat & precision){
    NonlinearModeValue result;
    const arma::uword n=response.n_elem, p=fixed.n_cols;
    const arma::uword k=nonlinearParameterCount(curve);
    const double weight=family ? 1 : std::exp(-parameters(parameters.n_elem-1));
    if(!mode.is_finite() || !std::isfinite(weight) || weight <= 0){
      result.reason="non-finite mode or residual precision"; return result;
    }
    result.mode=mode;
    result.mean=offset+additive*mode;
    if(p) result.mean+=fixed*parameters.head(p);
    const arma::vec deviations=incidence*mode;
    result.first.set_size(n); result.second.set_size(n); result.third.set_size(n);
    result.fixedDerivative.set_size(n,k);
    result.mixedFirst.set_size(n,k); result.mixedSecond.set_size(n,k);
    arma::mat jointDeviations;
    if(joint){
      jointDeviations.set_size(n,k);
      for(arma::uword j=0;j<k;++j) jointDeviations.col(j)=bases[j]*mode;
      result.derivatives.resize(n);
    }
    for(arma::uword i=0;i<n;++i){
      if(i%1024==0) Rcpp::checkUserInterrupt();
      arma::vec theta=parameters.subvec(p,p+k-1);
      if(joint) theta+=jointDeviations.row(i).t();
      else theta(coordinate)+=deviations(i);
      CurveDerivatives row;
      if(!nonlinearCurveDerivatives(curve,input(i),theta,row)){
        result.reason="non-finite random curve evaluation"; return result;
      }
      result.mean(i)+=row.mean;
      if(joint) result.derivatives[i]=row;
      result.first(i)=row.d[coordinate];
      result.second(i)=row.d2[coordinate][coordinate];
      result.third(i)=row.d3[coordinate][coordinate][coordinate];
      for(arma::uword j=0;j<k;++j){
        result.fixedDerivative(i,j)=row.d[j];
        result.mixedFirst(i,j)=row.d2[coordinate][j];
        result.mixedSecond(i,j)=row.d3[coordinate][coordinate][j];
      }
    }
    result.residual=response-result.mean;
    result.jacobian=additive;
    if(joint){
      for(arma::uword j=0;j<k;++j)
        for(auto cell=bases[j].begin();cell!=bases[j].end();++cell)
          result.jacobian(cell.row(),cell.col())+=
            (*cell)*result.derivatives[cell.row()].d[j];
    }else{
      for(auto cell=incidence.begin();cell!=incidence.end();++cell)
        result.jacobian(cell.row(),cell.col())+=(*cell)*result.first(cell.row());
    }
    arma::sp_mat correction(n,n);
    if(family){
      result.likelihood=family->evaluate(response,result.mean);
      if(!std::isfinite(result.likelihood.logDensity)){
        result.reason="non-finite family likelihood"; return result;
      }
      result.score=result.jacobian.t()*result.likelihood.score-precision*mode;
      correction.diag()=result.likelihood.score%result.second;
      const arma::sp_mat curvature=family->observationCurvature(result.likelihood);
      nonlinearSparseAdmission(arma::abs(additive)+arma::abs(incidence),
        precision,curvature,memoryMB);
      result.hessian=precision+result.jacobian.t()*curvature*result.jacobian;
      if(joint) result.hessian-=secondCorrection(result,result.likelihood.score);
      else result.hessian-=incidence.t()*correction*incidence;
      result.objective=-result.likelihood.logDensity+0.5*arma::dot(mode,precision*mode);
    }else{
      result.score=weight*(result.jacobian.t()*result.residual)-precision*mode;
      arma::sp_mat curvature(n,n);
      curvature.diag().fill(weight);
      nonlinearSparseAdmission(arma::abs(additive)+arma::abs(incidence),
        precision,curvature,memoryMB);
      correction.diag()=weight*(result.residual % result.second);
      result.hessian=precision+result.jacobian.t()*curvature*result.jacobian;
      if(joint) result.hessian-=secondCorrection(result,weight*result.residual);
      else result.hessian-=incidence.t()*correction*incidence;
      result.objective=0.5*(weight*arma::dot(result.residual,result.residual)+
        arma::dot(mode,precision*mode));
    }
    result.scoreNorm=arma::norm(result.score,"inf");
    result.valid=std::isfinite(result.objective) && result.score.is_finite() &&
      result.hessian.is_finite();
    result.reason=result.valid ? "evaluated" : "non-finite conditional objective";
    return result;
  }

  arma::sp_mat secondCorrection(const NonlinearModeValue & value,
      const arma::vec & score) const {
    const arma::uword k=bases.size(), n=response.n_elem, q=additive.n_cols;
    arma::sp_mat result(q,q);
    for(arma::uword a=0;a<k;++a) for(arma::uword b=0;b<k;++b){
      if(!bases[a].n_nonzero || !bases[b].n_nonzero) continue;
      arma::sp_mat diagonal(n,n);
      for(arma::uword i=0;i<n;++i)
        diagonal(i,i)=score(i)*value.derivatives[i].d2[a][b];
      result+=bases[a].t()*diagonal*bases[b];
    }
    return result;
  }

  void jointPartial(const NonlinearModeValue & value,
      const arma::sp_mat & additiveChange,
      const std::vector<arma::sp_mat> & baseChanges, const arma::mat & direction,
      arma::vec & eta, arma::sp_mat & jacobianChange) const {
    eta=additiveChange*value.mode;
    jacobianChange=additiveChange;
    const arma::uword k=bases.size(), n=response.n_elem;
    for(arma::uword a=0;a<k;++a){
      arma::vec first(n),change(n,arma::fill::zeros);
      for(arma::uword i=0;i<n;++i){
        first(i)=value.derivatives[i].d[a];
        for(arma::uword b=0;b<k;++b)
          change(i)+=value.derivatives[i].d2[a][b]*direction(i,b);
      }
      eta+=first%direction.col(a);
      for(auto cell=baseChanges[a].begin();cell!=baseChanges[a].end();++cell)
        jacobianChange(cell.row(),cell.col())+=(*cell)*first(cell.row());
      for(auto cell=bases[a].begin();cell!=bases[a].end();++cell)
        jacobianChange(cell.row(),cell.col())+=(*cell)*change(cell.row());
    }
  }

  arma::sp_mat jointHessianDirection(const NonlinearModeValue & value,
      const arma::vec & score, const arma::vec & scoreDirection,
      const arma::mat & direction,
      const std::vector<arma::sp_mat> & baseChanges) const {
    const arma::uword k=bases.size(), n=response.n_elem;
    arma::sp_mat result(additive.n_cols,additive.n_cols);
    for(arma::uword a=0;a<k;++a) for(arma::uword b=0;b<k;++b){
      arma::sp_mat diagonal(n,n), change(n,n);
      for(arma::uword i=0;i<n;++i){
        const auto & row=value.derivatives[i];
        double third=0;
        for(arma::uword c=0;c<k;++c) third+=row.d3[a][b][c]*direction(i,c);
        diagonal(i,i)=score(i)*row.d2[a][b];
        change(i,i)=scoreDirection(i)*row.d2[a][b]+score(i)*third;
      }
      result+=baseChanges[a].t()*diagonal*bases[b]+
        bases[a].t()*diagonal*baseChanges[b]+bases[a].t()*change*bases[b];
    }
    return result;
  }

  NonlinearModeValue optimize(const arma::vec & parameters, const arma::vec & start,
      const arma::sp_mat & precision, NonlinearModeSystem & system,
      int maxIters, double scoreTol, double stepTol){
    auto value=evaluate(parameters,start,precision);
    for(int iteration=0;iteration<maxIters && value.valid;++iteration){
      Rcpp::checkUserInterrupt();
      value.iterations=iteration+1;
      if(value.scoreNorm <= scoreTol*(1+arma::norm(value.mode,"inf"))){
        if(!system.factorize(value.hessian,value.logDet)){
          value.valid=false; value.reason="non-positive-definite final mode Hessian";
        }else value.reason="converged";
        return value;
      }
      double logDet=0, damping=0;
      bool factored=false;
      const double scale=std::max(1.0,arma::abs(value.hessian.diag()).max());
      for(int attempt=0;attempt<14;++attempt){
        arma::sp_mat hessian=value.hessian;
        if(damping) hessian.diag()+=damping;
        if(system.factorize(hessian,logDet)){ factored=true; break; }
        damping=attempt==0 ? 1e-6*scale : 10*damping;
      }
      if(!factored){ value.valid=false; value.reason="mode damping failed"; return value; }
      const arma::vec step=system.solve(arma::mat(value.score));
      double fraction=1;
      bool accepted=false;
      while(fraction >= 1.0/65536){
        auto candidate=evaluate(parameters,value.mode+fraction*step,precision);
        const double roundoff=64*std::numeric_limits<double>::epsilon()*
          (1+std::abs(value.objective));
        if(candidate.valid && candidate.objective <= value.objective+roundoff &&
           (candidate.objective < value.objective-roundoff ||
            candidate.scoreNorm < value.scoreNorm)){
          candidate.iterations=iteration+1;
          value=std::move(candidate); accepted=true; break;
        }
        fraction*=0.5;
      }
      if(!accepted){ value.valid=false; value.reason="mode line search failed"; return value; }
      if(arma::norm(fraction*step,"inf") <= stepTol*(1+arma::norm(value.mode,"inf")) &&
         value.scoreNorm > scoreTol*(1+arma::norm(value.mode,"inf"))){
        value.valid=false; value.reason="mode stalled above score tolerance"; return value;
      }
    }
    if(value.valid){ value.valid=false; value.reason="mode iteration limit reached"; }
    return value;
  }
};

} // namespace sommer
#endif
