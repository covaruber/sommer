#ifndef SOMMER_MME_NONLINEAR_FAMILY_H
#define SOMMER_MME_NONLINEAR_FAMILY_H

#include "mme_nonlinear_mode.h"
#include "mme_nonlinear_covariance.h"

namespace sommer {

class NonlinearFamilyObjective {
  NonlinearModeModel model;
  NonlinearCovariance covariance;
  GLMMFamilyModel family;
  NonlinearModeSystem system;
  arma::uword fullStart, q, p, k, c;
  int maxIters, evaluations=0;
  double scoreTol, stepTol, memoryMB;

  double trace(const arma::sp_mat & change){
    if(!change.n_nonzero) return 0;
    const double available=memoryMB*1048576-
      8.0*(150.0*model.response.n_elem+60.0*q+8.0*change.n_nonzero);
    const double capacity=std::floor(available/(24.0*q));
    if(capacity<1) Rcpp::stop("Nonlinear family trace workspace exceeds execution$denseMemoryMB.");
    const arma::uword panel=static_cast<arma::uword>(std::min(32.0,capacity));
    double result=0;
    for(arma::uword first=0;first<q;first+=panel){
      Rcpp::checkUserInterrupt();
      const arma::uword last=std::min(q-1,first+panel-1);
      const arma::mat solved=system.solve(arma::mat(change.cols(first,last)));
      for(arma::uword j=first;j<=last;++j) result+=solved(j,j-first);
    }
    return result;
  }

  arma::sp_mat curvatureDirection(const NonlinearModeValue & value,
      const arma::vec & direction){
    arma::sp_mat result=arma::sp_mat(arma::diagmat(value.likelihood.curvatureD1%direction));
    for(const auto & specification:family.specifications){
      const auto & groups=specification.nominalGroups;
      for(arma::uword i=0;i<groups.n_rows;++i){
        const arma::uvec rows=groups.row(i).t();
        const auto block=multinomialEvaluation(value.mean.elem(rows),
          static_cast<arma::uword>(model.response(rows(0)))-1,direction.elem(rows),1);
        for(arma::uword j=0;j<rows.n_elem;++j)
          for(arma::uword l=0;l<rows.n_elem;++l)
            result(rows(j),rows(l))=block.curvatureD1(j,l);
      }
    }
    return result;
  }

  double score(const NonlinearModeValue & value, const arma::vec & etaPartial,
      const arma::sp_mat & jPartial, const arma::sp_mat & bPartial,
      const arma::vec & secondPartial, const arma::sp_mat & qPartial,
      double priorEnvelope, double likelihoodEnvelope,
      const arma::vec & scorePartial, const arma::sp_mat & curvaturePartial){
    const auto & j=value.jacobian;
    const auto & b=model.incidence;
    const auto & s=value.likelihood.score;
    const arma::sp_mat curvature=family.observationCurvature(value.likelihood);
    const arma::vec rhs=jPartial.t()*s+j.t()*(scorePartial-curvature*etaPartial)-
      qPartial*value.mode;
    const arma::vec motion=system.solve(arma::mat(rhs));
    const arma::vec bd=b*motion;
    const arma::vec eta=etaPartial+j*motion;
    arma::sp_mat jt=jPartial;
    for(auto cell=b.begin();cell!=b.end();++cell)
      jt(cell.row(),cell.col())+=(*cell)*value.second(cell.row())*bd(cell.row());
    const arma::vec ds=scorePartial-curvature*eta;
    const arma::sp_mat dc=curvaturePartial+curvatureDirection(value,eta);
    const arma::sp_mat d=arma::sp_mat(arma::diagmat(s%value.second));
    const arma::vec dd=ds%value.second+s%(secondPartial+value.third%bd);
    const arma::sp_mat change=qPartial+jt.t()*curvature*j+j.t()*curvature*jt+
      j.t()*dc*j-bPartial.t()*d*b-b.t()*d*bPartial-
      b.t()*arma::sp_mat(arma::diagmat(dd))*b;
    return priorEnvelope+likelihoodEnvelope+arma::dot(s,etaPartial)-0.5*trace(change);
  }

  double jointScore(const NonlinearModeValue & value, arma::vec etaPartial,
      const arma::sp_mat & jPartial, arma::mat coordinatePartial,
      const std::vector<arma::sp_mat> & baseChanges,
      const arma::sp_mat & qPartial, double priorEnvelope,
      double likelihoodEnvelope, const arma::vec & scorePartial,
      const arma::sp_mat & curvaturePartial){
    const auto & j=value.jacobian;
    const auto & s=value.likelihood.score;
    const arma::sp_mat curvature=family.observationCurvature(value.likelihood);
    const arma::vec motion=system.solve(arma::mat(jPartial.t()*s+
      j.t()*(scorePartial-curvature*etaPartial)-qPartial*value.mode));
    arma::mat motionCoordinates(model.response.n_elem,k);
    for(arma::uword a=0;a<k;++a) motionCoordinates.col(a)=model.bases[a]*motion;
    arma::vec unused;
    arma::sp_mat jMotion;
    const arma::sp_mat empty(model.response.n_elem,q);
    const std::vector<arma::sp_mat> emptyBases(k,empty);
    model.jointPartial(value,empty,emptyBases,motionCoordinates,unused,jMotion);
    const arma::sp_mat jt=jPartial+jMotion;
    const arma::vec eta=etaPartial+j*motion;
    const arma::vec ds=scorePartial-curvature*eta;
    const arma::sp_mat dc=curvaturePartial+curvatureDirection(value,eta);
    const arma::sp_mat change=qPartial+jt.t()*curvature*j+j.t()*curvature*jt+
      j.t()*dc*j-model.jointHessianDirection(value,s,ds,
        coordinatePartial+motionCoordinates,baseChanges);
    return priorEnvelope+likelihoodEnvelope+arma::dot(s,etaPartial)-0.5*trace(change);
  }

public:
  NonlinearFamilyObjective(const arma::sp_mat & x, const arma::sp_mat & z,
      const arma::vec & y, const arma::vec & offset, const arma::vec & input,
      const std::string & curve, int coordinate, arma::uword start,
      const Rcpp::List & descriptors, const Rcpp::List & relationships,
      const Rcpp::List & specifications, const arma::uvec & index,
      const arma::vec & trials, double memory, int iterations,
      double scoreTolerance, double stepTolerance,
      const arma::mat & coordinateMap=arma::mat())
    : model(x,z,y,offset,input,curve,coordinate>0 ? coordinate-1 : 0,start,
        coordinate>0 || coordinateMap.n_elem,memory,coordinateMap),
      covariance(z,descriptors,relationships,memory),
      family("",y.n_elem,specifications,index,trials),
      fullStart(start),q(covariance.latentStart(z.n_cols)),p(x.n_cols),
      k(nonlinearParameterCount(curve)),c(covariance.parameterCount()),
      maxIters(iterations),scoreTol(scoreTolerance),stepTol(stepTolerance),memoryMB(memory){
    if(!std::isfinite(memory) || memory<=0 || iterations<1 ||
       !std::isfinite(scoreTol) || scoreTol<=0 || !std::isfinite(stepTol) || stepTol<=0)
      Rcpp::stop("Invalid nonlinear family resource or mode controls.");
    if(8.0*(150.0*y.n_elem+60.0*q+8.0*z.n_nonzero)>memory*1048576)
      Rcpp::stop("Nonlinear family workspace exceeds execution$denseMemoryMB.");
    family.validateResponse(y);
    model.family=&family;
  }

  NonlinearModeValue conditional(const arma::vec & parameters,const arma::vec & mode){
    if(parameters.n_elem!=p+k+c+family.parameterCount() || mode.n_elem!=q)
      Rcpp::stop("Nonlinear family parameter or latent mode length is incompatible.");
    NonlinearModeValue invalid;
    if(!parameters.is_finite()){ invalid.reason="non-finite outer parameters"; return invalid; }
    std::string reason;
    const arma::vec cp=c ? arma::vec(parameters.subvec(p+k,p+k+c-1)) : arma::vec();
    if(!covariance.update(cp,reason)){ invalid.reason=reason; return invalid; }
    if(!family.update(parameters.tail(family.parameterCount()))){
      invalid.reason="invalid family parameters"; return invalid;
    }
    if(model.joint){
      model.curveStart=covariance.latentStart(fullStart);
      model.additive=covariance.project(model.originalAdditive);
      model.incidence=arma::sp_mat(model.response.n_elem,q);
      for(arma::uword a=0;a<k;++a){
        model.bases[a]=covariance.project(model.originalBases[a]);
        model.incidence+=arma::abs(model.bases[a]);
      }
    }else model.setDesign(covariance.design,covariance.latentStart(fullStart));
    return model.evaluate(parameters,mode,covariance.prior().matrix);
  }

  Rcpp::List evaluate(const arma::vec & parameters,bool gradient){
    ++evaluations;
    auto invalid=[](const std::string & reason){
      return Rcpp::List::create(Rcpp::_["valid"]=false,Rcpp::_["reason"]=reason);
    };
    const auto initial=conditional(parameters,arma::vec(q,arma::fill::zeros));
    if(!initial.valid) return invalid(initial.reason);
    NonlinearModeValue best;
    int failed=0,competing=0;
    std::vector<std::string> failures;
    for(int start=0;start<3;++start){
      arma::vec mode(q,arma::fill::zeros);
      if(start){
        if(model.curveStart<q) mode.tail(q-model.curveStart).fill(start==1 ? 0.5 : -0.5);
        else mode.fill(start==1 ? 0.5 : -0.5);
      }
      auto candidate=model.optimize(parameters,mode,covariance.prior().matrix,system,
        maxIters,scoreTol,stepTol);
      if(!candidate.valid){ ++failed; failures.push_back(candidate.reason); continue; }
      if(best.valid && arma::norm(best.mode-candidate.mode,"inf")>1e-4) ++competing;
      if(!best.valid || candidate.objective<best.objective) best=std::move(candidate);
    }
    if(!best.valid){
      std::string reason="all deterministic family mode starts failed";
      for(const auto & failure:failures) reason+="; "+failure;
      return invalid(reason);
    }
    if(!system.factorize(best.hessian,best.logDet)) return invalid("invalid accepted family Hessian");
    const double likelihood=-best.objective+0.5*(covariance.prior().logDet-best.logDet);
    Rcpp::List result=Rcpp::List::create(Rcpp::_["valid"]=true,Rcpp::_["reason"]="evaluated",
      Rcpp::_["logLik"]=likelihood,Rcpp::_["mode"]=covariance.restore(best.mode),
      Rcpp::_["latentMode"]=best.mode,Rcpp::_["predictor"]=best.mean,
      Rcpp::_["fitted"]=best.likelihood.mean,Rcpp::_["modeScore"]=best.scoreNorm,
      Rcpp::_["modeIterations"]=best.iterations,Rcpp::_["failedStarts"]=failed,
      Rcpp::_["failedStartReasons"]=failures,Rcpp::_["competingModes"]=competing,
      Rcpp::_["logDetHessian"]=best.logDet);
    if(gradient){
      const arma::uword n=model.response.n_elem;
      const arma::vec zero(n,arma::fill::zeros);
      const arma::sp_mat empty(n,q), emptyQ(q,q), emptyC(n,n);
      arma::vec gradientValue(parameters.n_elem,arma::fill::zeros);
      const std::vector<arma::sp_mat> emptyBases(k,empty);
      for(arma::uword t=0;t<p+k;++t){
        arma::vec eta=zero,second=zero;
        arma::sp_mat jt=empty;
        if(model.joint){
          arma::mat direction(n,k,arma::fill::zeros);
          if(t<p) eta=arma::vec(model.fixed.col(t));
          else direction.col(t-p).ones();
          arma::vec curveEta;
          model.jointPartial(best,empty,emptyBases,direction,curveEta,jt);
          eta+=curveEta;
          gradientValue(t)=jointScore(best,eta,jt,direction,emptyBases,emptyQ,
            0,0,zero,emptyC);
          continue;
        }
        if(t<p) eta=arma::vec(model.fixed.col(t));
        else{
          eta=best.fixedDerivative.col(t-p);
          second=best.mixedSecond.col(t-p);
          for(auto cell=model.incidence.begin();cell!=model.incidence.end();++cell)
            jt(cell.row(),cell.col())=(*cell)*best.mixedFirst(cell.row(),t-p);
        }
        gradientValue(t)=score(best,eta,jt,empty,second,emptyQ,0,0,zero,emptyC);
      }
      const auto derivatives=covariance.derivatives();
      for(arma::uword t=0;t<c;++t){
        const auto & derivative=derivatives[t];
        arma::sp_mat qt(q,q);
        const arma::uword end=derivative.start+derivative.block.n_rows-1;
        qt.submat(derivative.start,derivative.start,end,end)=derivative.block;
        const double envelope=0.5*(derivative.logDet-arma::dot(best.mode,qt*best.mode));
        if(model.joint){
          const arma::sp_mat a=covariance.project(model.originalAdditive,t);
          std::vector<arma::sp_mat> changes;
          arma::mat direction(n,k);
          for(arma::uword b=0;b<k;++b){
            changes.push_back(covariance.project(model.originalBases[b],t));
            direction.col(b)=changes.back()*best.mode;
          }
          arma::vec eta;
          arma::sp_mat jt;
          model.jointPartial(best,a,changes,direction,eta,jt);
          gradientValue(p+k+t)=jointScore(best,eta,jt,direction,changes,qt,
            envelope,0,zero,emptyC);
          continue;
        }
        const arma::sp_mat moving=covariance.designDerivative(t);
        arma::sp_mat a=moving,b(n,q);
        for(auto cell=moving.begin();cell!=moving.end();++cell){
          if(cell.col()<model.curveStart) continue;
          b(cell.row(),cell.col())=*cell;
          a(cell.row(),cell.col())=0;
        }
        const arma::vec bd=b*best.mode;
        const arma::vec eta=a*best.mode+best.first%bd;
        arma::sp_mat jt=a;
        for(auto cell=b.begin();cell!=b.end();++cell)
          jt(cell.row(),cell.col())+=(*cell)*best.first(cell.row());
        for(auto cell=model.incidence.begin();cell!=model.incidence.end();++cell)
          jt(cell.row(),cell.col())+=(*cell)*best.second(cell.row())*bd(cell.row());
        gradientValue(p+k+t)=score(best,eta,jt,b,best.third%bd,qt,envelope,0,zero,emptyC);
      }
      arma::uword output=p+k+c;
      for(arma::uword f=0;f<family.specifications.size();++f){
        const auto & specification=family.specifications[f];
        if(!specification.free) continue;
        const arma::uvec rows=arma::find(family.index==f);
        if(specification.family=="ordinal"){
          const arma::uword count=specification.thresholds.n_elem;
          std::vector<arma::vec> ls(count,zero),sp(count,zero),cp(count,zero);
          for(const auto row:rows){
            const arma::uword category=static_cast<arma::uword>(model.response(row))-1;
            const auto density=glmmOrdinalLogDensity(model.response(row),best.mean(row),
              specification.thresholds,specification.link);
            for(int endpoint=0;endpoint<2;++endpoint){
              if((endpoint==0 && category==0) || (endpoint==1 && category==count)) continue;
              const arma::uword cut=endpoint==0 ? category-1 : category;
              const auto d=glmmOrdinalCutpointScore(specification.thresholds(cut),
                best.mean(row),specification.link,density);
              const double sign=endpoint==0 ? -1 : 1;
              ls[cut](row)=sign*d.coefficient[0][0];
              sp[cut](row)=sign*d.coefficient[1][0];
              cp[cut](row)=-2*sign*d.coefficient[2][0];
            }
          }
          for(arma::uword t=0;t<count;++t){
            arma::vec l=zero,s=zero,h=zero;
            if(t==0){
              for(const auto row:rows){
                l(row)=-best.likelihood.score(row);
                s(row)=best.likelihood.curvature(row);
                h(row)=-best.likelihood.curvatureD1(row);
              }
            }else{
              for(arma::uword cut=t;cut<count;++cut){ l+=ls[cut]; s+=sp[cut]; h+=cp[cut]; }
            }
            const double scale=t ? specification.thresholds(t)-specification.thresholds(t-1) : 1;
            if(model.joint)
              gradientValue(output++)=jointScore(best,zero,empty,
                arma::mat(n,k,arma::fill::zeros),emptyBases,emptyQ,0,
                scale*arma::accu(l),scale*s,arma::sp_mat(arma::diagmat(scale*h)));
            else gradientValue(output++)=score(best,zero,empty,empty,zero,emptyQ,0,
                scale*arma::accu(l),scale*s,arma::sp_mat(arma::diagmat(scale*h)));
          }
        }else{
          const arma::uword count=specification.zeroDesign.n_cols ?
            specification.zeroDesign.n_cols : 1;
          for(arma::uword t=0;t<count;++t){
            arma::vec l=zero,s=zero,h=zero;
            for(const auto row:rows){
              const double scale=specification.zeroDesign.n_cols ? specification.zeroDesign(row,t) : 1;
              l(row)=scale*best.likelihood.parameterScore(row);
              s(row)=scale*best.likelihood.scoreParameter(row);
              h(row)=scale*best.likelihood.curvatureParameter(row);
            }
            if(model.joint)
              gradientValue(output++)=jointScore(best,zero,empty,
                arma::mat(n,k,arma::fill::zeros),emptyBases,emptyQ,0,
                arma::accu(l),s,arma::sp_mat(arma::diagmat(h)));
            else gradientValue(output++)=score(best,zero,empty,empty,zero,emptyQ,0,
                arma::accu(l),s,arma::sp_mat(arma::diagmat(h)));
          }
        }
      }
      if(!gradientValue.is_finite()) return invalid("non-finite complete nonlinear family score");
      result["gradient"]=gradientValue;
    }
    result["work"]=Rcpp::List::create(Rcpp::_["evaluations"]=evaluations,
      Rcpp::_["factorizations"]=system.factorizations,Rcpp::_["solves"]=system.solves,
      Rcpp::_["rejectedCurvatures"]=system.rejectedCurvatures);
    return result;
  }
};

} // namespace sommer
#endif
