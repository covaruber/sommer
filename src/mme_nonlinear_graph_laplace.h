#ifndef SOMMER_MME_NONLINEAR_GRAPH_LAPLACE_H
#define SOMMER_MME_NONLINEAR_GRAPH_LAPLACE_H

#include "mme_nonlinear_mode.h"
#include "mme_nonlinear_covariance.h"
#include "mme_nonlinear_graph.h"

namespace sommer {

class NonlinearGraphLaplaceObjective {
  struct ModeValue {
    bool valid=false;
    std::string reason="not evaluated";
    double objective=0, scoreNorm=0, logDet=0;
    int iterations=0;
    arma::vec mode, mean, residual, score;
    arma::sp_mat jacobian, hessian;
  };

  arma::sp_mat fixed, originalRandom, originalAdditive;
  arma::vec response, offset;
  arma::mat graphData;
  std::unique_ptr<NonlinearGraph> graph;
  NonlinearCovariance covariance;
  NonlinearModeSystem system;
  std::vector<arma::sp_mat> originalBases, bases;
  arma::sp_mat additive;
  arma::uvec randomCoordinates;
  arma::uword p, k, c, q;
  double memoryMB, scoreTol, stepTol;
  int maxIters, evaluations=0;

  arma::mat parameterRows(const arma::vec & parameters,
      const arma::vec & mode) const {
    const arma::uword n=response.n_elem;
    arma::mat result(n,k);
    result.each_row()=parameters.subvec(p,p+k-1).t();
    for(arma::uword a=0;a<randomCoordinates.n_elem;++a)
      result.col(randomCoordinates(a))+=bases[randomCoordinates(a)]*mode;
    return result;
  }

  bool graphCoordinatesValid(const arma::mat & parameters) const {
    for(arma::uword row=0;row<parameters.n_rows;++row)
      if(!graph->validParameters(parameters.row(row).t())) return false;
    return true;
  }

  void graphDerivatives(const arma::mat & parameters, arma::vec & output,
      arma::mat & gradient, arma::cube & randomHessian) const {
    const arma::uword n=response.n_elem, r=randomCoordinates.n_elem;
    arma::mat directions(k,3*n,arma::fill::zeros);
    const Rcpp::List first=graph->evaluate(graphData,parameters,directions,memoryMB);
    output=Rcpp::as<arma::vec>(first["output"]);
    gradient=Rcpp::as<arma::mat>(first["gradient"]);
    randomHessian.zeros(n,r,r);
    for(arma::uword a=0;a<r;++a) for(arma::uword b=a;b<r;++b){
      directions.zeros();
      const arma::uword coordinateA=randomCoordinates(a);
      const arma::uword coordinateB=randomCoordinates(b);
      for(arma::uword row=0;row<n;++row){
        directions(coordinateA,3*row)=1;
        directions(coordinateB,3*row+1)=1;
      }
      const Rcpp::List value=graph->evaluate(graphData,parameters,directions,memoryMB);
      const Rcpp::NumericVector second=value["directionalSecond"];
      for(arma::uword row=0;row<n;++row){
        const double derivative=second[row+n*3];
        randomHessian(row,a,b)=derivative;
        randomHessian(row,b,a)=derivative;
      }
    }
    if(!output.is_finite() || !gradient.is_finite() ||
       !randomHessian.is_finite())
      Rcpp::stop("Non-finite graph value or conditional derivative.");
  }

  ModeValue conditional(const arma::vec & parameters,const arma::vec & mode){
    ModeValue result;
    result.mode=mode;
    const arma::uword n=response.n_elem;
    const double logVariance=parameters(parameters.n_elem-1);
    const double weight=std::exp(-logVariance);
    if(!std::isfinite(weight) || weight<=0 || !mode.is_finite()){
      result.reason="non-finite mode or residual precision";
      return result;
    }
    const arma::mat graphParameters=parameterRows(parameters,mode);
    if(!graphCoordinatesValid(graphParameters)){
      result.reason="random graph parameter outside transform range";
      return result;
    }
    arma::vec output;
    arma::mat gradient;
    arma::cube second;
    graphDerivatives(graphParameters,output,gradient,second);
    result.mean=offset+additive*mode+output;
    if(p) result.mean+=fixed*parameters.head(p);
    result.residual=response-result.mean;
    result.jacobian=additive;
    for(arma::uword a=0;a<randomCoordinates.n_elem;++a){
      const arma::uword coordinate=randomCoordinates(a);
      const auto & basis=bases[coordinate];
      for(auto cell=basis.begin();cell!=basis.end();++cell)
        result.jacobian(cell.row(),cell.col())+=
          (*cell)*gradient(cell.row(),coordinate);
    }
    const arma::vec scoreValue=weight*result.residual;
    const arma::sp_mat precision=covariance.prior().matrix;
    if(result.jacobian.n_rows!=n || result.jacobian.n_cols!=mode.n_elem ||
       precision.n_rows!=mode.n_elem || precision.n_cols!=mode.n_elem)
      Rcpp::stop("Graph conditional random design and precision dimensions differ.");
    const arma::vec observationScore=result.jacobian.t()*scoreValue;
    const arma::vec priorScore=precision*mode;
    if(observationScore.n_elem!=priorScore.n_elem)
      Rcpp::stop("Graph conditional observation and prior scores differ in dimension.");
    result.score=observationScore-priorScore;
    arma::sp_mat observationCurvature(n,n);
    observationCurvature.diag().fill(weight);
    arma::sp_mat designMagnitude(n,q);
    for(auto cell=additive.begin();cell!=additive.end();++cell)
      designMagnitude(cell.row(),cell.col())=std::abs(*cell);
    for(const auto & basis:bases)
      for(auto cell=basis.begin();cell!=basis.end();++cell)
        designMagnitude(cell.row(),cell.col())+=std::abs(*cell);
    nonlinearSparseAdmission(designMagnitude,precision,observationCurvature,
      memoryMB);
    result.hessian=precision+
      result.jacobian.t()*observationCurvature*result.jacobian;
    for(arma::uword a=0;a<randomCoordinates.n_elem;++a)
      for(arma::uword b=0;b<randomCoordinates.n_elem;++b){
        const arma::vec diagonal=weight*result.residual%second.slice(a).col(b);
        if(arma::any(diagonal!=0))
          result.hessian-=bases[randomCoordinates(a)].t()*
            arma::sp_mat(arma::diagmat(diagonal))*bases[randomCoordinates(b)];
      }
    result.objective=0.5*(weight*arma::dot(result.residual,result.residual)+
      arma::dot(mode,precision*mode));
    result.scoreNorm=arma::norm(result.score,"inf");
    result.valid=result.mean.is_finite() && result.score.is_finite() &&
      result.hessian.is_finite() && std::isfinite(result.objective);
    result.reason=result.valid ? "evaluated" : "non-finite conditional objective";
    return result;
  }

  ModeValue optimize(const arma::vec & parameters,const arma::vec & start){
    ModeValue value=conditional(parameters,start);
    for(int iteration=0;iteration<maxIters && value.valid;++iteration){
      Rcpp::checkUserInterrupt();
      value.iterations=iteration+1;
      if(value.scoreNorm<=scoreTol*(1+arma::norm(value.mode,"inf"))){
        if(!system.factorize(value.hessian,value.logDet)){
          value.valid=false;
          value.reason="non-positive-definite final mode Hessian";
        }
        return value;
      }
      double ignored=0;
      if(!system.factorize(value.hessian,ignored)){
        value.valid=false;
        value.reason="non-positive-definite trial mode Hessian";
        return value;
      }
      const arma::vec step=system.solve(arma::mat(value.score));
      if(!step.is_finite()){
        value.valid=false;
        value.reason="non-finite mode Newton step";
        return value;
      }
      bool accepted=false;
      double fraction=1;
      while(fraction>=1.0/65536){
        ModeValue candidate=conditional(parameters,value.mode+fraction*step);
        const double roundoff=64*std::numeric_limits<double>::epsilon()*
          (1+std::abs(value.objective));
        if(candidate.valid && candidate.objective<=value.objective+roundoff &&
           (candidate.objective<value.objective-roundoff ||
            candidate.scoreNorm<value.scoreNorm)){
          candidate.iterations=iteration+1;
          value=std::move(candidate);
          accepted=true;
          break;
        }
        fraction*=0.5;
      }
      if(!accepted){
        value.valid=false;
        value.reason="mode line search failed";
        return value;
      }
      if(arma::norm(fraction*step,"inf")<=stepTol*
           (1+arma::norm(value.mode,"inf")) &&
         value.scoreNorm>scoreTol*(1+arma::norm(value.mode,"inf"))){
        value.valid=false;
        value.reason="mode stalled above score tolerance";
        return value;
      }
    }
    if(value.valid){
      value.valid=false;
      value.reason="mode iteration limit reached";
    }
    return value;
  }

  Rcpp::List evaluateCore(const arma::vec & parameters){
    ++evaluations;
    auto invalid=[](const std::string & reason){
      return Rcpp::List::create(Rcpp::_["valid"]=false,Rcpp::_["reason"]=reason);
    };
    if(parameters.n_elem!=p+k+c+1 || !parameters.is_finite())
      return invalid("invalid or non-finite graph Laplace parameters");
    std::string reason;
    const arma::vec covarianceParameters=c ?
      arma::vec(parameters.subvec(p+k,p+k+c-1)) : arma::vec();
    if(!covariance.update(covarianceParameters,reason))
      return invalid(reason);
    additive=covariance.project(originalAdditive);
    for(arma::uword a=0;a<bases.size();++a)
      bases[a]=covariance.project(originalBases[a]);
    if(covariance.prior().matrix.n_rows!=q)
      Rcpp::stop("Graph random design and covariance dimensions differ.");
    ModeValue best;
    int failed=0,competing=0;
    std::vector<std::string> failures;
    for(int start=0;start<3;++start){
      arma::vec initial(q,arma::fill::zeros);
      if(start) initial.fill(start==1 ? 0.5 : -0.5);
      ModeValue candidate=optimize(parameters,initial);
      if(!candidate.valid){
        ++failed;
        failures.push_back(candidate.reason);
        continue;
      }
      if(best.valid && arma::norm(best.mode-candidate.mode,"inf")>1e-4)
        ++competing;
      if(!best.valid || candidate.objective<best.objective)
        best=std::move(candidate);
    }
    if(!best.valid){
      std::string message="all deterministic graph mode starts failed";
      for(const auto & failure:failures) message+="; "+failure;
      return invalid(message);
    }
    if(!system.factorize(best.hessian,best.logDet))
      return invalid("invalid accepted graph mode Hessian");
    const double logVariance=parameters(parameters.n_elem-1);
    const double likelihood=-best.objective-0.5*response.n_elem*
      (std::log(2*arma::datum::pi)+logVariance)+
      0.5*covariance.prior().logDet-0.5*best.logDet;
    if(!std::isfinite(likelihood)) return invalid("non-finite graph Laplace likelihood");
    arma::mat fixedParameters(response.n_elem,k);
    fixedParameters.each_row()=parameters.subvec(p,p+k-1).t();
    arma::mat noDirections(k,3*response.n_elem,arma::fill::zeros);
    const Rcpp::List fixedGraph=graph->evaluate(graphData,fixedParameters,
      noDirections,memoryMB);
    arma::vec fixedMean=offset+Rcpp::as<arma::vec>(fixedGraph["output"]);
    if(p) fixedMean+=fixed*parameters.head(p);
    Rcpp::List result=Rcpp::List::create(Rcpp::_["valid"]=true,
      Rcpp::_["reason"]="evaluated",Rcpp::_["logLik"]=likelihood,
      Rcpp::_["mode"]=covariance.restore(best.mode),Rcpp::_["latentMode"]=best.mode,
      Rcpp::_["fitted"]=best.mean,Rcpp::_["fixedMean"]=fixedMean,
      Rcpp::_["modeScore"]=best.scoreNorm,Rcpp::_["modeIterations"]=best.iterations,
      Rcpp::_["failedStarts"]=failed,Rcpp::_["failedStartReasons"]=failures,
      Rcpp::_["competingModes"]=competing,Rcpp::_["logDetHessian"]=best.logDet);
    return result;
  }

public:
  NonlinearGraphLaplaceObjective(const arma::sp_mat & x,
      const arma::sp_mat & z,const arma::vec & y,const arma::vec & o,
      const arma::mat & input,const Rcpp::List & specification,int output,
      const Rcpp::IntegerVector & transformed,const arma::mat & selected,
      const Rcpp::List & descriptors,const Rcpp::List & relationships,
      double memory,int iterations,double scoreTolerance,double stepTolerance)
    : fixed(x),originalRandom(z),response(y),offset(o),graphData(input),
      covariance(z,descriptors,relationships,memory),
      p(x.n_cols),k(transformed.size()),c(covariance.parameterCount()),
      q(covariance.latentStart(z.n_cols)),memoryMB(memory),
      scoreTol(scoreTolerance),stepTol(stepTolerance),maxIters(iterations) {
    if(output<0 || !y.n_elem || x.n_rows!=y.n_elem || z.n_rows!=y.n_elem ||
       o.n_elem!=y.n_elem || input.n_rows!=y.n_elem ||
       selected.n_rows!=z.n_cols || selected.n_cols!=k ||
       !x.is_finite() || !z.is_finite() || !y.is_finite() ||
       !o.is_finite() || !input.is_finite() ||
       !std::isfinite(memory) || memory<=0 || iterations<1 ||
       !std::isfinite(scoreTol) || scoreTol<=0 ||
       !std::isfinite(stepTol) || stepTol<=0)
      Rcpp::stop("Invalid nonlinear graph Laplace inputs or controls.");
    graph.reset(new NonlinearGraph(specification,static_cast<arma::uword>(output),
      input.n_cols,k,transformed));
    if(8.0*(100.0*y.n_elem+40.0*z.n_cols+8.0*z.n_nonzero+
       graph->workspaceBytes(y.n_elem,input.n_cols))>memory*1048576.0)
      Rcpp::stop("Nonlinear graph Laplace workspace exceeds execution$denseMemoryMB.");
    if(k==0 || !selected.is_finite() ||
       arma::any(arma::vectorise(selected)<0) ||
       arma::any(arma::vectorise(selected)>1) ||
       arma::any(arma::vectorise(selected)!=arma::floor(arma::vectorise(selected))))
      Rcpp::stop("Invalid nonlinear graph random-coordinate map.");
    std::vector<arma::uword> mapped;
    for(arma::uword j=0;j<k;++j){
      if(arma::accu(selected.col(j))>0) mapped.push_back(j);
    }
    if(mapped.empty())
      Rcpp::stop("Graph Laplace fitting requires at least one random graph coordinate.");
    randomCoordinates=arma::uvec(mapped);
    if(arma::accu(arma::sum(selected,1)>1)>0)
      Rcpp::stop("Random design columns cannot map to multiple graph coordinates.");
    originalAdditive=z;
    originalBases.assign(k,arma::sp_mat(y.n_elem,z.n_cols));
    for(auto cell=z.begin();cell!=z.end();++cell){
      for(arma::uword j=0;j<k;++j) if(selected(cell.col(),j)){
        if(*cell!=1)
          Rcpp::stop("Random graph parameter loadings must be unit incidence.");
        originalBases[j](cell.row(),cell.col())=*cell;
        originalAdditive(cell.row(),cell.col())=0;
      }
    }
    for(const auto j:randomCoordinates){
      arma::uvec counts(y.n_elem,arma::fill::zeros);
      for(auto cell=originalBases[j].begin();cell!=originalBases[j].end();++cell)
        ++counts(cell.row());
      if(arma::any(counts!=1))
        Rcpp::stop("Each random graph coordinate requires one group incidence per observation.");
    }
    bases=originalBases;
    const arma::uword r=randomCoordinates.n_elem;
    const double derivativeBytes=8.0*y.n_elem*r*r+
      24.0*y.n_elem*k;
    if(derivativeBytes+graph->workspaceBytes(y.n_elem,input.n_cols)>
       memory*1048576.0)
      Rcpp::stop("Nonlinear graph derivative workspace exceeds execution$denseMemoryMB.");
  }

  Rcpp::List evaluate(const arma::vec & parameters,bool gradient=true){
    Rcpp::List result=evaluateCore(parameters);
    if(!Rcpp::as<bool>(result["valid"]) || !gradient) return result;
    arma::vec derivative(parameters.n_elem);
    for(arma::uword j=0;j<parameters.n_elem;++j){
      Rcpp::checkUserInterrupt();
      const double step=1e-5*(1+std::abs(parameters(j)));
      arma::vec upper=parameters,lower=parameters;
      upper(j)+=step;
      lower(j)-=step;
      const Rcpp::List up=evaluateCore(upper);
      const Rcpp::List down=evaluateCore(lower);
      if(Rcpp::as<bool>(up["valid"]) && Rcpp::as<bool>(down["valid"])){
        derivative(j)=(Rcpp::as<double>(up["logLik"])-
          Rcpp::as<double>(down["logLik"]))/(2*step);
      }else if(Rcpp::as<bool>(up["valid"])){
        derivative(j)=(Rcpp::as<double>(up["logLik"])-
          Rcpp::as<double>(result["logLik"]))/step;
      }else if(Rcpp::as<bool>(down["valid"])){
        derivative(j)=(Rcpp::as<double>(result["logLik"])-
          Rcpp::as<double>(down["logLik"]))/step;
      }else{
        return Rcpp::List::create(Rcpp::_["valid"]=false,
          Rcpp::_["reason"]="numerical Laplace score could not find a valid perturbation");
      }
    }
    if(!derivative.is_finite())
      return Rcpp::List::create(Rcpp::_["valid"]=false,
        Rcpp::_["reason"]="non-finite numerical Laplace score");
    result["gradient"]=derivative;
    result["work"]=Rcpp::List::create(Rcpp::_["evaluations"]=evaluations,
      Rcpp::_["factorizations"]=system.factorizations,
      Rcpp::_["solves"]=system.solves,
      Rcpp::_["rejectedCurvatures"]=system.rejectedCurvatures,
      Rcpp::_["outerScore"]="central finite differences of exact Laplace objective");
    return result;
  }

  Rcpp::List conditionalAt(const arma::vec & parameters,const arma::vec & mode){
    if(parameters.n_elem!=p+k+c+1 || mode.n_elem!=q)
      Rcpp::stop("Graph Laplace parameter or conditional mode length is incompatible.");
    std::string reason;
    const arma::vec covarianceParameters=c ?
      arma::vec(parameters.subvec(p+k,p+k+c-1)) : arma::vec();
    if(!covariance.update(covarianceParameters,reason))
      return Rcpp::List::create(Rcpp::_["valid"]=false,Rcpp::_["reason"]=reason);
    additive=covariance.project(originalAdditive);
    for(arma::uword a=0;a<bases.size();++a)
      bases[a]=covariance.project(originalBases[a]);
    const ModeValue value=conditional(parameters,mode);
    return Rcpp::List::create(Rcpp::_["valid"]=value.valid,
      Rcpp::_["reason"]=value.reason,Rcpp::_["objective"]=value.valid ?
        value.objective : R_PosInf,Rcpp::_["score"]=value.score,
      Rcpp::_["hessian"]=value.hessian,Rcpp::_["mean"]=value.mean);
  }
};

class NonlinearGraphFamilyLaplaceObjective {
  struct ModeValue {
    bool valid=false;
    std::string reason="not evaluated";
    double objective=0, scoreNorm=0, logDet=0;
    int iterations=0;
    arma::vec mode, predictor, score;
    GLMMFamilyEvaluation likelihood;
    arma::sp_mat jacobian, hessian;
  };

  arma::sp_mat fixed, originalAdditive, additive;
  arma::vec response, offset;
  arma::mat graphData;
  std::unique_ptr<NonlinearGraph> graph;
  NonlinearCovariance covariance;
  GLMMFamilyModel family;
  NonlinearModeSystem system;
  std::vector<arma::sp_mat> originalBases, bases;
  arma::uvec randomCoordinates;
  arma::uword p, k, c, q;
  double memoryMB, scoreTol, stepTol;
  int maxIters, evaluations=0;

  arma::mat parameterRows(const arma::vec & parameters,
      const arma::vec & mode) const {
    arma::mat result(response.n_elem,k);
    result.each_row()=parameters.subvec(p,p+k-1).t();
    for(const auto coordinate:randomCoordinates)
      result.col(coordinate)+=bases[coordinate]*mode;
    return result;
  }

  void graphDerivatives(const arma::mat & parameters, arma::vec & output,
      arma::mat & gradient, arma::cube & randomHessian) const {
    const arma::uword n=response.n_elem, r=randomCoordinates.n_elem;
    arma::mat directions(k,3*n,arma::fill::zeros);
    const Rcpp::List value=graph->evaluate(graphData,parameters,directions,memoryMB);
    output=Rcpp::as<arma::vec>(value["output"]);
    gradient=Rcpp::as<arma::mat>(value["gradient"]);
    randomHessian.zeros(n,r,r);
    for(arma::uword a=0;a<r;++a) for(arma::uword b=a;b<r;++b){
      directions.zeros();
      const arma::uword coordinateA=randomCoordinates(a);
      const arma::uword coordinateB=randomCoordinates(b);
      for(arma::uword row=0;row<n;++row){
        directions(coordinateA,3*row)=1;
        directions(coordinateB,3*row+1)=1;
      }
      const Rcpp::List second=graph->evaluate(graphData,parameters,directions,memoryMB);
      const Rcpp::NumericVector directional=second["directionalSecond"];
      for(arma::uword row=0;row<n;++row){
        const double derivative=directional[row+n*3];
        randomHessian(row,a,b)=derivative;
        randomHessian(row,b,a)=derivative;
      }
    }
    if(!output.is_finite() || !gradient.is_finite() ||
       !randomHessian.is_finite())
      Rcpp::stop("Non-finite graph family value or conditional derivative.");
  }

  ModeValue conditional(const arma::vec & parameters,const arma::vec & mode){
    ModeValue result;
    result.mode=mode;
    if(parameters.n_elem!=p+k+c+family.parameterCount() ||
       !parameters.is_finite() || mode.n_elem!=q || !mode.is_finite()){
      result.reason="invalid or non-finite graph family parameters or mode";
      return result;
    }
    std::string reason;
    const arma::vec covarianceParameters=c ?
      arma::vec(parameters.subvec(p+k,p+k+c-1)) : arma::vec();
    if(!covariance.update(covarianceParameters,reason)){
      result.reason=reason;
      return result;
    }
    if(!family.update(parameters.tail(family.parameterCount()))){
      result.reason="invalid family parameters";
      return result;
    }
    additive=covariance.project(originalAdditive);
    for(arma::uword coordinate=0;coordinate<bases.size();++coordinate)
      bases[coordinate]=covariance.project(originalBases[coordinate]);
    const arma::sp_mat precision=covariance.prior().matrix;
    if(precision.n_rows!=q || precision.n_cols!=q){
      result.reason="graph family random design and precision dimensions differ";
      return result;
    }
    const arma::mat graphParameters=parameterRows(parameters,mode);
    for(arma::uword row=0;row<graphParameters.n_rows;++row)
      if(!graph->validParameters(graphParameters.row(row).t())){
        result.reason="random graph parameter outside transform range";
        return result;
      }
    arma::vec output;
    arma::mat gradient;
    arma::cube second;
    graphDerivatives(graphParameters,output,gradient,second);
    result.predictor=offset+additive*mode+output;
    if(p) result.predictor+=fixed*parameters.head(p);
    result.jacobian=additive;
    for(const auto coordinate:randomCoordinates){
      const arma::sp_mat & basis=bases[coordinate];
      for(auto cell=basis.begin();cell!=basis.end();++cell)
        result.jacobian(cell.row(),cell.col())+=
          (*cell)*gradient(cell.row(),coordinate);
    }
    result.likelihood=family.evaluate(response,result.predictor);
    if(!std::isfinite(result.likelihood.logDensity)){
      result.reason="non-finite family likelihood";
      return result;
    }
    const arma::sp_mat precisionMatrix=covariance.prior().matrix;
    const arma::sp_mat curvature=family.observationCurvature(result.likelihood);
    if(curvature.n_rows!=response.n_elem || curvature.n_cols!=response.n_elem ||
       result.jacobian.n_rows!=response.n_elem ||
       result.jacobian.n_cols!=q){
      result.reason="graph family observation curvature dimensions differ";
      return result;
    }
    arma::sp_mat designMagnitude=arma::abs(additive);
    for(const auto coordinate:randomCoordinates)
      designMagnitude+=arma::abs(bases[coordinate]);
    nonlinearSparseAdmission(designMagnitude,precisionMatrix,curvature,memoryMB);
    result.score=result.jacobian.t()*result.likelihood.score-precisionMatrix*mode;
    result.hessian=precisionMatrix+result.jacobian.t()*curvature*result.jacobian;
    for(arma::uword a=0;a<randomCoordinates.n_elem;++a)
      for(arma::uword b=0;b<randomCoordinates.n_elem;++b){
        const arma::vec diagonal=result.likelihood.score%second.slice(a).col(b);
        if(arma::any(diagonal!=0))
          result.hessian-=bases[randomCoordinates(a)].t()*
            arma::sp_mat(arma::diagmat(diagonal))*
            bases[randomCoordinates(b)];
      }
    result.objective=-result.likelihood.logDensity+
      0.5*arma::dot(mode,precisionMatrix*mode);
    result.scoreNorm=arma::norm(result.score,"inf");
    result.valid=result.predictor.is_finite() && result.score.is_finite() &&
      result.hessian.is_finite() && std::isfinite(result.objective);
    result.reason=result.valid ? "evaluated" : "non-finite conditional family objective";
    return result;
  }

  ModeValue optimize(const arma::vec & parameters,const arma::vec & start){
    ModeValue value=conditional(parameters,start);
    for(int iteration=0;iteration<maxIters && value.valid;++iteration){
      Rcpp::checkUserInterrupt();
      value.iterations=iteration+1;
      if(value.scoreNorm<=scoreTol*(1+arma::norm(value.mode,"inf"))){
        if(!system.factorize(value.hessian,value.logDet)){
          value.valid=false;
          value.reason="non-positive-definite final mode Hessian";
        }
        return value;
      }
      double ignored=0;
      if(!system.factorize(value.hessian,ignored)){
        value.valid=false;
        value.reason="non-positive-definite trial mode Hessian";
        return value;
      }
      const arma::vec step=system.solve(arma::mat(value.score));
      if(!step.is_finite()){
        value.valid=false;
        value.reason="non-finite mode Newton step";
        return value;
      }
      bool accepted=false;
      double fraction=1;
      while(fraction>=1.0/65536){
        ModeValue candidate=conditional(parameters,value.mode+fraction*step);
        const double roundoff=64*std::numeric_limits<double>::epsilon()*
          (1+std::abs(value.objective));
        if(candidate.valid && candidate.objective<=value.objective+roundoff &&
           (candidate.objective<value.objective-roundoff ||
            candidate.scoreNorm<value.scoreNorm)){
          candidate.iterations=iteration+1;
          value=std::move(candidate);
          accepted=true;
          break;
        }
        fraction*=0.5;
      }
      if(!accepted){
        value.valid=false;
        value.reason="mode line search failed";
        return value;
      }
      if(arma::norm(fraction*step,"inf")<=stepTol*
           (1+arma::norm(value.mode,"inf")) &&
         value.scoreNorm>scoreTol*(1+arma::norm(value.mode,"inf"))){
        value.valid=false;
        value.reason="mode stalled above score tolerance";
        return value;
      }
    }
    if(value.valid){
      value.valid=false;
      value.reason="mode iteration limit reached";
    }
    return value;
  }

  Rcpp::List evaluateCore(const arma::vec & parameters){
    ++evaluations;
    auto invalid=[](const std::string & reason){
      return Rcpp::List::create(Rcpp::_["valid"]=false,Rcpp::_["reason"]=reason);
    };
    if(parameters.n_elem!=p+k+c+family.parameterCount() ||
       !parameters.is_finite())
      return invalid("invalid or non-finite graph family parameters");
    std::string reason;
    const arma::vec covarianceParameters=c ?
      arma::vec(parameters.subvec(p+k,p+k+c-1)) : arma::vec();
    if(!covariance.update(covarianceParameters,reason)) return invalid(reason);
    if(!family.update(parameters.tail(family.parameterCount())))
      return invalid("invalid family parameters");
    additive=covariance.project(originalAdditive);
    for(arma::uword coordinate=0;coordinate<bases.size();++coordinate)
      bases[coordinate]=covariance.project(originalBases[coordinate]);
    ModeValue best;
    int failed=0, competing=0;
    std::vector<std::string> failures;
    for(int start=0;start<3;++start){
      arma::vec initial(q,arma::fill::zeros);
      if(start) initial.fill(start==1 ? 0.5 : -0.5);
      ModeValue candidate=optimize(parameters,initial);
      if(!candidate.valid){
        ++failed;
        failures.push_back(candidate.reason);
        continue;
      }
      if(best.valid && arma::norm(best.mode-candidate.mode,"inf")>1e-4)
        ++competing;
      if(!best.valid || candidate.objective<best.objective)
        best=std::move(candidate);
    }
    if(!best.valid){
      std::string message="all deterministic graph family mode starts failed";
      for(const auto & failure:failures) message+="; "+failure;
      return invalid(message);
    }
    if(!system.factorize(best.hessian,best.logDet))
      return invalid("invalid accepted graph family Hessian");
    const double likelihood=-best.objective+
      0.5*covariance.prior().logDet-0.5*best.logDet;
    if(!std::isfinite(likelihood)) return invalid("non-finite graph family Laplace likelihood");
    arma::mat fixedParameters(response.n_elem,k);
    fixedParameters.each_row()=parameters.subvec(p,p+k-1).t();
    const Rcpp::List fixedGraph=graph->evaluate(graphData,fixedParameters,
      arma::mat(k,3*response.n_elem,arma::fill::zeros),memoryMB);
    arma::vec fixedPredictor=offset+Rcpp::as<arma::vec>(fixedGraph["output"]);
    if(p) fixedPredictor+=fixed*parameters.head(p);
    Rcpp::List result=Rcpp::List::create(Rcpp::_["valid"]=true,
      Rcpp::_["reason"]="evaluated",Rcpp::_["logLik"]=likelihood,
      Rcpp::_["mode"]=covariance.restore(best.mode),Rcpp::_["latentMode"]=best.mode,
      Rcpp::_["predictor"]=best.predictor,Rcpp::_["fitted"]=best.likelihood.mean,
      Rcpp::_["fixedMean"]=fixedPredictor,Rcpp::_["modeScore"]=best.scoreNorm,
      Rcpp::_["modeIterations"]=best.iterations,Rcpp::_["failedStarts"]=failed,
      Rcpp::_["failedStartReasons"]=failures,Rcpp::_["competingModes"]=competing,
      Rcpp::_["logDetHessian"]=best.logDet);
    return result;
  }

public:
  NonlinearGraphFamilyLaplaceObjective(const arma::sp_mat & x,
      const arma::sp_mat & z,const arma::vec & y,const arma::vec & o,
      const arma::mat & input,const Rcpp::List & specification,int output,
      const Rcpp::IntegerVector & transformed,const arma::mat & selected,
      const Rcpp::List & descriptors,const Rcpp::List & relationships,
      const Rcpp::List & familySpecifications,const arma::uvec & familyIndex,
      const arma::vec & familyTrials,double memory,int iterations,
      double scoreTolerance,double stepTolerance)
    : fixed(x),originalAdditive(z),response(y),offset(o),
      graphData(input),covariance(z,descriptors,relationships,memory),
      family("",y.n_elem,familySpecifications,familyIndex,familyTrials),
      p(x.n_cols),k(transformed.size()),c(covariance.parameterCount()),
      q(covariance.latentStart(z.n_cols)),memoryMB(memory),
      scoreTol(scoreTolerance),stepTol(stepTolerance),maxIters(iterations) {
    if(output<0 || !y.n_elem || x.n_rows!=y.n_elem || z.n_rows!=y.n_elem ||
       o.n_elem!=y.n_elem || input.n_rows!=y.n_elem ||
       selected.n_rows!=z.n_cols || selected.n_cols!=k ||
       !x.is_finite() || !z.is_finite() || !y.is_finite() ||
       !o.is_finite() || !input.is_finite() || !std::isfinite(memory) ||
       memory<=0 || iterations<1 || !std::isfinite(scoreTol) ||
       scoreTol<=0 || !std::isfinite(stepTol) || stepTol<=0)
      Rcpp::stop("Invalid nonlinear graph family inputs or controls.");
    graph.reset(new NonlinearGraph(specification,static_cast<arma::uword>(output),
      input.n_cols,k,transformed));
    if(8.0*(120.0*y.n_elem+50.0*z.n_cols+8.0*z.n_nonzero+
       graph->workspaceBytes(y.n_elem,input.n_cols))>memory*1048576.0)
      Rcpp::stop("Nonlinear graph family workspace exceeds execution$denseMemoryMB.");
    if(k==0 || !selected.is_finite() ||
       arma::any(arma::vectorise(selected)<0) ||
       arma::any(arma::vectorise(selected)>1) ||
       arma::any(arma::vectorise(selected)!=arma::floor(arma::vectorise(selected))))
      Rcpp::stop("Invalid nonlinear graph family random-coordinate map.");
    std::vector<arma::uword> mapped;
    for(arma::uword j=0;j<k;++j)
      if(arma::accu(selected.col(j))>0) mapped.push_back(j);
    randomCoordinates=arma::uvec(mapped);
    if(arma::accu(arma::sum(selected,1)>1)>0)
      Rcpp::stop("Random design columns cannot map to multiple graph coordinates.");
    originalBases.assign(k,arma::sp_mat(y.n_elem,z.n_cols));
    for(auto cell=z.begin();cell!=z.end();++cell){
      for(arma::uword j=0;j<k;++j) if(selected(cell.col(),j)){
        if(*cell!=1)
          Rcpp::stop("Random graph parameter loadings must be unit incidence.");
        originalBases[j](cell.row(),cell.col())=*cell;
        originalAdditive(cell.row(),cell.col())=0;
      }
    }
    for(const auto coordinate:randomCoordinates){
      arma::uvec counts(y.n_elem,arma::fill::zeros);
      for(auto cell=originalBases[coordinate].begin();
          cell!=originalBases[coordinate].end();++cell)
        ++counts(cell.row());
      if(arma::any(counts!=1))
        Rcpp::stop("Each random graph coordinate requires one group incidence per observation.");
    }
    additive=originalAdditive;
    bases=originalBases;
    family.validateResponse(response);
    const arma::uword r=randomCoordinates.n_elem;
    const double bytes=8.0*y.n_elem*(24.0*k+static_cast<double>(r)*r)+
      graph->workspaceBytes(y.n_elem,input.n_cols);
    if(bytes>memory*1048576.0)
      Rcpp::stop("Nonlinear graph family derivative workspace exceeds execution$denseMemoryMB.");
  }

  Rcpp::List evaluate(const arma::vec & parameters,bool gradient=true){
    const Rcpp::List result=evaluateCore(parameters);
    if(!Rcpp::as<bool>(result["valid"]) || !gradient) return result;
    arma::vec derivative(parameters.n_elem);
    for(arma::uword j=0;j<parameters.n_elem;++j){
      Rcpp::checkUserInterrupt();
      const double step=1e-5*(1+std::abs(parameters(j)));
      arma::vec upper=parameters,lower=parameters;
      upper(j)+=step;
      lower(j)-=step;
      const Rcpp::List up=evaluateCore(upper);
      const Rcpp::List down=evaluateCore(lower);
      if(!Rcpp::as<bool>(up["valid"]) || !Rcpp::as<bool>(down["valid"]))
        return Rcpp::List::create(Rcpp::_["valid"]=false,
          Rcpp::_["reason"]="central finite-difference Laplace score could not evaluate both perturbations");
      derivative(j)=(Rcpp::as<double>(up["logLik"])-
        Rcpp::as<double>(down["logLik"]))/(2*step);
    }
    if(!derivative.is_finite())
      return Rcpp::List::create(Rcpp::_["valid"]=false,
        Rcpp::_["reason"]="non-finite central finite-difference Laplace score");
    Rcpp::List output=result;
    output["gradient"]=derivative;
    output["work"]=Rcpp::List::create(Rcpp::_["evaluations"]=evaluations,
      Rcpp::_["factorizations"]=system.factorizations,
      Rcpp::_["solves"]=system.solves,
      Rcpp::_["rejectedCurvatures"]=system.rejectedCurvatures,
      Rcpp::_["outerScore"]="central finite differences of exact Laplace objective");
    return output;
  }

  Rcpp::List conditionalAt(const arma::vec & parameters,const arma::vec & mode){
    if(parameters.n_elem!=p+k+c+family.parameterCount() || mode.n_elem!=q)
      Rcpp::stop("Graph family parameters or conditional mode length is incompatible.");
    std::string reason;
    const arma::vec covarianceParameters=c ?
      arma::vec(parameters.subvec(p+k,p+k+c-1)) : arma::vec();
    if(!covariance.update(covarianceParameters,reason))
      return Rcpp::List::create(Rcpp::_["valid"]=false,Rcpp::_["reason"]=reason);
    if(!family.update(parameters.tail(family.parameterCount())))
      return Rcpp::List::create(Rcpp::_["valid"]=false,
        Rcpp::_["reason"]="invalid family parameters");
    additive=covariance.project(originalAdditive);
    for(arma::uword coordinate=0;coordinate<bases.size();++coordinate)
      bases[coordinate]=covariance.project(originalBases[coordinate]);
    const ModeValue value=conditional(parameters,mode);
    return Rcpp::List::create(Rcpp::_["valid"]=value.valid,
      Rcpp::_["reason"]=value.reason,Rcpp::_["objective"]=value.valid ?
        value.objective : R_PosInf,Rcpp::_["score"]=value.score,
      Rcpp::_["hessian"]=value.hessian,Rcpp::_["predictor"]=value.predictor,
      Rcpp::_["mean"]=value.predictor);
  }
};

} // namespace sommer
#endif
