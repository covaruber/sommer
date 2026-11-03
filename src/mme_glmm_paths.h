#ifndef SOMMER_MME_GLMM_PATHS_H
#define SOMMER_MME_GLMM_PATHS_H

namespace sommer {

enum class GLMMEstimator { PQL, Laplace };
enum class GLMMPath { None, PQL, ScalarLaplace, CategoryBlockLaplace };
enum class GLMMPathReason {
  Eligible, CategoricalPQL, MixturePQL, LaplaceRequiresML,
  ResidualDependence
};

struct GLMMPathRequest {
  GLMMEstimator estimator;
  bool maximumLikelihood;
  bool ordinal;
  bool nominal;
  bool mixture;
  bool independentResiduals;
};

struct GLMMPathPlan {
  GLMMPath requested;
  bool eligible;
  GLMMPath selected;
  GLMMPathReason reason;
};

inline GLMMPathPlan planGLMMPath(const GLMMPathRequest & request){
  const GLMMPath requested = request.estimator == GLMMEstimator::PQL ?
    GLMMPath::PQL : (request.nominal ? GLMMPath::CategoryBlockLaplace :
      GLMMPath::ScalarLaplace);
  GLMMPathReason reason = GLMMPathReason::Eligible;
  if(request.estimator == GLMMEstimator::PQL){
    if(request.ordinal || request.nominal) reason = GLMMPathReason::CategoricalPQL;
    else if(request.mixture) reason = GLMMPathReason::MixturePQL;
  }else if(!request.maximumLikelihood){
    reason = GLMMPathReason::LaplaceRequiresML;
  }else if(!request.independentResiduals){
    reason = GLMMPathReason::ResidualDependence;
  }
  const bool eligible = reason == GLMMPathReason::Eligible;
  return {requested, eligible, eligible ? requested : GLMMPath::None, reason};
}

inline const char * glmmPathName(GLMMPath path){
  switch(path){
    case GLMMPath::PQL: return "pql";
    case GLMMPath::ScalarLaplace: return "native-laplace-scalar";
    case GLMMPath::CategoryBlockLaplace: return "native-laplace-category-block";
    default: return "none";
  }
}

inline const char * glmmPathReasonCode(GLMMPathReason reason){
  switch(reason){
    case GLMMPathReason::CategoricalPQL: return "categorical-pql-unsupported";
    case GLMMPathReason::MixturePQL: return "mixture-pql-unsupported";
    case GLMMPathReason::LaplaceRequiresML: return "laplace-requires-ml";
    case GLMMPathReason::ResidualDependence: return "laplace-residual-dependence";
    default: return "eligible";
  }
}

inline const char * glmmPathReasonMessage(GLMMPathReason reason){
  switch(reason){
    case GLMMPathReason::CategoricalPQL:
      return "ordinal() and multinomial() are not supported by PQL; categorical fitting is supported only by native Laplace ML.";
    case GLMMPathReason::MixturePQL:
      return "zipoisson() is currently supported only by native Laplace fitting.";
    case GLMMPathReason::LaplaceRequiresML:
      return "Laplace GLMMs use maximum likelihood; set task$objective='ml'.";
    case GLMMPathReason::ResidualDependence:
      return "Native Laplace requires independent observation likelihoods; structured rcov terms are not yet supported.";
    default: return "Requested GLMM likelihood path is eligible; no estimator fallback.";
  }
}

} // namespace sommer
#endif
