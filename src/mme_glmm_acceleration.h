#ifndef SOMMER_MME_GLMM_ACCELERATION_H
#define SOMMER_MME_GLMM_ACCELERATION_H

#include "mme_glmm_precision.h"

namespace sommer {

inline void laplaceChoice(const std::string & value,
    const std::vector<std::string> & choices, const std::string & label){
  if(std::find(choices.begin(),choices.end(),value)==choices.end())
    Rcpp::stop("Invalid Laplace "+label+" choice.");
}

struct LaplaceWarmStart {
  std::string requested;
  bool enabled;
  bool verifyCold;
  arma::vec mode;
  int attempts = 0;
  int coldRetries = 0;
  int coldChecks = 0;
  int iterations = 0;

  LaplaceWarmStart(const std::string & choice, bool concave)
    : requested(choice),enabled(choice=="on" || (choice=="auto" && concave)),
      verifyCold(choice=="on" && !concave) {
    laplaceChoice(choice,{"auto","off","on"},"warmStart");
  }
  arma::vec start(arma::uword dimension){
    if(enabled && mode.n_elem==dimension){ ++attempts; return mode; }
    return arma::vec(dimension,arma::fill::zeros);
  }
  void remember(const Rcpp::List & value){
    if(enabled && Rcpp::as<bool>(value["valid"])) mode=Rcpp::as<arma::vec>(value["mode"]);
  }
  void count(const Rcpp::List & value){
    if(value.containsElementNamed("iterations"))
      iterations+=Rcpp::as<int>(value["iterations"]);
  }
  Rcpp::List diagnostics() const {
    return Rcpp::List::create(Rcpp::_["requested"]=requested,
      Rcpp::_["enabled"]=enabled,Rcpp::_["coldVerification"]=verifyCold,
      Rcpp::_["attempts"]=attempts,Rcpp::_["coldRetries"]=coldRetries,
      Rcpp::_["coldChecks"]=coldChecks,Rcpp::_["modeIterations"]=iterations);
  }
};

class LaplaceSelectedInverse {
  std::vector<std::vector<int>> rows;
  std::vector<std::vector<double>> values;
  std::vector<int> permutation;
public:
  static double workspace(const CholmodState & state){
    return 48.0*(state.factor->is_super ? state.factor->xsize : state.factor->nzmax)+
      80.0*state.factor->n;
  }
  bool getPermuted(int a,int b,double & value) const {
    const int column=std::min(a,b),row=std::max(a,b);
    const auto & indices=rows[column];
    const auto found=std::lower_bound(indices.begin(),indices.end(),row);
    if(found==indices.end() || *found!=row) return false;
    value=values[column][found-indices.begin()];
    return true;
  }
  bool get(int a,int b,double & value) const {
    return getPermuted(permutation[a],permutation[b],value);
  }
  void build(CholmodState & state,double memoryMB){
    if(!state.factor || !state.factor->is_ll || state.factor->itype!=CHOLMOD_INT ||
       state.factor->xtype!=CHOLMOD_REAL || state.factor->dtype!=CHOLMOD_DOUBLE)
      Rcpp::stop("Laplace selected inverse requires a real, double LL' factor with integer indices.");
    if(workspace(state)>memoryMB*1048576.0)
      Rcpp::stop("Laplace selected-inverse workspace exceeds execution$denseMemoryMB.");
    auto releaseFactor=[&state](cholmod_factor * pointer){
      M_cholmod_free_factor(&pointer,&state.common);
    };
    auto releaseSparse=[&state](cholmod_sparse * pointer){
      M_cholmod_free_sparse(&pointer,&state.common);
    };
    std::unique_ptr<cholmod_factor,decltype(releaseFactor)> copy(
      M_cholmod_copy_factor(state.factor,&state.common),releaseFactor);
    if(!copy) Rcpp::stop("Laplace selected inverse could not copy the mode factor.");
    std::unique_ptr<cholmod_sparse,decltype(releaseSparse)> lower(
      M_cholmod_factor_to_sparse(copy.get(),&state.common),releaseSparse);
    copy.reset();
    if(!lower) Rcpp::stop("Laplace selected inverse could not extract the mode factor.");
    const int n=static_cast<int>(lower->ncol);
    const int * p=static_cast<const int *>(lower->p);
    const int * i=static_cast<const int *>(lower->i);
    const double * x=static_cast<const double *>(lower->x);
    const double bytes=48.0*p[n]+80.0*n;
    if(bytes>memoryMB*1048576.0){
      Rcpp::stop("Laplace selected-inverse workspace exceeds execution$denseMemoryMB.");
    }
    rows.assign(n,{});
    values.assign(n,{});
    permutation.assign(n,-1);
    const int * perm=static_cast<const int *>(state.factor->Perm);
    for(int column=0; column<n; ++column){
      permutation[perm ? perm[column] : column]=column;
      for(int k=p[column]; k<p[column+1]; ++k) rows[column].push_back(i[k]);
      values[column].assign(rows[column].size(),0.0);
    }
    // LL' Takahashi recurrence in the filled Cholesky pattern.
    for(int column=n-1; column>=0; --column){
      if(column%128==0) Rcpp::checkUserInterrupt();
      const double diagonal=x[p[column]];
      if(rows[column].empty() || rows[column][0]!=column ||
         !std::isfinite(diagonal) || diagonal<=0){
        Rcpp::stop("Laplace selected inverse encountered an invalid Cholesky pivot.");
      }
      for(std::size_t j=1; j<rows[column].size(); ++j){
        double sum=0;
        for(std::size_t k=1; k<rows[column].size(); ++k){
          double value;
          if(!getPermuted(rows[column][k],rows[column][j],value)){
            Rcpp::stop("Laplace Takahashi recurrence requested an unavailable entry.");
          }
          sum+=(x[p[column]+k]/diagonal)*value;
        }
        values[column][j]=-sum;
      }
      double correction=0;
      for(std::size_t k=1; k<rows[column].size(); ++k)
        correction+=(x[p[column]+k]/diagonal)*values[column][k];
      values[column][0]=1.0/(diagonal*diagonal)-correction;
    }
  }
};

} // namespace sommer
#endif
