#ifndef SOMMER_MME_NONLINEAR_COVARIANCE_H
#define SOMMER_MME_NONLINEAR_COVARIANCE_H

#include "mme_glmm_precision.h"
#include <memory>

namespace sommer {

// Keep zero-specific loading axes in the design and a proper density in latent space.
class NonlinearCovariance {
  Rcpp::List descriptors, latentDescriptors;
  arma::sp_mat original;
  std::unique_ptr<LaplacePrecision> precision;
  std::vector<arma::vec> parameters;
  std::vector<arma::uword> fullStarts, latentStarts, levels;
  std::vector<int> freeTerms, freePositions, precisionPositions;
  std::vector<arma::sp_mat> maps;
  bool reduced = false;
  double memoryMB;

  arma::sp_mat coordinateMap(int term, int derivative=-1) const {
    const Rcpp::List descriptor=descriptors[term], factors=descriptor["factors"];
    arma::sp_mat product(1,1); product(0,0)=1;
    for(int axis=0;axis<factors.size();++axis){
      const Rcpp::List factor=factors[axis];
      const arma::uword dimension=Rcpp::as<int>(factor["dim"]);
      arma::sp_mat value;
      if(zeroSpecificFactor(factor)){
        const int first=Rcpp::as<int>(factor["par_start"])-1;
        const int last=Rcpp::as<int>(factor["par_end"])-1;
        const arma::mat loading=normalizedFactorLoadings(factor,
          parameters[term].subvec(first,last));
        if(derivative>=first && derivative<=last){
          const int local=derivative-first;
          const Rcpp::IntegerVector rows=factor["rr_row"], columns=factor["rr_col"];
          const Rcpp::LogicalVector diagonal=factor["rr_diag"];
          arma::mat change(loading.n_rows,loading.n_cols,arma::fill::zeros);
          const double reference=parameters[term](first);
          change(rows[local]-1,columns[local]-1)=diagonal[local]
            ? loading(rows[local]-1,columns[local]-1) : std::exp(-reference);
          if(rows[local]==1 && columns[local]==1) change-=loading;
          value=arma::sp_mat(change);
        }else value=arma::sp_mat(loading);
      }else value=arma::speye<arma::sp_mat>(dimension,dimension);
      product=arma::kron(product,value);
    }
    return product;
  }

public:
  arma::sp_mat design;

  NonlinearCovariance(const arma::sp_mat & random, const Rcpp::List & supplied,
      const Rcpp::List & relationships, double memory)
    : descriptors(Rcpp::clone(supplied)), latentDescriptors(Rcpp::clone(supplied)),
      original(random), memoryMB(memory) {
    if(supplied.size()!=relationships.size() || supplied.size()==0)
      Rcpp::stop("Nonlinear covariance descriptors and relationships must match.");
    fullStarts.push_back(0); latentStarts.push_back(0);
    int precisionIndex=0;
    for(int term=0;term<descriptors.size();++term){
      const Rcpp::List descriptor=descriptors[term];
      Rcpp::List latent=latentDescriptors[term];
      Rcpp::List factors=Rcpp::clone(Rcpp::as<Rcpp::List>(descriptor["factors"]));
      Rcpp::LogicalVector free=Rcpp::clone(Rcpp::as<Rcpp::LogicalVector>(descriptor["free"]));
      const arma::vec originalPar=Rcpp::as<arma::vec>(descriptor["par"]);
      std::vector<double> latentPar(1,originalPar(0));
      std::vector<int> latentFree(1,free[0]);
      std::vector<bool> loading(free.size(),false);
      arma::uword dimension=1;
      for(int axis=0;axis<factors.size();++axis){
        Rcpp::List factor=factors[axis];
        const int first=Rcpp::as<int>(factor["par_start"])-1;
        const int last=Rcpp::as<int>(factor["par_end"])-1;
        factor["par_start"]=static_cast<int>(latentPar.size()+1);
        if(zeroSpecificFactor(factor)){
          reduced=true;
          const int rank=Rcpp::as<int>(factor["order"]);
          for(int j=first;j<=last;++j) loading[j]=true;
          factor["dim"]=rank;
          factor["zero_specific"]=false;
          factor["representation"]=Rcpp::List::create(Rcpp::_["kind"]="diagonal");
          factor["evaluator"]=Rcpp::List::create(Rcpp::_["backend"]="fixed",
            Rcpp::_["matrix"]=arma::eye<arma::mat>(rank,rank));
          factor["derivative"]=Rcpp::List::create(Rcpp::_["backend"]="none");
          factor["report"]=Rcpp::List::create(Rcpp::_["backend"]="builtin",
            Rcpp::_["transform"]=Rcpp::CharacterVector());
          factor["structurally_diagonal"]=true;
        }else{
          for(int j=first;j<=last;++j){
            latentPar.push_back(originalPar(j));
            latentFree.push_back(free[j]);
          }
        }
        factor["par_end"]=static_cast<int>(latentPar.size());
        factors[axis]=factor;
        dimension*=Rcpp::as<int>(factor["dim"]);
      }
      for(int j=0;j<free.size();++j){
        if(free[j]){
          freeTerms.push_back(term); freePositions.push_back(j);
          precisionPositions.push_back(loading[j] ? -1 : precisionIndex++);
        }
        if(loading[j]) free[j]=false;
      }
      latent["factors"]=factors;
      latent["par"]=arma::vec(latentPar);
      Rcpp::LogicalVector compactFree(latentFree.size());
      for(std::size_t j=0;j<latentFree.size();++j) compactFree[j]=latentFree[j];
      latent["free"]=compactFree;
      latent["dim"]=static_cast<int>(dimension);
      latentDescriptors[term]=latent;
      const arma::sp_mat relationship=Rcpp::as<arma::sp_mat>(relationships[term]);
      levels.push_back(relationship.n_rows);
      fullStarts.push_back(fullStarts.back()+
        Rcpp::as<int>(descriptor["dim"])*relationship.n_rows);
      latentStarts.push_back(latentStarts.back()+dimension*relationship.n_rows);
      parameters.push_back(Rcpp::as<arma::vec>(descriptor["par"]));
      maps.emplace_back();
    }
    if(fullStarts.back()!=original.n_cols)
      Rcpp::stop("Nonlinear covariance coordinate widths differ from the random design.");
    precision.reset(new LaplacePrecision(latentDescriptors,relationships));
    if(!reduced) design=original;
  }

  arma::uword parameterCount() const { return freeTerms.size(); }
  bool hasReducedRank() const { return reduced; }
  LaplacePrecision & prior() { return *precision; }

  arma::uword latentStart(arma::uword fullStart) const {
    for(std::size_t term=0;term<fullStarts.size();++term)
      if(fullStarts[term]==fullStart) return latentStarts[term];
    Rcpp::stop("Random curve must begin at a covariance term boundary.");
  }

  bool update(const arma::vec & free, std::string & reason){
    if(free.n_elem!=parameterCount())
      Rcpp::stop("Nonlinear covariance parameter vector has incompatible length.");
    arma::vec latent(precision->parameterCount());
    for(arma::uword j=0;j<free.n_elem;++j){
      parameters[freeTerms[j]](freePositions[j])=free(j);
      if(precisionPositions[j]>=0) latent(precisionPositions[j])=free(j);
    }
    if(!free.is_finite()){ reason="non-finite covariance coordinates"; return false; }
    if(!precision->update(latent,reason)) return false;
    if(!reduced) return true;
    design=arma::sp_mat(original.n_rows,latentStarts.back());
    for(int term=0;term<descriptors.size();++term){
      const Rcpp::List descriptor=descriptors[term];
      const double dimension=Rcpp::as<int>(descriptor["dim"]);
      const double width=latentStarts[term+1]-latentStarts[term];
      if(8*(4*dimension*width+
          4*original.n_nonzero+4*design.n_nonzero)>memoryMB*1048576){
        Rcpp::stop("Nonlinear reduced-rank map workspace exceeds execution$denseMemoryMB.");
      }
      const Rcpp::List factors=descriptor["factors"];
      for(int axis=0;axis<factors.size();++axis){
        const Rcpp::List factor=factors[axis];
        if(!zeroSpecificFactor(factor)) continue;
        const int first=Rcpp::as<int>(factor["par_start"])-1;
        const Rcpp::LogicalVector diagonal=factor["rr_diag"];
        for(int j=0;j<diagonal.size();++j){
          if(!diagonal[j]) continue;
          const double pivot=std::exp(parameters[term](first+j)-parameters[term](first));
          if(!std::isfinite(pivot) || pivot<=0){
            reason="reduced-rank loading pivot outside numerical range";
            return false;
          }
        }
      }
      const arma::sp_mat map=coordinateMap(term);
      if(!map.is_finite()){ reason="non-finite reduced-rank loadings"; return false; }
      maps[term]=arma::kron(map,arma::speye<arma::sp_mat>(levels[term],levels[term]));
      const arma::sp_mat block=original.cols(fullStarts[term],fullStarts[term+1]-1)*maps[term];
      design.cols(latentStarts[term],latentStarts[term+1]-1)=block;
    }
    return true;
  }

  arma::vec restore(const arma::vec & mode) const {
    if(!reduced) return mode;
    arma::vec full(fullStarts.back());
    for(int term=0;term<descriptors.size();++term)
      full.subvec(fullStarts[term],fullStarts[term+1]-1)=maps[term]*
        mode.subvec(latentStarts[term],latentStarts[term+1]-1);
    return full;
  }

  arma::sp_mat designDerivative(arma::uword coordinate) const {
    arma::sp_mat change(original.n_rows,latentStarts.back());
    if(precisionPositions[coordinate]>=0) return change;
    const int term=freeTerms[coordinate];
    const arma::sp_mat map=arma::kron(coordinateMap(term,freePositions[coordinate]),
      arma::speye<arma::sp_mat>(levels[term],levels[term]));
    change.cols(latentStarts[term],latentStarts[term+1]-1)=
      original.cols(fullStarts[term],fullStarts[term+1]-1)*map;
    return change;
  }

  arma::sp_mat project(const arma::sp_mat & source, int derivative=-1) const {
    if(derivative<0 && !reduced) return source;
    arma::sp_mat result(source.n_rows,latentStarts.back());
    if(derivative>=0){
      if(precisionPositions[derivative]>=0) return result;
      const int term=freeTerms[derivative];
      const arma::sp_mat map=arma::kron(coordinateMap(term,freePositions[derivative]),
        arma::speye<arma::sp_mat>(levels[term],levels[term]));
      result.cols(latentStarts[term],latentStarts[term+1]-1)=
        source.cols(fullStarts[term],fullStarts[term+1]-1)*map;
    }else{
      for(int term=0;term<descriptors.size();++term)
        result.cols(latentStarts[term],latentStarts[term+1]-1)=
          source.cols(fullStarts[term],fullStarts[term+1]-1)*maps[term];
    }
    return result;
  }

  std::vector<LaplacePrecisionDerivative> derivatives(){
    const auto latent=precision->derivatives();
    std::vector<LaplacePrecisionDerivative> result;
    for(arma::uword j=0;j<parameterCount();++j){
      if(precisionPositions[j]>=0) result.push_back(latent[precisionPositions[j]]);
      else{
        const int term=freeTerms[j];
        const arma::uword width=latentStarts[term+1]-latentStarts[term];
        result.push_back({latentStarts[term],arma::sp_mat(width,width),0});
      }
    }
    return result;
  }
};

} // namespace sommer
#endif
