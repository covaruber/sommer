#ifndef SOMMER_MME_DIRECT_RESIDUAL_H
#define SOMMER_MME_DIRECT_RESIDUAL_H

#include "mme_direct_covariance.h"

class DirectResidualCovariance {
  Rcpp::List descriptor, factors;
  bool compact, weighted;
  arma::mat weightInverse;
  std::vector<std::vector<std::pair<int,int>>> rows;
  std::vector<std::vector<arma::uword>> digits;
  std::vector<int> owners;
  arma::uword n;

public:
  DirectResidualCovariance(const Rcpp::List & cs, const arma::uvec & block,
                           const arma::uvec & index, const arma::sp_mat & H,
                           const bool useH)
    : descriptor(cs), factors(Rcpp::as<Rcpp::List>(cs["factors"])), compact(directResidualCompact(cs)),
      weighted(useH), n(block.n_elem) {
    if(index.n_elem != n){ Rcpp::stop("Residual mappings have inconsistent lengths."); }
    const int dimension = Rcpp::as<int>(cs["dim"]);
    if(dimension < 1){ Rcpp::stop("Residual covariance dimension must be positive."); }
    arma::uword blocks = 0;
    for(arma::uword r = 0; r < n; ++r){
      if(block(r) < 1 || index(r) < 1 || index(r) > static_cast<arma::uword>(dimension)){
        Rcpp::stop("Residual block/local indices are out of range.");
      }
      blocks = std::max(blocks, block(r));
    }
    rows.resize(blocks);
    for(arma::uword r = 0; r < n; ++r){
      rows[block(r)-1].push_back({static_cast<int>(r), static_cast<int>(index(r)-1)});
    }
    digits.resize(factors.size(), std::vector<arma::uword>(dimension));
    owners.assign(factors.size(), -1);
    arma::uword trailing = dimension;
    for(int f = 0; f < factors.size(); ++f){
      Rcpp::List factor = factors[f];
      const int dim = Rcpp::as<int>(factor["dim"]);
      if(dim < 1 || trailing % dim){ Rcpp::stop("Invalid residual factor dimensions."); }
      trailing /= dim;
      for(int local = 0; local < dimension; ++local){
        digits[f][local] = (local / trailing) % dim;
      }
      if(factor.containsElementNamed("section_factor")){
        owners[f] = Rcpp::as<int>(factor["section_factor"]) - 1;
        if(owners[f] < 0 || owners[f] >= factors.size()){
          Rcpp::stop("Invalid residual section owner.");
        }
      }
    }
    if(weighted){
      arma::mat upper;
      if(H.n_rows != n || H.n_cols != n || !arma::chol(upper, arma::mat(H))){
        Rcpp::stop("Residual precision must be positive definite and match the observations.");
      }
      weightInverse = arma::inv(arma::trimatu(upper));
    }
  }

  arma::mat evaluate(const arma::vec & par, const int derivative = -1) const {
    arma::mat out(n, n, arma::fill::zeros);
    if(!compact){
      const arma::mat covariance = derivative < 0
        ? directEvaluateDescriptor(descriptor, par)
        : directDescriptorD1(descriptor, par, static_cast<arma::uword>(derivative));
      for(const auto & block : rows){
        for(const auto & first : block){
          for(const auto & second : block){
            out(first.first, second.first) = covariance(first.second, second.second);
          }
        }
      }
    }else{
      std::vector<std::vector<arma::mat>> pieces(factors.size());
      int derivativeFactor = -1, derivativeSection = -1;
      const int k1 = derivative + 1;
      for(int f = 0; f < factors.size(); ++f){
        Rcpp::List factor = factors[f];
        const int start = factor.containsElementNamed("par_start") ? Rcpp::as<int>(factor["par_start"]) : 1;
        const int end = factor.containsElementNamed("par_end") ? Rcpp::as<int>(factor["par_end"]) : 0;
        arma::vec local;
        if(end >= start){ local = par.subvec(start-1, end-1); }
        const bool mine = derivative > 0 && end >= start && k1 >= start && k1 <= end;
        const int sections = factor.containsElementNamed("section_levels") ? Rf_length(factor["section_levels"]) : 0;
        if(!sections){
          pieces[f].assign(1, mine ? directFactorD1(factor, local, k1-start)
                                    : directEvalFactor(factor, local));
          if(mine){ derivativeFactor = f; }
        }else{
          if(owners[f] < 0){ Rcpp::stop("A sectioned residual factor has no section owner."); }
          const int count = Rcpp::as<int>(factor["section_npar"]);
          const int dim = Rcpp::as<int>(factor["dim"]);
          pieces[f].resize(sections);
          for(int section = 0; section < sections; ++section){
            const int listed = directSectionListed(factor, section);
            if(listed < 0){
              pieces[f][section] = arma::eye<arma::mat>(dim, dim);
            }else{
              const arma::vec slice = directSectionSlice(factor, local, listed);
              if(mine && (k1-start)/count == listed){
                pieces[f][section] = directFactorD1(factor, slice, (k1-start)%count);
                derivativeFactor = f;
                derivativeSection = section;
              }else{
                pieces[f][section] = directEvalFactor(factor, slice);
              }
            }
          }
        }
      }
      if(derivative > 0 && derivativeFactor < 0){
        Rcpp::stop("Residual covariance derivative parameter does not belong to any factor.");
      }
      const double scale = std::exp(par(0));
      for(const auto & block : rows){
        if(block.empty()){ continue; }
        std::vector<const arma::mat *> selected(factors.size());
        bool zero = false;
        for(int f = 0; f < factors.size(); ++f){
          arma::uword section = 0;
          if(pieces[f].size() > 1 || owners[f] >= 0){
            section = digits[owners[f]][block[0].second];
            if(f == derivativeFactor && derivativeSection >= 0 &&
               section != static_cast<arma::uword>(derivativeSection)){ zero = true; }
          }
          selected[f] = &pieces[f][section];
        }
        if(zero){ continue; }
        for(const auto & first : block){
          for(const auto & second : block){
            double value = scale;
            for(int f = 0; f < factors.size() && value != 0.0; ++f){
              value *= (*selected[f])(digits[f][first.second], digits[f][second.second]);
            }
            out(first.first, second.first) = value;
          }
        }
      }
    }
    if(weighted){ return weightInverse.t() * out * weightInverse; }
    return out;
  }
};

#endif
