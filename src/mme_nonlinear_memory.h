#ifndef SOMMER_MME_NONLINEAR_MEMORY_H
#define SOMMER_MME_NONLINEAR_MEMORY_H

#include <RcppArmadillo.h>
#include <numeric>

namespace sommer {

// Component-wise dense fill bounds retain linear admission for independent groups.
inline double nonlinearSparseAdmission(const arma::sp_mat & design,
    const arma::sp_mat & precision, const arma::sp_mat & curvature,
    double memoryMB){
  const arma::uword q=design.n_cols, n=design.n_rows;
  if(precision.n_rows!=q || precision.n_cols!=q ||
     curvature.n_rows!=n || curvature.n_cols!=n ||
     !std::isfinite(memoryMB) || memoryMB<=0)
    Rcpp::stop("Invalid nonlinear sparse workspace inputs.");
  const double baseline=8.0*(180.0*n+80.0*q+
    12.0*design.n_nonzero+8.0*precision.n_nonzero+8.0*curvature.n_nonzero);
  if(baseline>memoryMB*1048576)
    Rcpp::stop("Nonlinear sparse workspace exceeds execution$denseMemoryMB before Hessian assembly.");
  std::vector<arma::uword> parent(q), size(q,1), row(n,q);
  std::iota(parent.begin(),parent.end(),0);
  auto root=[&](arma::uword index){
    while(parent[index]!=index){
      parent[index]=parent[parent[index]];
      index=parent[index];
    }
    return index;
  };
  auto join=[&](arma::uword a,arma::uword b){
    a=root(a); b=root(b);
    if(a==b) return;
    if(size[a]<size[b]) std::swap(a,b);
    parent[b]=a; size[a]+=size[b];
  };
  for(auto cell=precision.begin();cell!=precision.end();++cell)
    join(cell.row(),cell.col());
  for(auto cell=design.begin();cell!=design.end();++cell){
    if(row[cell.row()]==q) row[cell.row()]=cell.col();
    else join(row[cell.row()],cell.col());
  }
  for(auto cell=curvature.begin();cell!=curvature.end();++cell){
    const arma::uword a=row[cell.row()], b=row[cell.col()];
    if(a<q && b<q) join(a,b);
  }
  double fill=0;
  for(arma::uword j=0;j<q;++j)
    if(parent[j]==j) fill+=static_cast<double>(size[j])*size[j];
  const double bytes=baseline+64.0*fill;
  if(!std::isfinite(bytes) || bytes>memoryMB*1048576)
    Rcpp::stop("Nonlinear component-wise factor-fill bound exceeds execution$denseMemoryMB before Hessian assembly; simplify connected random structure or increase the budget.");
  return bytes;
}

}
#endif
