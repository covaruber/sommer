#include <RcppEigen.h>
#include "mme_paths.h"

namespace {
using Sparse = Eigen::SparseMatrix<double>;
using Dense = Eigen::MatrixXd;
using Vector = Eigen::VectorXd;

struct ConditionalWorkspace {
  Eigen::SimplicialLDLT<Sparse> factor;
  std::vector<int> outer, inner;
  int analyses = 0;
  int factorizations = 0;
};
}

// Exact posterior moments on a requested sparse pattern. The workspace owns
// symbolic analysis only; every call performs a fresh numeric factorization.
// [[Rcpp::export]]
Rcpp::List mme_conditional_moments_cpp(
    const Eigen::MappedSparseMatrix<double>& H,
    const Eigen::MappedSparseMatrix<double>& pattern,
    const Eigen::VectorXd& rhs,
    SEXP workspace,
    bool moments,
    int panelSize,
    double memoryMB) {
  const int size = H.rows();
  if(size < 1 || H.cols() != size || pattern.rows() != size ||
     pattern.cols() != size || rhs.size() != size || panelSize < 1 ||
     !std::isfinite(memoryMB) || memoryMB <= 0 || !rhs.allFinite()){
    Rcpp::stop("Invalid conditional posterior system or workspace controls.");
  }
  const int width = std::min(size, panelSize);
  const double panelBytes = sommer::WorkspacePolicy::inversePanelBytes(size, width);
  if(moments && panelBytes > memoryMB * 1048576.0){
    Rcpp::stop("Conditional inverse panels exceed execution$denseMemoryMB; reduce output$panelSize.");
  }
  Rcpp::XPtr<ConditionalWorkspace> cache = Rf_isNull(workspace)
    ? Rcpp::XPtr<ConditionalWorkspace>(new ConditionalWorkspace(), true)
    : Rcpp::XPtr<ConditionalWorkspace>(workspace);
  Sparse system = pattern;
  for(int column = 0; column < size; ++column){
    for(Sparse::InnerIterator entry(system, column); entry; ++entry){
      entry.valueRef() = H.coeff(entry.row(), column);
      if(!std::isfinite(entry.value())) Rcpp::stop("Non-finite conditional precision.");
    }
  }
  if((system - H).norm() > 1e-12 * (H.norm() + 1.0)){
    Rcpp::stop("Conditional moment pattern omits precision entries.");
  }
  std::vector<int> outer(system.outerIndexPtr(), system.outerIndexPtr() + size + 1);
  std::vector<int> inner(system.innerIndexPtr(), system.innerIndexPtr() + system.nonZeros());
  if(cache->outer != outer || cache->inner != inner){
    cache->factor.analyzePattern(system);
    cache->outer = std::move(outer);
    cache->inner = std::move(inner);
    ++cache->analyses;
  }
  cache->factor.factorize(system);
  ++cache->factorizations;
  if(cache->factor.info() != Eigen::Success ||
     !cache->factor.vectorD().allFinite() ||
     (cache->factor.vectorD().array() <= 0).any()){
    Rcpp::stop("Conditional posterior precision is not positive definite; check fixed-effect estimability and covariance starting values.");
  }
  Vector mean = cache->factor.solve(rhs);
  if(!mean.allFinite() || (system * mean - rhs).norm() >
     1e-7 * (rhs.norm() + 1.0)){
    Rcpp::stop("Conditional posterior solve failed its residual check.");
  }
  Sparse selected;
  if(moments){
    selected = pattern;
    for(int start = 0; start < size; start += width){
      Rcpp::checkUserInterrupt();
      const int columns = std::min(width, size - start);
      Dense unit = Dense::Zero(size, columns);
      for(int j = 0; j < columns; ++j) unit(start + j, j) = 1.0;
      const Dense inverse = cache->factor.solve(unit);
      if(!inverse.allFinite()) Rcpp::stop("Non-finite conditional inverse panel.");
      for(int j = 0; j < columns; ++j){
        for(Sparse::InnerIterator entry(selected, start + j); entry; ++entry){
          entry.valueRef() = inverse(entry.row(), j);
        }
      }
    }
  }
  return Rcpp::List::create(
    Rcpp::Named("mean")=mean,
    Rcpp::Named("covariance")=moments ? Rcpp::wrap(selected) : R_NilValue,
    Rcpp::Named("logDet")=cache->factor.vectorD().array().log().sum(),
    Rcpp::Named("workspace")=cache,
    Rcpp::Named("symbolicAnalyses")=cache->analyses,
    Rcpp::Named("numericFactorizations")=cache->factorizations,
    Rcpp::Named("panelColumns")=moments ? width : 0,
    Rcpp::Named("panelBytes")=moments ? panelBytes : 0.0,
    Rcpp::Named("selectedEntries")=moments ? selected.nonZeros() : 0);
}

// [[Rcpp::export]]
Eigen::MatrixXd mme_conditional_solve_cpp(SEXP workspace,
    const Eigen::MatrixXd& rhs, double memoryMB) {
  if(Rf_isNull(workspace) || !rhs.allFinite() ||
     !std::isfinite(memoryMB) || memoryMB <= 0){
    Rcpp::stop("Invalid conditional factor solve.");
  }
  Rcpp::XPtr<ConditionalWorkspace> cache(workspace);
  if(cache->factorizations == 0 || rhs.rows() != cache->factor.rows() ||
     sommer::WorkspacePolicy::inversePanelBytes(rhs.rows(), rhs.cols()) >
       memoryMB * 1048576.0){
    Rcpp::stop("Conditional factor solve dimensions or dense workspace budget exceeded.");
  }
  Dense solution = cache->factor.solve(rhs);
  if(cache->factor.info() != Eigen::Success || !solution.allFinite()){
    Rcpp::stop("Conditional factor solve failed.");
  }
  return solution;
}
