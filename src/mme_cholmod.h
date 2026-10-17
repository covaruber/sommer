#ifndef SOMMER_MME_CHOLMOD_H
#define SOMMER_MME_CHOLMOD_H

#include <Rcpp.h>
#include <Eigen/Core>
#include <Matrix/Matrix.h>
#include <cstring>
#include <string>
#include "mme_sparse.h"

namespace sommer {

// A cholmod_sparse VIEW over an Eigen sparse matrix's existing compressed-
// storage arrays (zero-copy): CHOLMOD's analyze/factorize only read A, so
// a const_cast to CHOLMOD's void* fields is safe.
inline cholmod_sparse eigenToCholmodSparseView(const Eigen::SparseMatrix<double> & matrix) {
  cholmod_sparse view;
  std::memset(&view, 0, sizeof(view));
  view.nrow = static_cast<size_t>(matrix.rows());
  view.ncol = static_cast<size_t>(matrix.cols());
  view.nzmax = static_cast<size_t>(matrix.nonZeros());
  view.p = const_cast<int *>(matrix.outerIndexPtr());
  view.i = const_cast<int *>(matrix.innerIndexPtr());
  view.x = const_cast<double *>(matrix.valuePtr());
  // Both matrices are symmetric with both triangles populated; stype=-1
  // tells CHOLMOD to read only the lower triangle and assume symmetry
  // for the rest. stype=0 ("general") would instead make CHOLMOD
  // factorise A*A'.
  view.stype = -1;
  view.itype = CHOLMOD_INT;
  view.xtype = CHOLMOD_REAL;
  view.dtype = CHOLMOD_DOUBLE;
  view.sorted = 1;
  view.packed = 1;
  return view;
}

// CholmodState owns the CHOLMOD context and can survive successive PQL fits.
struct CholmodState {
  cholmod_common common;
  cholmod_factor * factor = nullptr;
  bool started = false;
  bool symbolicReady = false;
  std::vector<int> outerPattern;
  std::vector<int> innerPattern;

  CholmodState() = default;
  CholmodState(const CholmodState &) = delete;
  CholmodState & operator=(const CholmodState &) = delete;

  void ensureStarted() {
    if(!started) {
      M_cholmod_start(&common);
      common.supernodal = CHOLMOD_SUPERNODAL;
      started = true;
    }
  }

  bool analyzePattern(const Eigen::SparseMatrix<double> & matrix, int * analyses = nullptr) {
    ensureStarted();
    const bool samePattern = symbolicReady &&
      eigenSparsePatternMatches(matrix, outerPattern, innerPattern);
    if(!samePattern || factor == nullptr) {
      if(factor != nullptr) { M_cholmod_free_factor(&factor, &common); }
      cholmod_sparse view = eigenToCholmodSparseView(matrix);
      factor = M_cholmod_analyze(&view, &common);
      if(analyses != nullptr) { ++*analyses; }
      if(factor == nullptr || common.status != CHOLMOD_OK) { return false; }
      cacheEigenSparsePattern(matrix, outerPattern, innerPattern);
      symbolicReady = true;
    }
    return true;
  }

  bool factorize(const Eigen::SparseMatrix<double> & matrix) {
    cholmod_sparse view = eigenToCholmodSparseView(matrix);
    return M_cholmod_factorize(&view, factor, &common) != 0 &&
      common.status == CHOLMOD_OK;
  }

  Eigen::MatrixXd solve(const Eigen::Ref<const Eigen::MatrixXd> & rhs,
                        const std::string & failureMessage) {
    cholmod_dense rhsView;
    std::memset(&rhsView, 0, sizeof(rhsView));
    rhsView.nrow = static_cast<size_t>(rhs.rows());
    rhsView.ncol = static_cast<size_t>(rhs.cols());
    rhsView.nzmax = static_cast<size_t>(rhs.size());
    rhsView.d = static_cast<size_t>(rhs.rows());
    // CHOLMOD dense matrices are column-major, matching Eigen's default
    // storage order, so rhs's buffer can be handed over without copying.
    rhsView.x = const_cast<double *>(rhs.data());
    rhsView.xtype = CHOLMOD_REAL;
    rhsView.dtype = CHOLMOD_DOUBLE;
    cholmod_dense * solution = M_cholmod_solve(CHOLMOD_A, factor, &rhsView, &common);
    if(solution == nullptr) { Rcpp::stop(failureMessage); }
    Eigen::MatrixXd result = Eigen::Map<const Eigen::MatrixXd>(
      static_cast<double *>(solution->x), rhs.rows(), rhs.cols());
    M_cholmod_free_dense(&solution, &common);
    return result;
  }

  ~CholmodState() {
    if(factor != nullptr) { M_cholmod_free_factor(&factor, &common); }
    if(started) { M_cholmod_finish(&common); }
  }
};

struct CholmodFitCache {
  CholmodState joint;
  CholmodState random;
};

} // namespace sommer

#endif
