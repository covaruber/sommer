#ifndef SOMMER_MME_SPARSE_H
#define SOMMER_MME_SPARSE_H

#include <Eigen/SparseCore>
#include <algorithm>
#include <vector>

namespace sommer {

inline bool eigenSparsePatternMatches(
    const Eigen::SparseMatrix<double> & matrix,
    const std::vector<int> & outer, const std::vector<int> & inner) {
  if(outer.size() != static_cast<std::size_t>(matrix.outerSize() + 1) ||
     inner.size() != static_cast<std::size_t>(matrix.nonZeros())) {
    return false;
  }
  return std::equal(outer.begin(), outer.end(), matrix.outerIndexPtr()) &&
    std::equal(inner.begin(), inner.end(), matrix.innerIndexPtr());
}

inline void cacheEigenSparsePattern(
    const Eigen::SparseMatrix<double> & matrix,
    std::vector<int> & outer, std::vector<int> & inner) {
  outer.assign(matrix.outerIndexPtr(),
               matrix.outerIndexPtr() + matrix.outerSize() + 1);
  inner.assign(matrix.innerIndexPtr(),
               matrix.innerIndexPtr() + matrix.nonZeros());
}

enum class SparseFactorizationStatus { Success, SymbolicFailure, NumericFailure };

struct SparseFactorizationResult {
  SparseFactorizationStatus status;
  bool samePattern;
};

// Callbacks keep timing and selected-inverse invalidation with their owners.
template<class Factor, class OnAnalysis, class NumericFactorize>
SparseFactorizationResult factorizeLDLTWithCachedPattern(
    const Eigen::SparseMatrix<double> & matrix, Factor & factor,
    bool & symbolicReady, std::vector<int> & outer, std::vector<int> & inner,
    const OnAnalysis & onAnalysis, const NumericFactorize & numericFactorize) {
  const bool samePattern = symbolicReady &&
    eigenSparsePatternMatches(matrix, outer, inner);
  if(!samePattern) {
    factor.analyzePattern(matrix);
    onAnalysis();
    if(factor.info() != Eigen::Success) {
      return {SparseFactorizationStatus::SymbolicFailure, samePattern};
    }
    cacheEigenSparsePattern(matrix, outer, inner);
    symbolicReady = true;
  }

  numericFactorize();
  if(factor.info() != Eigen::Success && samePattern) {
    // Defensive retry with a fresh symbolic analysis.
    factor.analyzePattern(matrix);
    onAnalysis();
    if(factor.info() == Eigen::Success) {
      numericFactorize();
    }
  }
  return {factor.info() == Eigen::Success
    ? SparseFactorizationStatus::Success
    : SparseFactorizationStatus::NumericFailure, samePattern};
}

template<class Factor>
bool factorizeLDLTWithCachedPattern(
    const Eigen::SparseMatrix<double> & matrix, Factor & factor,
    bool & symbolicReady, std::vector<int> & outer, std::vector<int> & inner) {
  return factorizeLDLTWithCachedPattern(matrix, factor, symbolicReady, outer, inner,
    [](){}, [&](){ factor.factorize(matrix); }).status ==
      SparseFactorizationStatus::Success;
}

} // namespace sommer

#endif
