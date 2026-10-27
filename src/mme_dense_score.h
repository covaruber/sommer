#ifndef SOMMER_MME_DENSE_SCORE_H
#define SOMMER_MME_DENSE_SCORE_H

#include <Eigen/SparseCore>
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include "mme_control.h"

namespace sommer {

struct DenseScorePolicy {
  static bool select(control::Mode mode, double rowProducts,
                     std::size_t widestRow, std::size_t effects) {
    if(mode == control::Mode::Off || effects == 0 ||
       !std::isfinite(rowProducts) || rowProducts < 0.0) { return false; }
    if(mode == control::Mode::Force) { return true; }
    return widestRow >= 64 && rowProducts >= 1000000.0 &&
      rowProducts >= 8.0 * static_cast<double>(effects) * effects;
  }

  static double workspaceDoubles(std::size_t nonzeros, std::size_t effects) {
    return 3.0 * nonzeros + 6.0 * static_cast<double>(effects) * effects;
  }
};

inline Eigen::SparseMatrix<double> denseDesignScoreCrossproduct(
    const Eigen::SparseMatrix<double> & design, const Eigen::VectorXd & weights) {
  if(weights.size() != design.rows() || !weights.allFinite()) {
    throw std::invalid_argument("Invalid diagonal residual score cross-product weights.");
  }
  const Eigen::SparseMatrix<double> weighted = weights.asDiagonal() * design;
  Eigen::SparseMatrix<double> result =
    Eigen::SparseMatrix<double>(design.transpose()) * weighted;
  result.prune(0.0);
  result.makeCompressed();
  return result;
}

} // namespace sommer
#endif
