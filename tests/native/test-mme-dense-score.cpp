#include "../../src/mme_dense_score.h"
#include <cassert>

int main() {
  using sommer::DenseScorePolicy;
  using sommer::control::Mode;
  assert(DenseScorePolicy::select(Mode::Auto, 54930675, 409, 548));
  assert(!DenseScorePolicy::select(Mode::Off, 54930675, 409, 548));
  assert(!DenseScorePolicy::select(Mode::Auto, 54930675, 63, 548));
  assert(!DenseScorePolicy::select(Mode::Auto, 999999, 64, 100));
  assert(DenseScorePolicy::select(Mode::Auto, 1000000, 64, 100));
  assert(!DenseScorePolicy::select(Mode::Auto, 7999999, 64, 1000));
  assert(DenseScorePolicy::select(Mode::Auto, 8000000, 64, 1000));
  assert(DenseScorePolicy::select(Mode::Force, 1, 1, 1000));
  assert(!DenseScorePolicy::select(Mode::Force, 1, 1, 0));
  assert(!DenseScorePolicy::select(Mode::Force, -1, 1, 100));
  assert(DenseScorePolicy::workspaceDoubles(10, 4) == 126);

  Eigen::MatrixXd design(5, 3);
  design << 1, 0, 2, 0, 3, 0, 4, 5, 0, 1, 2, 3, 0, 0, 0;
  Eigen::VectorXd weights(5);
  weights << 1, 0, -2, 0.3, 2;
  const Eigen::SparseMatrix<double> sparse = design.sparseView();
  const Eigen::MatrixXd reference = design.transpose() * weights.asDiagonal() * design;
  const Eigen::MatrixXd candidate = sommer::denseDesignScoreCrossproduct(sparse, weights);
  assert((candidate - reference).norm() < 1e-12);
  const Eigen::SparseMatrix<double> zero =
    sommer::denseDesignScoreCrossproduct(sparse, Eigen::VectorXd::Zero(5));
  assert(zero.rows() == 3 && zero.cols() == 3 && zero.nonZeros() == 0);
  bool invalid = false;
  try {
    sommer::denseDesignScoreCrossproduct(sparse, Eigen::VectorXd::Zero(4));
  } catch(const std::invalid_argument &) {
    invalid = true;
  }
  assert(invalid);
}
