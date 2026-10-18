#include "../../src/mme_pcg.h"
#include <Eigen/LU>

#include <cassert>
#include <cmath>
#include <stdexcept>

int main() {
  Eigen::MatrixXd coefficient(6, 6);
  coefficient <<
    5.0, 0.3, 0.0, 0.1, 0.0, 0.0,
    0.3, 4.0, 0.2, 0.0, 0.0, 0.0,
    0.0, 0.2, 3.0, 0.4, 0.0, 0.0,
    0.1, 0.0, 0.4, 4.5, 0.2, 0.0,
    0.0, 0.0, 0.0, 0.2, 3.5, 0.3,
    0.0, 0.0, 0.0, 0.0, 0.3, 2.5;
  const Eigen::MatrixXd inverse = coefficient.inverse();
  const sommer::pcg::VectorOperator applyA =
    [&](const Eigen::VectorXd & value) -> Eigen::VectorXd {
      return coefficient * value;
    };
  const sommer::pcg::MatrixOperator applyAMatrix =
    [&](const Eigen::MatrixXd & value) -> Eigen::MatrixXd {
      return coefficient * value;
    };
  const sommer::pcg::VectorOperator applyM =
    [&](const Eigen::VectorXd & value) -> Eigen::VectorXd {
      return value.cwiseQuotient(coefficient.diagonal());
    };
  const sommer::pcg::MatrixOperator applyMMatrix =
    [&](const Eigen::MatrixXd & value) -> Eigen::MatrixXd {
      Eigen::MatrixXd result = value;
      for(Eigen::Index col = 0; col < result.cols(); ++col) {
        result.col(col).array() /= coefficient.diagonal().array();
      }
      return result;
    };

  sommer::pcg::Options options;
  options.tolerance = 1.0e-11;
  options.residualRefresh = 2;
  Eigen::MatrixXd rhs(6, 2);
  rhs << 1.0, -2.0, 2.0, 0.5, -1.0, 1.0,
         0.0, 3.0, 4.0, -2.0, 0.2, 1.0;
  Eigen::MatrixXd guess = Eigen::MatrixXd::Constant(6, 2, 0.25);
  const sommer::pcg::Result independent =
    sommer::pcg::solveIndependent(rhs, applyA, applyM, options, &guess);
  assert((independent.solution - inverse * rhs).norm() < 1.0e-9);
  const sommer::pcg::Result packed = sommer::pcg::solvePackedIndependent(
    rhs, applyAMatrix, applyMMatrix, options, &guess);
  assert((packed.solution - inverse * rhs).norm() < 1.0e-9);
  assert(packed.metrics.matrixPanelApplications > 0);
  assert(packed.metrics.vectorApplications >= packed.metrics.matrixPanelApplications);
  assert(independent.metrics.trueResidualChecks >= 2);
  assert(independent.metrics.operatorApplications >
         independent.metrics.iterations);
  assert(independent.metrics.residualReplacements > 0);
  assert((coefficient * independent.solution - rhs).norm() <
         options.tolerance * rhs.norm());

  Eigen::MatrixXd blockRhs(6, 4);
  blockRhs.col(0) = rhs.col(0);
  blockRhs.col(1) = rhs.col(1);
  blockRhs.col(2) = 2.0 * rhs.col(0) - 3.0 * rhs.col(1);
  blockRhs.col(3) = rhs.col(0) + rhs.col(1);
  options.residualRefresh = 0;
  options.blockSize = 4;
  const sommer::pcg::Result block = sommer::pcg::solveBlock(
    blockRhs, applyAMatrix, applyMMatrix, applyA, applyM, options);
  assert((block.solution - inverse * blockRhs).norm() < 1.0e-8);
  assert((coefficient * block.solution - blockRhs).norm() <
         options.tolerance * blockRhs.norm());
  assert(block.metrics.deflatedColumns == 2);
  assert(block.metrics.lockedColumns >= 2);
  assert(block.metrics.blockFallbacks == 0);

  Eigen::MatrixXd partlySolvedRhs = rhs;
  partlySolvedRhs.col(0) = coefficient.col(0);
  Eigen::MatrixXd partlySolvedGuess = Eigen::MatrixXd::Zero(6, 2);
  partlySolvedGuess(0,0) = 1.0;
  const sommer::pcg::Result partlySolved = sommer::pcg::solveBlock(
    partlySolvedRhs, applyAMatrix, applyMMatrix, applyA, applyM,
    options, &partlySolvedGuess);
  assert((partlySolved.solution - inverse * partlySolvedRhs).norm() < 1.0e-8);
  assert(partlySolved.metrics.lockedColumns >= 2);
  assert(partlySolved.metrics.blockFallbacks == 0);

  Eigen::MatrixXd unequalRhs = rhs;
  unequalRhs.col(1) *= 1.0e-14;
  std::size_t measuredCalls = 0;
  const sommer::pcg::Result unequal = sommer::pcg::solveBlock(
    unequalRhs,
    [&](const Eigen::MatrixXd & values) -> Eigen::MatrixXd {
      ++measuredCalls;
      return applyAMatrix(values);
    }, applyMMatrix,
    [&](const Eigen::VectorXd & value) -> Eigen::VectorXd {
      ++measuredCalls;
      return applyA(value);
    }, applyM, options);
  assert((coefficient * unequal.solution.col(1) - unequalRhs.col(1)).norm() <=
    options.tolerance * unequalRhs.col(1).norm());
  assert(unequal.metrics.blockFallbacks == 1);
  assert(unequal.metrics.operatorApplications == measuredCalls);

  options.maxWorkspaceBytes = 1;
  const sommer::pcg::Result bounded = sommer::pcg::solveBlock(
    rhs, applyAMatrix, applyMMatrix, applyA, applyM, options);
  assert((bounded.solution - inverse * rhs).norm() < 1.0e-8);
  assert(bounded.metrics.blockFallbacks > 0);
  assert(bounded.metrics.fallbackReason == "workspace-cap");

  sommer::pcg::Options limited;
  limited.tolerance = 1.0e-14;
  limited.maxIterations = 1;
  bool stoppedAtLimit = false;
  try {
    sommer::pcg::solveIndependent(rhs, applyA, applyM, limited);
  } catch(const std::runtime_error &) {
    stoppedAtLimit = true;
  }
  assert(stoppedAtLimit);

  Eigen::MatrixXd zeroRhs = Eigen::MatrixXd::Zero(6, 1);
  const sommer::pcg::Result zero = sommer::pcg::solveIndependent(
    zeroRhs, applyA, applyM, options);
  assert(zero.solution.norm() == 0.0);

  options.maxWorkspaceBytes = 64u * 1024u * 1024u;
  bool operatorFailurePropagated = false;
  int scalarFallbackCalls = 0;
  try {
    sommer::pcg::solveBlock(rhs,
      [](const Eigen::MatrixXd &) -> Eigen::MatrixXd {
        throw std::runtime_error("operator-failure-must-propagate");
      }, applyMMatrix,
      [&](const Eigen::VectorXd & value) -> Eigen::VectorXd {
        ++scalarFallbackCalls;
        return applyA(value);
      }, applyM, options);
  } catch(const std::runtime_error & error) {
    operatorFailurePropagated = std::string(error.what()) == "operator-failure-must-propagate";
  }
  assert(operatorFailurePropagated);
  assert(scalarFallbackCalls == 0);
}
