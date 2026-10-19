#include "../../src/mme_pcg_estimators.h"
#include <cassert>
#include <cmath>
#include <stdexcept>

using namespace sommer::pcg_estimators;

static VectorOperator matrixOperator(const Eigen::MatrixXd & matrix) {
  return [matrix](const Eigen::VectorXd & value) -> Eigen::VectorXd {
    return matrix * value;
  };
}

int main() {
  Eigen::MatrixXd identity = Eigen::MatrixXd::Identity(5, 5);
  const VectorOperator identityOperator = matrixOperator(identity);
  LanczosLogDetOptions logOptions;
  logOptions.probes = 3;
  logOptions.maxSteps = 5;
  LogDetResult logdet = estimateLogDet(5, identityOperator, Eigen::VectorXd(),
                                       logOptions);
  assert(std::abs(logdet.estimate) < 1e-12);
  assert(logdet.operatorApplications == 3);
  assert(logdet.probeEstimates.size() == 3);
  assert(logdet.probeSE < 1e-12);

  Eigen::MatrixXd diagonalMatrix = Eigen::MatrixXd::Zero(4, 4);
  diagonalMatrix.diagonal() << 0.03, 0.4, 3.0, 70.0;
  const VectorOperator diagonalOperator = matrixOperator(diagonalMatrix);
  Eigen::VectorXd diagonal = diagonalMatrix.diagonal();
  logOptions.seed = 0;
  logOptions.probes = 4;
  logOptions.maxSteps = 4;
  logOptions.scaling = LogDetScaling::Diagonal;
  logdet = estimateLogDet(4, diagonalOperator, diagonal, logOptions);
  assert(std::abs(logdet.estimate - diagonal.array().log().sum()) < 1e-10);
  const VectorOperator scaled = diagonallyScaledOperator(4, diagonalOperator, diagonal);
  assert((scaled(Eigen::VectorXd::Ones(4)) - Eigen::VectorXd::Ones(4)).norm() < 1e-12);
  Eigen::VectorXd originalGuess(4);
  originalGuess << 2.0, -1.0, 0.5, 4.0;
  const Eigen::VectorXd scaledGuess =
    mapInverseGuessToScaled(4, originalGuess, diagonal);
  assert((mapScaledSolutionToOriginal(4, scaledGuess, diagonal) -
          originalGuess).norm() < 1e-12);
  logOptions.scaling = LogDetScaling::None;
  logdet = estimateLogDet(4, diagonalOperator, Eigen::VectorXd(), logOptions);
  assert(std::abs(logdet.estimate - diagonal.array().log().sum()) < 1e-10);

  Eigen::MatrixXd illConditioned = Eigen::MatrixXd::Zero(3, 3);
  illConditioned.diagonal() << 1e-4, 2.0, 1e4;
  LanczosLogDetOptions illOptions;
  illOptions.probes = 2;
  illOptions.maxSteps = 3;
  logdet = estimateLogDet(3, matrixOperator(illConditioned),
                          Eigen::VectorXd(), illOptions);
  assert(std::abs(logdet.estimate -
                  illConditioned.diagonal().array().log().sum()) < 1e-7);
  const SensitivityAssessment depth =
    assessDepthSensitivity(3, matrixOperator(illConditioned),
                           Eigen::VectorXd(), illOptions, 2);
  assert(depth.absoluteDifference >= 0.0);
  const SensitivityAssessment probes = assessProbeSensitivity(logdet, 1);
  assert(probes.absoluteDifference >= 0.0);

  bool rejectedNonPositive = false;
  try {
    Eigen::MatrixXd indefinite(2, 2);
    indefinite << 1.0, 0.0, 0.0, -0.1;
    LanczosLogDetOptions invalidOptions;
    invalidOptions.probes = 1;
    invalidOptions.maxSteps = 2;
    estimateLogDet(2, matrixOperator(indefinite), Eigen::VectorXd(),
                   invalidOptions);
  } catch(const std::runtime_error &) {
    rejectedNonPositive = true;
  }
  assert(rejectedNonPositive);

  Eigen::MatrixXd coefficient(4, 4);
  coefficient << 4.0, 0.4, 0.1, 0.0,
                 0.4, 2.5, 0.2, 0.1,
                 0.1, 0.2, 1.7, 0.3,
                 0.0, 0.1, 0.3, 3.2;
  Eigen::MatrixXd derivative(4, 4);
  derivative << 1.0, 0.2, 0.0, 0.0,
                0.2, -2.0, 0.1, 0.0,
                0.0, 0.1, 0.5, 0.3,
                0.0, 0.0, 0.3, -0.7;
  const Eigen::MatrixXd inverse = coefficient.inverse();
  const double exactTrace = (inverse * derivative).trace();
  const VectorOperator applyInverse = matrixOperator(inverse);
  const VectorOperator applyDerivative = matrixOperator(derivative);
  TraceOptions traceOptions;
  traceOptions.seed = 19;
  traceOptions.probes = 6000;
  const TraceResult hutchinson = estimateInverseTrace(
    4, applyInverse, applyDerivative, traceOptions);
  assert(std::abs(hutchinson.estimate - exactTrace) < 0.12);
  assert(hutchinson.inverseApplications == traceOptions.probes);
  assert(hutchinson.derivativeApplications == traceOptions.probes);
  assert(hutchinson.probeSE > 0.0);

  const LandmarkControlVariate control =
    makeLandmarkControlVariate(4, matrixOperator(coefficient), 2, 7);
  const TraceResult corrected = estimateInverseTrace(
    4, applyInverse, applyDerivative, traceOptions, &control);
  assert(std::abs(corrected.estimate - exactTrace) < 0.12);
  assert(corrected.derivativeApplications == traceOptions.probes + 2);

  const LandmarkControlVariate fullControl =
    makeLandmarkControlVariate(4, matrixOperator(coefficient), 4, 7);
  const TraceResult exactControl = estimateInverseTrace(
    4, applyInverse, applyDerivative, traceOptions, &fullControl);
  assert(std::abs(exactControl.estimate - exactTrace) < 1e-10);
  assert(exactControl.probeSE < 1e-10);

  bool rejectedMemoryBudget = false;
  try {
    makeLandmarkControlVariate(4, matrixOperator(coefficient), 2, 0, 8);
  } catch(const std::runtime_error &) {
    rejectedMemoryBudget = true;
  }
  assert(rejectedMemoryBudget);

  bool rejectedBadDiagonal = false;
  try {
    Eigen::VectorXd invalidDiagonal = Eigen::VectorXd::Ones(4);
    invalidDiagonal(2) = 0.0;
    diagonallyScaledOperator(4, diagonalOperator, invalidDiagonal);
  } catch(const std::runtime_error &) {
    rejectedBadDiagonal = true;
  }
  assert(rejectedBadDiagonal);
}
