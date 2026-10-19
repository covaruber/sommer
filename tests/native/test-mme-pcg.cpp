#include "../../src/mme_pcg.h"
#include <Eigen/LU>

#include <cassert>
#include <cmath>
#include <iostream>
#include <limits>
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

  const Eigen::Index recycleDimension = 72;
  Eigen::VectorXd spectrum(recycleDimension);
  spectrum(0) = 1.0e-4;
  spectrum(1) = 2.0e-3;
  for(Eigen::Index i = 2; i < recycleDimension; ++i) {
    spectrum(i) = 1.0 + 4.0 * (i - 2);
  }
  const Eigen::MatrixXd recycleCoefficient = spectrum.asDiagonal();
  const sommer::pcg::VectorOperator recycleA =
    [&](const Eigen::VectorXd & value) -> Eigen::VectorXd {
      return recycleCoefficient * value;
    };
  const sommer::pcg::VectorOperator identityM =
    [](const Eigen::VectorXd & value) -> Eigen::VectorXd { return value; };
  sommer::pcg::RecycleSpace recycle(12, 4);
  recycle.capture(Eigen::VectorXd::Unit(recycleDimension, 0));
  recycle.capture(1.0e-250 * Eigen::VectorXd::Unit(recycleDimension, 1));
  recycle.capture(Eigen::VectorXd::Unit(recycleDimension, 0));
  recycle.capture(Eigen::VectorXd::Zero(recycleDimension));
  assert(recycle.candidateCount() == 3);
  sommer::pcg::RecycleSpace boundedCandidates(2, 1);
  boundedCandidates.capture(Eigen::VectorXd::Unit(recycleDimension, 0));
  boundedCandidates.capture(Eigen::VectorXd::Unit(recycleDimension, 1));
  boundedCandidates.capture(Eigen::VectorXd::Unit(recycleDimension, 2));
  const Eigen::MatrixXd boundedCandidateSnapshot =
    boundedCandidates.candidates();
  assert(boundedCandidates.candidateCount() == 2);
  assert(boundedCandidateSnapshot.cols() == 2);
  assert(boundedCandidateSnapshot.col(0).isApprox(
    Eigen::VectorXd::Unit(recycleDimension, 1)));
  sommer::pcg::RecycleSpace independentBounds(2, 4);
  assert(independentBounds.maxStoredDirections() == 10);
  Eigen::VectorXd metric = Eigen::VectorXd::Ones(recycleDimension);
  sommer::pcg::RecycleSpace::Selection selected =
    recycle.selectCandidates(recycleA, metric, 4);
  assert(selected.basis.cols() == 2);
  assert(std::abs(selected.ritzValues(0) - spectrum(0)) < 1.0e-12);
  assert(std::abs(selected.ritzValues(1) - spectrum(1)) < 1.0e-12);
  recycle.accept(selected);
  assert(recycle.candidateCount() == 0);
  recycle.prepare(recycleA, identityM);
  assert(recycle.prepared());
  assert(recycle.coarseMatrix().rows() == 2);
  Eigen::MatrixXd coarseAction = recycle.applyPanel(
    recycleCoefficient * selected.basis);
  assert((coarseAction - selected.basis).norm() < 1.0e-10);
  Eigen::VectorXd testVector(recycleDimension);
  for(Eigen::Index i = 0; i < recycleDimension; ++i) {
    testVector(i) = std::sin(0.37 * (i + 1));
  }
  Eigen::VectorXd secondVector(recycleDimension);
  for(Eigen::Index i = 0; i < recycleDimension; ++i) {
    secondVector(i) = std::cos(0.19 * (i + 1));
  }
  const Eigen::VectorXd recycledTest = recycle.apply(testVector);
  const Eigen::VectorXd recycledSecond = recycle.apply(secondVector);
  assert(testVector.dot(recycledTest) > 0.0);
  assert(std::abs(testVector.dot(recycledSecond) -
                  secondVector.dot(recycledTest)) < 1.0e-10);
  Eigen::LLT<Eigen::MatrixXd> recycleCoarseFactor(recycle.coarseMatrix());
  assert((recycle.coarseCorrection(testVector) -
          selected.basis * recycleCoarseFactor.solve(
            selected.basis.transpose() * testVector)).norm() < 1.0e-10);
  assert(recycle.storageBytesUpperBound() >=
         static_cast<std::size_t>(recycleDimension) *
           (recycle.maxCandidates() + 2 * recycle.maxRank()) *
           sizeof(double));

  sommer::pcg::RecycleSpace mergedRecycle(8, 3);
  mergedRecycle.capturePanel(selected.basis);
  sommer::pcg::RecycleSpace::Selection mergedFirst =
    mergedRecycle.selectCandidates(recycleA, metric, 2);
  mergedRecycle.accept(mergedFirst);
  const Eigen::MatrixXd retainedBasis = mergedRecycle.basis();
  mergedRecycle.capture(Eigen::VectorXd::Unit(recycleDimension, 2));
  sommer::pcg::RecycleSpace::Selection mergedSelection =
    mergedRecycle.selectCandidates(recycleA, metric, 2);
  assert(mergedSelection.basis.cols() == 2);
  assert((retainedBasis - mergedSelection.basis *
          (mergedSelection.basis.transpose() * retainedBasis)).norm() <
         1.0e-10);
  assert(std::abs(mergedSelection.ritzValues(0) - spectrum(0)) < 1.0e-12);
  assert(std::abs(mergedSelection.ritzValues(1) - spectrum(1)) < 1.0e-12);

  Eigen::MatrixXd denseSeed = Eigen::MatrixXd::Random(recycleDimension,
                                                       recycleDimension);
  Eigen::MatrixXd denseA = denseSeed.transpose() * denseSeed +
    Eigen::MatrixXd::Identity(recycleDimension, recycleDimension);
  Eigen::MatrixXd preconditionerSeed = Eigen::MatrixXd::Random(
    recycleDimension, recycleDimension);
  Eigen::MatrixXd denseM = preconditionerSeed * preconditionerSeed.transpose() +
    Eigen::MatrixXd::Identity(recycleDimension, recycleDimension);
  const sommer::pcg::VectorOperator denseApplyA =
    [&](const Eigen::VectorXd & value) -> Eigen::VectorXd {
      return denseA * value;
    };
  const sommer::pcg::VectorOperator denseApplyM =
    [&](const Eigen::VectorXd & value) -> Eigen::VectorXd {
      return denseM * value;
    };
  sommer::pcg::RecycleSpace denseRecycle(8, 3);
  denseRecycle.capture(Eigen::VectorXd::Unit(recycleDimension, 0));
  denseRecycle.capture(Eigen::VectorXd::Unit(recycleDimension, 1));
  denseRecycle.capture(Eigen::VectorXd::Unit(recycleDimension, 2));
  sommer::pcg::RecycleSpace::Selection denseSelection =
    denseRecycle.selectCandidates(denseApplyA, metric, 2);
  denseRecycle.accept(denseSelection);
  denseRecycle.prepare(denseApplyA, denseApplyM);
  const Eigen::MatrixXd & denseBasis = denseRecycle.basis();
  const Eigen::MatrixXd denseCoarse =
    denseBasis.transpose() * denseA * denseBasis;
  const Eigen::MatrixXd denseCoarseInverse = denseCoarse.inverse();
  const Eigen::MatrixXd denseP =
    Eigen::MatrixXd::Identity(recycleDimension, recycleDimension) -
    denseBasis * denseCoarseInverse * denseBasis.transpose() * denseA;
  const Eigen::MatrixXd densePt =
    Eigen::MatrixXd::Identity(recycleDimension, recycleDimension) -
    denseA * denseBasis * denseCoarseInverse * denseBasis.transpose();
  const Eigen::MatrixXd denseExpected =
    denseP * denseM * densePt +
    denseBasis * denseCoarseInverse * denseBasis.transpose();
  Eigen::VectorXd denseVector(recycleDimension);
  for(Eigen::Index i = 0; i < recycleDimension; ++i) {
    denseVector(i) = std::sin(0.11 * (i + 2));
  }
  const Eigen::VectorXd denseApplied = denseRecycle.apply(denseVector);
  assert((denseApplied - denseExpected * denseVector).norm() < 1.0e-8);
  Eigen::VectorXd denseSecondVector(recycleDimension);
  for(Eigen::Index i = 0; i < recycleDimension; ++i) {
    denseSecondVector(i) = std::cos(0.23 * (i + 3));
  }
  const Eigen::VectorXd denseSecondApplied =
    denseRecycle.apply(denseSecondVector);
  assert(denseVector.dot(denseApplied) > 0.0);
  assert(std::abs(denseVector.dot(denseSecondApplied) -
                  denseSecondVector.dot(denseApplied)) < 1.0e-8);
  assert((denseRecycle.coarseCorrection(denseVector) -
          denseBasis * denseCoarseInverse * denseBasis.transpose() *
            denseVector).norm() < 1.0e-10);

  sommer::pcg::RecycleSpace::Selection invalidCaptureSelection =
    denseRecycle.selectCandidates(denseApplyA, metric, 2);
  const Eigen::MatrixXd basisBeforeInvalidCapture = denseRecycle.basis();
  bool nonfiniteCaptureRejected = false;
  try {
    Eigen::VectorXd nonfiniteDirection =
      Eigen::VectorXd::Unit(recycleDimension, 4);
    nonfiniteDirection(0) = std::numeric_limits<double>::infinity();
    denseRecycle.capture(nonfiniteDirection);
  } catch(const std::invalid_argument &) {
    nonfiniteCaptureRejected = true;
  }
  assert(nonfiniteCaptureRejected);
  assert(denseRecycle.candidateCount() == 0);
  bool dimensionMismatchRejected = false;
  try {
    denseRecycle.capture(Eigen::VectorXd::Ones(recycleDimension + 1));
  } catch(const std::invalid_argument &) {
    dimensionMismatchRejected = true;
  }
  assert(dimensionMismatchRejected);
  assert((denseRecycle.basis() - basisBeforeInvalidCapture).norm() == 0.0);
  assert(denseRecycle.prepared());
  assert(invalidCaptureSelection.basis.cols() == denseRecycle.basis().cols());

  sommer::pcg::RecycleSpace changedOperatorSpace(12, 4);
  changedOperatorSpace.capturePanel(selected.basis);
  sommer::pcg::RecycleSpace::Selection changedSelection =
    changedOperatorSpace.selectCandidates(recycleA, metric, 2);
  changedOperatorSpace.accept(changedSelection);
  Eigen::VectorXd changedSpectrum = spectrum.array() + 3.0;
  const Eigen::MatrixXd changedCoefficient = changedSpectrum.asDiagonal();
  const sommer::pcg::VectorOperator changedA =
    [&](const Eigen::VectorXd & value) -> Eigen::VectorXd {
      return changedCoefficient * value;
    };
  changedOperatorSpace.prepare(changedA, identityM);
  assert((changedOperatorSpace.coarseMatrix() -
          changedSelection.basis.transpose() * changedCoefficient *
            changedSelection.basis).norm() < 1.0e-10);
  assert((changedOperatorSpace.applyPanel(changedCoefficient *
           changedSelection.basis) - changedSelection.basis).norm() < 1.0e-10);
  bool prepareOperatorFailurePropagated = false;
  try {
    changedOperatorSpace.prepare(
      [](const Eigen::VectorXd &) -> Eigen::VectorXd {
        throw std::runtime_error("recycle-prepare-failure-must-propagate");
      }, identityM);
  } catch(const std::runtime_error & error) {
    prepareOperatorFailurePropagated =
      std::string(error.what()) == "recycle-prepare-failure-must-propagate";
  }
  assert(prepareOperatorFailurePropagated);
  assert(!changedOperatorSpace.prepared());

  sommer::pcg::RecycleSpace capturedScalar(8, 3);
  sommer::pcg::Options captureOptions;
  captureOptions.tolerance = 1.0e-10;
  captureOptions.recycleSpace = &capturedScalar;
  Eigen::MatrixXd scalarRhs = Eigen::MatrixXd::Zero(recycleDimension, 1);
  scalarRhs(0, 0) = 1.0;
  scalarRhs(1, 0) = 2.0;
  sommer::pcg::solveIndependent(scalarRhs, recycleA, identityM, captureOptions);
  assert(capturedScalar.candidateCount() > 0);

  sommer::pcg::RecycleSpace capturedPacked(8, 3);
  captureOptions.recycleSpace = &capturedPacked;
  captureOptions.blockSize = 2;
  const sommer::pcg::MatrixOperator recycleAMatrix =
    [&](const Eigen::MatrixXd & values) -> Eigen::MatrixXd {
      return recycleCoefficient * values;
    };
  const sommer::pcg::MatrixOperator identityMMatrix =
    [](const Eigen::MatrixXd & values) -> Eigen::MatrixXd { return values; };
  Eigen::MatrixXd packedRhs = Eigen::MatrixXd::Zero(recycleDimension, 2);
  packedRhs(0, 0) = 1.0;
  packedRhs(1, 0) = 1.0;
  packedRhs(2, 1) = 1.0;
  sommer::pcg::solvePackedIndependent(packedRhs, recycleAMatrix,
    identityMMatrix, captureOptions);
  assert(capturedPacked.candidateCount() > 0);

  sommer::pcg::RecycleSpace capturedBlock(8, 3);
  captureOptions.recycleSpace = &capturedBlock;
  captureOptions.blockSize = 2;
  sommer::pcg::solveBlock(packedRhs, recycleAMatrix, identityMMatrix,
    recycleA, identityM, captureOptions);
  assert(capturedBlock.candidateCount() > 0);

  captureOptions.recycleSpace = &recycle;
  sommer::pcg::Result recycledPacked = sommer::pcg::solvePackedIndependent(
    packedRhs, recycleAMatrix, identityMMatrix, captureOptions);
  assert((recycleCoefficient * recycledPacked.solution - packedRhs).norm() <
         captureOptions.tolerance * packedRhs.norm());
  sommer::pcg::Result recycledBlock = sommer::pcg::solveBlock(
    packedRhs, recycleAMatrix, identityMMatrix, recycleA, identityM,
    captureOptions);
  assert((recycleCoefficient * recycledBlock.solution - packedRhs).norm() <
         captureOptions.tolerance * packedRhs.norm());

  Eigen::MatrixXd difficultRhs(recycleDimension, 1);
  for(Eigen::Index i = 0; i < recycleDimension; ++i) {
    difficultRhs(i, 0) = 1.0 / std::sqrt(static_cast<double>(i + 1));
  }
  sommer::pcg::Options unrecycledOptions;
  unrecycledOptions.tolerance = 1.0e-8;
  sommer::pcg::Result unrecycled =
    sommer::pcg::solveIndependent(difficultRhs, recycleA, identityM,
                                  unrecycledOptions);
  sommer::pcg::RecycleSpace measuredRecycle(8, 3);
  measuredRecycle.capture(Eigen::VectorXd::Unit(recycleDimension, 0));
  measuredRecycle.capture(Eigen::VectorXd::Unit(recycleDimension, 1));
  sommer::pcg::RecycleSpace::Selection measuredSelection =
    measuredRecycle.selectCandidates(recycleA, metric, 2);
  measuredRecycle.accept(measuredSelection);
  measuredRecycle.prepare(recycleA, identityM);
  sommer::pcg::Options recycledOptions = unrecycledOptions;
  recycledOptions.recycleSpace = &measuredRecycle;
  sommer::pcg::Result recycled =
    sommer::pcg::solveIndependent(difficultRhs, recycleA, identityM,
                                  recycledOptions);
  const std::size_t setupCalls =
    measuredRecycle.selectionOperatorApplications() +
    measuredRecycle.preparationOperatorApplications();
  const double baselineTrueRelativeResidual =
    (recycleCoefficient * unrecycled.solution - difficultRhs).norm() /
    difficultRhs.norm();
  const double recycledTrueRelativeResidual =
    (recycleCoefficient * recycled.solution - difficultRhs).norm() /
    difficultRhs.norm();
  assert((recycleCoefficient * recycled.solution - difficultRhs).norm() <
         unrecycledOptions.tolerance * difficultRhs.norm());
  assert(baselineTrueRelativeResidual < unrecycledOptions.tolerance);
  assert(recycledTrueRelativeResidual < unrecycledOptions.tolerance);
  std::cout << "PCG recycle benchmark: n=" << recycleDimension
            << ", rhs=1, rank=" << measuredRecycle.basis().cols()
            << ", baseline iterations=" << unrecycled.metrics.iterations
            << ", baseline operator applications="
            << unrecycled.metrics.operatorApplications
            << ", recycled iterations=" << recycled.metrics.iterations
            << ", recycled solve operator applications="
            << recycled.metrics.operatorApplications
            << ", selection operator applications="
            << measuredRecycle.selectionOperatorApplications()
            << ", preparation operator applications="
            << measuredRecycle.preparationOperatorApplications()
            << ", recycled total operator applications="
            << recycled.metrics.operatorApplications + setupCalls
            << ", baseline max true relative residual="
            << baselineTrueRelativeResidual
            << ", recycled max true relative residual="
            << recycledTrueRelativeResidual << std::endl;
  assert(recycled.metrics.operatorApplications + setupCalls <
         unrecycled.metrics.operatorApplications);

  bool invalidMetricRejected = false;
  try {
    Eigen::VectorXd invalidMetric = metric;
    invalidMetric(0) = 0.0;
    measuredRecycle.selectCandidates(recycleA, invalidMetric, 2);
  } catch(const std::invalid_argument &) {
    invalidMetricRejected = true;
  }
  assert(invalidMetricRejected);

  bool indefiniteOperatorRejected = false;
  try {
    sommer::pcg::RecycleSpace invalidRecycle(4, 2);
    invalidRecycle.capture(Eigen::VectorXd::Unit(recycleDimension, 0));
    const sommer::pcg::VectorOperator negativeA =
      [](const Eigen::VectorXd & value) -> Eigen::VectorXd { return -value; };
    invalidRecycle.selectCandidates(negativeA, metric, 1);
  } catch(const std::invalid_argument &) {
    indefiniteOperatorRejected = true;
  }
  assert(indefiniteOperatorRejected);
}
