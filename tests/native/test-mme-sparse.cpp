#include "../../src/mme_sparse.h"
#include "../../src/mme_recovery.h"
#include "../../src/mme_posterior.h"
#include <Eigen/SparseCholesky>
#include <Eigen/LU>
#include <cassert>
#include <cmath>
#include <string>

using sommer::SparseFactorizationStatus;
typedef Eigen::SparseMatrix<double> SparseMatrix;

struct ScriptedFactor {
  Eigen::ComputationInfo status = Eigen::Success;
  int analyses = 0;
  int failAnalysis = 0;
  void analyzePattern(const SparseMatrix &) {
    ++analyses;
    status = analyses == failAnalysis ? Eigen::InvalidInput : Eigen::Success;
  }
  Eigen::ComputationInfo info() const { return status; }
};

int main() {
  Eigen::MatrixXd coefficient(3, 3);
  coefficient << 4, 1, 0.5, 1, 3, 0.2, 0.5, 0.2, 2;
  const Eigen::MatrixXd jointInverse = coefficient.inverse();
  const Eigen::MatrixXd randomInverse = coefficient.bottomRightCorner(2, 2).inverse();
  Eigen::MatrixXd rhs(3, 2);
  rhs << 1, 2, 3, 4, 5, 6;
  int jointCalls = 0, randomCalls = 0;
  const auto jointSolve = [&](const Eigen::Ref<const Eigen::MatrixXd> & values) -> Eigen::MatrixXd {
    ++jointCalls;
    return jointInverse * values;
  };
  const auto randomSolve = [&](const Eigen::Ref<const Eigen::MatrixXd> & values) -> Eigen::MatrixXd {
    ++randomCalls;
    assert(values.outerStride() == values.rows());
    return randomInverse * values;
  };
  const auto remlMoment = sommer::posteriorCoefficientSolve(rhs, 1, true, jointSolve, randomSolve);
  assert((remlMoment - jointInverse * rhs).norm() < 1e-12);
  const auto mlMoment = sommer::posteriorCoefficientSolve(rhs, 1, false, jointSolve, randomSolve);
  assert(mlMoment.topRows(1).norm() == 0);
  assert((mlMoment.bottomRows(2) - randomInverse * rhs.bottomRows(2)).norm() < 1e-12);
  assert((mlMoment - remlMoment).norm() > 1e-3);
  const auto noRandom = sommer::posteriorCoefficientSolve(rhs, 3, false, jointSolve, randomSolve);
  assert(noRandom.norm() == 0 && jointCalls == 1 && randomCalls == 1);
  bool invalidDimension = false;
  try { sommer::posteriorCoefficientSolve(rhs, 4, false, jointSolve, randomSolve); }
  catch(const std::runtime_error &) { invalidDimension = true; }
  assert(invalidDimension);
  sommer::AIRecoveryPolicy policy;
  assert(policy.observe(false, 99, 0.0, 1e-4, false) == nullptr);
  assert(policy.observe(true, 3, 1.0, 1e-4, false) == nullptr);
  assert(policy.observe(true, 4, 1.0, 1e-4, false) == nullptr);
  assert(policy.observe(true, 4, 1.0, 1e-4, false) == nullptr);
  assert(std::string(policy.observe(true, 4, 1.0, 1e-4, false)) == "repeated-backtracking");
  policy.reset();
  assert(policy.observe(true, 0, 0.0, 1e-4, false) == nullptr);
  assert(policy.observe(true, 0, 0.0, 1e-4, false) == nullptr);
  assert(policy.observe(true, 0, 1e-4, 1e-4, false) == nullptr);
  assert(policy.negligibleGains == 0);
  assert(policy.observe(true, 0, 0.0, 1e-4, false) == nullptr);
  assert(policy.observe(true, 0, 0.0, 1e-4, false) == nullptr);
  assert(std::string(policy.observe(true, 0, 0.0, 1e-4, false)) == "negligible-gains");
  assert(policy.observe(true, 99, 0.0, 1e-4, true) == nullptr);
  assert(policy.negligibleGains == 0 && policy.difficultSteps == 0);
  SparseMatrix matrix(4, 4);
  matrix.setIdentity();
  matrix.coeffRef(1, 0) = 0.2;
  matrix.coeffRef(0, 1) = 0.2;
  matrix.makeCompressed();
  std::vector<int> outer, inner;
  assert(!sommer::eigenSparsePatternMatches(matrix, outer, inner));
  sommer::cacheEigenSparsePattern(matrix, outer, inner);
  assert(sommer::eigenSparsePatternMatches(matrix, outer, inner));

  Eigen::SimplicialLDLT<SparseMatrix> factor;
  bool ready = false;
  int analyses = 0, numericCalls = 0;
  bool topologyReady = true;
  const auto factorize = [&](const SparseMatrix & value) {
    return sommer::factorizeLDLTWithCachedPattern(
      value, factor, ready, outer, inner,
      [&](){ ++analyses; topologyReady = false; },
      [&](){ ++numericCalls; factor.factorize(value); });
  };
  auto result = factorize(matrix);
  assert(result.status == SparseFactorizationStatus::Success && !result.samePattern);
  assert(ready && analyses == 1 && numericCalls == 1 && !topologyReady);
  const Eigen::MatrixXd identity = Eigen::MatrixXd::Identity(4, 4);
  Eigen::MatrixXd inverse = factor.solve(identity);
  assert((Eigen::MatrixXd(matrix) * inverse - identity).norm() < 1e-12);
  assert(std::abs(factor.vectorD().array().log().sum() - std::log(0.96)) < 1e-12);

  matrix.coeffRef(0, 0) = 2.0;
  topologyReady = true;
  result = factorize(matrix);
  assert(result.status == SparseFactorizationStatus::Success && result.samePattern);
  assert(analyses == 1 && numericCalls == 2 && topologyReady);
  inverse = factor.solve(identity);
  assert((Eigen::MatrixXd(matrix) * inverse - identity).norm() < 1e-12);

  SparseMatrix changed(4, 4);
  changed.setIdentity();
  changed.coeffRef(2, 0) = 0.1;
  changed.coeffRef(0, 2) = 0.1;
  changed.makeCompressed();
  assert(changed.nonZeros() == matrix.nonZeros());
  assert(!sommer::eigenSparsePatternMatches(changed, outer, inner));
  result = factorize(changed);
  assert(result.status == SparseFactorizationStatus::Success && !result.samePattern);
  assert(analyses == 2 && numericCalls == 3 && !topologyReady);

  SparseMatrix diagonal(4, 4);
  diagonal.setIdentity();
  result = factorize(diagonal);
  assert(result.status == SparseFactorizationStatus::Success && !result.samePattern);
  assert(analyses == 3 && numericCalls == 4);
  ready = false;
  assert(sommer::factorizeLDLTWithCachedPattern(diagonal, factor, ready, outer, inner));

  for(int scenario = 0; scenario < 6; ++scenario) {
    ScriptedFactor scripted;
    bool symbolicReady = scenario >= 2;
    std::vector<int> cachedOuter, cachedInner;
    if(symbolicReady) { sommer::cacheEigenSparsePattern(diagonal, cachedOuter, cachedInner); }
    if(scenario == 0 || scenario == 3) { scripted.failAnalysis = 1; }
    int callbacks = 0, factorizations = 0;
    const auto outcome = sommer::factorizeLDLTWithCachedPattern(
      diagonal, scripted, symbolicReady, cachedOuter, cachedInner,
      [&](){ ++callbacks; },
      [&](){
        ++factorizations;
        scripted.status = scenario == 4 || (factorizations == 1 && scenario != 5)
          ? Eigen::NumericalIssue : Eigen::Success;
      });
    assert(callbacks == scripted.analyses);
    if(scenario == 0) {
      assert(outcome.status == SparseFactorizationStatus::SymbolicFailure);
      assert(!symbolicReady && cachedOuter.empty() && factorizations == 0);
    }else if(scenario == 1) {
      assert(outcome.status == SparseFactorizationStatus::NumericFailure);
      assert(!outcome.samePattern && callbacks == 1 && factorizations == 1);
    }else if(scenario == 2) {
      assert(outcome.status == SparseFactorizationStatus::Success);
      assert(outcome.samePattern && callbacks == 1 && factorizations == 2);
    }else if(scenario == 3) {
      assert(outcome.status == SparseFactorizationStatus::NumericFailure);
      assert(callbacks == 1 && factorizations == 1);
    }else if(scenario == 4) {
      assert(outcome.status == SparseFactorizationStatus::NumericFailure);
      assert(callbacks == 1 && factorizations == 2);
    }else {
      assert(outcome.status == SparseFactorizationStatus::Success);
      assert(outcome.samePattern && callbacks == 0 && factorizations == 1);
    }
  }
}
