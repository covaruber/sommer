#include "../../src/mme_saddle.h"
#include "../../src/mme_exact.h"
#include <Eigen/Cholesky>
#include <Eigen/LU>
#include <cassert>
#include <cmath>
#include <functional>
#include <limits>
#include <stdexcept>

namespace {

using sommer::SaddleSystem;
using Sparse = SaddleSystem::Sparse;
using Dense = Eigen::MatrixXd;

Sparse sparse(const Dense & matrix) {
  return Sparse(matrix.sparseView());
}

void expectFailure(const std::function<void()> & operation) {
  bool failed = false;
  try {
    operation();
  } catch(const std::runtime_error &) {
    failed = true;
  }
  assert(failed);
}

Dense makeK(const Dense & r, const Dense & a, const Dense & q) {
  Dense k = Dense::Zero(r.rows() + q.rows(), r.cols() + q.cols());
  k.topLeftCorner(r.rows(), r.cols()) = r;
  k.topRightCorner(r.rows(), a.cols()) = a;
  k.bottomLeftCorner(a.cols(), r.cols()) = a.transpose();
  k.bottomRightCorner(q.rows(), q.cols()) = -q;
  return k;
}

}  // namespace

int main() {
  {
    Dense inverse(3, 3);
    inverse << 2.0, 0.2, -0.1, 0.2, 1.5, 0.3, -0.1, 0.3, 1.0;
    Dense design(6, 3);
    design << 1, 0, 2, 0, 0, 0, 0, -1, 0,
              1, 1, 0, 0, 2, 1, 1, 0, -1;
    for(int count : {2, 6}) for(int availability : {0, 1, 2}) {
      const Sparse a = sparse(design.topRows(count));
      int solvedColumns = 0;
      const Eigen::VectorXd actual = sommer::exactInverseRowQuadratics(a, 2,
        [&](const Dense & rhs) -> Dense {
          solvedColumns += rhs.cols();
          return inverse * rhs;
        },
        [&](Eigen::Index row, Eigen::Index column, double & value) {
          value = inverse(row, column);
          return availability == 2 || (availability == 1 && row == column);
        });
      const Dense expected = Dense(a) * inverse * Dense(a.transpose());
      assert((actual - expected.diagonal()).norm() < 1e-12);
      assert(solvedColumns <= std::min(count, 3));
      if(availability == 2) { assert(solvedColumns == 0); }
    }
    const Sparse empty(4, 0);
    assert(sommer::exactInverseRowQuadratics(empty, 2, nullptr, nullptr).isZero());
    expectFailure([&]() {
      sommer::exactInverseRowQuadratics(sparse(design), 0, nullptr, nullptr);
    });
    expectFailure([&]() {
      sommer::exactInverseRowQuadratics(sparse(design), 2,
        [](const Dense &) -> Dense { return Dense::Zero(1, 1); }, nullptr);
    });
  }
  Dense r(3, 3);
  r << 4.0, 0.5, 0.0,
       0.5, 3.0, 0.25,
       0.0, 0.25, 2.0;
  Dense a(3, 2);
  a << 1.0, 0.0,
       0.5, 1.0,
       0.0, 2.0;
  Dense q(2, 2);
  q << 2.0, 0.25,
       0.25, 1.5;
  const Dense denseK = makeK(r, a, q);
  const Dense c = a.transpose() * r.ldlt().solve(a) + q;
  const double denseDeterminant = denseK.determinant();
  assert(denseDeterminant > 0.0);

  SaddleSystem system;
  system.factorize(sparse(r), sparse(a), sparse(q));
  assert(system.determinantSign() == 1);
  assert(std::abs(system.logAbsDeterminant() -
                  std::log(std::abs(denseDeterminant))) < 1e-10);
  assert(system.symbolicAnalyses() == 1);
  assert(system.numericFactorizations() == 1);
  assert(system.nonzeros() == static_cast<std::uint64_t>(sparse(denseK).nonZeros()));
  assert(system.factorNonzeros() > 0);

  Dense rhs(2, 35);
  for(Eigen::Index row = 0; row < rhs.rows(); ++row) {
    for(Eigen::Index col = 0; col < rhs.cols(); ++col) {
      rhs(row, col) = 0.25 + row + 0.1 * col;
    }
  }
  const Dense coefficients = system.solveCoefficients(rhs);
  const Dense expectedCoefficients = c.ldlt().solve(rhs);
  assert((coefficients - expectedCoefficients).norm() < 1e-10);
  assert(system.solveColumns() == 35);
  assert(system.maxRelativeResidual() <= 1e-10);

  Eigen::VectorXd y(3);
  y << 1.0, -2.0, 0.5;
  const Eigen::VectorXd responseSolution = system.solveResponse(y);
  const Eigen::VectorXd expectedResponse =
    c.ldlt().solve(a.transpose() * r.ldlt().solve(y));
  assert((responseSolution - expectedResponse).norm() < 1e-10);
  assert(system.solveColumns() == 36);

  Dense observationRhs(3, 35);
  for(Eigen::Index row = 0; row < observationRhs.rows(); ++row) {
    for(Eigen::Index col = 0; col < observationRhs.cols(); ++col) {
      observationRhs(row, col) = 0.5 + row - 0.05 * col;
    }
  }
  Dense augmentedRhs = Dense::Zero(5, observationRhs.cols());
  augmentedRhs.topRows(3) = observationRhs;
  const Dense expectedResidual =
    denseK.fullPivLu().solve(augmentedRhs).topRows(3);
  const Dense residualSolution = system.solveResidual(observationRhs);
  assert((residualSolution - expectedResidual).norm() < 1e-10);
  assert(system.solveColumns() == 71);
  const Sparse assembled = sommer::assembleNormalEquations(
    sparse(a), sparse(q), [&](const Dense & values) -> Dense { return r.ldlt().solve(values); },
    1, 100);
  assert((Dense(assembled) - c).norm() < 1e-10);
  const Eigen::VectorXd traces = sommer::exactInverseTraces(2, 2, 1,
    [&](const Dense & values) -> Dense { return system.solveCoefficients(values); },
    [&](Eigen::Index parameter, const Dense & values) -> Dense {
      return parameter == 0 ? Dense(q * values) : values;
    });
  const Dense inverse = c.inverse();
  assert(std::abs(traces(0) - (q * inverse).trace()) < 1e-10);
  assert(std::abs(traces(1) - inverse.trace()) < 1e-10);
  const Eigen::MatrixXd blocks = sommer::exactInverseBlockTraces(2, {0, 1},
    sparse(Dense::Constant(1, 1, 2.0)), 1,
    [&](const Dense & values) -> Dense { return system.solveCoefficients(values); },
    [](std::size_t, std::size_t) { return true; });
  assert((blocks - 2.0 * inverse).norm() < 1e-10);
  int cachedSolves = 0;
  for(bool partial : {false, true}) {
    const Eigen::MatrixXd cached = sommer::exactInverseBlockTraces(2, {0, 1},
      sparse(Dense::Constant(1, 1, 2.0)), 1,
      [&](const Dense & values) -> Dense {
        ++cachedSolves;
        return system.solveCoefficients(values);
      },
      [](std::size_t, std::size_t) { return true; },
      [&](Eigen::Index row, Eigen::Index column, double & value) {
        if(partial && row != column) { return false; }
        value = inverse(row, column);
        return true;
      });
    assert((cached - 2.0 * inverse).norm() < 1e-10);
    assert(cachedSolves == (partial ? 1 : 0));
  }
  expectFailure([&]() {
    sommer::exactInverseBlockTraces(2, {0}, Sparse(0, 0), 1,
      [&](const Dense & values) -> Dense { return values; },
      [](std::size_t, std::size_t) { return true; });
  });
  expectFailure([&]() {
    sommer::assembleNormalEquations(sparse(a), sparse(q),
      [&](const Dense & values) -> Dense { return values; }, 1, 0);
  });

  Dense rSamePattern = r;
  rSamePattern(0, 0) += 0.75;
  system.factorize(sparse(rSamePattern), sparse(a), sparse(q));
  assert(system.symbolicAnalyses() == 1);
  assert(system.numericFactorizations() == 2);

  Dense aNewPattern = a;
  aNewPattern(0, 1) = 0.125;
  system.factorize(sparse(rSamePattern), sparse(aNewPattern), sparse(q));
  assert(system.symbolicAnalyses() == 2);
  assert(system.numericFactorizations() == 3);

  Dense emptyDesign(3, 0);
  Dense emptyPrior(0, 0);
  system.factorize(sparse(r), sparse(emptyDesign), sparse(emptyPrior));
  assert(system.determinantSign() == 1);
  assert(system.solveCoefficients(Dense(0, 2)).rows() == 0);
  assert(system.solveResponse(y).size() == 0);

  Dense oneCoefficientDesign(3, 1);
  oneCoefficientDesign << 1.0, 0.0, 0.5;
  Dense oneCoefficientPrior(1, 1);
  oneCoefficientPrior << 1.25;
  SaddleSystem oddCoefficientSystem;
  oddCoefficientSystem.factorize(
    sparse(r), sparse(oneCoefficientDesign), sparse(oneCoefficientPrior));
  assert(oddCoefficientSystem.determinantSign() == -1);
  assert(std::abs(oddCoefficientSystem.logAbsDeterminant() -
    std::log(std::abs(makeK(r, oneCoefficientDesign,
                            oneCoefficientPrior).determinant()))) < 1e-10);

  expectFailure([&]() {
    system.factorize(sparse(r), sparse(a), sparse(Dense::Identity(1, 1)));
  });
  Dense nonFinite = r;
  nonFinite(0, 0) = std::numeric_limits<double>::infinity();
  expectFailure([&]() {
    system.factorize(sparse(nonFinite), sparse(a), sparse(q));
  });
  Dense asymmetric = r;
  asymmetric(0, 1) += 0.5;
  expectFailure([&]() {
    system.factorize(sparse(asymmetric), sparse(a), sparse(q));
  });
  const Dense fixedPrior = Dense::Zero(2, 2);
  SaddleSystem fixedSystem;
  fixedSystem.factorize(sparse(r), sparse(a), sparse(fixedPrior));
  const Dense fixedC = a.transpose() * r.ldlt().solve(a);
  assert((fixedSystem.solveCoefficients(rhs) - fixedC.ldlt().solve(rhs)).norm() < 1e-10);
  const Dense deficient = Dense::Ones(3, 2);
  expectFailure([&]() {
    fixedSystem.factorize(sparse(r), sparse(deficient), sparse(fixedPrior));
    fixedSystem.solveCoefficients(Dense::Identity(2, 2));
  });
  expectFailure([&]() {
    system.factorize(sparse(Dense::Zero(2, 2)),
                     sparse(Dense::Zero(2, 1)),
                     sparse(Dense::Zero(1, 1)));
  });
  expectFailure([&]() {
    system.solveCoefficients(Dense(1, 1));
  });

  return 0;
}
