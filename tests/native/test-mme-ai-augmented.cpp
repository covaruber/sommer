#include "../../src/mme_ai_augmented.h"
#include "../../src/mme_control.h"
#include <Eigen/Cholesky>
#include <cassert>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace {

double scaledError(const Eigen::MatrixXd & actual,
                   const Eigen::MatrixXd & expected) {
  assert(actual.rows() == expected.rows() && actual.cols() == expected.cols());
  double error = 0.0;
  for(Eigen::Index i = 0; i < actual.size(); ++i) {
    error = std::max(error, std::abs(actual.data()[i] - expected.data()[i]) /
      (1.0 + std::abs(expected.data()[i])));
  }
  return error;
}

}

int main() {
  using sommer::ai::absorbAugmentedBorder;
  sommer::control::AIAssembly parsedAssembly = sommer::control::AIAssembly::Sensitivity;
  assert(sommer::control::parseAIAssembly("sensitivity", parsedAssembly));
  assert(parsedAssembly == sommer::control::AIAssembly::Sensitivity);
  assert(sommer::control::parseAIAssembly("augmented-working", parsedAssembly));
  assert(parsedAssembly == sommer::control::AIAssembly::AugmentedWorking);
  assert(!sommer::control::parseAIAssembly("unknown", parsedAssembly));

  const Eigen::Index observations = 13;
  const Eigen::Index coefficients = 7;
  const Eigen::Index fixed = 3;
  const Eigen::Index random = coefficients - fixed;
  const Eigen::Index randomDirections = 3;
  const Eigen::Index residualDirections = 2;
  const Eigen::Index directions = randomDirections + residualDirections;
  Eigen::MatrixXd W(observations, coefficients);
  for(Eigen::Index row = 0; row < W.rows(); ++row) {
    for(Eigen::Index column = 0; column < W.cols(); ++column) {
      W(row, column) = std::sin(0.31 * (row + 1) * (column + 2)) +
        0.07 * (row - 2 * column);
    }
  }
  Eigen::MatrixXd relationship = Eigen::MatrixXd::Identity(random, random);
  for(Eigen::Index i = 0; i < random; ++i) {
    relationship(i, i) += 0.2 * (i + 1);
    if(i + 1 < random && i % 2 == 0) {
      relationship(i, i + 1) = relationship(i + 1, i) = -0.13;
    }
  }
  Eigen::MatrixXd J = Eigen::MatrixXd::Zero(coefficients, coefficients);
  J.bottomRightCorner(random, random) = relationship;
  Eigen::MatrixXd Ri = Eigen::MatrixXd::Identity(observations, observations);
  for(Eigen::Index i = 0; i < observations; ++i) {
    Ri(i, i) = 0.8 + 0.04 * i;
    if(i + 1 < observations) {
      Ri(i, i + 1) = Ri(i + 1, i) = 0.035;
    }
  }
  Eigen::MatrixXd C = W.transpose() * Ri * W + J;
  const Eigen::LDLT<Eigen::MatrixXd> factor(C);
  assert(factor.info() == Eigen::Success);

  Eigen::MatrixXd T(coefficients, randomDirections);
  for(Eigen::Index row = 0; row < coefficients; ++row) {
    for(Eigen::Index column = 0; column < randomDirections; ++column) {
      const bool inTerm = column < 2 ? row >= fixed && row < fixed + 2
        : row >= fixed + 2;
      T(row, column) = inTerm
        ? std::cos(0.21 * (row + 2) * (column + 1)) + 0.11 * (row - column)
        : 0.0;
    }
  }
  Eigen::MatrixXd H(observations, residualDirections);
  for(Eigen::Index row = 0; row < observations; ++row) {
    H(row, 0) = std::sin(0.17 * (row + 1));
    H(row, 1) = std::cos(0.23 * (row + 3)) - 0.2 * H(row, 0);
  }

  Eigen::MatrixXd F(observations, directions);
  F.leftCols(randomDirections) = W * T;
  F.rightCols(residualDirections) = H;
  const Eigen::MatrixXd P = Ri - Ri * W * factor.solve(W.transpose() * Ri);
  const Eigen::MatrixXd ordinary = F.transpose() * P * F;

  Eigen::MatrixXd border = Eigen::MatrixXd::Zero(coefficients, directions);
  border.leftCols(randomDirections) = -J * T;
  border.rightCols(residualDirections) = W.transpose() * Ri * H;
  Eigen::MatrixXd direct = Eigen::MatrixXd::Zero(directions, directions);
  direct.topLeftCorner(randomDirections, randomDirections) = T.transpose() * J * T;
  direct.bottomRightCorner(residualDirections, residualDirections) =
    H.transpose() * Ri * H;

  const auto solve = [&](const Eigen::MatrixXd & rhs) { return factor.solve(rhs); };
  const auto full = absorbAugmentedBorder(border, direct, solve);
  assert(!full.tiled && full.solveCalls == 1);
  assert(scaledError(full.information, ordinary) <= 1.0e-12);
  assert(std::abs(ordinary(0, 2)) > 1.0e-8);
  const auto randomOnly = absorbAugmentedBorder(
    border.leftCols(randomDirections),
    direct.topLeftCorner(randomDirections, randomDirections), solve);
  assert(scaledError(randomOnly.information,
    ordinary.topLeftCorner(randomDirections, randomDirections)) <= 1.0e-12);
  const auto residualOnly = absorbAugmentedBorder(
    border.rightCols(residualDirections),
    direct.bottomRightCorner(residualDirections, residualDirections), solve);
  assert(scaledError(residualOnly.information,
    ordinary.bottomRightCorner(residualDirections, residualDirections)) <= 1.0e-12);

  const Eigen::MatrixXd reducedMap = (Eigen::MatrixXd(5, 3) <<
    1.0, 0.0, 0.2,
    0.0, 1.0, -0.1,
    0.0, 0.0, 0.7,
    0.5, 0.2, 0.0,
    -0.3, 0.0, 1.0).finished();
  const Eigen::MatrixXd reducedBorder = border * reducedMap;
  const Eigen::MatrixXd reducedDirect = reducedMap.transpose() * direct * reducedMap;
  const auto reduced = absorbAugmentedBorder(reducedBorder, reducedDirect, solve);
  assert(scaledError(reduced.information,
    reducedMap.transpose() * ordinary * reducedMap) <= 1.0e-12);

  const std::size_t q = static_cast<std::size_t>(directions);
  const std::size_t minimumBudget =
    (2 * q * q + 8 * static_cast<std::size_t>(coefficients) + q) * sizeof(double);
  const auto tiled = absorbAugmentedBorder(
    border, direct, solve, minimumBudget);
  assert(tiled.tiled && tiled.panelColumns == 1);
  assert(tiled.solveCalls == static_cast<std::size_t>(directions));
  assert(scaledError(tiled.information, ordinary) <= 1.0e-12);

  Eigen::MatrixXd emptyBorder(coefficients, 0);
  Eigen::MatrixXd emptyGram(0, 0);
  const auto empty = absorbAugmentedBorder(emptyBorder, emptyGram, solve, 0);
  assert(empty.information.rows() == 0 && empty.solveCalls == 0);

  Eigen::MatrixXd nearBorder = Eigen::MatrixXd::Zero(coefficients, 2);
  nearBorder.col(0) = border.col(0);
  nearBorder.col(1) = border.col(0) * (1.0 + 1.0e-12);
  Eigen::MatrixXd nearDirect = Eigen::MatrixXd::Zero(2, 2);
  nearDirect(0, 0) = direct(0, 0);
  nearDirect(1, 1) = direct(0, 0) * (1.0 + 2.0e-12);
  const auto nearlyDependent = absorbAugmentedBorder(
    nearBorder, nearDirect, solve);
  assert(nearlyDependent.information.allFinite());

  Eigen::MatrixXd termLoad0(3, 2), termLoad1(3, 2);
  termLoad0 << 1.0, -0.4, 0.2, 1.3, -0.7, 0.6;
  termLoad1 << -0.2, 0.8, 1.1, -0.3, 0.5, 1.4;
  Eigen::MatrixXd relPrecision(3, 3);
  relPrecision << 1.8, 0.2, -0.1, 0.2, 1.5, 0.3, -0.1, 0.3, 1.7;
  Eigen::MatrixXd covPrecision(2, 2);
  covPrecision << 1.4, -0.25, -0.25, 1.9;
  const Eigen::MatrixXd prior = sommer::ai::priorMetricGram(
    {termLoad0, termLoad1}, relPrecision, covPrecision);
  const Eigen::MatrixXd action0 = relPrecision * termLoad0 * covPrecision;
  const Eigen::MatrixXd action1 = relPrecision * termLoad1 * covPrecision;
  assert(std::abs(prior(0, 1) -
    (termLoad0.array() * action1.array()).sum()) < 1.0e-13);
  assert(std::abs(prior(1, 0) -
    (termLoad1.array() * action0.array()).sum()) < 1.0e-13);

  bool invalid = false;
  try {
    Eigen::MatrixXd wrong(2, 2);
    absorbAugmentedBorder(border, wrong, solve);
  } catch(const std::invalid_argument &) {
    invalid = true;
  }
  assert(invalid);
  invalid = false;
  try {
    Eigen::MatrixXd nonfinite = direct;
    nonfinite(0, 0) = std::numeric_limits<double>::quiet_NaN();
    absorbAugmentedBorder(border, nonfinite, solve);
  } catch(const std::invalid_argument &) {
    invalid = true;
  }
  assert(invalid);
  invalid = false;
  try {
    sommer::ai::checkedProduct(std::numeric_limits<std::size_t>::max(), 2,
      "test overflow");
  } catch(const std::overflow_error &) {
    invalid = true;
  }
  assert(invalid);
  invalid = false;
  try {
    absorbAugmentedBorder(border, direct, solve, minimumBudget - 1);
  } catch(const std::length_error &) {
    invalid = true;
  }
  assert(invalid);
  return 0;
}
