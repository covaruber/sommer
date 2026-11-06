#ifndef SOMMER_MME_AI_AUGMENTED_H
#define SOMMER_MME_AI_AUGMENTED_H

#include <Eigen/Core>
#include <algorithm>
#include <cstddef>
#include <functional>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace sommer {
namespace ai {

struct AugmentedAssemblyResult {
  Eigen::MatrixXd information;
  std::size_t directions = 0;
  std::size_t solveCalls = 0;
  std::size_t panelColumns = 0;
  std::size_t workspaceBytes = 0;
  bool tiled = false;
};

struct AugmentedWorkspacePlan {
  std::size_t directions = 0;
  std::size_t panelColumns = 0;
  std::size_t workspaceBytes = 0;
  bool tiled = false;
};

inline std::size_t checkedProduct(std::size_t left, std::size_t right,
                                  const char * context) {
  if(right != 0 && left > std::numeric_limits<std::size_t>::max() / right) {
    throw std::overflow_error(std::string("Workspace size overflow in ") + context + ".");
  }
  return left * right;
}

inline std::size_t checkedSum(std::size_t left, std::size_t right,
                              const char * context) {
  if(left > std::numeric_limits<std::size_t>::max() - right) {
    throw std::overflow_error(std::string("Workspace size overflow in ") + context + ".");
  }
  return left + right;
}

inline AugmentedWorkspacePlan planAugmentedWorkspace(
    std::size_t coefficientCount, std::size_t directionCount,
    std::size_t residentBytes,
    std::size_t maxWorkspaceBytes = std::numeric_limits<std::size_t>::max()) {
  const std::size_t qSquared = checkedProduct(
    directionCount, directionCount, "augmented information workspace");
  const std::size_t resultAndSymmetryBytes = checkedProduct(
    checkedProduct(qSquared, 2, "augmented result and symmetry scratch"),
    sizeof(double), "augmented result and symmetry scratch");
  const std::size_t panelDoubles = checkedSum(
    checkedProduct(coefficientCount, 8, "augmented backend solve panels"),
    directionCount, "augmented solve panels");
  const std::size_t panelBytes = checkedProduct(
    panelDoubles, sizeof(double), "augmented solve panels");
  const std::size_t fixedBytes = checkedSum(
    residentBytes, resultAndSymmetryBytes, "augmented AI resident workspace");
  const std::size_t minimumBytes = checkedSum(
    fixedBytes, panelBytes, "augmented AI minimum workspace");
  if(directionCount != 0 && minimumBytes > maxWorkspaceBytes) {
    throw std::length_error(
      "The memory budget cannot admit one exact augmented AI solve column.");
  }
  AugmentedWorkspacePlan plan;
  plan.directions = directionCount;
  if(directionCount == 0) {
    plan.workspaceBytes = residentBytes;
    return plan;
  }
  std::size_t panelWidth = directionCount;
  if(directionCount > (maxWorkspaceBytes - fixedBytes) / panelBytes) {
    plan.tiled = true;
    panelWidth = std::min(directionCount,
      (maxWorkspaceBytes - fixedBytes) / panelBytes);
  }
  plan.panelColumns = panelWidth;
  plan.workspaceBytes = checkedSum(
    fixedBytes, checkedProduct(panelBytes, panelWidth,
      "augmented solve panel plan"), "augmented solve workspace");
  return plan;
}

inline Eigen::MatrixXd priorMetricGram(
    const std::vector<Eigen::MatrixXd> & loadings,
    const Eigen::MatrixXd & relationshipPrecision,
    const Eigen::MatrixXd & covariancePrecision) {
  const Eigen::Index levels = relationshipPrecision.rows();
  const Eigen::Index traits = covariancePrecision.rows();
  if(relationshipPrecision.cols() != levels || covariancePrecision.cols() != traits) {
    throw std::invalid_argument("Prior precision factors must be square.");
  }
  if(!relationshipPrecision.allFinite() || !covariancePrecision.allFinite()) {
    throw std::invalid_argument("Prior precision factors must be finite.");
  }
  for(const Eigen::MatrixXd & loading : loadings) {
    if(loading.rows() != levels || loading.cols() != traits) {
      throw std::invalid_argument("Random covariance loading has incompatible dimensions.");
    }
    if(!loading.allFinite()) {
      throw std::invalid_argument("Random covariance loading must be finite.");
    }
  }
  Eigen::MatrixXd result = Eigen::MatrixXd::Zero(
    static_cast<Eigen::Index>(loadings.size()),
    static_cast<Eigen::Index>(loadings.size()));
  std::vector<Eigen::MatrixXd> precisionActions;
  precisionActions.reserve(loadings.size());
  for(const Eigen::MatrixXd & loading : loadings) {
    precisionActions.emplace_back(
      relationshipPrecision * loading * covariancePrecision);
  }
  for(std::size_t column = 0; column < loadings.size(); ++column) {
    for(std::size_t row = 0; row <= column; ++row) {
      const double value = (loadings[row].array() *
        precisionActions[column].array()).sum();
      result(static_cast<Eigen::Index>(row), static_cast<Eigen::Index>(column)) = value;
      result(static_cast<Eigen::Index>(column), static_cast<Eigen::Index>(row)) = value;
    }
  }
  if(!result.allFinite()) {
    throw std::runtime_error("Non-finite random prior Gram in augmented AI assembly.");
  }
  return result;
}

template<class SolvePanel>
AugmentedAssemblyResult absorbAugmentedBorder(
    const Eigen::Ref<const Eigen::MatrixXd> & border,
    const Eigen::Ref<const Eigen::MatrixXd> & directGram,
    SolvePanel solvePanel,
    const AugmentedWorkspacePlan & plan) {
  if(border.cols() != directGram.rows() ||
     directGram.rows() != directGram.cols()) {
    throw std::invalid_argument("Augmented border and direct Gram dimensions do not agree.");
  }
  if(!border.allFinite() || !directGram.allFinite()) {
    throw std::invalid_argument("Augmented border and direct Gram must be finite.");
  }
  const std::size_t coefficientCount = static_cast<std::size_t>(border.rows());
  const std::size_t directionCount = static_cast<std::size_t>(border.cols());
  if(plan.directions != directionCount ||
     (directionCount > 0 &&
      (plan.panelColumns == 0 || plan.panelColumns > directionCount))) {
    throw std::invalid_argument("Augmented AI workspace plan does not match the border.");
  }

  AugmentedAssemblyResult result;
  result.directions = directionCount;
  result.panelColumns = plan.panelColumns;
  result.workspaceBytes = plan.workspaceBytes;
  result.tiled = plan.tiled;
  result.information = directGram;
  if(directionCount == 0) { return result; }

  for(std::size_t first = 0; first < directionCount; first += plan.panelColumns) {
    const std::size_t count = std::min(plan.panelColumns, directionCount - first);
    const Eigen::MatrixXd rhs = border.middleCols(
      static_cast<Eigen::Index>(first), static_cast<Eigen::Index>(count));
    const Eigen::MatrixXd solved = solvePanel(rhs);
    if(solved.rows() != border.rows() ||
       solved.cols() != static_cast<Eigen::Index>(count) ||
       !solved.allFinite()) {
      throw std::runtime_error("Coefficient solve returned an invalid augmented AI panel.");
    }
    const Eigen::MatrixXd cross = border.transpose() * solved;
    result.information.middleCols(
      static_cast<Eigen::Index>(first), static_cast<Eigen::Index>(count)) -= cross;
    ++result.solveCalls;
  }
  result.information =
    0.5 * (result.information + result.information.transpose()).eval();
  if(!result.information.allFinite()) {
    throw std::runtime_error("Non-finite augmented working-variate AI information.");
  }
  return result;
}

template<class SolvePanel>
AugmentedAssemblyResult absorbAugmentedBorder(
    const Eigen::Ref<const Eigen::MatrixXd> & border,
    const Eigen::Ref<const Eigen::MatrixXd> & directGram,
    SolvePanel solvePanel,
    std::size_t maxWorkspaceBytes = std::numeric_limits<std::size_t>::max(),
    std::size_t residentBytes = 0) {
  const AugmentedWorkspacePlan plan = planAugmentedWorkspace(
    static_cast<std::size_t>(border.rows()),
    static_cast<std::size_t>(border.cols()), residentBytes, maxWorkspaceBytes);
  return absorbAugmentedBorder(border, directGram, solvePanel, plan);
}

} // namespace ai
} // namespace sommer

#endif
