#ifndef SOMMER_MME_SADDLE_H
#define SOMMER_MME_SADDLE_H

#include <Eigen/SparseLU>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace sommer {

class SaddleSystem {
public:
  using Sparse = Eigen::SparseMatrix<double>;
  using Ordering = Eigen::COLAMDOrdering<int>;
  using Factor = Eigen::SparseLU<Sparse, Ordering>;

  SaddleSystem()
    : analyzed_(false), factorized_(false), emptySystem_(false),
      coefficientCount_(0), logAbsDeterminant_(0.0), determinantSign_(1),
      symbolicAnalyses_(0), numericFactorizations_(0), solveColumns_(0),
      nonzeros_(0), factorNonzeros_(0), maxRelativeResidual_(0.0) {}

  void factorize(const Sparse & residualCovariance,
                 const Sparse & design,
                 const Sparse & priorPrecision,
                 const std::string & context = "SaddleSystem::factorize") {
    factorized_ = false;
    const Eigen::Index n = residualCovariance.rows();
    const Eigen::Index p = design.cols();
    if(residualCovariance.cols() != n) {
      fail(context, "residual covariance R must be square");
    }
    if(design.rows() != n) {
      fail(context, "design A row count must match R");
    }
    if(priorPrecision.rows() != p || priorPrecision.cols() != p) {
      fail(context, "prior precision Q dimensions must match A columns");
    }

    const Eigen::Index indexLimit =
      static_cast<Eigen::Index>((std::numeric_limits<int>::max)());
    if(n > indexLimit || p > indexLimit || n > indexLimit - p) {
      fail(context, "saddle system dimension exceeds SparseLU index range");
    }
    validateFinite(residualCovariance, context, "R");
    validateFinite(design, context, "A");
    validateFinite(priorPrecision, context, "Q");
    if((n > 0 && !residualCovariance.isApprox(Sparse(residualCovariance.transpose()), 1e-12)) ||
       (p > 0 && !priorPrecision.isApprox(Sparse(priorPrecision.transpose()), 1e-12))) {
      fail(context, "R and Q must be symmetric");
    }

    const Eigen::Index rNnz = residualCovariance.nonZeros();
    const Eigen::Index aNnz = design.nonZeros();
    const Eigen::Index qNnz = priorPrecision.nonZeros();
    if(rNnz < 0 || aNnz < 0 || qNnz < 0 ||
       aNnz > indexLimit / 2 ||
       rNnz > indexLimit - 2 * aNnz ||
       qNnz > indexLimit - rNnz - 2 * aNnz) {
      fail(context, "saddle system nonzero count exceeds SparseLU index range");
    }
    const Eigen::Index total = n + p;
    std::vector<Eigen::Triplet<double, int>> triplets;
    triplets.reserve(static_cast<std::size_t>(rNnz + 2 * aNnz + qNnz));

    for(Eigen::Index outer = 0; outer < residualCovariance.outerSize(); ++outer) {
      for(Sparse::InnerIterator it(residualCovariance, outer); it; ++it) {
        triplets.emplace_back(static_cast<int>(it.row()),
                              static_cast<int>(it.col()), it.value());
      }
    }
    for(Eigen::Index outer = 0; outer < design.outerSize(); ++outer) {
      for(Sparse::InnerIterator it(design, outer); it; ++it) {
        const int row = static_cast<int>(it.row());
        const int col = static_cast<int>(it.col());
        triplets.emplace_back(row, static_cast<int>(n) + col, it.value());
        triplets.emplace_back(static_cast<int>(n) + col, row, it.value());
      }
    }
    for(Eigen::Index outer = 0; outer < priorPrecision.outerSize(); ++outer) {
      for(Sparse::InnerIterator it(priorPrecision, outer); it; ++it) {
        triplets.emplace_back(static_cast<int>(n) + static_cast<int>(it.row()),
                              static_cast<int>(n) + static_cast<int>(it.col()),
                              -it.value());
      }
    }

    Sparse assembled(static_cast<int>(total), static_cast<int>(total));
    assembled.setFromTriplets(triplets.begin(), triplets.end());
    assembled.makeCompressed();
    if(assembled.nonZeros() > indexLimit) {
      fail(context, "assembled saddle system exceeds SparseLU index range");
    }

    coefficientCount_ = p;
    analyzedRows_ = total;
    nonzeros_ = static_cast<std::uint64_t>(assembled.nonZeros());
    assembledMatrix_ = assembled;
    factorized_ = false;
    emptySystem_ = total == 0;
    if(emptySystem_) {
      analyzed_ = true;
      logAbsDeterminant_ = 0.0;
      determinantSign_ = 1;
      factorNonzeros_ = 0;
      factorized_ = true;
      return;
    }

    if(!hasSamePattern(assembled)) {
      factor_.analyzePattern(assembled);
      ++symbolicAnalyses_;
      cachePattern(assembled);
      analyzed_ = true;
    }
    if(!analyzed_) {
      fail(context, "SparseLU symbolic analysis is unavailable");
    }
    factor_.factorize(assembledMatrix_);
    ++numericFactorizations_;
    if(factor_.info() != Eigen::Success) {
      fail(context, "SparseLU numeric factorization failed: " +
           factor_.lastErrorMessage());
    }

    const int expectedSign = (p % 2 == 0) ? 1 : -1;
    const double sign = factor_.signDeterminant();
    const double logDeterminant = factor_.logAbsDeterminant();
    if(!std::isfinite(sign) || sign != static_cast<double>(expectedSign) ||
       !std::isfinite(logDeterminant)) {
      fail(context, "saddle determinant sign or log-absolute-determinant is invalid");
    }
    determinantSign_ = expectedSign;
    logAbsDeterminant_ = logDeterminant;

    factorNonzeros_ = static_cast<std::uint64_t>(factor_.nnzL()) +
                      static_cast<std::uint64_t>(factor_.nnzU());
    factorized_ = true;
  }

  Eigen::MatrixXd solveCoefficients(
      const Eigen::MatrixXd & rhs,
      const std::string & context = "SaddleSystem::solveCoefficients") {
    requireFactorized(context);
    if(rhs.rows() != coefficientCount_) {
      fail(context, "right-hand side row count must match coefficient count");
    }
    if(!rhs.allFinite()) {
      fail(context, "right-hand side contains non-finite values");
    }
    Eigen::MatrixXd result(coefficientCount_, rhs.cols());
    if(coefficientCount_ == 0 || rhs.cols() == 0) {
      return result;
    }
    for(Eigen::Index first = 0; first < rhs.cols(); first += panelColumns()) {
      const Eigen::Index columns =
        (std::min)(panelColumns(), rhs.cols() - first);
      Eigen::MatrixXd saddleRhs =
        Eigen::MatrixXd::Zero(coefficientCount_ + residualCount(),
                              columns);
      saddleRhs.bottomRows(coefficientCount_).noalias() =
        -rhs.middleCols(first, columns);
      const Eigen::MatrixXd solved =
        solvePanel(saddleRhs, context, columns);
      result.middleCols(first, columns) =
        solved.bottomRows(coefficientCount_);
    }
    return result;
  }

  Eigen::VectorXd solveResponse(
      const Eigen::VectorXd & response,
      const std::string & context = "SaddleSystem::solveResponse") {
    requireFactorized(context);
    if(response.size() != residualCount()) {
      fail(context, "response length must match residual dimension");
    }
    if(!response.allFinite()) {
      fail(context, "response contains non-finite values");
    }
    if(coefficientCount_ == 0) {
      return Eigen::VectorXd(0);
    }
    Eigen::MatrixXd saddleRhs =
      Eigen::MatrixXd::Zero(coefficientCount_ + residualCount(), 1);
    saddleRhs.topRows(residualCount()) = response;
    const Eigen::MatrixXd solved = solvePanel(saddleRhs, context, 1);
    return solved.bottomRows(coefficientCount_).col(0);
  }

  Eigen::MatrixXd solveResidual(
      const Eigen::MatrixXd & rhs,
      const std::string & context = "SaddleSystem::solveResidual") {
    requireFactorized(context);
    if(rhs.rows() != residualCount()) {
      fail(context, "right-hand side row count must match residual dimension");
    }
    if(!rhs.allFinite()) {
      fail(context, "right-hand side contains non-finite values");
    }
    Eigen::MatrixXd result(residualCount(), rhs.cols());
    if(rhs.cols() == 0 || residualCount() == 0) {
      return result;
    }
    for(Eigen::Index first = 0; first < rhs.cols(); first += panelColumns()) {
      const Eigen::Index columns =
        (std::min)(panelColumns(), rhs.cols() - first);
      Eigen::MatrixXd saddleRhs =
        Eigen::MatrixXd::Zero(coefficientCount_ + residualCount(), columns);
      saddleRhs.topRows(residualCount()) = rhs.middleCols(first, columns);
      const Eigen::MatrixXd solved =
        solvePanel(saddleRhs, context, columns);
      result.middleCols(first, columns) = solved.topRows(residualCount());
    }
    return result;
  }

  double logAbsDeterminant() const {
    requireFactorized("SaddleSystem::logAbsDeterminant");
    return logAbsDeterminant_;
  }

  int determinantSign() const {
    requireFactorized("SaddleSystem::determinantSign");
    return determinantSign_;
  }

  std::uint64_t symbolicAnalyses() const { return symbolicAnalyses_; }
  std::uint64_t numericFactorizations() const { return numericFactorizations_; }
  std::uint64_t solveColumns() const { return solveColumns_; }
  std::uint64_t nonzeros() const { return nonzeros_; }
  std::uint64_t factorNonzeros() const { return factorNonzeros_; }
  double maxRelativeResidual() const { return maxRelativeResidual_; }

private:
  static Eigen::Index panelColumns() { return 32; }

  static void fail(const std::string & context, const std::string & message) {
    throw std::runtime_error(context + ": " + message);
  }

  static void validateFinite(const Sparse & matrix,
                             const std::string & context,
                             const char * name) {
    for(Eigen::Index outer = 0; outer < matrix.outerSize(); ++outer) {
      for(Sparse::InnerIterator it(matrix, outer); it; ++it) {
        if(!std::isfinite(it.value())) {
          fail(context, std::string(name) + " contains non-finite values");
        }
      }
    }
  }

  bool hasSamePattern(const Sparse & matrix) const {
    if(!analyzed_ || matrix.rows() != cachedRows_ ||
       matrix.cols() != cachedCols_ ||
       matrix.nonZeros() != static_cast<Eigen::Index>(cachedInner_.size())) {
      return false;
    }
    for(Eigen::Index i = 0; i <= matrix.cols(); ++i) {
      if(matrix.outerIndexPtr()[i] != cachedOuter_[static_cast<std::size_t>(i)]) {
        return false;
      }
    }
    for(Eigen::Index i = 0; i < matrix.nonZeros(); ++i) {
      if(matrix.innerIndexPtr()[i] != cachedInner_[static_cast<std::size_t>(i)]) {
        return false;
      }
    }
    return true;
  }

  void cachePattern(const Sparse & matrix) {
    cachedRows_ = matrix.rows();
    cachedCols_ = matrix.cols();
    cachedOuter_.assign(matrix.outerIndexPtr(),
                        matrix.outerIndexPtr() + matrix.cols() + 1);
    if(matrix.nonZeros() == 0) {
      cachedInner_.clear();
    } else {
      cachedInner_.assign(matrix.innerIndexPtr(),
                          matrix.innerIndexPtr() + matrix.nonZeros());
    }
  }

  Eigen::Index residualCount() const {
    return analyzedRows_ - coefficientCount_;
  }

  void requireFactorized(const std::string & context) const {
    if(!factorized_) {
      fail(context, "factorize must succeed before this operation");
    }
  }

  Eigen::MatrixXd solvePanel(const Eigen::MatrixXd & rhs,
                             const std::string & context,
                             Eigen::Index requestedColumns) {
    Eigen::MatrixXd solution = factor_.solve(rhs);
    if(factor_.info() != Eigen::Success || !solution.allFinite()) {
      fail(context, "SparseLU solve failed");
    }
    double relativeResidual = 0.0;
    for(int refinement = 0; refinement <= 2; ++refinement) {
      const Eigen::MatrixXd residual = rhs - assembledMatrix_ * solution;
      if(!residual.allFinite()) {
        fail(context, "solve residual contains non-finite values");
      }
      relativeResidual = 0.0;
      for(Eigen::Index column = 0; column < rhs.cols(); ++column) {
        const double rhsNorm = rhs.col(column).stableNorm();
        const double residualNorm = residual.col(column).stableNorm();
        const double columnResidual = rhsNorm == 0.0 ? residualNorm : residualNorm / rhsNorm;
        if(!std::isfinite(columnResidual)) {
          fail(context, "solve residual is non-finite");
        }
        relativeResidual = (std::max)(relativeResidual, columnResidual);
      }
      if(!std::isfinite(relativeResidual)) {
        fail(context, "solve residual is non-finite");
      }
      if(relativeResidual <= 1e-10) {
        break;
      }
      if(refinement == 2) {
        fail(context, "solve relative residual exceeds 1e-10");
      }
      const Eigen::MatrixXd correction = factor_.solve(residual);
      if(factor_.info() != Eigen::Success || !correction.allFinite()) {
        fail(context, "SparseLU iterative-refinement solve failed");
      }
      solution += correction;
      if(!solution.allFinite()) {
        fail(context, "iteratively refined solution contains non-finite values");
      }
    }
    solveColumns_ += static_cast<std::uint64_t>(requestedColumns);
    maxRelativeResidual_ = (std::max)(maxRelativeResidual_, relativeResidual);
    return solution;
  }

  Factor factor_;
  bool analyzed_;
  bool factorized_;
  bool emptySystem_;
  Eigen::Index coefficientCount_;
  Eigen::Index cachedRows_ = 0;
  Eigen::Index cachedCols_ = 0;
  Eigen::Index analyzedRows_ = 0;
  double logAbsDeterminant_;
  int determinantSign_;
  std::uint64_t symbolicAnalyses_;
  std::uint64_t numericFactorizations_;
  std::uint64_t solveColumns_;
  std::uint64_t nonzeros_;
  std::uint64_t factorNonzeros_;
  double maxRelativeResidual_;
  std::vector<int> cachedOuter_;
  std::vector<int> cachedInner_;
  Sparse assembledMatrix_;
};

}  // namespace sommer

#endif
