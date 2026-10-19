#ifndef SOMMER_MME_PCG_RECYCLE_H
#define SOMMER_MME_PCG_RECYCLE_H

#include <Eigen/Cholesky>
#include <Eigen/Core>
#include <Eigen/Eigenvalues>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
#include <stdexcept>
#include <vector>

namespace sommer {
namespace pcg {

class RecycleSpace {
public:
  typedef std::function<Eigen::VectorXd(const Eigen::VectorXd &)> Operator;

  struct Selection {
    Eigen::MatrixXd basis;
    Eigen::VectorXd ritzValues;
  };

  explicit RecycleSpace(std::size_t maxCandidates = 32,
                        std::size_t maxRank = 8)
    : maxCandidates_(maxCandidates), maxRank_(maxRank), dimension_(0),
      candidateCount_(0), candidateStart_(0), prepared_(false),
      selectionOperatorApplications_(0), preparationOperatorApplications_(0) {
    if(maxCandidates_ == 0 || maxRank_ == 0 ||
       maxCandidates_ > static_cast<std::size_t>(
         std::numeric_limits<Eigen::Index>::max()) ||
       maxRank_ > static_cast<std::size_t>(
         std::numeric_limits<Eigen::Index>::max()) ||
       maxCandidates_ > static_cast<std::size_t>(
         std::numeric_limits<Eigen::Index>::max()) - maxRank_) {
      throw std::invalid_argument("PCG recycle-space bounds are invalid.");
    }
  }

  std::size_t maxCandidates() const { return maxCandidates_; }
  std::size_t maxRank() const { return maxRank_; }
  std::size_t candidateCount() const { return candidateCount_; }
  std::size_t maxStoredDirections() const {
    const std::size_t max = static_cast<std::size_t>(-1);
    return maxRank_ > (max - maxCandidates_) / 2
      ? max : maxCandidates_ + 2 * maxRank_;
  }
  // Persistent capacity only; callers must budget selection scratch separately.
  std::size_t storageBytesUpperBound() const {
    if(dimension_ == 0) { return 0; }
    const std::size_t directions = maxStoredDirections();
    const std::size_t dimension = static_cast<std::size_t>(dimension_);
    const std::size_t max = static_cast<std::size_t>(-1);
    if(directions == max || dimension > max / directions) { return max; }
    const std::size_t directionCapacity = dimension * directions;
    if(maxRank_ != 0 && maxRank_ > max / maxRank_) { return max; }
    const std::size_t square = maxRank_ * maxRank_;
    if(square > (max - directionCapacity) / 2) { return max; }
    const std::size_t scalarCapacity = directionCapacity + 2 * square;
    return scalarCapacity > max / sizeof(double)
      ? max : scalarCapacity * sizeof(double);
  }
  std::size_t selectionOperatorApplications() const {
    return selectionOperatorApplications_;
  }
  std::size_t preparationOperatorApplications() const {
    return preparationOperatorApplications_;
  }
  bool prepared() const { return prepared_; }
  Eigen::MatrixXd candidates() const {
    Eigen::MatrixXd result(dimension_,
      static_cast<Eigen::Index>(candidateCount_));
    for(std::size_t j = 0; j < candidateCount_; ++j) {
      result.col(static_cast<Eigen::Index>(j)) =
        candidates_.col(static_cast<Eigen::Index>(
          (candidateStart_ + j) % maxCandidates_));
    }
    return result;
  }
  const Eigen::MatrixXd & basis() const { return basis_; }
  const Eigen::MatrixXd & coarseMatrix() const { return coarse_; }

  void clearCandidates() {
    candidateCount_ = 0;
    candidateStart_ = 0;
  }

  void capture(const Eigen::VectorXd & direction) {
    validateDirectionDimension(direction.size());
    if(!direction.allFinite()) {
      throw std::invalid_argument("PCG recycle direction must be finite.");
    }
    if(direction.size() == 0) {
      return;
    }
    const double norm = direction.stableNorm();
    if(!std::isfinite(norm) || norm == 0.0) {
      return;
    }
    if(dimension_ == 0) {
      dimension_ = direction.size();
      candidates_.resize(dimension_, static_cast<Eigen::Index>(maxCandidates_));
    }
    Eigen::VectorXd normalized = direction / norm;
    if(candidateCount_ < maxCandidates_) {
      const std::size_t slot = (candidateStart_ + candidateCount_) % maxCandidates_;
      candidates_.col(static_cast<Eigen::Index>(slot)) = normalized;
      ++candidateCount_;
    } else {
      candidates_.col(static_cast<Eigen::Index>(candidateStart_)) = normalized;
      candidateStart_ = (candidateStart_ + 1) % maxCandidates_;
    }
  }

  void capturePanel(const Eigen::MatrixXd & directions) {
    validateDirectionDimension(directions.rows());
    if(!directions.allFinite()) {
      throw std::invalid_argument("PCG recycle direction panel must be finite.");
    }
    for(Eigen::Index col = 0; col < directions.cols(); ++col) {
      capture(directions.col(col));
    }
  }

  Selection selectCandidates(const Operator & applyA,
                             const Eigen::VectorXd & diagonalMetric,
                             std::size_t requestedRank) {
    validateMetric(diagonalMetric);
    if(diagonalMetric.size() != dimension_) {
      throw std::invalid_argument("PCG recycle metric dimension is invalid.");
    }
    Selection selected;
    const std::size_t totalDirections = candidateCount_ +
      static_cast<std::size_t>(basis_.cols());
    if(totalDirections == 0 || requestedRank == 0) {
      selected.basis.resize(dimension_, 0);
      selected.ritzValues.resize(0);
      return selected;
    }

    const Eigen::VectorXd sqrtMetric = diagonalMetric.array().sqrt();
    const Eigen::VectorXd inverseSqrtMetric = sqrtMetric.cwiseInverse();
    Eigen::MatrixXd orthonormal(dimension_,
      static_cast<Eigen::Index>(totalDirections));
    Eigen::Index rank = 0;
    for(std::size_t j = 0; j < totalDirections; ++j) {
      Eigen::VectorXd source;
      if(j < static_cast<std::size_t>(basis_.cols())) {
        source = basis_.col(static_cast<Eigen::Index>(j));
      } else {
        const std::size_t candidate = j -
          static_cast<std::size_t>(basis_.cols());
        source = candidates_.col(static_cast<Eigen::Index>(
          (candidateStart_ + candidate) % maxCandidates_));
      }
      Eigen::VectorXd vector = sqrtMetric.cwiseProduct(source);
      const double originalNorm = vector.stableNorm();
      if(!std::isfinite(originalNorm) || originalNorm == 0.0) { continue; }
      vector /= originalNorm;
      for(int pass = 0; pass < 2; ++pass) {
        for(Eigen::Index i = 0; i < rank; ++i) {
          vector.noalias() -= orthonormal.col(i).dot(vector) *
            orthonormal.col(i);
        }
      }
      const double remainder = vector.stableNorm();
      if(!std::isfinite(remainder) || remainder <= 1.0e-10) { continue; }
      orthonormal.col(rank++) = vector / remainder;
    }
    if(rank == 0) {
      selected.basis.resize(dimension_, 0);
      selected.ritzValues.resize(0);
      return selected;
    }
    orthonormal.conservativeResize(Eigen::NoChange, rank);
    const Eigen::MatrixXd vectors =
      inverseSqrtMetric.asDiagonal() * orthonormal;
    Eigen::MatrixXd operatorVectors(dimension_, rank);
    for(Eigen::Index col = 0; col < rank; ++col) {
      ++selectionOperatorApplications_;
      const Eigen::VectorXd applied = applyA(vectors.col(col));
      if(applied.size() != dimension_ || !applied.allFinite()) {
        throw std::runtime_error("PCG recycle operator returned an invalid vector.");
      }
      operatorVectors.col(col) = applied;
    }
    Eigen::MatrixXd projected = vectors.transpose() * operatorVectors;
    projected = 0.5 * (projected + projected.transpose());
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(projected);
    if(eigensolver.info() != Eigen::Success ||
       !eigensolver.eigenvalues().allFinite() ||
       !eigensolver.eigenvectors().allFinite()) {
      throw std::runtime_error("PCG recycle Ritz decomposition failed.");
    }
    if((eigensolver.eigenvalues().array() <= 0.0).any()) {
      throw std::invalid_argument("PCG recycle operator is not positive definite.");
    }
    const std::size_t keepCount = std::min<std::size_t>(maxRank_,
      std::max<std::size_t>(requestedRank,
        static_cast<std::size_t>(basis_.cols())));
    const Eigen::Index keep = static_cast<Eigen::Index>(
      std::min<std::size_t>(keepCount, static_cast<std::size_t>(rank)));
    selected.basis = vectors * eigensolver.eigenvectors().leftCols(keep);
    selected.ritzValues = eigensolver.eigenvalues().head(keep);
    return selected;
  }

  void accept(const Selection & selection) {
    if(selection.basis.rows() != dimension_ ||
       selection.basis.cols() != selection.ritzValues.size() ||
       selection.basis.cols() > static_cast<Eigen::Index>(maxRank_) ||
       !selection.basis.allFinite() || !selection.ritzValues.allFinite()) {
      throw std::invalid_argument("PCG recycle basis is invalid.");
    }
    basis_ = selection.basis;
    clearCandidates();
    prepared_ = false;
    coarse_.resize(0, 0);
    coarseVectors_.resize(dimension_, 0);
  }

  void prepare(const Operator & applyA, const Operator & applyM) {
    prepared_ = false;
    coarse_.resize(0, 0);
    coarseVectors_.resize(dimension_, 0);
    if(basis_.cols() == 0) {
      return;
    }
    coarseVectors_.resize(dimension_, basis_.cols());
    for(Eigen::Index col = 0; col < basis_.cols(); ++col) {
      ++preparationOperatorApplications_;
      const Eigen::VectorXd applied = applyA(basis_.col(col));
      if(applied.size() != dimension_ || !applied.allFinite()) {
        throw std::runtime_error("PCG recycle operator returned an invalid vector.");
      }
      coarseVectors_.col(col) = applied;
    }
    coarse_ = basis_.transpose() * coarseVectors_;
    coarse_ = 0.5 * (coarse_ + coarse_.transpose());
    if(!factorPositiveDefinite(coarse_, coarseFactor_)) {
      prepared_ = false;
      throw std::invalid_argument("PCG recycle coarse matrix is not positive definite.");
    }
    applyM_ = applyM;
    prepared_ = true;
  }

  Eigen::VectorXd apply(const Eigen::VectorXd & value) const {
    if(!prepared_) {
      throw std::logic_error("PCG recycle space is not prepared.");
    }
    if(value.size() != dimension_ || !value.allFinite()) {
      throw std::invalid_argument("PCG recycle input is invalid.");
    }
    const Eigen::VectorXd coarseRhs = basis_.transpose() * value;
    const Eigen::VectorXd coarseSolution = coarseFactor_.solve(coarseRhs);
    const Eigen::VectorXd projected = value - coarseVectors_ * coarseSolution;
    const Eigen::VectorXd base = applyM_(projected);
    if(base.size() != dimension_ || !base.allFinite()) {
      throw std::runtime_error("PCG preconditioner returned an invalid vector.");
    }
    Eigen::VectorXd result = base - basis_ *
      coarseFactor_.solve(coarseVectors_.transpose() * base) +
      basis_ * coarseSolution;
    const double positivityCheck = value.dot(result);
    if(!result.allFinite() ||
       (value.stableNorm() > 0.0 &&
        (!std::isfinite(positivityCheck) || positivityCheck <= 0.0))) {
      throw std::runtime_error("PCG recycle preconditioner is not positive definite.");
    }
    return result;
  }

  Eigen::VectorXd coarseCorrection(const Eigen::VectorXd & residual) const {
    if(!prepared_) {
      throw std::logic_error("PCG recycle space is not prepared.");
    }
    if(residual.size() != dimension_ || !residual.allFinite()) {
      throw std::invalid_argument("PCG recycle residual is invalid.");
    }
    return basis_ * coarseFactor_.solve(basis_.transpose() * residual);
  }

  Eigen::MatrixXd coarseCorrectionPanel(
      const Eigen::MatrixXd & residuals) const {
    if(residuals.rows() != dimension_ || !residuals.allFinite()) {
      throw std::invalid_argument("PCG recycle residual panel is invalid.");
    }
    if(!prepared_) {
      throw std::logic_error("PCG recycle space is not prepared.");
    }
    return basis_ * coarseFactor_.solve(basis_.transpose() * residuals);
  }

  Eigen::MatrixXd applyPanel(const Eigen::MatrixXd & values) const {
    if(values.rows() != dimension_ || !values.allFinite()) {
      throw std::invalid_argument("PCG recycle input panel is invalid.");
    }
    Eigen::MatrixXd result(dimension_, values.cols());
    for(Eigen::Index col = 0; col < values.cols(); ++col) {
      result.col(col) = apply(Eigen::VectorXd(values.col(col)));
    }
    return result;
  }

  void reset() {
    dimension_ = 0;
    candidateCount_ = 0;
    candidateStart_ = 0;
    candidates_.resize(0, 0);
    basis_.resize(0, 0);
    coarse_.resize(0, 0);
    coarseVectors_.resize(0, 0);
    prepared_ = false;
  }

private:
  void validateDirectionDimension(Eigen::Index size) const {
    if(dimension_ != 0 && size != dimension_) {
      throw std::invalid_argument("PCG recycle direction dimension is invalid.");
    }
  }

  static bool factorPositiveDefinite(const Eigen::MatrixXd & matrix,
                                    Eigen::LLT<Eigen::MatrixXd> & factor) {
    if(matrix.rows() == 0 || matrix.rows() != matrix.cols() ||
       !matrix.allFinite()) {
      return false;
    }
    factor.compute(matrix);
    return factor.info() == Eigen::Success;
  }

  static void validateMetric(const Eigen::VectorXd & metric) {
    if(metric.size() == 0 || !metric.allFinite() ||
       (metric.array() <= 0.0).any()) {
      throw std::invalid_argument("PCG recycle diagonal metric must be positive and finite.");
    }
  }

  const std::size_t maxCandidates_;
  const std::size_t maxRank_;
  Eigen::Index dimension_;
  std::size_t candidateCount_;
  std::size_t candidateStart_;
  Eigen::MatrixXd candidates_;
  Eigen::MatrixXd basis_;
  Eigen::MatrixXd coarse_;
  Eigen::MatrixXd coarseVectors_;
  Eigen::LLT<Eigen::MatrixXd> coarseFactor_;
  Operator applyM_;
  bool prepared_;
  std::size_t selectionOperatorApplications_;
  std::size_t preparationOperatorApplications_;
};

}  // namespace pcg
}  // namespace sommer

#endif
