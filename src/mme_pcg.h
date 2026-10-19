#ifndef SOMMER_MME_PCG_H
#define SOMMER_MME_PCG_H

#include <Eigen/Cholesky>
#include <Eigen/Core>
#include <Eigen/QR>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include "mme_pcg_recycle.h"

namespace sommer {
namespace pcg {

typedef std::function<Eigen::VectorXd(const Eigen::VectorXd &)> VectorOperator;
typedef std::function<Eigen::MatrixXd(const Eigen::MatrixXd &)> MatrixOperator;

struct Options {
  double tolerance;
  int maxIterations;
  int residualRefresh;
  Eigen::Index blockSize;
  std::size_t maxWorkspaceBytes;
  RecycleSpace * recycleSpace;

  Options()
    : tolerance(1.0e-8), maxIterations(0), residualRefresh(0),
      blockSize(16), maxWorkspaceBytes(64u * 1024u * 1024u),
      recycleSpace(0) {}
};

struct Metrics {
  std::size_t operatorApplications;
  std::size_t vectorApplications;
  std::size_t matrixPanelApplications;
  std::size_t preconditionerApplications;
  std::size_t trueResidualChecks;
  std::size_t residualReplacements;
  std::size_t iterations;
  std::size_t lockedColumns;
  std::size_t deflatedColumns;
  std::size_t blockFallbacks;
  std::string fallbackReason;

  Metrics()
    : operatorApplications(0), vectorApplications(0), matrixPanelApplications(0),
      preconditionerApplications(0), trueResidualChecks(0), residualReplacements(0),
      iterations(0), lockedColumns(0), deflatedColumns(0), blockFallbacks(0) {}

  void add(const Metrics & other) {
    operatorApplications += other.operatorApplications;
    vectorApplications += other.vectorApplications;
    matrixPanelApplications += other.matrixPanelApplications;
    preconditionerApplications += other.preconditionerApplications;
    trueResidualChecks += other.trueResidualChecks;
    residualReplacements += other.residualReplacements;
    iterations += other.iterations;
    lockedColumns += other.lockedColumns;
    deflatedColumns += other.deflatedColumns;
    blockFallbacks += other.blockFallbacks;
    if(fallbackReason.empty() && !other.fallbackReason.empty()) {
      fallbackReason = other.fallbackReason;
    }
  }
};

struct BlockBreakdown : std::runtime_error {
  Metrics metrics;
  BlockBreakdown(const std::string & message, const Metrics & completed)
    : std::runtime_error(message), metrics(completed) {}
};

struct Result {
  Eigen::MatrixXd solution;
  Metrics metrics;
};

inline void validateOptions(const Options & options) {
  if(!std::isfinite(options.tolerance) || options.tolerance <= 0.0) {
    throw std::invalid_argument("PCG tolerance must be positive and finite.");
  }
  if(options.maxIterations < 0 || options.residualRefresh < 0 ||
     options.blockSize < 1 || options.maxWorkspaceBytes == 0) {
    throw std::invalid_argument("PCG iteration, refresh, block, and workspace controls are invalid.");
  }
}

inline int iterationLimit(Eigen::Index n, const Options & options) {
  if(options.maxIterations > 0) { return options.maxIterations; }
  const Eigen::Index bounded = std::min<Eigen::Index>(
    std::max<Eigen::Index>(1000, n),
    static_cast<Eigen::Index>(std::numeric_limits<int>::max()));
  return static_cast<int>(bounded);
}

inline Result solvePackedIndependent(const Eigen::MatrixXd & rhs,
                                     const MatrixOperator & applyA,
                                     const MatrixOperator & applyM,
                                     const Options & options,
                                     const Eigen::MatrixXd * initialGuess = 0) {
  validateOptions(options);
  if(!rhs.allFinite()){ throw std::invalid_argument("PCG RHS must be finite."); }
  const long double needed = static_cast<long double>(rhs.rows()) *
    std::min<Eigen::Index>(rhs.cols(), options.blockSize) * 12 * sizeof(double);
  if(needed > options.maxWorkspaceBytes){
    throw std::runtime_error("PCG packed panel exceeds its workspace budget.");
  }
  Result result;
  result.solution = Eigen::MatrixXd::Zero(rhs.rows(), rhs.cols());
  const Eigen::Index panelSize = options.blockSize;
  for(Eigen::Index first = 0; first < rhs.cols(); first += panelSize){
    const Eigen::Index count = std::min(panelSize, rhs.cols() - first);
    const Eigen::MatrixXd panel = rhs.middleCols(first, count);
    Eigen::MatrixXd x = Eigen::MatrixXd::Zero(rhs.rows(), count);
    if(initialGuess && initialGuess->rows() == rhs.rows() && initialGuess->cols() == rhs.cols()){
      x = initialGuess->middleCols(first, count);
    }
    if(!x.allFinite()){ throw std::invalid_argument("PCG initial guess must be finite."); }
    const Eigen::RowVectorXd norms = panel.colwise().norm();
    for(Eigen::Index col = 0; col < count; ++col){
      if(norms(col) == 0.0){ x.col(col).setZero(); }
    }
    auto checkedApply = [&](const Eigen::MatrixXd & values) -> Eigen::MatrixXd {
      ++result.metrics.operatorApplications;
      ++result.metrics.matrixPanelApplications;
      result.metrics.vectorApplications += static_cast<std::size_t>(values.cols());
      Eigen::MatrixXd output = applyA(values);
      if(output.rows() != values.rows() || output.cols() != values.cols() || !output.allFinite()){
        throw std::runtime_error("PCG operator returned an invalid panel.");
      }
      return output;
    };
    auto precondition = [&](const Eigen::MatrixXd & values) -> Eigen::MatrixXd {
      ++result.metrics.preconditionerApplications;
      Eigen::MatrixXd output =
        options.recycleSpace != 0 && options.recycleSpace->prepared()
          ? options.recycleSpace->applyPanel(values) : applyM(values);
      if(output.rows() != values.rows() || output.cols() != values.cols() || !output.allFinite()){
        throw std::runtime_error("PCG preconditioner returned an invalid panel.");
      }
      return output;
    };
    Eigen::MatrixXd r = panel - checkedApply(x);
    result.metrics.trueResidualChecks += count;
    Eigen::MatrixXd z = precondition(r), p = z;
    if(options.recycleSpace != 0) { options.recycleSpace->capturePanel(p); }
    Eigen::RowVectorXd rho = (r.array() * z.array()).colwise().sum();
    std::vector<bool> active(static_cast<std::size_t>(count), true);
    for(Eigen::Index col = 0; col < count; ++col){
      if(norms(col) == 0.0 || r.col(col).norm() <= options.tolerance * norms(col)){
        active[col] = false;
        p.col(col).setZero();
        ++result.metrics.lockedColumns;
      }else if(!std::isfinite(rho(col)) || rho(col) <= 0.0){
        throw std::runtime_error("PCG encountered a non-positive preconditioned residual.");
      }
    }
    int iteration = 0;
    for(; iteration < iterationLimit(rhs.rows(), options) &&
        std::find(active.begin(), active.end(), true) != active.end(); ++iteration){
      const Eigen::MatrixXd ap = checkedApply(p);
      if(options.recycleSpace != 0) { options.recycleSpace->capturePanel(p); }
      for(Eigen::Index col = 0; col < count; ++col){
        if(!active[col]) continue;
        const double curvature = p.col(col).dot(ap.col(col));
        if(!std::isfinite(curvature) || curvature <= 0.0){
          throw std::runtime_error("PCG encountered non-positive curvature.");
        }
        const double alpha = rho(col) / curvature;
        x.col(col).noalias() += alpha * p.col(col);
        r.col(col).noalias() -= alpha * ap.col(col);
      }
      const bool periodic = options.residualRefresh > 0 && (iteration + 1) % options.residualRefresh == 0;
      bool verify = periodic;
      for(Eigen::Index col = 0; col < count; ++col){
        verify = verify || (active[col] && r.col(col).norm() <= options.tolerance * norms(col));
      }
      std::vector<bool> restarted(static_cast<std::size_t>(count), false);
      if(verify){
        const Eigen::MatrixXd trueResidual = panel - checkedApply(x);
        result.metrics.trueResidualChecks += count;
        for(Eigen::Index col = 0; col < count; ++col){
          if(!active[col]) continue;
          const bool candidate = r.col(col).norm() <= options.tolerance * norms(col);
          if(trueResidual.col(col).norm() <= options.tolerance * norms(col)){
            active[col] = false;
            p.col(col).setZero();
            ++result.metrics.lockedColumns;
          }else if(periodic || candidate){
            r.col(col) = trueResidual.col(col);
            restarted[col] = true;
            ++result.metrics.residualReplacements;
          }
        }
      }
      if(std::find(active.begin(), active.end(), true) == active.end()) { ++iteration; break; }
      z = precondition(r);
      for(Eigen::Index col = 0; col < count; ++col){
        if(!active[col]) continue;
        const double next = r.col(col).dot(z.col(col));
        if(!std::isfinite(next) || next <= 0.0){
          throw std::runtime_error("PCG encountered a non-positive residual update.");
        }
        if(restarted[col]) { p.col(col) = z.col(col); }
        else { p.col(col) = z.col(col) + (next / rho(col)) * p.col(col); }
        rho(col) = next;
      }
    }
    result.metrics.iterations += iteration;
    if(std::find(active.begin(), active.end(), true) != active.end()){
      throw std::runtime_error("PCG failed to converge within its iteration limit.");
    }
    result.solution.middleCols(first, count) = x;
  }
  return result;
}

inline Eigen::VectorXd applyVector(const VectorOperator & applyA,
                                   const Eigen::VectorXd & value,
                                   Metrics & metrics) {
  ++metrics.operatorApplications;
  ++metrics.vectorApplications;
  Eigen::VectorXd result = applyA(value);
  if(result.size() != value.size() || !result.allFinite()) {
    throw std::runtime_error("PCG operator returned an invalid vector.");
  }
  return result;
}

inline Eigen::MatrixXd applyMatrix(const MatrixOperator & applyA,
                                   const Eigen::MatrixXd & value,
                                   Metrics & metrics) {
  ++metrics.operatorApplications;
  ++metrics.matrixPanelApplications;
  Eigen::MatrixXd result = applyA(value);
  if(result.rows() != value.rows() || result.cols() != value.cols() ||
     !result.allFinite()) {
    throw std::runtime_error("PCG operator returned an invalid panel.");
  }
  return result;
}

inline Eigen::VectorXd applyPreconditioner(
    const VectorOperator & applyM, const Eigen::VectorXd & value,
    Metrics & metrics, RecycleSpace * recycleSpace = 0) {
  ++metrics.preconditionerApplications;
  Eigen::VectorXd result = recycleSpace != 0 && recycleSpace->prepared()
    ? recycleSpace->apply(value) : applyM(value);
  if(result.size() != value.size() || !result.allFinite()) {
    throw std::runtime_error("PCG preconditioner returned an invalid vector.");
  }
  return result;
}

inline Eigen::MatrixXd applyPreconditioner(
    const MatrixOperator & applyM, const Eigen::MatrixXd & value,
    Metrics & metrics, RecycleSpace * recycleSpace = 0) {
  ++metrics.preconditionerApplications;
  Eigen::MatrixXd result = recycleSpace != 0 && recycleSpace->prepared()
    ? recycleSpace->applyPanel(value) : applyM(value);
  if(result.rows() != value.rows() || result.cols() != value.cols() ||
     !result.allFinite()) {
    throw std::runtime_error("PCG preconditioner returned an invalid panel.");
  }
  return result;
}

inline Eigen::VectorXd solveOne(const Eigen::VectorXd & rhs,
                                const VectorOperator & applyA,
                                const VectorOperator & applyM,
                                const Options & options,
                                Metrics & metrics,
                                const Eigen::VectorXd * initialGuess = 0) {
  const Eigen::Index n = rhs.size();
  if(!rhs.allFinite()) { throw std::invalid_argument("PCG RHS must be finite."); }
  const double rhsNorm = rhs.norm();
  if(rhsNorm == 0.0) { return Eigen::VectorXd::Zero(n); }

  Eigen::VectorXd solution = initialGuess != 0 && initialGuess->size() == n &&
      initialGuess->allFinite()
    ? *initialGuess : Eigen::VectorXd::Zero(n);
  Eigen::VectorXd residual = rhs - applyVector(applyA, solution, metrics);
  ++metrics.trueResidualChecks;
  double relativeResidual = residual.norm() / rhsNorm;
  if(relativeResidual <= options.tolerance) { return solution; }

  Eigen::VectorXd preconditioned = applyPreconditioner(
    applyM, residual, metrics, options.recycleSpace);
  double residualProduct = residual.dot(preconditioned);
  if(!std::isfinite(residualProduct) || residualProduct <= 0.0) {
    throw std::runtime_error("PCG encountered a non-positive preconditioned residual.");
  }
  Eigen::VectorXd direction = preconditioned;
  if(options.recycleSpace != 0) { options.recycleSpace->capture(direction); }
  const int limit = iterationLimit(n, options);

  for(int iteration = 0; iteration < limit; ++iteration) {
    const Eigen::VectorXd operatorDirection =
      applyVector(applyA, direction, metrics);
    const double curvature = direction.dot(operatorDirection);
    if(!std::isfinite(curvature) || curvature <= 0.0) {
      throw std::runtime_error("PCG encountered non-positive curvature.");
    }
    const double alpha = residualProduct / curvature;
    if(!std::isfinite(alpha)) { throw std::runtime_error("PCG step is non-finite."); }
    solution.noalias() += alpha * direction;
    residual.noalias() -= alpha * operatorDirection;
    ++metrics.iterations;

    const bool refresh = options.residualRefresh > 0 &&
      (iteration + 1) % options.residualRefresh == 0;
    relativeResidual = residual.norm() / rhsNorm;
    if(refresh || relativeResidual <= options.tolerance) {
      Eigen::VectorXd trueResidual =
        rhs - applyVector(applyA, solution, metrics);
      ++metrics.trueResidualChecks;
      const double trueRelativeResidual = trueResidual.norm() / rhsNorm;
      if(trueRelativeResidual <= options.tolerance) { return solution; }
      residual.swap(trueResidual);
      ++metrics.residualReplacements;
      if(refresh || relativeResidual <= options.tolerance) {
        preconditioned = applyPreconditioner(
          applyM, residual, metrics, options.recycleSpace);
        residualProduct = residual.dot(preconditioned);
        if(!std::isfinite(residualProduct) || residualProduct <= 0.0) {
          throw std::runtime_error("PCG residual replacement produced a non-positive update.");
        }
        direction = preconditioned;
        if(options.recycleSpace != 0) { options.recycleSpace->capture(direction); }
        continue;
      }
    }

    preconditioned = applyPreconditioner(
      applyM, residual, metrics, options.recycleSpace);
    const double nextProduct = residual.dot(preconditioned);
    if(!std::isfinite(nextProduct) || nextProduct <= 0.0) {
      throw std::runtime_error("PCG encountered a non-positive residual update.");
    }
    const double beta = nextProduct / residualProduct;
    if(!std::isfinite(beta)) { throw std::runtime_error("PCG recurrence is non-finite."); }
    direction = preconditioned + beta * direction;
    if(options.recycleSpace != 0) { options.recycleSpace->capture(direction); }
    residualProduct = nextProduct;
  }

  throw std::runtime_error("PCG failed to converge within its iteration limit.");
}

inline Result solveIndependent(const Eigen::MatrixXd & rhs,
                               const VectorOperator & applyA,
                               const VectorOperator & applyM,
                               const Options & options,
                               const Eigen::MatrixXd * initialGuess = 0) {
  validateOptions(options);
  Result result;
  result.solution.resize(rhs.rows(), rhs.cols());
  for(Eigen::Index col = 0; col < rhs.cols(); ++col) {
    Eigen::VectorXd localGuess;
    if(initialGuess != 0 && initialGuess->rows() == rhs.rows() &&
       initialGuess->cols() == rhs.cols()) {
      localGuess = initialGuess->col(col);
      result.solution.col(col) = solveOne(rhs.col(col), applyA, applyM,
        options, result.metrics, &localGuess);
    } else {
      result.solution.col(col) = solveOne(rhs.col(col), applyA, applyM,
        options, result.metrics, 0);
    }
  }
  return result;
}

inline bool factorPositiveDefinite(const Eigen::MatrixXd & matrix,
                                   Eigen::LLT<Eigen::MatrixXd> & factor) {
  if(matrix.rows() == 0 || matrix.rows() != matrix.cols() ||
     !matrix.allFinite()) { return false; }
  factor.compute(0.5 * (matrix + matrix.transpose()));
  return factor.info() == Eigen::Success;
}

inline Result solveBlockPanel(const Eigen::MatrixXd & rhs,
                              const MatrixOperator & applyA,
                              const MatrixOperator & applyM,
                              const Options & options,
                              const Eigen::MatrixXd * initialGuess) {
  Result result;
  const Eigen::Index n = rhs.rows();
  const Eigen::Index columns = rhs.cols();
  if(columns == 0) { result.solution.resize(n, 0); return result; }
  if(!rhs.allFinite()) { throw std::invalid_argument("PCG RHS must be finite."); }

  Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(rhs);
  qr.setThreshold(std::max(options.tolerance * 0.01, 1.0e-13));
  const Eigen::Index rank = qr.rank();
  if(rank == 0) {
    result.solution = Eigen::MatrixXd::Zero(n, columns);
    return result;
  }

  std::vector<Eigen::Index> pivots(static_cast<std::size_t>(rank));
  for(Eigen::Index i = 0; i < rank; ++i) {
    pivots[static_cast<std::size_t>(i)] =
      qr.colsPermutation().indices()(i);
  }
  Eigen::MatrixXd independentRhs(n, rank);
  for(Eigen::Index i = 0; i < rank; ++i) {
    independentRhs.col(i) = rhs.col(pivots[static_cast<std::size_t>(i)]);
  }
  const Eigen::MatrixXd qrCoefficients = qr.solve(rhs);
  Eigen::MatrixXd coefficients(rank, columns);
  for(Eigen::Index i = 0; i < rank; ++i) {
    coefficients.row(i) =
      qrCoefficients.row(pivots[static_cast<std::size_t>(i)]);
  }
  const Eigen::MatrixXd reconstruction = independentRhs * coefficients;
  const double rhsScale = std::max(1.0, rhs.norm());
  if(!reconstruction.allFinite() ||
     (reconstruction - rhs).norm() > options.tolerance * rhsScale) {
    throw BlockBreakdown("PCG block RHS dependency could not be safely deflated.", result.metrics);
  }
  result.metrics.deflatedColumns = static_cast<std::size_t>(columns - rank);

  Eigen::MatrixXd solution = Eigen::MatrixXd::Zero(n, rank);
  if(initialGuess != 0 && initialGuess->rows() == n &&
     initialGuess->cols() == columns) {
    for(Eigen::Index i = 0; i < rank; ++i) {
      solution.col(i) = initialGuess->col(pivots[static_cast<std::size_t>(i)]);
    }
  }
  Eigen::MatrixXd residual = independentRhs -
    applyMatrix(applyA, solution, result.metrics);
  ++result.metrics.trueResidualChecks;
  const Eigen::RowVectorXd rhsNorms = independentRhs.colwise().norm();
  std::vector<bool> active(static_cast<std::size_t>(rank), true);
  Eigen::MatrixXd preconditioned = applyPreconditioner(
    applyM, residual, result.metrics, options.recycleSpace);
  Eigen::MatrixXd direction = preconditioned;
  if(options.recycleSpace != 0) { options.recycleSpace->capturePanel(direction); }
  Eigen::MatrixXd gamma = residual.transpose() * preconditioned;
  Eigen::Index activeCount = rank;
  const int limit = iterationLimit(n, options);

  for(Eigen::Index col = 0; col < rank; ++col) {
    if(rhsNorms(col) == 0.0 ||
       residual.col(col).norm() <= options.tolerance * rhsNorms(col)) {
      active[static_cast<std::size_t>(col)] = false;
      --activeCount;
      ++result.metrics.lockedColumns;
    }
  }

  int iteration = 0;
  while(activeCount > 0 && iteration < limit) {
    std::vector<Eigen::Index> activeColumns;
    activeColumns.reserve(static_cast<std::size_t>(activeCount));
    for(Eigen::Index col = 0; col < rank; ++col) {
      if(active[static_cast<std::size_t>(col)]) { activeColumns.push_back(col); }
    }
    Eigen::MatrixXd r(residual.rows(), activeCount);
    Eigen::MatrixXd p(direction.rows(), activeCount);
    Eigen::MatrixXd g(activeCount, activeCount);
    for(Eigen::Index i = 0; i < activeCount; ++i) {
      const Eigen::Index source = activeColumns[static_cast<std::size_t>(i)];
      r.col(i) = residual.col(source);
      p.col(i) = direction.col(source);
      for(Eigen::Index j = 0; j < activeCount; ++j) {
        g(i, j) = gamma(source, activeColumns[static_cast<std::size_t>(j)]);
      }
    }
    Eigen::MatrixXd ap = applyMatrix(applyA, p, result.metrics);
    if(options.recycleSpace != 0) { options.recycleSpace->capturePanel(p); }
    Eigen::MatrixXd curvature = p.transpose() * ap;
    Eigen::LLT<Eigen::MatrixXd> curvatureFactor;
    Eigen::LLT<Eigen::MatrixXd> gammaFactor;
    if(!factorPositiveDefinite(curvature, curvatureFactor) ||
       !factorPositiveDefinite(g, gammaFactor)) {
      throw BlockBreakdown("PCG block recurrence lost positive definiteness.", result.metrics);
    }
    const Eigen::MatrixXd alpha = curvatureFactor.solve(g);
    if(!alpha.allFinite()) { throw std::runtime_error("PCG block step is non-finite."); }
    for(Eigen::Index i = 0; i < activeCount; ++i) {
      const Eigen::Index target = activeColumns[static_cast<std::size_t>(i)];
      solution.col(target).noalias() += p * alpha.col(i);
      residual.col(target).noalias() -= ap * alpha.col(i);
    }
    ++iteration;
    ++result.metrics.iterations;

    const bool refresh = options.residualRefresh > 0 &&
      iteration % options.residualRefresh == 0;
    bool candidate = refresh;
    for(Eigen::Index i = 0; i < activeCount; ++i) {
      const Eigen::Index col = activeColumns[static_cast<std::size_t>(i)];
      candidate = candidate ||
        residual.col(col).norm() <= options.tolerance * rhsNorms(col);
    }
    bool restart = refresh;
    bool lockedAny = false;
    if(candidate) {
      Eigen::MatrixXd trueResidual = independentRhs -
        applyMatrix(applyA, solution, result.metrics);
      ++result.metrics.trueResidualChecks;
      for(Eigen::Index i = 0; i < activeCount; ++i) {
        const Eigen::Index col = activeColumns[static_cast<std::size_t>(i)];
        const double norm = rhsNorms(col);
        const double relative = norm == 0.0 ? 0.0 : trueResidual.col(col).norm() / norm;
        if(relative <= options.tolerance) {
          active[static_cast<std::size_t>(col)] = false;
          --activeCount;
          ++result.metrics.lockedColumns;
          lockedAny = true;
        } else if(refresh || residual.col(col).norm() <= options.tolerance * norm) {
          restart = true;
        }
      }
      restart = restart || lockedAny;
      if(restart) {
        residual.swap(trueResidual);
        ++result.metrics.residualReplacements;
      }
    }

    if(activeCount == 0) { break; }
    if(restart) {
      preconditioned = applyPreconditioner(
        applyM, residual, result.metrics, options.recycleSpace);
      direction = preconditioned;
      if(options.recycleSpace != 0) { options.recycleSpace->capturePanel(direction); }
      gamma = residual.transpose() * preconditioned;
      continue;
    }

    preconditioned = applyPreconditioner(
      applyM, residual, result.metrics, options.recycleSpace);
    const Eigen::MatrixXd nextGamma = residual.transpose() * preconditioned;
    Eigen::MatrixXd nextActiveGamma(activeCount, activeCount);
    for(Eigen::Index i = 0; i < activeCount; ++i){
      for(Eigen::Index j = 0; j < activeCount; ++j){
        nextActiveGamma(i,j) = nextGamma(activeColumns[static_cast<std::size_t>(i)],
          activeColumns[static_cast<std::size_t>(j)]);
      }
    }
    const Eigen::MatrixXd beta = gammaFactor.solve(nextActiveGamma);
    if(!beta.allFinite()) { throw std::runtime_error("PCG block recurrence is non-finite."); }
    for(Eigen::Index i = 0; i < activeCount; ++i){
      const Eigen::Index target = activeColumns[static_cast<std::size_t>(i)];
      direction.col(target) = preconditioned.col(target) + p * beta.col(i);
    }
    if(options.recycleSpace != 0) { options.recycleSpace->capturePanel(direction); }
    gamma = nextGamma;
  }

  if(activeCount > 0) {
    throw BlockBreakdown("PCG block solve failed to converge within its iteration limit.", result.metrics);
  }
  result.solution = solution * coefficients;

  const Eigen::MatrixXd trueResidual = rhs -
    applyMatrix(applyA, result.solution, result.metrics);
  ++result.metrics.trueResidualChecks;
  for(Eigen::Index col = 0; col < columns; ++col) {
    const double norm = rhs.col(col).norm();
    if(norm > 0.0 && trueResidual.col(col).norm() > options.tolerance * norm) {
      throw BlockBreakdown("PCG block solve failed final true-residual verification.", result.metrics);
    }
  }
  return result;
}

inline Result solveBlock(const Eigen::MatrixXd & rhs,
                         const MatrixOperator & applyA,
                         const MatrixOperator & applyM,
                         const VectorOperator & applyAVector,
                         const VectorOperator & applyMVector,
                         const Options & options,
                         const Eigen::MatrixXd * initialGuess = 0) {
  validateOptions(options);
  if(rhs.rows() < 0 || rhs.cols() < 0) {
    throw std::invalid_argument("PCG block dimensions are invalid.");
  }
  Result result;
  result.solution.resize(rhs.rows(), rhs.cols());
  if(rhs.cols() == 0) { return result; }

  const std::size_t maxBytes = options.maxWorkspaceBytes;
  const std::size_t rows = static_cast<std::size_t>(rhs.rows());
  const std::size_t bytesPerColumn =
    rows > 0 && rows > std::numeric_limits<std::size_t>::max() /
        (12u * sizeof(double))
      ? std::numeric_limits<std::size_t>::max()
      : rows * 12u * sizeof(double);
  const std::size_t capByMemory = bytesPerColumn == 0
    ? static_cast<std::size_t>(options.blockSize)
    : maxBytes / bytesPerColumn;
  if(capByMemory == 0) {
    result.metrics.blockFallbacks = 1;
    result.metrics.fallbackReason = "workspace-cap";
    Result fallback = solveIndependent(rhs, applyAVector, applyMVector,
                                       options, initialGuess);
    fallback.metrics.add(result.metrics);
    return fallback;
  }
  const Eigen::Index panelLimit = std::max<Eigen::Index>(1,
    std::min<Eigen::Index>(options.blockSize,
      static_cast<Eigen::Index>(std::min<std::size_t>(
        capByMemory, static_cast<std::size_t>(
          std::numeric_limits<Eigen::Index>::max())))));

  for(Eigen::Index first = 0; first < rhs.cols(); first += panelLimit) {
    const Eigen::Index count = std::min(panelLimit, rhs.cols() - first);
    const Eigen::MatrixXd panel = rhs.middleCols(first, count);
    Eigen::MatrixXd guess;
    const bool haveGuess = initialGuess != 0 &&
      initialGuess->rows() == rhs.rows() && initialGuess->cols() == rhs.cols();
    if(haveGuess) { guess = initialGuess->middleCols(first, count); }
    try {
      Result panelResult = solveBlockPanel(
        panel, applyA, applyM, options, haveGuess ? &guess : 0);
      result.solution.middleCols(first, count) = panelResult.solution;
      result.metrics.add(panelResult.metrics);
    } catch(const BlockBreakdown & error) {
      result.metrics.add(error.metrics);
      result.metrics.blockFallbacks++;
      if(result.metrics.fallbackReason.empty()) {
        result.metrics.fallbackReason = error.what();
      }
      Result fallback = solveIndependent(
        panel, applyAVector, applyMVector, options, haveGuess ? &guess : 0);
      result.solution.middleCols(first, count) = fallback.solution;
      result.metrics.add(fallback.metrics);
    }
  }
  return result;
}

}  // namespace pcg
}  // namespace sommer

#endif
