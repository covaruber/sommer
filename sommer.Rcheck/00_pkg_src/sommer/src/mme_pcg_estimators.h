#ifndef SOMMER_MME_PCG_ESTIMATORS_H
#define SOMMER_MME_PCG_ESTIMATORS_H

#include <Eigen/Cholesky>
#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include <Eigen/QR>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <functional>
#include <limits>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

namespace sommer {
namespace pcg_estimators {

typedef std::function<Eigen::VectorXd(const Eigen::VectorXd &)> VectorOperator;

enum class LogDetScaling { None, Diagonal };

struct LanczosLogDetOptions {
  LanczosLogDetOptions()
      : seed(0), probes(16), maxSteps(32), scaling(LogDetScaling::None) {}
  unsigned seed;
  Eigen::Index probes;
  Eigen::Index maxSteps;
  LogDetScaling scaling;
};

struct LogDetResult {
  LogDetResult()
      : estimate(0.0), probeSE(std::numeric_limits<double>::quiet_NaN()),
        operatorApplications(0) {}
  double estimate;
  double probeSE;
  std::vector<double> probeEstimates;
  Eigen::Index operatorApplications;
};

struct SensitivityAssessment {
  SensitivityAssessment()
      : referenceEstimate(0.0), comparisonEstimate(0.0),
        absoluteDifference(0.0), empiricalSE(std::numeric_limits<double>::quiet_NaN()) {}
  double referenceEstimate;
  double comparisonEstimate;
  double absoluteDifference;
  double empiricalSE;
};

struct LandmarkControlVariate {
  LandmarkControlVariate() : operatorApplications(0) {}
  Eigen::MatrixXd basis;
  Eigen::MatrixXd projectedInverse;
  Eigen::Index operatorApplications;
};

struct TraceOptions {
  TraceOptions() : seed(0), probes(32) {}
  unsigned seed;
  Eigen::Index probes;
};

struct TraceResult {
  TraceResult()
      : estimate(0.0), probeSE(std::numeric_limits<double>::quiet_NaN()),
        inverseApplications(0), derivativeApplications(0),
        controlOperatorApplications(0) {}
  double estimate;
  double probeSE;
  std::vector<double> probeEstimates;
  Eigen::Index inverseApplications;
  Eigen::Index derivativeApplications;
  Eigen::Index controlOperatorApplications;
};

inline std::uint32_t domainSeed(unsigned seed, std::uint32_t domain) {
  std::uint32_t value = static_cast<std::uint32_t>(seed) + domain;
  value = (value ^ (value >> 16)) * 0x7feb352dU;
  value = (value ^ (value >> 15)) * 0x846ca68bU;
  return value ^ (value >> 16);
}

inline void validateOperatorOutput(const Eigen::VectorXd & output,
                                   Eigen::Index dimension,
                                   const char * message) {
  if(output.size() != dimension || !output.allFinite()) {
    throw std::runtime_error(message);
  }
}

inline Eigen::VectorXd diagonalRoot(const Eigen::VectorXd & diagonal,
                                   Eigen::Index dimension) {
  if(diagonal.size() != dimension || !diagonal.allFinite() ||
     (diagonal.array() <= 0.0).any()) {
    throw std::runtime_error("Diagonal scaling requires a finite positive diagonal.");
  }
  return diagonal.array().sqrt().matrix();
}

inline VectorOperator diagonallyScaledOperator(
    Eigen::Index dimension, const VectorOperator & applyA,
    const Eigen::VectorXd & diagonal) {
  const Eigen::VectorXd root = diagonalRoot(diagonal, dimension);
  if(dimension < 0 || !applyA) {
    throw std::runtime_error("Invalid diagonal-scaled operator dimensions or callback.");
  }
  return [dimension, root, applyA](const Eigen::VectorXd & value) -> Eigen::VectorXd {
    if(value.size() != dimension || !value.allFinite()) {
      throw std::runtime_error("Invalid input to diagonal-scaled operator.");
    }
    const Eigen::VectorXd applied = applyA(value.cwiseQuotient(root));
    validateOperatorOutput(applied, dimension, "Invalid matrix action in diagonal scaling.");
    return applied.cwiseQuotient(root);
  };
}

inline Eigen::VectorXd mapInverseGuessToScaled(
    Eigen::Index dimension, const Eigen::VectorXd & originalGuess,
    const Eigen::VectorXd & diagonal) {
  const Eigen::VectorXd root = diagonalRoot(diagonal, dimension);
  if(originalGuess.size() != dimension || !originalGuess.allFinite()) {
    throw std::runtime_error("Invalid inverse guess for diagonal scaling.");
  }
  return originalGuess.cwiseProduct(root);
}

inline Eigen::VectorXd mapScaledSolutionToOriginal(
    Eigen::Index dimension, const Eigen::VectorXd & scaledSolution,
    const Eigen::VectorXd & diagonal) {
  const Eigen::VectorXd root = diagonalRoot(diagonal, dimension);
  if(scaledSolution.size() != dimension || !scaledSolution.allFinite()) {
    throw std::runtime_error("Invalid scaled solution for diagonal scaling.");
  }
  return scaledSolution.cwiseQuotient(root);
}

inline double logQuadraticForm(const Eigen::MatrixXd & tridiagonal) {
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(tridiagonal);
  if(solver.info() != Eigen::Success ||
     !solver.eigenvalues().allFinite() ||
     (solver.eigenvalues().array() <= 0.0).any()) {
    throw std::runtime_error("Lanczos log-determinant requires a positive-definite operator; minimum Ritz value is " +
      std::to_string(solver.eigenvalues().size() == 0 ? 0.0 : solver.eigenvalues().minCoeff()) + ".");
  }
  const Eigen::VectorXd weights = solver.eigenvectors().row(0).array().square().matrix().transpose();
  return (weights.array() * solver.eigenvalues().array().log()).sum();
}

inline double sampleSE(const std::vector<double> & values) {
  if(values.size() < 2) { return std::numeric_limits<double>::quiet_NaN(); }
  double mean = 0.0;
  for(std::size_t i = 0; i < values.size(); ++i) { mean += values[i]; }
  mean /= static_cast<double>(values.size());
  double sumSquares = 0.0;
  for(std::size_t i = 0; i < values.size(); ++i) {
    const double difference = values[i] - mean;
    sumSquares += difference * difference;
  }
  return std::sqrt(sumSquares /
                   (static_cast<double>(values.size()) *
                    static_cast<double>(values.size() - 1)));
}

inline double lanczosProbe(
    Eigen::Index dimension, const VectorOperator & applyA,
    Eigen::VectorXd residual, Eigen::Index steps,
    Eigen::Index & operatorApplications) {
  const double probeNormSquared = residual.squaredNorm();
  residual.normalize();
  std::vector<Eigen::VectorXd> basis;
  basis.reserve(static_cast<std::size_t>(steps));
  Eigen::VectorXd previous = Eigen::VectorXd::Zero(dimension);
  Eigen::VectorXd diagonal = Eigen::VectorXd::Zero(steps);
  Eigen::VectorXd offDiagonal = Eigen::VectorXd::Zero(std::max<Eigen::Index>(0, steps - 1));
  double previousBeta = 0.0;
  Eigen::Index actualSteps = 0;
  for(Eigen::Index step = 0; step < steps; ++step) {
    basis.push_back(residual);
    Eigen::VectorXd applied = applyA(residual);
    ++operatorApplications;
    validateOperatorOutput(applied, dimension, "Invalid matrix action in Lanczos log-determinant.");
    const double appliedNorm = applied.norm();
    if(step > 0) { applied.noalias() -= previousBeta * previous; }
    const double alpha = residual.dot(applied);
    applied.noalias() -= alpha * residual;
    for(std::size_t i = 0; i < basis.size(); ++i) {
      applied.noalias() -= basis[i].dot(applied) * basis[i];
    }
    diagonal(step) = alpha;
    actualSteps = step + 1;
    if(step + 1 == steps) { break; }
    const double beta = applied.norm();
    if(!std::isfinite(beta)) {
      throw std::runtime_error("Non-finite Lanczos recurrence in log-determinant.");
    }
    if(beta <= std::sqrt(std::numeric_limits<double>::epsilon()) *
                 std::max(1.0, appliedNorm)) {
      break;
    }
    offDiagonal(step) = beta;
    previous = residual;
    residual = applied / beta;
    previousBeta = beta;
  }
  Eigen::MatrixXd tridiagonal = Eigen::MatrixXd::Zero(actualSteps, actualSteps);
  tridiagonal.diagonal() = diagonal.head(actualSteps);
  for(Eigen::Index i = 0; i + 1 < actualSteps; ++i) {
    tridiagonal(i, i + 1) = offDiagonal(i);
    tridiagonal(i + 1, i) = offDiagonal(i);
  }
  return probeNormSquared * logQuadraticForm(tridiagonal);
}

// SLQ assumes a symmetric positive-definite operator; sampled Ritz checks do not certify SPD.
inline LogDetResult estimateLogDet(
    Eigen::Index dimension, const VectorOperator & applyA,
    const Eigen::VectorXd & diagonal, const LanczosLogDetOptions & options) {
  if(dimension < 0 || !applyA || options.probes < 1 ||
     (dimension > 0 && options.maxSteps < 1)) {
    throw std::runtime_error("Invalid Lanczos log-determinant dimensions or settings.");
  }
  LogDetResult result;
  if(dimension == 0) {
    if(diagonal.size() != 0) {
      throw std::runtime_error("Invalid diagonal for empty log-determinant operator.");
    }
    return result;
  }
  VectorOperator scaled = applyA;
  double diagonalLogDet = 0.0;
  if(options.scaling == LogDetScaling::Diagonal) {
    diagonalRoot(diagonal, dimension);
    scaled = diagonallyScaledOperator(dimension, applyA, diagonal);
    diagonalLogDet = diagonal.array().log().sum();
  } else if(diagonal.size() != 0 && diagonal.size() != dimension) {
    throw std::runtime_error("Invalid diagonal for log-determinant operator.");
  }
  const Eigen::Index steps = std::min(dimension, options.maxSteps);
  std::mt19937 generator(domainSeed(options.seed, 0x53a9e2d1U));
  std::uniform_int_distribution<int> sign(0, 1);
  result.probeEstimates.reserve(static_cast<std::size_t>(options.probes));
  for(Eigen::Index probe = 0; probe < options.probes; ++probe) {
    Eigen::VectorXd vector(dimension);
    for(Eigen::Index i = 0; i < dimension; ++i) {
      vector(i) = sign(generator) == 0 ? -1.0 : 1.0;
    }
    result.probeEstimates.push_back(
      lanczosProbe(dimension, scaled, vector, steps, result.operatorApplications) +
      diagonalLogDet);
  }
  for(std::size_t i = 0; i < result.probeEstimates.size(); ++i) {
    result.estimate += result.probeEstimates[i];
  }
  result.estimate /= static_cast<double>(result.probeEstimates.size());
  result.probeSE = sampleSE(result.probeEstimates);
  return result;
}

// These paired differences are empirical diagnostics, not total-error bounds.
inline SensitivityAssessment assessProbeSensitivity(
    const LogDetResult & estimate, Eigen::Index nestedProbeCount) {
  if(nestedProbeCount < 1 ||
     static_cast<std::size_t>(nestedProbeCount) > estimate.probeEstimates.size()) {
    throw std::runtime_error("Invalid nested probe count for log-determinant sensitivity.");
  }
  SensitivityAssessment result;
  result.referenceEstimate = estimate.estimate;
  for(Eigen::Index i = 0; i < nestedProbeCount; ++i) {
    result.comparisonEstimate += estimate.probeEstimates[static_cast<std::size_t>(i)];
  }
  result.comparisonEstimate /= static_cast<double>(nestedProbeCount);
  result.absoluteDifference = std::abs(result.referenceEstimate - result.comparisonEstimate);
  return result;
}

inline SensitivityAssessment assessDepthSensitivity(
    Eigen::Index dimension, const VectorOperator & applyA,
    const Eigen::VectorXd & diagonal, const LanczosLogDetOptions & options,
    Eigen::Index nestedSteps) {
  if(nestedSteps < 1 || nestedSteps >= options.maxSteps) {
    throw std::runtime_error("Nested Lanczos depth must be positive and below maximum depth.");
  }
  const LogDetResult full = estimateLogDet(dimension, applyA, diagonal, options);
  LanczosLogDetOptions nested = options;
  nested.maxSteps = nestedSteps;
  const LogDetResult shallow = estimateLogDet(dimension, applyA, diagonal, nested);
  SensitivityAssessment result;
  result.referenceEstimate = full.estimate;
  result.comparisonEstimate = shallow.estimate;
  result.absoluteDifference = std::abs(full.estimate - shallow.estimate);
  std::vector<double> paired(full.probeEstimates.size(), 0.0);
  for(std::size_t i = 0; i < paired.size(); ++i) {
    paired[i] = full.probeEstimates[i] - shallow.probeEstimates[i];
  }
  result.empiricalSE = sampleSE(paired);
  return result;
}

inline LandmarkControlVariate makeLandmarkControlVariate(
    Eigen::Index dimension, const VectorOperator & applyA,
    Eigen::Index rank, unsigned seed = 0,
    std::size_t maxBytes = 64U * 1024U * 1024U) {
  if(dimension < 1 || !applyA || rank < 1 || rank > dimension) {
    throw std::runtime_error("Invalid landmark control-variate dimensions or rank.");
  }
  const std::size_t n = static_cast<std::size_t>(dimension);
  const std::size_t r = static_cast<std::size_t>(rank);
  if(n > (std::numeric_limits<std::size_t>::max() / r) ||
     n * r > (std::numeric_limits<std::size_t>::max() - r * r) ||
     n * r + r * r > std::numeric_limits<std::size_t>::max() / 2 ||
     2 * (n * r + r * r) > maxBytes / sizeof(double)) {
    throw std::runtime_error("Landmark control-variate rank exceeds its memory budget.");
  }
  std::mt19937 generator(domainSeed(seed, 0x16f11fe9U));
  std::uniform_int_distribution<int> sign(0, 1);
  LandmarkControlVariate control;
  {
    Eigen::MatrixXd randomBasis(dimension, rank);
    for(Eigen::Index column = 0; column < rank; ++column) {
      for(Eigen::Index row = 0; row < dimension; ++row) {
        randomBasis(row, column) = sign(generator) == 0 ? -1.0 : 1.0;
      }
    }
    Eigen::HouseholderQR<Eigen::MatrixXd> qr(randomBasis);
    randomBasis.resize(0, 0);
    control.basis = qr.householderQ() * Eigen::MatrixXd::Identity(dimension, rank);
  }
  Eigen::MatrixXd projected = Eigen::MatrixXd::Zero(rank, rank);
  for(Eigen::Index column = 0; column < rank; ++column) {
    const Eigen::VectorXd basisColumn = control.basis.col(column);
    const Eigen::VectorXd applied = applyA(basisColumn);
    ++control.operatorApplications;
    validateOperatorOutput(applied, dimension, "Invalid matrix action in landmark construction.");
    projected.col(column) = control.basis.transpose() * applied;
  }
  projected = 0.5 * (projected + projected.transpose()).eval();
  Eigen::LLT<Eigen::MatrixXd> factor(projected);
  if(factor.info() != Eigen::Success) {
    throw std::runtime_error("Landmark projection is not positive definite.");
  }
  control.projectedInverse = factor.solve(Eigen::MatrixXd::Identity(rank, rank));
  if(factor.info() != Eigen::Success || !control.projectedInverse.allFinite()) {
    throw std::runtime_error("Failed to invert landmark projection.");
  }
  return control;
}

// The fixed K=Q(Q'AQ)^-1Q' correction is unbiased for symmetric or indefinite B.
inline TraceResult estimateInverseTrace(
    Eigen::Index dimension, const VectorOperator & applyAinv,
    const VectorOperator & applyB, const TraceOptions & options,
    const LandmarkControlVariate * control = nullptr) {
  if(dimension < 0 || !applyAinv || !applyB || options.probes < 1) {
    throw std::runtime_error("Invalid inverse trace dimensions or settings.");
  }
  if(control && (control->basis.rows() != dimension ||
                 control->projectedInverse.rows() != control->basis.cols() ||
                 control->projectedInverse.cols() != control->basis.cols() ||
                 !control->basis.allFinite() || !control->projectedInverse.allFinite())) {
    throw std::runtime_error("Invalid landmark control-variate dimensions.");
  }
  TraceResult result;
  if(control) { result.controlOperatorApplications = control->operatorApplications; }
  std::mt19937 generator(domainSeed(options.seed, 0xa4b6c7d3U));
  std::uniform_int_distribution<int> sign(0, 1);
  double controlTrace = 0.0;
  if(control) {
    Eigen::MatrixXd projectedB(control->basis.cols(), control->basis.cols());
    for(Eigen::Index column = 0; column < control->basis.cols(); ++column) {
      const Eigen::VectorXd basisColumn = control->basis.col(column);
      const Eigen::VectorXd applied = applyB(basisColumn);
      ++result.derivativeApplications;
      validateOperatorOutput(applied, dimension, "Invalid derivative action in control-variate trace.");
      projectedB.col(column) = control->basis.transpose() * applied;
    }
    controlTrace = (control->projectedInverse.transpose().cwiseProduct(projectedB)).sum();
  }
  result.probeEstimates.reserve(static_cast<std::size_t>(options.probes));
  for(Eigen::Index probe = 0; probe < options.probes; ++probe) {
    Eigen::VectorXd vector(dimension);
    for(Eigen::Index i = 0; i < dimension; ++i) {
      vector(i) = sign(generator) == 0 ? -1.0 : 1.0;
    }
    const Eigen::VectorXd appliedB = applyB(vector);
    ++result.derivativeApplications;
    validateOperatorOutput(appliedB, dimension, "Invalid derivative action in inverse trace.");
    const Eigen::VectorXd inverseAppliedB = applyAinv(appliedB);
    ++result.inverseApplications;
    validateOperatorOutput(inverseAppliedB, dimension, "Invalid inverse action in inverse trace.");
    double estimate = vector.dot(inverseAppliedB);
    if(control) {
      const Eigen::VectorXd projected = control->basis *
        (control->projectedInverse * (control->basis.transpose() * appliedB));
      estimate -= vector.dot(projected);
      estimate += controlTrace;
    }
    result.probeEstimates.push_back(estimate);
  }
  for(std::size_t i = 0; i < result.probeEstimates.size(); ++i) {
    result.estimate += result.probeEstimates[i];
  }
  result.estimate /= static_cast<double>(result.probeEstimates.size());
  result.probeSE = sampleSE(result.probeEstimates);
  return result;
}

} // namespace pcg_estimators
} // namespace sommer

#endif
