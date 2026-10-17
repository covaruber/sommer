#ifndef SOMMER_MME_POSTERIOR_H
#define SOMMER_MME_POSTERIOR_H

#include <Eigen/Core>
#include <stdexcept>

namespace sommer {

// ML conditions on fixed effects; only random coefficients have uncertainty.
// Keep joint coefficient coordinates so all moment contractions share layouts.
template<class JointSolve, class RandomSolve>
Eigen::MatrixXd posteriorCoefficientSolve(
    const Eigen::Ref<const Eigen::MatrixXd> & rhs, Eigen::Index fixed,
    bool reml, const JointSolve & jointSolve, const RandomSolve & randomSolve) {
  if(fixed < 0 || fixed > rhs.rows()) {
    throw std::runtime_error("Invalid fixed-effect dimension in posterior solve.");
  }
  if(reml) { return jointSolve(rhs); }
  Eigen::MatrixXd result = Eigen::MatrixXd::Zero(rhs.rows(), rhs.cols());
  const Eigen::Index random = rhs.rows() - fixed;
  if(random > 0) {
    // CHOLMOD dense views require packed panels, not strided coefficient slices.
    const Eigen::MatrixXd randomRhs = rhs.bottomRows(random);
    result.bottomRows(random) = randomSolve(randomRhs);
  }
  return result;
}

} // namespace sommer

#endif
