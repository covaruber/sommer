#ifndef SOMMER_MME_LAYOUT_H
#define SOMMER_MME_LAYOUT_H

#include "mme_sp2_support.h"

namespace sommer {

inline void assembleMMEDesign(Eigen::SparseMatrix<double> & design,
                              const arma::sp_mat & fixed,
                              const arma::field<arma::sp_mat> & random,
                              int randomCount) {
  std::vector<Eigen::Triplet<double>> entries;
  const Eigen::SparseMatrix<double> fixedEigen = armaSparseToEigenGlobal(fixed);
  entries.reserve(static_cast<std::size_t>(fixedEigen.nonZeros()));
  for(int col = 0; col < fixedEigen.outerSize(); ++col){
    for(Eigen::SparseMatrix<double>::InnerIterator it(fixedEigen, col); it; ++it){
      entries.emplace_back(static_cast<int>(it.row()),
                           static_cast<int>(it.col()), it.value());
    }
  }
  int columnOffset = static_cast<int>(fixedEigen.cols());
  for(int term = 0; term < randomCount; ++term){
    const Eigen::SparseMatrix<double> randomEigen = armaSparseToEigenGlobal(random(term));
    entries.reserve(entries.size() + static_cast<std::size_t>(randomEigen.nonZeros()));
    for(int col = 0; col < randomEigen.outerSize(); ++col){
      for(Eigen::SparseMatrix<double>::InnerIterator it(randomEigen, col); it; ++it){
        entries.emplace_back(static_cast<int>(it.row()),
                             columnOffset + static_cast<int>(it.col()), it.value());
      }
    }
    columnOffset += static_cast<int>(randomEigen.cols());
  }
  design.setFromTriplets(entries.begin(), entries.end());
  design.makeCompressed();
}

} // namespace sommer

#endif
