#ifndef SOMMER_MME_EXACT_H
#define SOMMER_MME_EXACT_H

#include <Eigen/SparseCore>
#include <algorithm>
#include <functional>
#include <stdexcept>
#include <vector>

namespace sommer {

typedef std::function<Eigen::MatrixXd(const Eigen::MatrixXd &)> DensePanelOperator;

// Contract inverse coefficient blocks against one shared relationship precision.
inline Eigen::MatrixXd exactInverseBlockTraces(
    Eigen::Index dimension, const std::vector<Eigen::Index> & starts,
    const Eigen::SparseMatrix<double> & relationship, Eigen::Index panelColumns,
    const DensePanelOperator & solve,
    const std::function<bool(std::size_t, std::size_t)> & needed,
    const std::function<bool(Eigen::Index, Eigen::Index, double &)> & lookup = nullptr) {
  const Eigen::Index width = relationship.rows();
  if(width != relationship.cols() || (!starts.empty() && width == 0) ||
     panelColumns < 1 || dimension < 0) {
    throw std::runtime_error("Invalid inverse block trace dimensions.");
  }
  for(Eigen::Index start : starts) {
    if(start < 0 || start > dimension || width > dimension - start) {
      throw std::runtime_error("Inverse block trace exceeds coefficient dimensions.");
    }
  }
  Eigen::MatrixXd traces = Eigen::MatrixXd::Zero(starts.size(), starts.size());
  for(std::size_t column = 0; column < starts.size(); ++column) {
    std::vector<bool> missing(column + 1, false);
    bool active = false;
    for(std::size_t row = 0; row <= column; ++row) {
      if(!needed(row, column)) { continue; }
      bool complete = static_cast<bool>(lookup);
      double trace = 0.0;
      if(complete) {
        for(Eigen::Index local = 0; local < width && complete; ++local) {
          for(Eigen::SparseMatrix<double>::InnerIterator entry(relationship, local); entry; ++entry) {
            double value = 0.0;
            if(!lookup(starts[row] + entry.row(), starts[column] + local, value)) {
              complete = false;
              break;
            }
            trace += entry.value() * value;
          }
        }
      }
      if(complete) { traces(row, column) = trace; }
      else { missing[row] = true; active = true; }
    }
    if(!active) { continue; }
    for(Eigen::Index first = 0; first < width; first += panelColumns) {
      const Eigen::Index count = std::min(panelColumns, width - first);
      Eigen::MatrixXd rhs = Eigen::MatrixXd::Zero(dimension, count);
      for(Eigen::Index local = 0; local < count; ++local) {
        rhs(starts[column] + first + local, local) = 1.0;
      }
      const Eigen::MatrixXd solved = solve(rhs);
      if(solved.rows() != dimension || solved.cols() != count || !solved.allFinite()) {
        throw std::runtime_error("Invalid inverse action in block trace panels.");
      }
      for(std::size_t row = 0; row <= column; ++row) {
        if(!missing[row]) { continue; }
        for(Eigen::Index local = 0; local < count; ++local) {
          for(Eigen::SparseMatrix<double>::InnerIterator entry(relationship, first + local); entry; ++entry) {
            traces(row, column) += entry.value() * solved(starts[row] + entry.row(), local);
          }
        }
      }
    }
  }
  return Eigen::MatrixXd(traces.selfadjointView<Eigen::Upper>());
}

// Share identity solves across derivative contractions without retaining an inverse.
inline Eigen::VectorXd exactInverseTraces(
    Eigen::Index dimension, Eigen::Index derivatives, Eigen::Index panelColumns,
    const DensePanelOperator & solve,
    const std::function<Eigen::MatrixXd(Eigen::Index, const Eigen::MatrixXd &)> & derivative) {
  if(dimension < 0 || derivatives < 0 || panelColumns < 1) {
    throw std::runtime_error("Invalid exact trace panel dimensions.");
  }
  Eigen::VectorXd traces = Eigen::VectorXd::Zero(derivatives);
  if(derivatives == 0) { return traces; }
  for(Eigen::Index first = 0; first < dimension; first += panelColumns) {
    const Eigen::Index count = std::min(panelColumns, dimension - first);
    Eigen::MatrixXd basis = Eigen::MatrixXd::Zero(dimension, count);
    for(Eigen::Index column = 0; column < count; ++column) {
      basis(first + column, column) = 1.0;
    }
    const Eigen::MatrixXd inverse = solve(basis);
    if(inverse.rows() != dimension || inverse.cols() != count || !inverse.allFinite()) {
      throw std::runtime_error("Invalid inverse action in exact trace panels.");
    }
    for(Eigen::Index parameter = 0; parameter < derivatives; ++parameter) {
      const Eigen::MatrixXd applied = derivative(parameter, inverse);
      if(applied.rows() != dimension || applied.cols() != count || !applied.allFinite()) {
        throw std::runtime_error("Invalid derivative action in exact trace panels.");
      }
      for(Eigen::Index column = 0; column < count; ++column) {
        traces(parameter) += applied(first + column, column);
      }
    }
  }
  return traces;
}

// Materialize compatibility output only; fitting backends can retain operators.
inline Eigen::SparseMatrix<double> assembleNormalEquations(
    const Eigen::SparseMatrix<double> & design,
    const Eigen::SparseMatrix<double> & prior,
    const DensePanelOperator & precision, Eigen::Index panelColumns,
    std::size_t maxEntries) {
  if(panelColumns < 1 || prior.rows() != design.cols() || prior.cols() != design.cols()) {
    throw std::runtime_error("Invalid normal-equation output dimensions.");
  }
  std::vector<Eigen::Triplet<double>> entries;
  for(Eigen::Index first = 0; first < design.cols(); first += panelColumns) {
    const Eigen::Index count = std::min(panelColumns, design.cols() - first);
    const Eigen::MatrixXd rhs = Eigen::MatrixXd(design.middleCols(first, count));
    const Eigen::MatrixXd applied = precision(rhs);
    if(applied.rows() != design.rows() || applied.cols() != count || !applied.allFinite()) {
      throw std::runtime_error("Invalid precision action in normal-equation output.");
    }
    const Eigen::MatrixXd cross = design.transpose() * applied;
    for(Eigen::Index column = 0; column < count; ++column) {
      for(Eigen::Index row = 0; row < design.cols(); ++row) {
        if(cross(row, column) != 0.0) {
          if(entries.size() >= maxEntries) {
            throw std::runtime_error("Final compatibility C exceeds the sparse entry limit.");
          }
          entries.emplace_back(row, first + column, cross(row, column));
        }
      }
    }
  }
  Eigen::SparseMatrix<double> result(design.cols(), design.cols());
  result.setFromTriplets(entries.begin(), entries.end());
  result += prior;
  if(static_cast<std::size_t>(result.nonZeros()) > maxEntries) {
    throw std::runtime_error("Final compatibility C exceeds the sparse entry limit.");
  }
  result.makeCompressed();
  return result;
}

} // namespace sommer

#endif
