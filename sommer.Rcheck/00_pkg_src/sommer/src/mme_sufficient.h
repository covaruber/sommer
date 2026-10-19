#ifndef SOMMER_MME_SUFFICIENT_H
#define SOMMER_MME_SUFFICIENT_H

#include <Eigen/SparseCore>
#include <Eigen/Core>
#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
#include <stdexcept>
#include <unordered_map>
#include <utility>
#include <vector>

namespace sommer {

class SufficientStatistics {
public:
  typedef Eigen::SparseMatrix<double, Eigen::RowMajor, int> RowMatrix;
  typedef Eigen::SparseMatrix<double, Eigen::ColMajor, int> SparseMatrix;
  struct Group {
    int residualClass;
    std::size_t count;
    double weight, mean, centeredSquares;
    std::vector<std::pair<int, double>> row;
    Group(int category)
      : residualClass(category), count(0), weight(0), mean(0), centeredSquares(0) {}
    void add(double response, double precision) {
      const double total = weight + precision;
      const double delta = response - mean;
      mean += delta * precision / total;
      centeredSquares += precision * delta * (response - mean);
      weight = total;
      ++count;
      if(!std::isfinite(weight) || !std::isfinite(mean) || !std::isfinite(centeredSquares)){
        throw std::runtime_error("Sufficient statistics overflowed.");
      }
    }
    double squares() const { return centeredSquares + weight * mean * mean; }
  };

  std::vector<Group> groups;
  RowMatrix design;
  std::vector<std::size_t> classCounts;
  std::vector<SparseMatrix> classCross;
  Eigen::MatrixXd classRhs;
  Eigen::VectorXd classSquares;
  std::size_t retainedBytes;
  double setupSeconds;

  SufficientStatistics() : retainedBytes(0), setupSeconds(0) {}

  void build(const SparseMatrix & input, const Eigen::VectorXd & response,
             const std::vector<int> & classes, const Eigen::VectorXd & precision,
             int classCount, bool cacheClasses, std::size_t maxBytes) {
    if(input.rows() != response.size() || input.rows() != precision.size() ||
       classes.size() != static_cast<std::size_t>(input.rows()) || classCount < 1 ||
       !response.allFinite() || !precision.allFinite() || (precision.array() <= 0).any()){
      throw std::invalid_argument("Invalid sufficient-statistic inputs.");
    }
    groups.clear();
    classCross.clear();
    classCounts.assign(classCount, 0);
    const RowMatrix rows = input;
    std::unordered_map<std::size_t, std::vector<std::size_t>> buckets;
    std::size_t entries = 0;
    for(int index = 0; index < rows.rows(); ++index){
      const int category = classes[index];
      if(category < 0 || category >= classCount){
        throw std::invalid_argument("Invalid sufficient-statistic residual class.");
      }
      std::vector<std::pair<int, double>> row;
      std::size_t hash = std::hash<int>()(category);
      auto combine = [&](std::size_t value){ hash ^= value + 0x9e3779b9U + (hash << 6) + (hash >> 2); };
      for(RowMatrix::InnerIterator entry(rows, index); entry; ++entry){
        if(!std::isfinite(entry.value())){ throw std::invalid_argument("Non-finite compressed design."); }
        if(entry.value() == 0.0){ continue; }
        row.emplace_back(entry.col(), entry.value());
        combine(std::hash<int>()(entry.col()));
        combine(std::hash<double>()(entry.value()));
      }
      std::vector<std::size_t> & bucket = buckets[hash];
      std::size_t found = groups.size();
      for(std::size_t candidate : bucket){
        if(groups[candidate].residualClass == category && groups[candidate].row == row){
          found = candidate;
          break;
        }
      }
      if(found == groups.size()){
        entries += row.size();
        const long double needed = static_cast<long double>(entries) * 64 +
          static_cast<long double>(groups.size() + 1) * 256 +
          static_cast<long double>(input.nonZeros()) * 16 + (input.rows() + 1.0L) * sizeof(int);
        if(needed > maxBytes){ throw std::runtime_error("Sufficient-statistic grouping exceeds denseMemoryMB."); }
        groups.emplace_back(category);
        groups.back().row = std::move(row);
        bucket.push_back(found);
      }
      groups[found].add(response(index), precision(index));
      ++classCounts[category];
    }
    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(entries);
    for(std::size_t group = 0; group < groups.size(); ++group){
      for(const auto & entry : groups[group].row){ triplets.emplace_back(group, entry.first, entry.second); }
    }
    design.resize(groups.size(), input.cols());
    design.setFromTriplets(triplets.begin(), triplets.end());
    retainedBytes = entries * (sizeof(double) + sizeof(int) + sizeof(std::pair<int,double>)) +
      groups.size() * (sizeof(Group) + sizeof(int)) + classCounts.size() * sizeof(std::size_t);
    if(cacheClasses){
      const long double denseNeeded = static_cast<long double>(input.cols()) * classCount * sizeof(double);
      long double productEntries = 0;
      for(const Group & group : groups){ productEntries += static_cast<long double>(group.row.size()) * group.row.size(); }
      if(retainedBytes + denseNeeded + productEntries * 64 > maxBytes){
        throw std::runtime_error("Sufficient-statistic class cross-products exceed denseMemoryMB.");
      }
      classRhs = Eigen::MatrixXd::Zero(input.cols(), classCount);
      classSquares = Eigen::VectorXd::Zero(classCount);
      classCross.resize(classCount);
      for(int category = 0; category < classCount; ++category){
        triplets.clear();
        for(const Group & group : groups){
          if(group.residualClass != category) continue;
          classSquares(category) += group.squares();
          for(const auto & a : group.row){
            classRhs(a.first, category) += a.second * group.weight * group.mean;
            for(const auto & b : group.row){
              triplets.emplace_back(a.first, b.first, group.weight * a.second * b.second);
            }
          }
        }
        classCross[category].resize(input.cols(), input.cols());
        classCross[category].setFromTriplets(triplets.begin(), triplets.end());
        retainedBytes += classCross[category].nonZeros() * (sizeof(double) + sizeof(int)) +
          (input.cols() + 1) * sizeof(int);
      }
      retainedBytes += static_cast<std::size_t>(denseNeeded) + classCount * sizeof(double);
    }
  }

  void assemble(const Eigen::VectorXd & variance, SparseMatrix & cross,
                Eigen::VectorXd & rhs, double & quadratic) const {
    if(variance.size() != static_cast<Eigen::Index>(classCounts.size()) ||
       !variance.allFinite() || (variance.array() <= 0).any()){
      throw std::invalid_argument("Invalid compressed residual variances.");
    }
    rhs = Eigen::VectorXd::Zero(design.cols());
    quadratic = 0.0;
    if(!classCross.empty()){
      cross.resize(design.cols(), design.cols());
      cross.setZero();
      for(int category = 0; category < variance.size(); ++category){
        cross += classCross[category] * (1.0 / variance(category));
        rhs.noalias() += classRhs.col(category) / variance(category);
        quadratic += classSquares(category) / variance(category);
      }
    }else{
      Eigen::VectorXd diagonal(groups.size()), weightedResponse(groups.size());
      for(std::size_t index = 0; index < groups.size(); ++index){
        const Group & group = groups[index];
        diagonal(index) = group.weight / variance(group.residualClass);
        weightedResponse(index) = diagonal(index) * group.mean;
        quadratic += group.squares() / variance(group.residualClass);
      }
      const SparseMatrix weighted = diagonal.asDiagonal() * design;
      cross = SparseMatrix(design.transpose()) * weighted;
      rhs.noalias() = design.transpose() * weightedResponse;
    }
    cross.makeCompressed();
  }

  Eigen::VectorXd residualSquares(const Eigen::VectorXd & fitted) const {
    if(fitted.size() != static_cast<Eigen::Index>(groups.size()) || !fitted.allFinite()){
      throw std::invalid_argument("Invalid compressed fitted values.");
    }
    Eigen::VectorXd result(groups.size());
    for(std::size_t index = 0; index < groups.size(); ++index){
      const Group & group = groups[index];
      const double difference = group.mean - fitted(index);
      result(index) = group.centeredSquares + group.weight * difference * difference;
    }
    return result;
  }

  SparseMatrix derivativeCross(const Eigen::VectorXd & scales) const {
    if(scales.size() != static_cast<Eigen::Index>(classCounts.size()) || !scales.allFinite()){
      throw std::invalid_argument("Invalid compressed derivative scales.");
    }
    SparseMatrix result(design.cols(), design.cols());
    if(!classCross.empty()){
      for(int category = 0; category < scales.size(); ++category){
        if(scales(category) != 0.0){ result += classCross[category] * scales(category); }
      }
    }else{
      Eigen::VectorXd weights(groups.size());
      for(std::size_t group = 0; group < groups.size(); ++group){
        weights(group) = groups[group].weight * scales(groups[group].residualClass);
      }
      const SparseMatrix weighted = weights.asDiagonal() * design;
      result = SparseMatrix(design.transpose()) * weighted;
    }
    result.makeCompressed();
    return result;
  }
};

} // namespace sommer
#endif
