#ifndef SOMMER_MME_PATHS_H
#define SOMMER_MME_PATHS_H

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
#include <utility>
#include <vector>

namespace sommer {

enum class ResidualPath { Diagonal, Grid, Kronecker, Sparse };
enum class AssemblyPath { Sparse, Blocks, Batched };
enum class SolverBackend { LDLT, Cholmod, PCG };
enum class TracePreparation { SelectedInverse, InverseBlocks, Probes, None };

struct TensorComponentPlan {
  bool eligible = false;
  std::size_t coordinates = 0;
  std::vector<std::vector<int>> groups;
};

class DisjointComponents {
  std::vector<int> parent;
  int root(int vertex) {
    while(parent[vertex] != vertex) {
      parent[vertex] = parent[parent[vertex]];
      vertex = parent[vertex];
    }
    return vertex;
  }
public:
  explicit DisjointComponents(std::size_t size) {
    if(size > static_cast<std::size_t>(std::numeric_limits<int>::max())) { return; }
    parent.resize(size);
    for(std::size_t vertex = 0; vertex < size; ++vertex) { parent[vertex] = static_cast<int>(vertex); }
  }
  bool unite(int first, int second) {
    if(first < 0 || second < 0 || static_cast<std::size_t>(first) >= parent.size() ||
       static_cast<std::size_t>(second) >= parent.size()) { return false; }
    const int firstRoot = root(first);
    const int secondRoot = root(second);
    parent[std::max(firstRoot, secondRoot)] = std::min(firstRoot, secondRoot);
    return true;
  }
  std::vector<int> labels() {
    std::vector<int> result(parent.size());
    std::vector<int> mapping(parent.size(), -1);
    int count = 0;
    for(std::size_t vertex = 0; vertex < parent.size(); ++vertex) {
      const int representative = root(static_cast<int>(vertex));
      if(mapping[representative] < 0) { mapping[representative] = count++; }
      result[vertex] = mapping[representative];
    }
    return result;
  }
};

inline TensorComponentPlan planTensorGraphComponents(
    const std::vector<std::size_t> & dimensions,
    const std::vector<std::vector<int>> & labels,
    std::size_t maxCoordinates) {
  TensorComponentPlan result;
  if(dimensions.size() != labels.size() || maxCoordinates == 0) { return result; }
  maxCoordinates = std::min(maxCoordinates,
    static_cast<std::size_t>(std::numeric_limits<int>::max()));
  std::size_t coordinates = 1;
  std::size_t components = 1;
  std::vector<std::vector<int>> canonical(labels.size());
  std::vector<std::size_t> counts(labels.size());
  for(std::size_t axis = 0; axis < dimensions.size(); ++axis) {
    const std::size_t dimension = dimensions[axis];
    if(dimension == 0 || coordinates > maxCoordinates / dimension ||
       labels[axis].size() != dimension) { return result; }
    coordinates *= dimension;
    std::vector<int> mapping(dimension, -1);
    canonical[axis].resize(dimension);
    for(std::size_t vertex = 0; vertex < dimension; ++vertex) {
      const int label = labels[axis][vertex];
      if(label < 0 || static_cast<std::size_t>(label) >= dimension) { return result; }
      if(mapping[label] < 0) { mapping[label] = static_cast<int>(counts[axis]++); }
      canonical[axis][vertex] = mapping[label];
    }
    components *= counts[axis];
  }
  result.coordinates = coordinates;
  result.groups.resize(components);
  for(std::size_t coordinate = 0; coordinate < coordinates; ++coordinate) {
    std::size_t trailing = coordinates;
    std::size_t component = 0;
    for(std::size_t axis = 0; axis < dimensions.size(); ++axis) {
      trailing /= dimensions[axis];
      component = component * counts[axis] +
        canonical[axis][(coordinate / trailing) % dimensions[axis]];
    }
    result.groups[component].push_back(static_cast<int>(coordinate));
  }
  result.eligible = true;
  return result;
}

inline TensorComponentPlan planTensorComponents(
    const std::vector<std::size_t> & dimensions,
    const std::vector<bool> & diagonalAxes,
    std::size_t maxCoordinates) {
  TensorComponentPlan result;
  if(dimensions.size() != diagonalAxes.size() || maxCoordinates == 0) {
    return result;
  }
  maxCoordinates = std::min(maxCoordinates,
    static_cast<std::size_t>(std::numeric_limits<int>::max()));
  std::size_t coordinates = 1;
  std::size_t components = 1;
  for(std::size_t axis = 0; axis < dimensions.size(); ++axis) {
    const std::size_t dimension = dimensions[axis];
    if(dimension == 0 || coordinates > maxCoordinates / dimension) {
      return result;
    }
    coordinates *= dimension;
    if(diagonalAxes[axis]) { components *= dimension; }
  }
  result.coordinates = coordinates;
  result.groups.resize(components);
  for(std::size_t coordinate = 0; coordinate < coordinates; ++coordinate) {
    std::size_t trailing = coordinates;
    std::size_t component = 0;
    for(std::size_t axis = 0; axis < dimensions.size(); ++axis) {
      trailing /= dimensions[axis];
      if(diagonalAxes[axis]) {
        component = component * dimensions[axis]
          + (coordinate / trailing) % dimensions[axis];
      }
    }
    result.groups[component].push_back(static_cast<int>(coordinate));
  }
  result.eligible = true;
  return result;
}

struct ResidualPolicy {
  static bool compactStorage(bool sectioned, bool nonDiagonalFactor,
                             bool completeHint, bool weighted) {
    return !(weighted && completeHint) && (sectioned || nonDiagonalFactor);
  }
  static bool gridCandidate(bool compact, bool complete, bool sectioned,
                            bool diagonal, bool weighted) {
    return compact && (!complete || sectioned) && (!diagonal || sectioned)
      && !weighted;
  }
  static bool gridFits(double rectangle, double observed) {
    return rectangle <= 3.0 * observed + 1000.0 && rectangle - observed <= 4000.0;
  }
};

struct SolverPlan {
  SolverBackend backend;
  bool reml;
  int inverseMode;
  bool diagonalPreconditioner;

  bool ldlt() const { return backend == SolverBackend::LDLT; }
  bool cholmod() const { return backend == SolverBackend::Cholmod; }
  bool pcg() const { return backend == SolverBackend::PCG; }
  bool matrixFreeEligible() const {
    return pcg() && diagonalPreconditioner && reml && inverseMode == 0;
  }
  bool matrixFree(bool hasRandomEffects, bool factorWise) const {
    return matrixFreeEligible() && hasRandomEffects && factorWise;
  }
  bool factorRandomOnly(bool hasRandomEffects) const {
    return !reml && hasRandomEffects;
  }
  bool dFirstML(std::size_t fixed, std::size_t random,
                double maxDoubles, bool retainJointCache) const {
    const double fixedSize = static_cast<double>(fixed);
    const double randomSize = static_cast<double>(random);
    return !pcg() && !reml && inverseMode != 1 && random > 0 && fixed <= 256
      && !retainJointCache && 3.0 * randomSize * fixedSize +
        4.0 * fixedSize * fixedSize <= maxDoubles;
  }
  TracePreparation traces() const {
    return ldlt() ? TracePreparation::SelectedInverse
      : pcg() ? TracePreparation::Probes
      : reml ? TracePreparation::InverseBlocks
      : TracePreparation::None;
  }
};

struct ResidualPlan {
  ResidualPath path;
  bool diagonalWeights;

  static ResidualPlan select(bool structurallyDiagonal, bool completeBlocks,
                             bool gridBlocks, bool useWeights,
                             bool weightsDiagonal) {
    const ResidualPath path = gridBlocks ? ResidualPath::Grid
      : structurallyDiagonal ? ResidualPath::Diagonal
      : completeBlocks ? ResidualPath::Kronecker
      : ResidualPath::Sparse;
    return {path, !useWeights || weightsDiagonal};
  }

  bool diagonal() const { return path == ResidualPath::Diagonal; }
  bool grid() const { return path == ResidualPath::Grid; }
  bool kronecker() const { return path == ResidualPath::Kronecker; }
  bool elementwisePrecision() const { return diagonal() && diagonalWeights; }
  bool blockAssembly() const {
    return (grid() || kronecker()) && diagonalWeights;
  }
  AssemblyPath assembly(bool sparseWeightFactors) const {
    return elementwisePrecision() || (diagonal() && sparseWeightFactors)
      ? AssemblyPath::Sparse
      : blockAssembly() ? AssemblyPath::Blocks
      : AssemblyPath::Batched;
  }
};

struct BlockSchurPolicy {
  std::size_t maxBorder = 2000;
  double maxDoubles = 1.25e8;
  double smallGroup = 64.0;
  double minDensity = 0.1;

  double storageDoubles(const std::vector<std::vector<int>> & groups,
                        std::size_t border) const {
    const double borderSize = static_cast<double>(border);
    double result = 2.0 * borderSize * borderSize;
    for(const auto & group : groups) {
      const double size = static_cast<double>(group.size());
      result += 2.0 * size * size + 3.0 * size * borderSize;
    }
    return result;
  }

  double cost(const std::vector<std::vector<int>> & groups,
              std::size_t border) const {
    if(groups.empty() || border > maxBorder ||
       storageDoubles(groups, border) > maxDoubles) {
      return std::numeric_limits<double>::infinity();
    }
    const double borderSize = static_cast<double>(border);
    double result = (5.0 / 3.0) * borderSize * borderSize * borderSize;
    for(const auto & group : groups) {
      const double size = static_cast<double>(group.size());
      result += (5.0 / 3.0) * size * size * size
        + 4.0 * borderSize * size * size
        + 2.0 * size * borderSize * borderSize;
    }
    return result;
  }

  bool needsDensityCheck(double effects, double group) const {
    return !(group <= smallGroup || effects * group <= maxDoubles);
  }
  bool inverseCacheFits(double effects, double group) const {
    return 2.0 * group * group + effects * group <= maxDoubles;
  }
  double chainStorageDoubles(const std::vector<std::vector<int>> & blocks,
                            const std::vector<int> & components,
                            std::size_t border) const {
    if(blocks.empty() || blocks.size() != components.size() || border > maxBorder) {
      return std::numeric_limits<double>::infinity();
    }
    const double borderSize = static_cast<double>(border);
    double result = 3.0 * borderSize * borderSize;
    for(std::size_t block = 0; block < blocks.size(); ++block) {
      const double size = static_cast<double>(blocks[block].size());
      result += 4.0 * size * size + 5.0 * size * borderSize;
      if(block > 0 && components[block] == components[block - 1]) {
        result += 3.0 * size * blocks[block - 1].size();
      }
    }
    return result;
  }
  double bandStorageDoubles(const std::vector<std::vector<int>> & blocks,
                            std::size_t bandwidth,
                            std::size_t border) const {
    if(blocks.empty() || bandwidth == 0 || bandwidth >= blocks.size() ||
       border > maxBorder) { return std::numeric_limits<double>::infinity(); }
    const double borderSize = static_cast<double>(border);
    double result = 2.0 * borderSize * borderSize;
    for(std::size_t block = 0; block < blocks.size(); ++block) {
      const double size = static_cast<double>(blocks[block].size());
      result += 3.0 * size * borderSize + size * size;
      for(std::size_t offset = 1; offset <= bandwidth && offset <= block; ++offset) {
        result += 2.0 * size * blocks[block - offset].size();
      }
    }
    return result;
  }
  double latentStorageDoubles(const std::vector<std::vector<int>> & blocks,
                             std::size_t border, std::size_t dimension) const {
    if(blocks.empty() || dimension == 0) {
      return std::numeric_limits<double>::infinity();
    }
    const double borderSize = static_cast<double>(border);
    const double mainSize = static_cast<double>(dimension);
    double result = 3.0 * borderSize * borderSize + mainSize * mainSize;
    for(const auto & block : blocks) {
      if(block.size() != dimension) { return std::numeric_limits<double>::infinity(); }
      result += 4.0 * mainSize * mainSize + 6.0 * mainSize * borderSize;
    }
    return result;
  }
  bool denseEnough(double nonzeros, double group) const {
    return nonzeros >= minDensity * group * group;
  }
};

struct BlockSchurPlan {
  bool eligible = false;
  std::vector<int> border;
  std::vector<std::vector<int>> groups;
};

struct BlockChainOrder {
  bool eligible = false;
  std::vector<int> blocks;
  std::vector<int> components;
};

struct BlockBandPlan {
  bool eligible = false;
  std::size_t bandwidth = 0;
};

struct BlockBandPolicy {
  std::size_t maxBandwidth = 4;
  std::size_t minBlocks = 40;
  static bool isAdvantageous(std::size_t blockCount, std::size_t bandwidth,
                             std::size_t borderSize, std::size_t totalEffects) {
    if(bandwidth > 4 || blockCount < 40 || totalEffects == 0) return false;
    const double borderRatio = static_cast<double>(borderSize) / static_cast<double>(totalEffects);
    return borderRatio <= 0.05;
  }
};

struct LatentSchurPolicy {
  static bool isAdvantageousMultiFactor(std::size_t lowRankDimension, std::size_t borderSize) {
    // Multi-factor latent Woodbury Schur is advantageous when the low-rank
    // dimension is sufficiently large (>= 8) and the virtual border is compact (<= 300).
    return (lowRankDimension >= 8 && borderSize <= 300);
  }
};

inline BlockBandPlan planBlockBand(
    std::size_t count, const std::vector<std::pair<int, int>> & edges,
    std::size_t maxBandwidth) {
  BlockBandPlan result;
  if(count == 0 || count > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
    return result;
  }
  for(const auto & edge : edges) {
    if(edge.first < 0 || edge.second < 0 || edge.first == edge.second ||
       static_cast<std::size_t>(edge.first) >= count ||
       static_cast<std::size_t>(edge.second) >= count) { return result; }
    const std::size_t distance = static_cast<std::size_t>(std::abs(edge.first - edge.second));
    result.bandwidth = std::max(result.bandwidth, distance);
    if(result.bandwidth > maxBandwidth) { return BlockBandPlan(); }
  }
  result.eligible = result.bandwidth > 0;
  return result;
}

inline BlockChainOrder orderBlockPaths(std::size_t count,
    const std::vector<std::pair<int, int>> & edges) {
  BlockChainOrder result;
  std::vector<std::vector<int>> neighbors(count);
  for(const auto & edge : edges) {
    if(edge.first < 0 || edge.second < 0 || edge.first == edge.second ||
       static_cast<std::size_t>(edge.first) >= count ||
       static_cast<std::size_t>(edge.second) >= count) { return result; }
    neighbors[edge.first].push_back(edge.second);
    neighbors[edge.second].push_back(edge.first);
  }
  for(auto & adjacent : neighbors) {
    std::sort(adjacent.begin(), adjacent.end());
    adjacent.erase(std::unique(adjacent.begin(), adjacent.end()), adjacent.end());
    if(adjacent.size() > 2) { return result; }
  }
  std::vector<bool> visited(count, false);
  int component = 0;
  for(std::size_t start = 0; start < count; ++start) {
    if(visited[start] || neighbors[start].size() > 1) { continue; }
    int current = static_cast<int>(start);
    while(current >= 0) {
      if(visited[current]) { return result; }
      visited[current] = true;
      result.blocks.push_back(current);
      result.components.push_back(component);
      int next = -1;
      for(const int adjacent : neighbors[current]) {
        if(!visited[adjacent]) { next = adjacent; }
      }
      current = next;
    }
    ++component;
  }
  result.eligible = count > 0 && result.blocks.size() == count;
  return result;
}
inline BlockSchurPlan planBlockSchur(
    const std::vector<std::vector<int>> & candidates,
    const std::vector<std::pair<int, int>> & couplings,
  int fixedEffects, const BlockSchurPolicy & policy,
  const std::function<bool(const BlockSchurPlan &)> & admissible = nullptr) {
  std::vector<char> inBorder(candidates.size(), 0);
  std::vector<int> degree(candidates.size(), 0);
  for(const auto & pair : couplings) {
    ++degree[pair.first];
    ++degree[pair.second];
  }
  std::vector<std::size_t> order(couplings.size());
  for(std::size_t index = 0; index < order.size(); ++index) {
    order[index] = index;
  }
  std::sort(order.begin(), order.end(), [&](std::size_t first, std::size_t second) {
    const auto & firstPair = couplings[first];
    const auto & secondPair = couplings[second];
    const std::size_t firstSize = std::min(candidates[firstPair.first].size(),
                                          candidates[firstPair.second].size());
    const std::size_t secondSize = std::min(candidates[secondPair.first].size(),
                                           candidates[secondPair.second].size());
    return firstSize != secondSize ? firstSize < secondSize : first < second;
  });
  for(const std::size_t index : order) {
    const int first = couplings[index].first;
    const int second = couplings[index].second;
    if(inBorder[first] || inBorder[second]) { continue; }
    const std::size_t firstSize = candidates[first].size();
    const std::size_t secondSize = candidates[second].size();
    const bool pickFirst = firstSize != secondSize ? firstSize < secondSize
      : degree[first] >= degree[second];
    inBorder[pickFirst ? first : second] = 1;
  }

  BlockSchurPlan bordered;
  for(int column = 0; column < fixedEffects; ++column) {
    bordered.border.push_back(column);
  }
  for(std::size_t group = 0; group < candidates.size(); ++group) {
    if(inBorder[group]) {
      bordered.border.insert(bordered.border.end(), candidates[group].begin(),
                              candidates[group].end());
    } else {
      bordered.groups.push_back(candidates[group]);
    }
  }

  std::vector<int> root(candidates.size());
  for(std::size_t group = 0; group < root.size(); ++group) {
    root[group] = static_cast<int>(group);
  }
  auto findRoot = [&](int group) -> int {
    while(root[static_cast<std::size_t>(group)] != group) {
      root[static_cast<std::size_t>(group)] =
        root[static_cast<std::size_t>(root[static_cast<std::size_t>(group)])];
      group = root[static_cast<std::size_t>(group)];
    }
    return group;
  };
  for(const auto & pair : couplings) {
    const int first = findRoot(pair.first);
    const int second = findRoot(pair.second);
    if(first != second) {
      root[static_cast<std::size_t>(std::max(first, second))] = std::min(first, second);
    }
  }
  BlockSchurPlan merged;
  for(int column = 0; column < fixedEffects; ++column) {
    merged.border.push_back(column);
  }
  std::vector<int> componentIndex(candidates.size(), -1);
  for(std::size_t group = 0; group < candidates.size(); ++group) {
    const std::size_t component = static_cast<std::size_t>(findRoot(static_cast<int>(group)));
    if(componentIndex[component] < 0) {
      componentIndex[component] = static_cast<int>(merged.groups.size());
      merged.groups.emplace_back();
    }
    auto & columns = merged.groups[static_cast<std::size_t>(componentIndex[component])];
    columns.insert(columns.end(), candidates[group].begin(), candidates[group].end());
  }

  double borderedCost = policy.cost(bordered.groups, bordered.border.size());
  double mergedCost = policy.cost(merged.groups, merged.border.size());
  if(admissible){
    if(std::isfinite(borderedCost) && !admissible(bordered)){
      borderedCost = std::numeric_limits<double>::infinity();
    }
    if(std::isfinite(mergedCost) && !admissible(merged)){
      mergedCost = std::numeric_limits<double>::infinity();
    }
  }
  BlockSchurPlan selected = mergedCost < borderedCost ? std::move(merged) : std::move(bordered);
  if(!std::isfinite(borderedCost) && !std::isfinite(mergedCost)) { return selected; }
  selected.eligible = true;
  return selected;
}

}

#endif