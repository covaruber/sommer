#ifndef SOMMER_MME_PATHS_H
#define SOMMER_MME_PATHS_H

#include <algorithm>
#include <array>
#include <chrono>
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

enum class PrecisionFamily {
  Dense, Diagonal, AR1, StationaryAR, ModifiedCholesky,
  LowRankDiagonal, CompoundSymmetry, NormalEquations, Affine
};

struct FactorPolicy {
  static bool shiftedCholesky(std::size_t factors, std::size_t dimension) {
    return factors == 1 && dimension >= 128;
  }
  static bool supportContraction(double support, double dimension, std::size_t capacity) {
    return support <= 0.5 * dimension * dimension && support <= static_cast<double>(capacity);
  }
  static bool spatial(PrecisionFamily family) {
    return family == PrecisionFamily::NormalEquations || family == PrecisionFamily::Affine;
  }
  static bool innovationDerivative(PrecisionFamily family, bool sparseEnabled) {
    return sparseEnabled && family == PrecisionFamily::ModifiedCholesky;
  }
  static double logDetMultiplicity(std::size_t total, std::size_t factor) {
    return factor > 0 && total % factor == 0
      ? static_cast<double>(total / factor) : std::numeric_limits<double>::quiet_NaN();
  }
};

struct WorkspacePolicy {
  static std::size_t cacheEntries() { return 1048576; }
  static std::size_t panelColumns() { return 32; }
  static std::size_t inversePanelColumns() { return 64; }
  static double inversePanelBytes(std::size_t rows, std::size_t columns) {
    return 3.0 * sizeof(double) * static_cast<double>(rows) *
      static_cast<double>(columns);
  }
  static std::size_t nystromRank() { return 32; }
  static std::size_t crossGroups() { return 128; }
  static std::size_t indexedColumnRows() { return 64; }
  static std::size_t sparseEntryLimit() { return 16000000; }
  static bool cacheFits(std::size_t capacity, std::size_t retained, std::size_t requested) {
    return retained <= capacity && requested <= capacity - retained;
  }
  static bool crossBlockWorthCaching(double dotFlops, double materializeFlops) {
    return dotFlops >= 2.0 * materializeFlops;
  }
  static std::size_t cacheCapacity(double doubles, std::size_t entryBytes) {
    return static_cast<std::size_t>(std::min(static_cast<double>(cacheEntries()),
      doubles * sizeof(double) / (48.0 * static_cast<double>(entryBytes))));
  }
};

struct InverseCachePolicy {
  static bool preferSelected(bool patternLocal, double selectedDoubles, double budgetDoubles) {
    return patternLocal && std::isfinite(selectedDoubles) && selectedDoubles > 0.0 &&
      std::isfinite(budgetDoubles) && selectedDoubles <= budgetDoubles;
  }
};

// Phases are consecutive laps of one fit; Factorization also covers the
// symbolic analysis and any PCG preparation, Final covers post-loop output.
enum class FitPhase {
  Setup, FactorEval, Assembly, Factorization, Acceptance, TracePrep,
  Sensitivities, AIAssembly, ScoreTraces, Update, Constraints, Proposal, Final
};

const std::size_t kFitPhaseCount = 13;

inline const char * fitPhaseName(std::size_t phase) {
  static const char * const names[kFitPhaseCount] = {
    "setup", "factorEval", "assembly", "factorization", "acceptance",
    "tracePrep", "sensitivities", "aiAssembly", "scoreTraces", "update",
    "constraints", "proposal", "final"
  };
  return phase < kFitPhaseCount ? names[phase] : "";
}

class PhaseProfile {
public:
  typedef std::array<double, kFitPhaseCount> Row;

private:
  bool on;
  std::size_t current;
  std::chrono::steady_clock::time_point mark;
  Row totals;
  std::vector<Row> rows;

  void lap() {
    const std::chrono::steady_clock::time_point now = std::chrono::steady_clock::now();
    const double seconds = std::chrono::duration<double>(now - mark).count();
    mark = now;
    totals[current] += seconds;
    if(!rows.empty() && current != static_cast<std::size_t>(FitPhase::Final)) {
      rows.back()[current] += seconds;
    }
  }

public:
  explicit PhaseProfile(bool enabled = false)
    : on(enabled), current(static_cast<std::size_t>(FitPhase::Setup)), mark(), totals() {
    if(on) { mark = std::chrono::steady_clock::now(); }
    totals.fill(0.0);
  }
  bool enabled() const { return on; }
  void enter(FitPhase phase) {
    if(!on) { return; }
    lap();
    current = static_cast<std::size_t>(phase);
  }
  void startIteration() {
    if(!on) { return; }
    enter(FitPhase::FactorEval);
    Row empty;
    empty.fill(0.0);
    rows.push_back(empty);
  }
  void finish() {
    if(!on) { return; }
    lap();
  }
  const Row & total() const { return totals; }
  const std::vector<Row> & iterations() const { return rows; }
};

struct SymbolicFactor {
  bool valid = false;
  std::size_t nonzeros = 0;
  double flops = 0.0;
};

// Sparse Cholesky fill and work for a symmetric pattern under permutation[new] = old.
// Counts follow the row-subtree method on the elimination tree; flops = sum of squared column counts.
inline SymbolicFactor symbolicCholeskyCounts(
    std::size_t size, const std::vector<int> & outer, const std::vector<int> & inner,
    const std::vector<int> & permutation, std::size_t maxNonzeros) {
  SymbolicFactor result;
  if(size == 0 || outer.size() != size + 1 || permutation.size() != size ||
     size > static_cast<std::size_t>(std::numeric_limits<int>::max())) { return result; }
  std::vector<int> position(size, -1);
  for(std::size_t index = 0; index < size; ++index) {
    const int original = permutation[index];
    if(original < 0 || static_cast<std::size_t>(original) >= size || position[original] >= 0) { return result; }
    position[original] = static_cast<int>(index);
  }
  std::vector<std::vector<int>> lowerNeighbours(size);
  for(std::size_t column = 0; column < size; ++column) {
    for(int entry = outer[column]; entry < outer[column + 1]; ++entry) {
      const int row = inner[entry];
      if(row < 0 || static_cast<std::size_t>(row) >= size) { return result; }
      const int first = position[column];
      const int second = position[row];
      if(first == second) { continue; }
      lowerNeighbours[std::max(first, second)].push_back(std::min(first, second));
    }
  }
  std::vector<int> parent(size, -1), ancestor(size, -1), mark(size, -1);
  std::vector<std::size_t> columnCount(size, 1);
  for(std::size_t node = 0; node < size; ++node) {
    const int current = static_cast<int>(node);
    for(const int start : lowerNeighbours[node]) {
      int walker = start;
      while(ancestor[walker] != -1 && ancestor[walker] != current) {
        const int next = ancestor[walker];
        ancestor[walker] = current;
        walker = next;
      }
      if(ancestor[walker] == -1) { ancestor[walker] = current; parent[walker] = current; }
    }
    mark[node] = current;
    for(const int start : lowerNeighbours[node]) {
      int walker = start;
      while(walker >= 0 && mark[walker] != current) {
        ++columnCount[walker];
        mark[walker] = current;
        walker = parent[walker];
      }
    }
  }
  std::size_t total = 0;
  double flops = 0.0;
  for(std::size_t column = 0; column < size; ++column) {
    total += columnCount[column];
    flops += static_cast<double>(columnCount[column]) * static_cast<double>(columnCount[column]);
    if(total > maxNonzeros) {
      result.nonzeros = total;
      result.flops = std::numeric_limits<double>::infinity();
      return result;
    }
  }
  result.valid = true;
  result.nonzeros = total;
  result.flops = flops;
  return result;
}

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
  double minFillRatio = 0.05;

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
      result += blockWork(size, borderSize);
    }
    return result;
  }

  static double blockWork(double size, double border) {
    return (5.0 / 3.0) * size * size * size +
      4.0 * border * size * size + 2.0 * size * border * border;
  }

  bool fillFavoursSchur(double sparseFlops, double schurCost) const {
    if(!std::isfinite(schurCost) || schurCost <= 0.0) { return false; }
    if(!std::isfinite(sparseFlops) || sparseFlops < 0.0) { return false; }
    return sparseFlops >= minFillRatio * schurCost;
  }
  bool needsDensityCheck(double /*effects*/, double group) const {
    return group > smallGroup;
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
    // Preserve the established compact-border region; work estimates admit
    // additional cases without changing those existing decisions.
    return (lowRankDimension >= 8 && borderSize <= 300);
  }
  static double latentWork(std::size_t blocks, std::size_t width, std::size_t border) {
    const double b = static_cast<double>(border);
    return (5.0 / 3.0) * b * b * b +
      static_cast<double>(blocks) * BlockSchurPolicy::blockWork(
        static_cast<double>(width), b);
  }
  static double marginalWork(const std::vector<std::size_t> & widths,
                             std::size_t border) {
    const double b = static_cast<double>(border);
    double work = (5.0 / 3.0) * b * b * b;
    if(widths.empty()) { return std::numeric_limits<double>::quiet_NaN(); }
    for(const std::size_t width : widths) {
      if(width == 0) { return std::numeric_limits<double>::quiet_NaN(); }
      work += BlockSchurPolicy::blockWork(static_cast<double>(width), b);
    }
    return work;
  }
  static bool workFavoursLatent(double latent, double marginal) {
    return std::isfinite(latent) && latent > 0.0 &&
      std::isfinite(marginal) && marginal > 0.0 && latent <= 0.8 * marginal;
  }
  static bool specificsWellScaled(double loadingEnergy, double smallestSpecific) {
    // The latent Schur subtraction and information contractions can compound
    // loss of precision as specific variances approach zero.
    return std::isfinite(loadingEnergy) && loadingEnergy >= 0.0 &&
      std::isfinite(smallestSpecific) && smallestSpecific > 0.0 &&
      loadingEnergy / smallestSpecific <=
        1.0 / std::sqrt(std::sqrt(std::numeric_limits<double>::epsilon()));
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