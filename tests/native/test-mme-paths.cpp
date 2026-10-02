#include "../../src/mme_paths.h"
#include <cassert>
#include <cmath>
#include <initializer_list>

int main() {
  using sommer::ResidualPath;
  using sommer::ResidualPlan;
  for(bool sectioned : {false, true}) {
    for(bool nonDiagonal : {false, true}) {
      for(bool complete : {false, true}) {
        for(bool weighted : {false, true}) {
          const bool compact = sommer::ResidualPolicy::compactStorage(
            sectioned, nonDiagonal, complete, weighted);
          assert(compact == (!(weighted && complete) && (sectioned || nonDiagonal)));
          for(bool storedCompact : {false, true}) {
            assert(sommer::ResidualPolicy::gridCandidate(storedCompact, complete,
              sectioned, !nonDiagonal, weighted) ==
              (storedCompact && (!complete || sectioned) && (nonDiagonal || sectioned)
               && !weighted));
          }
        }
      }
    }
  }
  assert(sommer::ResidualPolicy::gridFits(1300, 100));
  assert(!sommer::ResidualPolicy::gridFits(1301, 100));
  assert(sommer::ResidualPolicy::gridFits(6000, 2000));
  assert(!sommer::ResidualPolicy::gridFits(6001, 2000));
  for(bool diagonal : {false, true}) {
    for(bool complete : {false, true}) {
      for(bool grid : {false, true}) {
        for(bool weighted : {false, true}) {
          for(bool diagonalWeights : {false, true}) {
            const auto plan = ResidualPlan::select(
              diagonal, complete, grid, weighted, diagonalWeights);
            const bool oldDiagonal = diagonal && !grid;
            const bool oldKronecker = !oldDiagonal && complete && !grid;
            assert(plan.diagonal() == oldDiagonal);
            assert(plan.kronecker() == oldKronecker);
            assert(plan.grid() == grid);
            assert(plan.elementwisePrecision() ==
              (oldDiagonal && (!weighted || diagonalWeights)));
            assert(plan.blockAssembly() ==
              ((grid || oldKronecker) && (!weighted || diagonalWeights)));
            assert((plan.path == ResidualPath::Sparse) ==
              (!oldDiagonal && !grid && !oldKronecker));
            for(bool sparseWeights : {false, true}) {
              const bool sparse = (oldDiagonal && (!weighted || diagonalWeights))
                || (oldDiagonal && sparseWeights);
              const auto expected = sparse ? sommer::AssemblyPath::Sparse
                : ((grid || oldKronecker) && (!weighted || diagonalWeights))
                  ? sommer::AssemblyPath::Blocks : sommer::AssemblyPath::Batched;
              assert(plan.assembly(sparseWeights) == expected);
            }
          }
        }
      }
    }
  }
  for(auto backend : {sommer::SolverBackend::LDLT, sommer::SolverBackend::Cholmod,
                      sommer::SolverBackend::PCG}) {
    for(bool reml : {false, true}) {
      for(int inverse : {0, 1, 2}) {
        for(bool diagonal : {false, true}) {
          const sommer::SolverPlan plan{backend, reml, inverse, diagonal};
          for(bool random : {false, true}) {
            for(bool factorWise : {false, true}) {
              assert(plan.matrixFree(random, factorWise) ==
                (backend == sommer::SolverBackend::PCG && reml && inverse == 0
                 && diagonal && random && factorWise));
            }
            assert(plan.factorRandomOnly(random) == (!reml && random));
          }
          const auto expected = plan.ldlt() ? sommer::TracePreparation::SelectedInverse
            : plan.cholmod() && reml ? sommer::TracePreparation::InverseBlocks
            : plan.pcg() ? sommer::TracePreparation::Probes : sommer::TracePreparation::None;
          assert(plan.traces() == expected);
        }
      }
    }
  }
  const sommer::BlockSchurPolicy policy;
  assert(std::isinf(policy.cost({}, 1)));
  assert(std::isinf(policy.cost({std::vector<int>(2)}, 2001)));
  assert(std::isinf(policy.cost({std::vector<int>(8001)}, 1)));
  assert(std::isfinite(policy.cost({std::vector<int>(8000)}, 2000)));
  assert(policy.cost({std::vector<int>(3)}, 2) ==
    (5.0 / 3.0) * 8.0 + (5.0 / 3.0) * 27.0 + 72.0 + 24.0);
  assert(!policy.needsDensityCheck(1e9, 64));
  assert(!policy.needsDensityCheck(1.25e8 / 100, 100));
  assert(policy.needsDensityCheck(1.25e8 / 100 + 1, 100));
  assert(policy.denseEnough(1000, 100));
  assert(!policy.denseEnough(999, 100));

  const auto independent = sommer::planBlockSchur({{2, 3}, {4, 5}}, {}, 2, policy);
  assert(independent.eligible);
  assert(independent.border == std::vector<int>({0, 1}));
  assert(independent.groups == std::vector<std::vector<int>>({{2, 3}, {4, 5}}));
  const auto merged = sommer::planBlockSchur({{2}, {3}}, {{0, 1}}, 2, policy);
  assert(merged.eligible);
  assert(merged.border == std::vector<int>({0, 1}));
  assert(merged.groups == std::vector<std::vector<int>>({{2, 3}}));
  const auto bordered = sommer::planBlockSchur({{0, 1}, {2, 3}}, {{0, 1}}, 0, policy);
  assert(bordered.eligible);
  assert(bordered.border == std::vector<int>({0, 1}));
  assert(bordered.groups == std::vector<std::vector<int>>({{2, 3}}));
  assert(!sommer::planBlockSchur({}, {}, 0, policy).eligible);
  sommer::BlockSchurPolicy limited = policy;
  limited.maxDoubles = 28;
  assert(sommer::planBlockSchur({{2, 3}}, {}, 2, limited).eligible);
  limited.maxDoubles = 27;
  assert(!sommer::planBlockSchur({{2, 3}}, {}, 2, limited).eligible);
  limited = policy;
  limited.maxGroup = 1;
  limited.maxBorder = 0;
  assert(!sommer::planBlockSchur({{0}, {1}}, {{0, 1}}, 0, limited).eligible);

  for(int fixedEffects : {0, 2}) {
    std::vector<std::vector<int>> candidates(5);
    int effects = fixedEffects;
    for(int group = 0; group < 5; ++group) {
      for(int index = 0; index <= group; ++index) {
        candidates[group].push_back(effects++);
      }
    }
    std::vector<std::pair<int, int>> edges;
    for(int first = 0; first < 5; ++first) {
      for(int second = first + 1; second < 5; ++second) {
        edges.emplace_back(first, second);
      }
    }
    for(unsigned int graph = 0; graph < (1U << edges.size()); ++graph) {
      std::vector<std::pair<int, int>> couplings;
      for(std::size_t edge = 0; edge < edges.size(); ++edge) {
        if(graph & (1U << edge)) { couplings.push_back(edges[edge]); }
      }
      const auto plan = sommer::planBlockSchur(candidates, couplings, fixedEffects, policy);
      assert(plan.eligible);
      const auto repeated = sommer::planBlockSchur(candidates, couplings, fixedEffects, policy);
      assert(plan.border == repeated.border && plan.groups == repeated.groups);
      std::vector<int> membership(effects, -2);
      for(int column : plan.border) {
        assert(membership[column] == -2);
        membership[column] = -1;
      }
      for(std::size_t group = 0; group < plan.groups.size(); ++group) {
        for(int column : plan.groups[group]) {
          assert(membership[column] == -2);
          membership[column] = static_cast<int>(group);
        }
      }
      for(int column = 0; column < effects; ++column) {
        assert(membership[column] != -2);
        if(column < fixedEffects) { assert(membership[column] == -1); }
      }
      for(const auto & coupling : couplings) {
        for(int first : candidates[coupling.first]) {
          for(int second : candidates[coupling.second]) {
            assert(membership[first] == -1 || membership[second] == -1
              || membership[first] == membership[second]);
          }
        }
      }
    }
  }
}