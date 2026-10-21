#include "../../src/mme_paths.h"
#include <cassert>
#include <chrono>
#include <cmath>
#include <initializer_list>
#include <string>

int main() {
  assert(sommer::FactorPolicy::spatial(sommer::PrecisionFamily::NormalEquations));
  assert(sommer::FactorPolicy::spatial(sommer::PrecisionFamily::Affine));
  assert(!sommer::FactorPolicy::spatial(sommer::PrecisionFamily::ModifiedCholesky));
  assert(sommer::FactorPolicy::innovationDerivative(sommer::PrecisionFamily::ModifiedCholesky, true));
  assert(!sommer::FactorPolicy::innovationDerivative(sommer::PrecisionFamily::ModifiedCholesky, false));
  assert(sommer::FactorPolicy::logDetMultiplicity(12, 4) == 3.0);
  assert(std::isnan(sommer::FactorPolicy::logDetMultiplicity(12, 5)));
  assert(sommer::WorkspacePolicy::cacheFits(100, 80, 20));
  assert(!sommer::WorkspacePolicy::cacheFits(100, 80, 21));
  assert(!sommer::WorkspacePolicy::cacheFits(100, 101, 0));
  assert(sommer::WorkspacePolicy::crossBlockWorthCaching(200, 100));
  assert(!sommer::WorkspacePolicy::crossBlockWorthCaching(199, 100));
  assert(sommer::WorkspacePolicy::panelColumns() == 32);
  assert(sommer::InverseCachePolicy::preferSelected(true, 100.0, 100.0));
  assert(!sommer::InverseCachePolicy::preferSelected(false, 100.0, 100.0));
  assert(!sommer::InverseCachePolicy::preferSelected(true, 101.0, 100.0));
  for(double invalid : {0.0, -1.0, std::numeric_limits<double>::infinity(),
                        std::numeric_limits<double>::quiet_NaN()}){
    assert(!sommer::InverseCachePolicy::preferSelected(true, invalid, 100.0));
  }
  assert(!sommer::InverseCachePolicy::preferSelected(
    true, 100.0, std::numeric_limits<double>::quiet_NaN()));
  assert(sommer::WorkspacePolicy::inversePanelColumns() == 64);
  assert(sommer::FactorPolicy::shiftedCholesky(1, 128));
  assert(!sommer::FactorPolicy::shiftedCholesky(2, 128));
  assert(sommer::FactorPolicy::supportContraction(50, 10, 50));
  assert(!sommer::FactorPolicy::supportContraction(51, 10, 100));
  using sommer::ResidualPath;
  using sommer::ResidualPlan;
  const auto identityTensor = sommer::planTensorComponents({}, {}, 1);
  assert(identityTensor.eligible);
  assert(identityTensor.groups == std::vector<std::vector<int>>({{0}}));
  sommer::DisjointComponents graph(4);
  assert(graph.unite(0, 2));
  assert(graph.unite(1, 3));
  assert(!graph.unite(0, 4));
  assert(graph.labels() == std::vector<int>({0, 1, 0, 1}));
  const auto graphTensor = sommer::planTensorGraphComponents(
    {2, 3}, {{0, 0}, {0, 1, 0}}, 6);
  assert(graphTensor.eligible);
  assert(graphTensor.groups == std::vector<std::vector<int>>({{0, 2, 3, 5}, {1, 4}}));
  assert(!sommer::planTensorGraphComponents({2}, {{0, 2}}, 2).eligible);
  assert(!sommer::planTensorGraphComponents({2}, {{0}}, 2).eligible);
  assert(!sommer::planTensorGraphComponents({2, 3}, {{0, 1}, {0, 0, 0}}, 5).eligible);
  assert(!sommer::planTensorComponents({2}, {}, 2).eligible);
  assert(!sommer::planTensorComponents({0}, {true}, 1).eligible);
  assert(!sommer::planTensorComponents({2, 3}, {true, false}, 5).eligible);
  assert(!sommer::planTensorComponents(
    {std::numeric_limits<std::size_t>::max(), 2}, {true, false},
    std::numeric_limits<std::size_t>::max()).eligible);
  const auto partialTensor = sommer::planTensorComponents(
    {2, 3, 2}, {true, false, true}, 12);
  assert(partialTensor.eligible);
  assert(partialTensor.coordinates == 12);
  assert(partialTensor.groups == std::vector<std::vector<int>>(
    {{0, 2, 4}, {1, 3, 5}, {6, 8, 10}, {7, 9, 11}}));
  for(std::size_t axes = 1; axes <= 5; ++axes) {
    const std::vector<std::size_t> dimensions(axes, 2);
    const std::size_t coordinates = std::size_t(1) << axes;
    for(std::size_t mask = 0; mask < coordinates; ++mask) {
      std::vector<bool> diagonal(axes, false);
      std::size_t expectedGroups = 1;
      for(std::size_t axis = 0; axis < axes; ++axis) {
        diagonal[axis] = (mask & (std::size_t(1) << axis)) != 0;
        if(diagonal[axis]) { expectedGroups *= 2; }
      }
      const auto tensor = sommer::planTensorComponents(dimensions, diagonal, coordinates);
      assert(tensor.eligible && tensor.groups.size() == expectedGroups);
      std::vector<bool> visited(coordinates, false);
      for(const auto & group : tensor.groups) {
        assert(group.size() == coordinates / expectedGroups);
        for(int coordinate : group) {
          assert(!visited[coordinate]);
          visited[coordinate] = true;
          for(std::size_t axis = 0; axis < axes; ++axis) {
            if(!diagonal[axis]) { continue; }
            const std::size_t shift = axes - axis - 1;
            assert(((coordinate >> shift) & 1) == ((group.front() >> shift) & 1));
          }
        }
      }
      assert(std::all_of(visited.begin(), visited.end(), [](bool value) { return value; }));
    }
  }
  const auto chain = sommer::orderBlockPaths(4, {{0, 2}, {2, 1}, {1, 3}});
  assert(chain.eligible);
  assert(chain.blocks == std::vector<int>({0, 2, 1, 3}));
  assert(chain.components == std::vector<int>({0, 0, 0, 0}));
  const auto disconnected = sommer::orderBlockPaths(4, {{0, 2}});
  assert(disconnected.eligible);
  assert(disconnected.blocks == std::vector<int>({0, 2, 1, 3}));
  assert(disconnected.components == std::vector<int>({0, 0, 1, 2}));
  assert(!sommer::orderBlockPaths(3, {{0, 1}, {1, 2}, {2, 0}}).eligible);
  assert(!sommer::orderBlockPaths(4, {{0, 1}, {0, 2}, {0, 3}}).eligible);
  assert(!sommer::orderBlockPaths(2, {{0, 2}}).eligible);
  assert(!sommer::orderBlockPaths(0, {}).eligible);
  const auto bandwidthTwo = sommer::planBlockBand(6, {{0, 1}, {1, 2}, {0, 2}, {2, 4}}, 2);
  assert(bandwidthTwo.eligible && bandwidthTwo.bandwidth == 2);
  assert(!sommer::planBlockBand(5, {{0, 3}}, 2).eligible);
  assert(!sommer::planBlockBand(3, {{0, 3}}, 3).eligible);
  assert(!sommer::planBlockBand(4, {{1, 1}}, 2).eligible);
  assert(sommer::BlockBandPolicy::isAdvantageous(50, 2, 2, 100));
  assert(!sommer::BlockBandPolicy::isAdvantageous(30, 2, 2, 100));
  assert(!sommer::BlockBandPolicy::isAdvantageous(50, 5, 2, 100));
  assert(!sommer::BlockBandPolicy::isAdvantageous(50, 2, 10, 100));
  assert(!sommer::BlockBandPolicy::isAdvantageous(50, 2, 2, 0));
  assert(sommer::LatentSchurPolicy::isAdvantageousMultiFactor(8, 200));
  assert(sommer::LatentSchurPolicy::isAdvantageousMultiFactor(16, 300));
  assert(!sommer::LatentSchurPolicy::isAdvantageousMultiFactor(4, 200));
  assert(!sommer::LatentSchurPolicy::isAdvantageousMultiFactor(8, 350));
  using sommer::LatentSchurPolicy;
  const double lowRankWork = LatentSchurPolicy::latentWork(6, 128, 134);
  const double diagonalMarginalWork = LatentSchurPolicy::marginalWork({384, 384}, 6);
  assert(LatentSchurPolicy::workFavoursLatent(lowRankWork, diagonalMarginalWork));
  assert(!LatentSchurPolicy::workFavoursLatent(
    LatentSchurPolicy::latentWork(6, 128, 390), diagonalMarginalWork));
  const double sameLatentWork = LatentSchurPolicy::latentWork(4, 192, 196);
  assert(!LatentSchurPolicy::workFavoursLatent(
    sameLatentWork, LatentSchurPolicy::marginalWork({256, 256, 256}, 4)));
  assert(LatentSchurPolicy::workFavoursLatent(
    sameLatentWork, LatentSchurPolicy::marginalWork({768}, 4)));
  assert(LatentSchurPolicy::workFavoursLatent(80, 100));
  assert(!LatentSchurPolicy::workFavoursLatent(80.001, 100));
  for(double invalid : {0.0, -1.0, std::numeric_limits<double>::infinity(),
                        std::numeric_limits<double>::quiet_NaN()}){
    assert(!LatentSchurPolicy::workFavoursLatent(invalid, 100));
    assert(!LatentSchurPolicy::workFavoursLatent(100, invalid));
  }
  assert(std::isnan(LatentSchurPolicy::marginalWork({}, 1)));
  assert(std::isnan(LatentSchurPolicy::marginalWork({0}, 1)));
  assert(LatentSchurPolicy::specificsWellScaled(1.0, 0.5));
  assert(LatentSchurPolicy::specificsWellScaled(1e-8, 5e-9));
  assert(!LatentSchurPolicy::specificsWellScaled(1.0, 1e-6));
  const double numericalLimit =
    1.0 / std::sqrt(std::sqrt(std::numeric_limits<double>::epsilon()));
  assert(LatentSchurPolicy::specificsWellScaled(numericalLimit, 1.0));
  assert(!LatentSchurPolicy::specificsWellScaled(numericalLimit + 1.0, 1.0));
  for(double invalid : {0.0, -1.0, std::numeric_limits<double>::infinity(),
                        std::numeric_limits<double>::quiet_NaN()}){
    assert(!LatentSchurPolicy::specificsWellScaled(1.0, invalid));
  }
  assert(!LatentSchurPolicy::specificsWellScaled(
    std::numeric_limits<double>::infinity(), 1.0));
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
  assert(std::isfinite(policy.bandStorageDoubles({{0, 1}, {2, 3}, {4, 5}}, 2, 1)));
  assert(std::isinf(policy.bandStorageDoubles({{0}, {1}}, 2, 1)));
  const sommer::SolverPlan mlPlan{sommer::SolverBackend::LDLT, false, 0, true};
  assert(mlPlan.dFirstML(2, 10, 76, false));
  assert(!mlPlan.dFirstML(2, 10, 75, false));
  assert(!mlPlan.dFirstML(2, 10, 76, true));
  assert(!mlPlan.dFirstML(257, 10, 1e9, false));
  assert(!mlPlan.dFirstML(2, 0, 76, false));
  const sommer::SolverPlan subsetPlan{sommer::SolverBackend::LDLT, false, 1, true};
  assert(!subsetPlan.dFirstML(2, 10, 76, false));
  assert(std::isinf(policy.cost({}, 1)));
  assert(std::isinf(policy.cost({std::vector<int>(2)}, 2001)));
  assert(std::isinf(policy.cost({std::vector<int>(8001)}, 1)));
  assert(std::isinf(policy.cost({std::vector<int>(8000)}, 2000)));
  sommer::BlockSchurPolicy expanded = policy;
  expanded.maxDoubles = 2.1e8;
  assert(std::isfinite(expanded.cost({std::vector<int>(10000)}, 30)));
  assert(std::isinf(expanded.cost({std::vector<int>(11000)}, 30)));
  assert(std::isfinite(policy.cost({std::vector<int>(1000), std::vector<int>(1000)}, 1030)));
  assert(policy.storageDoubles({std::vector<int>(3)}, 2) == 44.0);
  assert(policy.cost({std::vector<int>(3)}, 2) ==
    (5.0 / 3.0) * 8.0 + (5.0 / 3.0) * 27.0 + 72.0 + 24.0);
  assert(!policy.needsDensityCheck(1e9, 64));
  assert(policy.fillFavoursSchur(1e7, 1e8));
  assert(!policy.fillFavoursSchur(1e4, 1e8));
  assert(!policy.fillFavoursSchur(0.049999 * 1e8, 1e8));
  assert(policy.fillFavoursSchur(0.05 * 1e8, 1e8));
  assert(!policy.fillFavoursSchur(-1.0, 1e8));
  assert(!policy.fillFavoursSchur(std::numeric_limits<double>::infinity(), 1e8));
  assert(!policy.fillFavoursSchur(std::numeric_limits<double>::quiet_NaN(), 1e8));
  assert(!policy.fillFavoursSchur(1e6, std::numeric_limits<double>::infinity()));
  assert(!policy.fillFavoursSchur(1e6, 0.0));
  assert(policy.needsDensityCheck(1.25e8 / 100, 100));
  assert(policy.needsDensityCheck(1.25e8 / 100 + 1, 100));
  assert(policy.inverseCacheFits(10000, 1000));
  assert(!policy.inverseCacheFits(10030, 10000));
  sommer::BlockSchurPolicy cacheBoundary = policy;
  cacheBoundary.maxDoubles = 28;
  assert(cacheBoundary.inverseCacheFits(10, 2));
  cacheBoundary.maxDoubles = 27;
  assert(!cacheBoundary.inverseCacheFits(10, 2));
  assert(policy.chainStorageDoubles({{0, 1}, {2, 3}}, {0, 0}, 1) == 67.0);
  assert(policy.chainStorageDoubles({{0, 1}, {2, 3}}, {0, 1}, 1) == 55.0);
  assert(std::isinf(policy.chainStorageDoubles({{0}}, {}, 1)));
  assert(std::isinf(policy.chainStorageDoubles({{0}}, {0}, 2001)));
  assert(policy.latentStorageDoubles({{0, 1}, {2, 3}}, 1, 2) == 63.0);
  assert(std::isinf(policy.latentStorageDoubles({{0}, {1, 2}}, 1, 2)));
  assert(std::isinf(policy.latentStorageDoubles({}, 1, 2)));
  assert(std::isfinite(policy.latentStorageDoubles({std::vector<int>(1000)}, 2030, 1000)));
  assert(policy.denseEnough(1000, 100));
  assert(!policy.denseEnough(999, 100));

  const auto independent = sommer::planBlockSchur({{2, 3}, {4, 5}}, {}, 2, policy);
  assert(independent.eligible);
  assert(independent.border == std::vector<int>({0, 1}));
  assert(independent.groups == std::vector<std::vector<int>>({{2, 3}, {4, 5}}));
  assert(!sommer::planBlockSchur({{2, 3}, {4, 5}}, {}, 2, policy,
    [&](const sommer::BlockSchurPlan & plan) {
      return policy.fillFavoursSchur(1.0, policy.cost(plan.groups, plan.border.size()));
    }).eligible);
  assert(sommer::planBlockSchur({{2, 3}, {4, 5}}, {}, 2, policy,
    [&](const sommer::BlockSchurPlan & plan) {
      return policy.fillFavoursSchur(1e6, policy.cost(plan.groups, plan.border.size()));
    }).eligible);
  const auto merged = sommer::planBlockSchur({{2}, {3}}, {{0, 1}}, 2, policy);
  assert(merged.eligible);
  assert(merged.border == std::vector<int>({0, 1}));
  assert(merged.groups == std::vector<std::vector<int>>({{2, 3}}));
  const auto bordered = sommer::planBlockSchur({{0, 1}, {2, 3}}, {{0, 1}}, 0, policy);
  assert(bordered.eligible);
  assert(bordered.border == std::vector<int>({0, 1}));
  assert(bordered.groups == std::vector<std::vector<int>>({{2, 3}}));
  const auto alternative = sommer::planBlockSchur({{0, 1}, {2, 3}}, {{0, 1}}, 0, policy,
    [](const sommer::BlockSchurPlan & plan) { return plan.border.empty(); });
  assert(alternative.eligible && alternative.border.empty());
  assert(alternative.groups == std::vector<std::vector<int>>({{0, 1, 2, 3}}));
  assert(!sommer::planBlockSchur({{0, 1}, {2, 3}}, {{0, 1}}, 0, policy,
    [](const sommer::BlockSchurPlan &) { return false; }).eligible);
  assert(!sommer::planBlockSchur({}, {}, 0, policy).eligible);
  sommer::BlockSchurPolicy limited = policy;
  limited.maxDoubles = 28;
  assert(sommer::planBlockSchur({{2, 3}}, {}, 2, limited).eligible);
  limited.maxDoubles = 27;
  assert(!sommer::planBlockSchur({{2, 3}}, {}, 2, limited).eligible);
  limited = policy;
  limited.maxDoubles = 2;
  limited.maxBorder = 0;
  assert(!sommer::planBlockSchur({{0}, {1}}, {{0, 1}}, 0, limited).eligible);

  limited = policy;
  limited.maxDoubles = 90;
  const auto affordable = sommer::planBlockSchur({{0, 1, 2}, {3, 4, 5, 6}}, {{0, 1}}, 0, limited);
  assert(affordable.eligible);
  assert(affordable.border == std::vector<int>({0, 1, 2}));
  assert(affordable.groups == std::vector<std::vector<int>>({{3, 4, 5, 6}}));

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
  {
    using sommer::FitPhase;
    using sommer::PhaseProfile;
    auto spin = [](double seconds) {
      const auto until = std::chrono::steady_clock::now() +
        std::chrono::duration<double>(seconds);
      while(std::chrono::steady_clock::now() < until) {}
    };
    PhaseProfile off(false);
    off.startIteration();
    off.enter(FitPhase::Factorization);
    off.finish();
    assert(!off.enabled() && off.iterations().empty());
    for(double seconds : off.total()) { assert(seconds == 0.0); }

    PhaseProfile profile(true);
    spin(0.01);
    for(int iteration = 0; iteration < 3; ++iteration) {
      profile.startIteration();
      spin(0.002);
      profile.enter(FitPhase::Factorization);
      spin(0.004);
      profile.enter(FitPhase::Update);
    }
    profile.enter(FitPhase::Final);
    spin(0.003);
    profile.finish();
    const auto factorization = static_cast<std::size_t>(FitPhase::Factorization);
    const auto factorEval = static_cast<std::size_t>(FitPhase::FactorEval);
    const auto setup = static_cast<std::size_t>(FitPhase::Setup);
    const auto final = static_cast<std::size_t>(FitPhase::Final);
    assert(profile.iterations().size() == 3);
    assert(profile.total()[setup] >= 0.01);
    assert(profile.total()[final] >= 0.003);
    assert(profile.total()[factorization] >= 3 * 0.004);
    assert(profile.total()[factorEval] >= 3 * 0.002);
    double rowSum = 0.0;
    for(const auto & row : profile.iterations()) {
      assert(row[setup] == 0.0 && row[final] == 0.0);
      for(double seconds : row) { rowSum += seconds; }
    }
    double totalSum = 0.0;
    for(double seconds : profile.total()) { totalSum += seconds; }
    assert(std::fabs(totalSum - rowSum - profile.total()[setup] - profile.total()[final]) < 1e-9);
    assert(std::string(sommer::fitPhaseName(sommer::kFitPhaseCount - 1)) == "final");
    assert(std::string(sommer::fitPhaseName(sommer::kFitPhaseCount)).empty());
  }
  {
    // Column-major pattern builder: neighbours[j] lists the rows of column j.
    auto pattern = [](const std::vector<std::vector<int>> & neighbours, std::vector<int> & outer, std::vector<int> & inner) {
      outer.assign(1, 0);
      inner.clear();
      for(const auto & rows : neighbours) {
        inner.insert(inner.end(), rows.begin(), rows.end());
        outer.push_back(static_cast<int>(inner.size()));
      }
    };
    std::vector<int> outer, inner;
    const std::size_t n = 10;
    std::vector<int> identity(n);
    for(std::size_t index = 0; index < n; ++index) { identity[index] = static_cast<int>(index); }
    std::vector<std::vector<int>> chain(n);
    for(std::size_t index = 0; index < n; ++index) {
      if(index > 0) { chain[index].push_back(static_cast<int>(index) - 1); }
      if(index + 1 < n) { chain[index].push_back(static_cast<int>(index) + 1); }
    }
    pattern(chain, outer, inner);
    const auto tridiagonal = sommer::symbolicCholeskyCounts(n, outer, inner, identity, 1000);
    assert(tridiagonal.valid && tridiagonal.nonzeros == 2 * n - 1);
    assert(std::fabs(tridiagonal.flops - (4.0 * (n - 1) + 1.0)) < 1e-12);
    std::vector<std::vector<int>> star(n);
    for(std::size_t index = 1; index < n; ++index) {
      star[0].push_back(static_cast<int>(index));
      star[index].push_back(0);
    }
    pattern(star, outer, inner);
    const auto hubFirst = sommer::symbolicCholeskyCounts(n, outer, inner, identity, 1000);
    assert(hubFirst.valid && hubFirst.nonzeros == n * (n + 1) / 2);
    std::vector<int> hubLast(n);
    for(std::size_t index = 0; index < n; ++index) { hubLast[index] = static_cast<int>((index + 1) % n); }
    const auto reordered = sommer::symbolicCholeskyCounts(n, outer, inner, hubLast, 1000);
    assert(reordered.valid && reordered.nonzeros == 2 * n - 1);
    assert(reordered.flops < hubFirst.flops);
    const auto capped = sommer::symbolicCholeskyCounts(n, outer, inner, identity, 20);
    assert(!capped.valid && std::isinf(capped.flops));
    assert(!sommer::symbolicCholeskyCounts(n, outer, inner, std::vector<int>(n, 0), 1000).valid);
  }
}