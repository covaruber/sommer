#include "../src/mme_pcg.h"
#include <Eigen/Dense>
#include <chrono>
#include <iostream>
#include <random>

int main() {
  const int dimension = 160, systems = 6, columns = 8, rank = 8;
  std::mt19937 generator(617);
  std::normal_distribution<double> normal;
  Eigen::MatrixXd seed(dimension, dimension), rhs(dimension, columns);
  for(int j = 0; j < dimension; ++j){
    for(int i = 0; i < dimension; ++i){ seed(i, j) = normal(generator); }
  }
  for(int j = 0; j < columns; ++j){
    for(int i = 0; i < dimension; ++i){ rhs(i, j) = normal(generator); }
  }
  const Eigen::MatrixXd rotation = seed.householderQr().householderQ() *
    Eigen::MatrixXd::Identity(dimension, dimension);
  sommer::pcg::Options options;
  options.tolerance = 1.0e-9;
  options.maxIterations = 2000;
  std::size_t baselineCalls = 0, baselineIterations = 0;
  std::cout << "variant,systems,dimension,rhs,rank,iterations,operatorVectors,elapsed,maxTrueRelativeResidual\n";
  for(bool recycling : {false, true}){
    sommer::pcg::RecycleSpace space(32, rank);
    sommer::pcg::Metrics total;
    double maxResidual = 0;
    std::size_t calls = 0;
    const auto started = std::chrono::steady_clock::now();
    for(int system = 0; system < systems; ++system){
      Eigen::VectorXd spectrum(dimension);
      for(int i = 0; i < dimension; ++i){
        const double eigenvalue = i < 6 ? 0.001 * (i + 1) : 1.0 + (i - 6) * 0.4;
        spectrum(i) = eigenvalue * (1.0 + 0.015 * system * std::sin(i + 1.0));
      }
      const Eigen::MatrixXd matrix = rotation * spectrum.asDiagonal() * rotation.transpose();
      const Eigen::VectorXd diagonal = matrix.diagonal();
      const sommer::pcg::VectorOperator applyA =
        [&](const Eigen::VectorXd & value) -> Eigen::VectorXd {
          ++calls;
          return matrix * value;
        };
      const sommer::pcg::VectorOperator applyM =
        [&](const Eigen::VectorXd & value) -> Eigen::VectorXd {
          return value.cwiseQuotient(diagonal);
        };
      if(recycling){
        space.prepare(applyA, applyM);
        options.recycleSpace = &space;
      }else{
        options.recycleSpace = nullptr;
      }
      const sommer::pcg::Result result = sommer::pcg::solveIndependent(rhs,
        applyA, applyM, options);
      total.add(result.metrics);
      for(int column = 0; column < columns; ++column){
        maxResidual = std::max(maxResidual,
          (matrix * result.solution.col(column) - rhs.col(column)).norm() / rhs.col(column).norm());
      }
      if(recycling && system + 1 < systems){
        const auto selection = space.selectCandidates(applyA, diagonal, rank);
        space.accept(selection);
      }
    }
    const double seconds = std::chrono::duration<double>(
      std::chrono::steady_clock::now() - started).count();
    std::cout << (recycling ? "krylov" : "baseline") << ',' << systems << ','
      << dimension << ',' << columns << ',' << (recycling ? rank : 0) << ','
      << total.iterations << ',' << calls << ',' << seconds << ',' << maxResidual << '\n';
    if(!std::isfinite(maxResidual) || maxResidual > options.tolerance){ return 1; }
    if(!recycling){
      baselineCalls = calls;
      baselineIterations = total.iterations;
    }else if(calls >= baselineCalls || total.iterations >= baselineIterations){
      std::cerr << "Learned recycling did not reduce work on the difficult-system fixture.\n";
      return 2;
    }
  }
}
