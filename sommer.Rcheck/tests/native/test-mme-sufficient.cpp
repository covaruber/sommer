#include "../../src/mme_sufficient.h"
#include <Eigen/Dense>
#include <cassert>
#include <iostream>

int main() {
  typedef sommer::SufficientStatistics Statistics;
  const int n = 24, p = 4, categories = 2, parameters = 3;
  Eigen::MatrixXd design = Eigen::MatrixXd::Zero(n, p);
  Eigen::VectorXd y(n), weights(n), variances(categories), coefficients(p);
  Eigen::MatrixXd derivatives(categories, parameters);
  std::vector<int> classes(n);
  for(int i = 0; i < n; ++i){
    design(i, 0) = 1.0;
    design(i, 1 + i % 3) = 0.5;
    classes[i] = i % 2;
    weights(i) = 0.3 + (i % 5);
    y(i) = std::sin(i * 0.3) + (i % 3);
  }
  variances << 0.7, 2.3;
  derivatives << 0.7, 0.3, -0.1, 2.3, -0.2, 0.8;
  coefficients << 0.4, -0.7, 0.8, 0.2;
  Eigen::VectorXd ri(n);
  for(int i = 0; i < n; ++i){ ri(i) = weights(i) / variances(classes[i]); }
  const Eigen::VectorXd residual = y - design * coefficients;
  const Eigen::MatrixXd expectedCross = design.transpose() * ri.asDiagonal() * design;
  const Eigen::VectorXd expectedRhs = design.transpose() * ri.asDiagonal() * y;
  for(bool cache : {false, true}){
    Statistics statistics;
    statistics.build(design.sparseView(), y, classes, weights, categories, cache, 1000000);
    assert(statistics.groups.size() == 6);
    assert(statistics.classCounts[0] == 12 && statistics.classCounts[1] == 12);
    Statistics::SparseMatrix cross;
    Eigen::VectorXd rhs;
    double quadratic = 0;
    statistics.assemble(variances, cross, rhs, quadratic);
    assert((Eigen::MatrixXd(cross) - expectedCross).norm() < 1e-12);
    assert((rhs - expectedRhs).norm() < 1e-12);
    assert(std::abs(quadratic - y.dot(ri.cwiseProduct(y))) < 1e-12);
    for(int parameter = 0; parameter < parameters; ++parameter){
      const Eigen::VectorXd scales = derivatives.col(parameter).array() / variances.array().square();
      Eigen::VectorXd recordWeights(n);
      for(int i = 0; i < n; ++i){ recordWeights(i) = weights(i) * scales(classes[i]); }
      assert((Eigen::MatrixXd(statistics.derivativeCross(scales)) -
        design.transpose() * recordWeights.asDiagonal() * design).norm() < 1e-12);
    }
    const Eigen::VectorXd fitted = statistics.design * coefficients;
    const Eigen::VectorXd squares = statistics.residualSquares(fitted);
    Eigen::MatrixXd groupedRhs = Eigen::MatrixXd::Zero(statistics.groups.size(), parameters);
    Eigen::MatrixXd gram = Eigen::MatrixXd::Zero(parameters, parameters);
    Eigen::VectorXd score = Eigen::VectorXd::Zero(parameters);
    double logdet = 0;
    for(std::size_t g = 0; g < statistics.groups.size(); ++g){
      const auto & group = statistics.groups[g];
      const int c = group.residualClass;
      const double v = variances(c);
      logdet += group.count * std::log(v);
      for(int a = 0; a < parameters; ++a){
        groupedRhs(g, a) = group.weight * (group.mean - fitted(g)) * derivatives(c, a) / (v * v);
        score(a) += group.count * derivatives(c, a) / v - squares(g) * derivatives(c, a) / (v * v);
        for(int b = 0; b < parameters; ++b){
          gram(a, b) += squares(g) * derivatives(c, a) * derivatives(c, b) / (v * v * v);
        }
      }
    }
    Eigen::MatrixXd working(n, parameters);
    Eigen::VectorXd expectedScore = Eigen::VectorXd::Zero(parameters);
    double expectedLogdet = 0;
    for(int i = 0; i < n; ++i){
      const int c = classes[i];
      expectedLogdet += std::log(variances(c));
      for(int a = 0; a < parameters; ++a){
        working(i, a) = derivatives(c, a) * residual(i) / variances(c);
        expectedScore(a) += derivatives(c, a) / variances(c) -
          weights(i) * residual(i) * residual(i) * derivatives(c, a) / std::pow(variances(c), 2);
      }
    }
    assert(std::abs(logdet - expectedLogdet) < 1e-12);
    assert((score - expectedScore).norm() < 1e-12);
    assert((gram - working.transpose() * ri.asDiagonal() * working).norm() < 1e-11);
    assert((statistics.design.transpose() * groupedRhs -
      design.transpose() * ri.asDiagonal() * working).norm() < 1e-12);
  }
  design(0, 1) = std::nextafter(design(0, 1), 1.0);
  Statistics exact;
  exact.build(design.sparseView(), y, classes, weights, categories, false, 1000000);
  assert(exact.groups.size() == 7);
  bool rejected = false;
  try { exact.build(design.sparseView(), y, classes, weights, categories, true, 1); }
  catch(const std::runtime_error &) { rejected = true; }
  assert(rejected);
  weights(0) = 0.0;
  rejected = false;
  try { exact.build(design.sparseView(), y, classes, weights, categories, false, 1000000); }
  catch(const std::invalid_argument &) { rejected = true; }
  assert(rejected);
  std::cout << "Sufficient-statistic cross-products, likelihood, score and AI checks passed.\n";
}
