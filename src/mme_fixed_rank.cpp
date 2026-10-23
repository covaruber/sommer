#include <RcppEigen.h>
#include <Eigen/SparseQR>
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

// [[Rcpp::export]]
Rcpp::List fixed_design_basis_cpp(const Eigen::MappedSparseMatrix<double> & fixedDesign){
  typedef Eigen::SparseMatrix<double> SparseDesign;
  typedef Eigen::Triplet<double> DesignEntry;
  const int observationTotal = fixedDesign.rows();
  const int columnTotal = fixedDesign.cols();
  const double rankTolerance = std::max(observationTotal, columnTotal) *
    std::numeric_limits<double>::epsilon();
  std::vector<int> nonzeroColumns;
  std::vector<double> columnScales(columnTotal, 0.0);
  std::vector<DesignEntry> scaledEntries;
  scaledEntries.reserve(fixedDesign.nonZeros());
  for(int columnIndex = 0; columnIndex < columnTotal; ++columnIndex){
    double columnNorm = 0.0;
    for(Eigen::MappedSparseMatrix<double>::InnerIterator entry(fixedDesign, columnIndex); entry; ++entry){
      if(!std::isfinite(entry.value())){
        Rcpp::stop("The fixed-effect design contains non-finite values.");
      }
      columnNorm = std::hypot(columnNorm, entry.value());
    }
    columnScales[columnIndex] = columnNorm;
    if(columnNorm == 0.0) continue;
    if(!std::isfinite(columnNorm)){
      Rcpp::stop("A fixed-effect column norm exceeds the numerical range; rescale the covariate.");
    }
    const int scaledColumn = nonzeroColumns.size();
    nonzeroColumns.push_back(columnIndex);
    for(Eigen::MappedSparseMatrix<double>::InnerIterator entry(fixedDesign, columnIndex); entry; ++entry){
      if(entry.value() != 0.0){
        scaledEntries.emplace_back(entry.row(), scaledColumn, entry.value() / columnNorm);
      }
    }
  }
  if(nonzeroColumns.empty()){
    Rcpp::stop("No nonzero fixed-effect columns remain after observation filtering; revise the fixed formula.");
  }
  SparseDesign scaledDesign(observationTotal, nonzeroColumns.size());
  scaledDesign.setFromTriplets(scaledEntries.begin(), scaledEntries.end());
  Eigen::SparseQR<SparseDesign, Eigen::COLAMDOrdering<int>> basisFactor;
  basisFactor.setPivotThreshold(rankTolerance);
  basisFactor.compute(scaledDesign);
  if(basisFactor.info() != Eigen::Success){
    Rcpp::stop("Sparse QR failed while selecting a full-rank fixed-effect basis.");
  }
  const int basisRank = basisFactor.rank();
  if(basisRank == 0){
    Rcpp::stop("The fixed-effect design has numerical rank zero; rescale covariates or revise the fixed formula.");
  }
  const Eigen::VectorXi basisPermutation = basisFactor.colsPermutation().indices();
  std::vector<int> retainedColumns;
  for(int basisIndex = 0; basisIndex < basisRank; ++basisIndex){
    retainedColumns.push_back(nonzeroColumns[basisPermutation[basisIndex]]);
  }
  std::sort(retainedColumns.begin(), retainedColumns.end());
  std::vector<int> retainedPositions(columnTotal, -1);
  Rcpp::IntegerVector retainedIndices(basisRank);
  for(int basisIndex = 0; basisIndex < basisRank; ++basisIndex){
    retainedPositions[retainedColumns[basisIndex]] = basisIndex;
    retainedIndices[basisIndex] = retainedColumns[basisIndex] + 1;
  }
  std::vector<int> omittedIndices;
  for(int columnIndex = 0; columnIndex < columnTotal; ++columnIndex){
    if(retainedPositions[columnIndex] < 0) omittedIndices.push_back(columnIndex + 1);
  }
  std::vector<DesignEntry> relationEntries;
  for(int basisIndex = 0; basisIndex < basisRank; ++basisIndex){
    relationEntries.emplace_back(basisIndex, retainedColumns[basisIndex], 1.0);
  }
  // X_original = X_retained * relation, without a dense rank-by-column workspace.
  const SparseDesign triangularFactor = basisFactor.matrixR().topLeftCorner(basisRank, basisRank);
  const SparseDesign fullFactor = basisFactor.matrixR();
  for(int pivotIndex = basisRank; pivotIndex < static_cast<int>(nonzeroColumns.size()); ++pivotIndex){
    if(pivotIndex % 32 == 0) Rcpp::checkUserInterrupt();
    const int originalColumn = nonzeroColumns[basisPermutation[pivotIndex]];
    Eigen::VectorXd columnRhs = Eigen::VectorXd::Zero(basisRank);
    for(SparseDesign::InnerIterator entry(fullFactor, pivotIndex); entry; ++entry){
      if(entry.row() < basisRank) columnRhs[entry.row()] = entry.value();
    }
    const Eigen::VectorXd dependence = triangularFactor.triangularView<Eigen::Upper>().solve(columnRhs);
    if(!dependence.allFinite()){
      Rcpp::stop("Non-finite dependence coefficients in the fixed-effect basis; rescale covariates.");
    }
    for(int basisIndex = 0; basisIndex < basisRank; ++basisIndex){
      const int originalBasis = nonzeroColumns[basisPermutation[basisIndex]];
      const double relationValue = dependence[basisIndex] *
        (columnScales[originalColumn] / columnScales[originalBasis]);
      if(!std::isfinite(relationValue)){
        Rcpp::stop("Fixed-effect dependence coefficients exceed the numerical range; rescale covariates.");
      }
      if(relationValue != 0.0){
        relationEntries.emplace_back(retainedPositions[originalBasis], originalColumn, relationValue);
      }
    }
  }
  SparseDesign columnRelation(basisRank, columnTotal);
  columnRelation.setFromTriplets(relationEntries.begin(), relationEntries.end());
  return Rcpp::List::create(
    Rcpp::Named("retained")=retainedIndices,
    Rcpp::Named("omitted")=Rcpp::wrap(omittedIndices),
    Rcpp::Named("relation")=columnRelation,
    Rcpp::Named("rank")=basisRank,
    Rcpp::Named("tolerance")=rankTolerance,
    Rcpp::Named("method")="column-scaled sparse QR (COLAMD)");
}
