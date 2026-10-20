#ifndef SOMMER_MME_SP2_SUPPORT_H
#define SOMMER_MME_SP2_SUPPORT_H

// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// RcppArmadillo pulls in Rcpp.h for us
#define ARMA_DONT_PRINT_ERRORS

#include <RcppArmadillo.h>
#include "stdlib.h"

#include <progress.hpp>

// Eigen sparse LDLT used by ai_mme_sp()
#include <Eigen/SparseCore>
#include <Eigen/SparseCholesky>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/Eigenvalues>
#include <cctype>

// Nested-dissection sparse ordering for ai_mme_sp2()'s LDLT factorisations,
// enabled only when configure detected a usable METIS installation; falls
// back to Eigen's built-in AMD ordering otherwise.
#ifdef SOMMER_HAVE_METIS
#include <Eigen/MetisSupport>
typedef Eigen::MetisOrdering<int> SommerSparseOrdering;
#else
typedef Eigen::AMDOrdering<int> SommerSparseOrdering;
#endif
typedef Eigen::SimplicialLDLT<
  Eigen::SparseMatrix<double>,
  Eigen::Lower,
  SommerSparseOrdering
> SommerSpatialLDLT;

// Supernodal (BLAS-3) sparse Cholesky backend for ai_mme_sp2()'s
// solver="cholmod" option, provided by R's Matrix package. Matrix.h declares
// the M_cholmod_*() wrappers; stubs.c defines them via lazy R_GetCCallable()
// lookups into the already-loaded Matrix package, so no extra link flags
// are required beyond LinkingTo: Matrix.
#include <Matrix/Matrix.h>


// Standard C++ headers used by the new implementation
#include <vector>
#include <limits>
#include <cmath>
#include <algorithm>
#include <atomic>
#include <cstring>
#include <functional>
#include <map>
#include <memory>
#include <unordered_map>
#include <cstdint>
#include <chrono>
#include <ctime>
#include <iomanip>
#include <sstream>
#include "mme_paths.h"
#include "mme_control_r.h"
#include "mme_factors.h"

#ifdef _WIN32
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#endif

#ifdef _OPENMP
#include <omp.h>
#endif


const std::string currentDateTime();
arma::vec seqCpp(const int & a,
                 const int & b);
arma::vec mat_to_vecCpp(const arma::mat & x,
                        const arma::mat & x2);
arma::mat vec_to_matCpp(const arma::vec & x,
                        const arma::mat & x2);
arma::cube vec_to_cubeCpp(const arma::vec & x,
                          const Rcpp::List & g);
arma::vec varCols(const arma::mat & x);
arma::mat scaleCpp(const arma::mat & x);
arma::mat makeFull(const arma::mat & X);
bool isIdentity_mat(const arma::mat x);
bool isIdentity_spmat(const arma::sp_mat x);
bool isDiagonal_mat(const arma::mat x);
bool isDiagonal_spmat(const arma::sp_mat x);
arma::mat amat(const arma::mat & Xo, const bool & vanraden, double minMAF);
arma::mat dmat(const arma::mat & Xo, const bool & nishio, double minMAF);
arma::mat emat(const arma::mat & X1, const arma::mat & X2);
arma::mat emat_genotypes(const arma::mat & Xo, bool nishio,
                         const std::string & type, double minMAF);
arma::mat hmat(const arma::mat & A, const arma::mat & G22,
               const arma::vec & index, double tolparinv,
               double tau, double omega);
arma::cube scorecalc(const arma::mat & Mimv,
                     const arma::mat & Ymv, // Y is provided as multitrait
                     const arma::mat & Zmv, // Z is provided as univariate
                     const arma::mat & Xmv, // X is provided as univariate
                     const arma::mat & Vinv, // multivariate inverse of V
                     int nt, double minMAF
);
arma::cube gwasForLoop(const arma::mat & M, // marker matrix
                       const arma::mat & Y, // Y is provided as multitrait
                       const arma::mat & Z, // Z is provided as univariate
                       const arma::mat & X, // X is provided as univariate
                       const arma::mat & Vinv, // multivariate inverse of V
                       double minMAF,
                       bool display_progress);
Rcpp::List newton_di_sp(const arma::sp_mat & Y, const Rcpp::List & X,
                        const Rcpp::List & Gx,
                        const Rcpp::List & Z, const Rcpp::List & K,
                        const Rcpp::List & R, 
                        const Rcpp::List & Ge, const Rcpp::List & GeI, // theta and thetaC
                        const arma::sp_mat & W, const bool & isInvW,
                        int iters, double tolpar, double tolparinv,
                        const bool & ai, const bool & pev,
                        const bool & verbose,const bool & retscaled,
                        const arma::vec & stepweight, // const arma::vec & emupdate,
                        const arma::vec & emweight, const Rcpp::List & thetaConstOri,
                        const arma::vec & thetaIndex);
arma::sp_mat convertSparse(Rcpp::S4 mat);
arma::vec mat_to_vecCpp2(const arma::mat & x,
                         const arma::mat & x2);
arma::mat nearPDcpp(const arma::mat X0,
                    const int & maxit,
                    const double & eig_tol,
                    const double & conv_tol);
Rcpp::List post_mme_Cinverse_cpp(Rcpp::List model, const int mode);
arma::mat predict_mmes_vcov_cpp(Rcpp::List model, const arma::sp_mat & Dmat);
Rcpp::List ai_mme_sp(const arma::sp_mat & X, const Rcpp::List & ZI,  const arma::vec & Zind,
                     const Rcpp::List & AiI, const arma::sp_mat & y0,
                     const Rcpp::List & SI, const Rcpp::List & partitionsS,
                     const arma::sp_mat & H, const bool & useH,
                     int nIters, double tolParConvLL, double tolParConvNorm,
                     double tolParInv, const Rcpp::List & thetaI,
                     const Rcpp::List & thetaCI, const arma::mat & thetaF,
                     const arma::vec & addScaleParam, const arma::vec & weightEmInf,
                     const arma::vec & weightInf, const bool & verbose
);
Rcpp::List MNR(const arma::mat & Y, const Rcpp::List & X,
               const Rcpp::List & Gx,
               const Rcpp::List & Z, const Rcpp::List & K,
               const Rcpp::List & R, const Rcpp::List & Ge,
               const Rcpp::List & GeI, const arma::mat & W, const bool & isInvW,
               int iters, double tolpar, double tolparinv,
               const bool & ai, const bool & pev,
               const bool & verbose,const bool & retscaled,
               const arma::vec & stepweight, // const arma::vec & emupdate,
               const arma::vec & emweight);

static bool sparseFactorKernelsEnabled(const Rcpp::List & factor){
  return factor.containsElementNamed("sparse_kernels") && Rcpp::as<bool>(factor["sparse_kernels"]);
}

static inline arma::sp_mat eigenSparseToArmaGlobal(const Eigen::SparseMatrix<double> & src);

struct SparseSpatialFactor {
  arma::mat covariance;
  arma::sp_mat precision;
  arma::mat covarianceD1;
  arma::sp_mat precisionD1;
  std::shared_ptr<SommerSpatialLDLT> factorization;
  double normalization = 1.0;
  double logDet = 0.0;
  double logDetD1 = 0.0;
};

static SparseSpatialFactor evaluateSparseSpatialFactor(const Rcpp::List & factor,
                                                       const arma::vec & parameters,
                                                       const bool derivatives,
                                                       const SparseSpatialFactor * cached = nullptr){
  const bool normalEquations = sommer::precisionFamily(factor) == sommer::PrecisionFamily::NormalEquations;
  const Rcpp::List capability = sommer::factorCapability(factor, "spatial");
  const arma::mat weights = Rcpp::as<arma::mat>(capability["weights"]);
  const int dimension = static_cast<int>(weights.n_rows);
  if(dimension < 2 || weights.n_cols != weights.n_rows || !weights.is_finite() || parameters.n_elem != 1){
    Rcpp::stop("Malformed sparse spatial factor.");
  }
  const double eta = parameters(0);
  const double logistic = eta >= 0.0 ? 1.0 / (1.0 + std::exp(-eta))
    : std::exp(eta) / (1.0 + std::exp(eta));
  const double lower = Rcpp::as<double>(capability["lower"]);
  const double upper = Rcpp::as<double>(capability["upper"]);
  const double rho = lower + (upper - lower) * logistic;
  const double rhoD1 = (upper - lower) * logistic * (1.0 - logistic);
  std::vector<Eigen::Triplet<double> > entries;
  for(int column = 0; column < dimension; ++column){
    for(int row = 0; row < dimension; ++row){
      if(weights(row, column) != 0.0){ entries.emplace_back(row, column, weights(row, column)); }
    }
  }
  Eigen::SparseMatrix<double> adjacency(dimension, dimension);
  adjacency.setFromTriplets(entries.begin(), entries.end());
  Eigen::SparseMatrix<double> diagonal(dimension, dimension);
  diagonal.setIdentity();
  if(!normalEquations){
    const arma::vec degree = Rcpp::as<arma::vec>(capability["degree"]);
    if(degree.n_elem != weights.n_rows || !degree.is_finite() || degree.min() <= 0.0){
      Rcpp::stop("Malformed sparse CAR degrees.");
    }
    for(int row = 0; row < dimension; ++row){ diagonal.coeffRef(row, row) = degree(row); }
  }
  Eigen::SparseMatrix<double> system = diagonal - rho * adjacency;
  Eigen::SparseMatrix<double> precision = normalEquations
    ? Eigen::SparseMatrix<double>(system.transpose() * system) : system;
  SparseSpatialFactor result;
  arma::mat covariance;
  if(cached && cached->factorization && cached->covariance.n_rows == static_cast<arma::uword>(dimension)){
    result = *cached;
    covariance = cached->covariance;
  }else{
    result.factorization = std::make_shared<SommerSpatialLDLT>();
    result.factorization->compute(precision);
    if(result.factorization->info() != Eigen::Success ||
       !result.factorization->vectorD().allFinite() ||
       result.factorization->vectorD().minCoeff() <= 0.0){
      Rcpp::stop("Sparse spatial precision factorisation failed.");
    }
    covariance.set_size(dimension, dimension);
    const int panelColumns = static_cast<int>(sommer::WorkspacePolicy::panelColumns());
    for(int first = 0; first < dimension; first += panelColumns){
      const int count = std::min(panelColumns, dimension - first);
      Eigen::MatrixXd rhs = Eigen::MatrixXd::Zero(dimension, count);
      for(int column = 0; column < count; ++column){ rhs(first + column, column) = 1.0; }
      const Eigen::MatrixXd solution = result.factorization->solve(rhs);
      for(int column = 0; column < count; ++column){
        std::copy(solution.col(column).data(), solution.col(column).data() + dimension,
          covariance.colptr(first + column));
      }
    }
    covariance = 0.5 * (covariance + covariance.t());
    result.normalization = covariance(0, 0);
    if(!covariance.is_finite() || !std::isfinite(result.normalization) || result.normalization <= 0.0){
      Rcpp::stop("Sparse spatial normalization solve failed.");
    }
    covariance /= result.normalization;
    result.precision = result.normalization * eigenSparseToArmaGlobal(precision);
    result.logDet = -result.factorization->vectorD().array().log().sum() -
      dimension * std::log(result.normalization);
    result.covariance = covariance;
  }
  if(derivatives){
    Eigen::SparseMatrix<double> systemD1 = -rhoD1 * adjacency;
    Eigen::SparseMatrix<double> precisionD1 = normalEquations
      ? Eigen::SparseMatrix<double>(systemD1.transpose() * system + system.transpose() * systemD1) : systemD1;
    const double normalization = result.normalization;
    Eigen::Map<const Eigen::VectorXd> normalizedFirstColumn(covariance.memptr(), dimension);
    const Eigen::VectorXd firstColumn = normalization * normalizedFirstColumn;
    const double normalizationD1 = -firstColumn.dot(precisionD1 * firstColumn);
    double determinantD1 = 0.0;
    for(int column = 0; column < dimension; ++column){
      for(Eigen::SparseMatrix<double>::InnerIterator entry(precisionD1, column); entry; ++entry){
        determinantD1 += entry.value() * covariance(entry.col(), entry.row());
      }
    }
    determinantD1 *= normalization;
    result.precisionD1 = normalizationD1 * eigenSparseToArmaGlobal(precision) +
      normalization * eigenSparseToArmaGlobal(precisionD1);
    result.logDetD1 = -determinantD1 - dimension * normalizationD1 / normalization;
    result.covarianceD1.set_size(dimension, dimension);
    const int panelColumns = static_cast<int>(sommer::WorkspacePolicy::panelColumns());
    for(int first = 0; first < dimension; first += panelColumns){
      const int count = std::min(panelColumns, dimension - first);
      Eigen::Map<const Eigen::MatrixXd> panel(covariance.colptr(first), dimension, count);
      const Eigen::MatrixXd rhs = -precisionD1 * panel;
      const Eigen::MatrixXd solution = result.factorization->solve(rhs);
      for(int column = 0; column < count; ++column){
        std::copy(solution.col(column).data(), solution.col(column).data() + dimension,
          result.covarianceD1.colptr(first + column));
      }
    }
    result.covarianceD1 -= covariance * (normalizationD1 / normalization);
    result.covarianceD1 = 0.5 * (result.covarianceD1 + result.covarianceD1.t());
  }
  return result;
}

static arma::mat anteInnovationCovariance(const arma::uword dimension,
                                         const Rcpp::IntegerVector & rows,
                                         const Rcpp::IntegerVector & columns,
                                         const arma::vec & parameters,
                                         const bool sparseKernels){
  const arma::uword coefficients = static_cast<arma::uword>(rows.size());
  if(!sparseKernels){
    arma::mat transition(dimension, dimension, arma::fill::eye);
    for(arma::uword coefficient = 0; coefficient < coefficients; ++coefficient){
      transition(rows[coefficient] - 1, columns[coefficient] - 1) = -parameters(coefficient);
    }
    arma::vec innovation(dimension, arma::fill::ones);
    for(arma::uword row = 1; row < dimension; ++row){ innovation(row) = std::exp(parameters(coefficients + row - 1)); }
    const arma::mat inverse = arma::inv(arma::trimatl(transition));
    arma::mat covariance = inverse * arma::diagmat(innovation) * inverse.t();
    return covariance / covariance(0, 0);
  }
  std::vector<std::vector<std::pair<arma::uword, double> > > innovationRows(dimension);
  for(arma::uword coefficient = 0; coefficient < coefficients; ++coefficient){
    const int row = rows[coefficient] - 1;
    const int column = columns[coefficient] - 1;
    if(row <= column || column < 0 || row >= static_cast<int>(dimension)){
      Rcpp::stop("Malformed antedependence innovation coordinates.");
    }
    innovationRows[static_cast<arma::uword>(row)].emplace_back(
      static_cast<arma::uword>(column), parameters(coefficient));
  }
  arma::mat covariance(dimension, dimension, arma::fill::zeros);
  for(arma::uword row = 0; row < dimension; ++row){
    auto & active = innovationRows[row];
    std::sort(active.begin(), active.end());
    for(arma::uword column = 0; column < row; ++column){
      double value = 0.0;
      for(const auto & coefficient : active){ value += coefficient.second * covariance(coefficient.first, column); }
      covariance(row, column) = value;
      covariance(column, row) = value;
    }
    double variance = row == 0 ? 1.0 : std::exp(parameters(coefficients + row - 1));
    for(const auto & coefficient : active){ variance += coefficient.second * covariance(row, coefficient.first); }
    covariance(row, row) = variance;
  }
  return covariance;
}

static arma::sp_mat antePrecisionDerivative(const arma::uword dimension,
                                           const Rcpp::IntegerVector & rows,
                                           const Rcpp::IntegerVector & columns,
                                           const arma::vec & parameters,
                                           const arma::uword derivative){
  const arma::uword coefficients = static_cast<arma::uword>(rows.size());
  if(derivative >= coefficients + dimension - 1){
    Rcpp::stop("Invalid ANTE precision derivative index.");
  }
  std::vector<std::vector<std::pair<arma::uword, double> > > innovationRows(dimension);
  for(arma::uword row = 0; row < dimension; ++row){ innovationRows[row].emplace_back(row, 1.0); }
  for(arma::uword coefficient = 0; coefficient < coefficients; ++coefficient){
    innovationRows[static_cast<arma::uword>(rows[coefficient] - 1)].emplace_back(
      static_cast<arma::uword>(columns[coefficient] - 1), -parameters(coefficient));
  }
  arma::sp_mat derivativePrecision(dimension, dimension);
  if(derivative < coefficients){
    const arma::uword row = static_cast<arma::uword>(rows[derivative] - 1);
    const arma::uword column = static_cast<arma::uword>(columns[derivative] - 1);
    const double weight = row == 0 ? 1.0 : std::exp(-parameters(coefficients + row - 1));
    for(const auto & entry : innovationRows[row]){
      derivativePrecision(column, entry.first) -= weight * entry.second;
      derivativePrecision(entry.first, column) -= weight * entry.second;
    }
  }else{
    const arma::uword row = derivative - coefficients + 1;
    const double weight = std::exp(-parameters(derivative));
    for(const auto & first : innovationRows[row]){
      for(const auto & second : innovationRows[row]){
        derivativePrecision(first.first, second.first) -=
          weight * first.second * second.second;
      }
    }
  }
  derivativePrecision.sync();
  return derivativePrecision;
}

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

// simple example of creating two matrices and
// returning the result of an operatioon on them
//
// via the exports attribute we tell Rcpp to make this function
// available from R
//



































// [[Rcpp::depends(RcppProgress)]]




















// Reusable arma::sp_mat <-> Eigen::SparseMatrix<double> bridges for the
// pieces of ai_mme_sp2() that are being migrated off Armadillo incrementally.
static inline Eigen::SparseMatrix<double> armaSparseToEigenGlobal(
    const arma::sp_mat & src){
  std::vector<Eigen::Triplet<double>> triplets;
  triplets.reserve(src.n_nonzero);
  for(arma::sp_mat::const_iterator it = src.begin(); it != src.end(); ++it){
    triplets.emplace_back(
      static_cast<int>(it.row()),
      static_cast<int>(it.col()),
      *it
    );
  }
  Eigen::SparseMatrix<double> out(
    static_cast<int>(src.n_rows),
    static_cast<int>(src.n_cols)
  );
  out.setFromTriplets(triplets.begin(), triplets.end());
  out.makeCompressed();
  return out;
}

static inline arma::sp_mat eigenSparseToArmaGlobal(
    const Eigen::SparseMatrix<double> & input){
  // Compressed Eigen storage is already CSC with sorted rows: copy it directly
  // instead of sorting coordinate triplets.
  Eigen::SparseMatrix<double> compressedCopy;
  if(!input.isCompressed()){
    compressedCopy = input;
    compressedCopy.makeCompressed();
  }
  const Eigen::SparseMatrix<double> & src = input.isCompressed() ? input : compressedCopy;
  const arma::uword nnz = static_cast<arma::uword>(src.nonZeros());
  const arma::uword columns = static_cast<arma::uword>(src.outerSize());
  arma::uvec rowIndices(nnz), columnPointers(columns + 1);
  arma::vec values(nnz);
  for(arma::uword column = 0; column <= columns; ++column){
    columnPointers(column) = static_cast<arma::uword>(src.outerIndexPtr()[column]);
  }
  for(arma::uword k = 0; k < nnz; ++k){
    rowIndices(k) = static_cast<arma::uword>(src.innerIndexPtr()[k]);
    values(k) = src.valuePtr()[k];
  }
  return arma::sp_mat(
    rowIndices,
    columnPointers,
    values,
    static_cast<arma::uword>(src.rows()),
    static_cast<arma::uword>(src.cols())
  );
}

// Sparse LDLT log-determinant of a sparse SPD precision/relationship matrix,
// avoiding the O(n^2) densification that arma::log_det() would require.
static inline double sparseSpdLogDet(const arma::sp_mat & A){
  const Eigen::SparseMatrix<double> Ae = armaSparseToEigenGlobal(A);
  Eigen::SimplicialLDLT<
    Eigen::SparseMatrix<double>,
    Eigen::Lower,
    SommerSparseOrdering
  > factor;
  factor.compute(Ae);
  if(factor.info() != Eigen::Success){
    Rcpp::stop("Sparse LDLT factorisation failed while computing a log-determinant.");
  }
  const Eigen::VectorXd pivots = factor.vectorD();
  double logDet = 0.0;
  for(Eigen::Index k = 0; k < pivots.size(); ++k){
    const double pivot = pivots(k);
    if(!std::isfinite(pivot) || pivot <= 0.0){
      Rcpp::stop("Non-positive LDLT pivot while computing a sparse log-determinant.");
    }
    logDet += std::log(pivot);
  }
  return logDet;
}

// Dense SPD inverse via Eigen::LLT, used in place of arma::inv_sympd() for
// the small per-random-effect covariance matrices in ai_mme_sp2().
static inline bool eigenSpdInverse(const arma::mat & A, arma::mat & out,
                                   double * logDet = nullptr){
  const arma::uword n = A.n_rows;
  // Armadillo and Eigen dense matrices are both column-major, so the
  // conversion is a straight buffer copy rather than an element loop.
  const Eigen::Map<const Eigen::MatrixXd> Ae(A.memptr(), n, n);
  Eigen::LLT<Eigen::MatrixXd> llt(Ae);
  if(llt.info() != Eigen::Success){
    return false;
  }
  if(logDet != nullptr){
    const Eigen::MatrixXd L = llt.matrixL();
    double value = 0.0;
    for(Eigen::Index k = 0; k < L.rows(); ++k){
      const double diagonal = L(k,k);
      if(!std::isfinite(diagonal) || diagonal <= 0.0){
        return false;
      }
      value += 2.0 * std::log(diagonal);
    }
    *logDet = value;
  }
  const Eigen::MatrixXd inv = llt.solve(Eigen::MatrixXd::Identity(n, n));
  out.set_size(n, n);
  std::copy(inv.data(), inv.data() + inv.size(), out.memptr());
  return true;
}

// Affine map between working covariance parameters and reduced parameters,
// w = o + T phi, supplied by mmes(vcc=) as the "vcmap" attribute of covStructI.
struct HybridAIState {
  arma::mat metric;
  arma::vec parameters;
  arma::vec gradient;
  arma::uvec free;
  bool ready = false;

  void initialize(const arma::mat & information, const arma::vec & point,
                  const arma::vec & score, const arma::uvec & indices){
    metric = information;
    parameters = point;
    gradient = score;
    free = indices;
    ready = free.n_elem > 0;
    if(!ready){ return; }
    arma::mat block = metric.submat(free, free);
    block = 0.5 * (block + block.t());
    arma::vec eigenvalues;
    if(!arma::eig_sym(eigenvalues, block) || !eigenvalues.is_finite()){
      Rcpp::stop("Unable to initialize hybrid AI curvature.");
    }
    const double floor = 1.0e-8 * std::max(1.0, arma::abs(eigenvalues).max());
    if(eigenvalues.min() < floor){
      block.diag() += floor - eigenvalues.min();
    }
    metric.submat(free, free) = block;
  }

  bool update(const arma::vec & point, const arma::vec & score,
              const arma::uvec & indices){
    if(!ready || indices.n_elem != free.n_elem || arma::any(indices != free)){
      return false;
    }
    const arma::vec step = point(free) - parameters(free);
    arma::vec change = score(free) - gradient(free);
    if(!step.is_finite() || !change.is_finite() ||
       arma::norm(step, 2) <= std::numeric_limits<double>::epsilon() *
         (1.0 + arma::norm(point(free), 2))){
      return false;
    }
    arma::mat block = metric.submat(free, free);
    const arma::vec product = block * step;
    const double predicted = arma::dot(step, product);
    const double observed = arma::dot(step, change);
    if(!std::isfinite(predicted) || !std::isfinite(observed) || predicted <= 0.0){
      return false;
    }
    // Powell damping preserves positive curvature even away from a local minimum.
    if(observed < 0.2 * predicted){
      const double weight = 0.8 * predicted / (predicted - observed);
      change = weight * change + (1.0 - weight) * product;
    }
    const double curvature = arma::dot(step, change);
    if(!std::isfinite(curvature) || curvature <= 0.0){ return false; }
    block += (change * change.t()) / curvature - (product * product.t()) / predicted;
    block = 0.5 * (block + block.t());
    if(!block.is_finite()){ return false; }
    metric.submat(free, free) = block;
    parameters = point;
    gradient = score;
    return true;
  }
};

struct VarianceParameterMap {
  bool active = false;
  arma::mat T;
  arma::vec o;
};

static VarianceParameterMap readVarianceParameterMap(const Rcpp::List & covStructI,
                                                     const arma::uword nVcTotal){
  VarianceParameterMap map;
  SEXP attribute = covStructI.attr("vcmap");
  if(Rf_isNull(attribute)){
    return map;
  }
  Rcpp::List spec(attribute);
  map.T = Rcpp::as<arma::mat>(spec["T"]);
  map.o = Rcpp::as<arma::vec>(spec["o"]);
  if(map.T.n_rows != nVcTotal || map.o.n_elem != nVcTotal){
    Rcpp::stop("The variance-parameter constraint map does not match the covariance parameters.");
  }
  map.active = true;
  return map;
}

// Project a proposal onto the constraint manifold over the rows in idx using
// the metric Imetric; with the information matrix this is the constrained
// Newton step, with the identity it is an exact Euclidean repair.
static arma::vec projectVarianceParameters(const VarianceParameterMap & map,
                                           const arma::vec & w,
                                           const arma::uvec & idx,
                                           const arma::mat & Imetric){
  if(!map.active || idx.n_elem == 0){
    return w;
  }
  arma::mat Ts = map.T.rows(idx);
  arma::uvec used = arma::find(arma::sum(arma::abs(Ts), 0).t() > 0);
  if(used.n_elem == 0){
    return w;
  }
  Ts = Ts.cols(used);
  const arma::vec r = w(idx) - map.o(idx);
  arma::mat Is = Imetric.submat(idx, idx);
  Is = 0.5 * (Is + Is.t());
  arma::mat lhs = Ts.t() * Is * Ts;
  arma::vec rhs = Ts.t() * Is * r;
  arma::vec phi;
  bool ok = arma::rcond(lhs) > 1.0e-12 && arma::solve(phi, lhs, rhs, arma::solve_opts::no_approx);
  if(!ok){
    phi = arma::solve(Ts.t() * Ts, Ts.t() * r);
  }
  arma::vec out = w;
  out(idx) = map.o(idx) + Ts * phi;
  return out;
}

static bool find_intercept_vector(const arma::sp_mat & X, arma::vec & cOut) {
  if(X.n_cols == 0 || X.n_rows == 0){
    return false;
  }
  cOut.zeros(X.n_cols);
  // Check if any single column is all 1s
  for(arma::uword j = 0; j < X.n_cols; ++j){
    if(arma::all(arma::vec(X.col(j)) == 1.0)){
      cOut(j) = 1.0;
      return true;
    }
  }
  // General: solve (X'X) c = X' 1
  arma::mat XtX = arma::mat(X.t() * X);
  arma::vec Xt1 = arma::vec(arma::sum(X, 0).t());
  arma::vec c;
  bool ok = arma::solve(c, XtX, Xt1, arma::solve_opts::fast);
  if(!ok){
    c = arma::pinv(XtX) * Xt1;
    ok = true;
  }
  if(ok){
    arma::vec fit = arma::vec(X * c);
    if(arma::max(arma::abs(fit - 1.0)) < 1.0e-4){
      cOut = c;
      return true;
    }
  }
  return false;
}

#include "mme_covariance.h"


#endif
