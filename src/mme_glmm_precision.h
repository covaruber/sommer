#ifndef SOMMER_MME_GLMM_PRECISION_H
#define SOMMER_MME_GLMM_PRECISION_H

#include "mme_cholmod.h"
#include "mme_sp2_support.h"

namespace sommer {

inline Eigen::SparseMatrix<double> laplaceEigenSparse(const arma::sp_mat & source){
  std::vector<Eigen::Triplet<double>> entries;
  entries.reserve(source.n_nonzero);
  for(auto cell = source.begin(); cell != source.end(); ++cell){
    entries.emplace_back(static_cast<int>(cell.row()),
      static_cast<int>(cell.col()), *cell);
  }
  Eigen::SparseMatrix<double> converted(
    static_cast<int>(source.n_rows), static_cast<int>(source.n_cols));
  converted.setFromTriplets(entries.begin(), entries.end());
  converted.makeCompressed();
  return converted;
}

inline bool laplaceFactorLogDeterminant(CholmodState & factor, double & value){
  // factor_to_sparse consumes numeric storage; preserve the reusable factor.
  cholmod_factor * factorCopy = M_cholmod_copy_factor(factor.factor, &factor.common);
  if(factorCopy == nullptr) return false;
  cholmod_sparse * sparseFactor =
    M_cholmod_factor_to_sparse(factorCopy, &factor.common);
  M_cholmod_free_factor(&factorCopy, &factor.common);
  if(sparseFactor == nullptr || factor.common.status != CHOLMOD_OK){
    if(sparseFactor != nullptr) M_cholmod_free_sparse(&sparseFactor, &factor.common);
    return false;
  }
  const int * columnStart = static_cast<const int *>(sparseFactor->p);
  const int * rowIndex = static_cast<const int *>(sparseFactor->i);
  const double * entries = static_cast<const double *>(sparseFactor->x);
  value = 0.0;
  bool valid = true;
  for(int columnIndex = 0; columnIndex < static_cast<int>(sparseFactor->ncol);
      ++columnIndex){
    double diagonalValue = 0.0;
    for(int entryIndex = columnStart[columnIndex];
        entryIndex < columnStart[columnIndex + 1]; ++entryIndex){
      if(rowIndex[entryIndex] == columnIndex){
        diagonalValue = entries[entryIndex];
        break;
      }
    }
    if(!std::isfinite(diagonalValue) || diagonalValue <= 0.0){
      valid = false;
      break;
    }
    value += 2.0 * std::log(diagonalValue);
  }
  M_cholmod_free_sparse(&sparseFactor, &factor.common);
  return valid && std::isfinite(value);
}

struct LaplacePrecisionDerivative {
  arma::uword start;
  arma::sp_mat block;
  double logDet;
};

class LaplacePrecision {
  CovarianceFactorEngine engine;
  Rcpp::List descriptors;
  std::vector<arma::sp_mat> relationships;
  std::vector<double> relationshipLogDet;
  std::vector<FactorPrecisionState> states;
  std::vector<arma::vec> cachedParameters;
  std::vector<arma::sp_mat> cachedBlocks;
  std::vector<bool> ready;
  std::vector<bool> derivativesReady;
  bool matrixReady = false;

  arma::sp_mat descriptorMatrix(const FactorPrecisionState & state,
      const int derivativeIndex = -1){
    arma::sp_mat product(1, 1);
    product(0, 0) = state.inverseScale;
    const int derivativeAxis = derivativeIndex < 0 ? -2 :
      state.derivativeFactor[static_cast<std::size_t>(derivativeIndex)];
    if(derivativeAxis == -1) product *= -1.0;
    for(std::size_t axisIndex = 0; axisIndex < state.factorPrecision.size();
        ++axisIndex){
      arma::sp_mat factor;
      if(static_cast<int>(axisIndex) == derivativeAxis){
        const std::size_t parameterIndex = static_cast<std::size_t>(derivativeIndex);
        factor = state.sparseLocalPrecisionD1[parameterIndex].is_empty()
          ? arma::sp_mat(state.localPrecisionD1[parameterIndex])
          : state.sparseLocalPrecisionD1[parameterIndex];
      }else{
        factor = state.sparseFactorPrecision[axisIndex].is_empty()
          ? arma::sp_mat(state.factorPrecision[axisIndex])
          : state.sparseFactorPrecision[axisIndex];
      }
      product = arma::kron(product, factor);
    }
    return product;
  }

public:
  arma::sp_mat matrix;
  double logDet = 0.0;
  int relationshipFactorizations = 0;
  int descriptorUpdates = 0;
  int descriptorReuses = 0;

  arma::uword parameterCount() const {
    arma::uword result = 0;
    for(int termPosition = 0; termPosition < descriptors.size(); ++termPosition){
      const Rcpp::List descriptor = descriptors[termPosition];
      const Rcpp::LogicalVector freeFlags = descriptor["free"];
      for(int parameterPosition = 0; parameterPosition < freeFlags.size(); ++parameterPosition)
        result += freeFlags[parameterPosition] == TRUE;
    }
    return result;
  }

  LaplacePrecision(const Rcpp::List & suppliedDescriptors,
      const Rcpp::List & suppliedRelationships)
      : descriptors(Rcpp::clone(suppliedDescriptors)){
    if(descriptors.size() == 0 || descriptors.size() != suppliedRelationships.size()){
      Rcpp::stop("Native Laplace descriptors and relationships must have matching terms.");
    }
    for(int termIndex = 0; termIndex < descriptors.size(); ++termIndex){
      arma::sp_mat relationship = Rcpp::as<arma::sp_mat>(suppliedRelationships[termIndex]);
      if(relationship.n_rows == 0 || relationship.n_rows != relationship.n_cols ||
         !relationship.is_finite() || !relationship.is_symmetric(1e-12)){
        Rcpp::stop("Native Laplace relationship precision must be finite, square and symmetric.");
      }
      CholmodState relationshipFactor;
      const auto eigenRelationship = laplaceEigenSparse(relationship);
      double determinantValue = NA_REAL;
      if(!relationshipFactor.analyzePattern(eigenRelationship) ||
         !relationshipFactor.factorize(eigenRelationship) ||
         !laplaceFactorLogDeterminant(relationshipFactor, determinantValue)){
        Rcpp::stop("Native Laplace relationship precision must be positive definite.");
      }
      ++relationshipFactorizations;
      relationships.push_back(std::move(relationship));
      relationshipLogDet.push_back(determinantValue);
    }
    states.resize(descriptors.size());
    cachedParameters.resize(descriptors.size());
    cachedBlocks.resize(descriptors.size());
    ready.assign(descriptors.size(), false);
    derivativesReady.assign(descriptors.size(), false);
  }

  bool update(const arma::vec & freeParameters, std::string & reason){
    arma::uword freeIndex = 0;
    arma::uword totalDimension = 0;
    arma::uword totalNonzeros = 0;
    logDet = 0.0;
    for(int termIndex = 0; termIndex < descriptors.size(); ++termIndex){
      const Rcpp::List descriptor = descriptors[termIndex];
      arma::vec parameters = Rcpp::as<arma::vec>(descriptor["par"]);
      const Rcpp::LogicalVector freeFlags = descriptor["free"];
      if(parameters.n_elem == 0 ||
         freeFlags.size() != static_cast<int>(parameters.n_elem) ||
         Rcpp::is_true(Rcpp::any(Rcpp::is_na(freeFlags)))){
        Rcpp::stop("Native Laplace descriptor free flags have incompatible length.");
      }
      for(arma::uword parameterIndex = 0; parameterIndex < parameters.n_elem;
          ++parameterIndex){
        if(freeFlags[parameterIndex] == TRUE){
          if(freeIndex >= freeParameters.n_elem){
            Rcpp::stop("Native Laplace covariance parameter vector has incompatible length.");
          }
          parameters(parameterIndex) = freeParameters(freeIndex++);
        }
      }
      if(!parameters.is_finite()){
        reason = "non-finite covariance parameters";
        return false;
      }
      const double exponentialLimit = std::log(std::numeric_limits<double>::max()) / 4.0;
      bool outsideRange = std::abs(parameters(0)) > exponentialLimit;
      const Rcpp::List factors = descriptor["factors"];
      for(int axisIndex = 0; axisIndex < factors.size(); ++axisIndex){
        const Rcpp::List factor = factors[axisIndex];
        const Rcpp::List report = factor["report"];
        if(Rcpp::as<std::string>(report["backend"]) != "builtin") continue;
        const Rcpp::CharacterVector transforms = report["transform"];
        const int firstParameter = Rcpp::as<int>(factor["par_start"]) - 1;
        for(int localIndex = 0; localIndex < transforms.size(); ++localIndex){
          if(Rcpp::as<std::string>(transforms[localIndex]) == "exp" &&
             std::abs(parameters(firstParameter + localIndex)) > exponentialLimit){
            outsideRange = true;
          }
        }
      }
      if(outsideRange){
        reason = "exponential covariance coordinate outside numerical range";
        return false;
      }
      if(ready[termIndex] && arma::all(parameters == cachedParameters[termIndex])){
        ++descriptorReuses;
      }else{
        matrixReady = false;
        ready[termIndex] = false;
        derivativesReady[termIndex] = false;
        if(!engine.buildDescriptorPrecision(descriptor, parameters, states[termIndex])){
          reason = "non-positive-definite covariance";
          return false;
        }
        cachedBlocks[termIndex] = arma::kron(
          descriptorMatrix(states[termIndex]), relationships[termIndex]);
        if(!cachedBlocks[termIndex].is_finite()){
          reason = "non-finite random precision";
          return false;
        }
        cachedParameters[termIndex] = parameters;
        ready[termIndex] = true;
        ++descriptorUpdates;
      }
      const auto & state = states[termIndex];
      logDet += state.totalDim * relationshipLogDet[termIndex] -
        relationships[termIndex].n_rows * state.logDet;
      totalDimension += cachedBlocks[termIndex].n_rows;
      totalNonzeros += cachedBlocks[termIndex].n_nonzero;
    }
    if(freeIndex != freeParameters.n_elem){
      Rcpp::stop("Native Laplace covariance parameter vector has incompatible length.");
    }
    if(!std::isfinite(logDet)){
      reason = "non-finite random precision determinant";
      return false;
    }
    if(matrixReady) return true;
    arma::umat locations(2, totalNonzeros);
    arma::vec values(totalNonzeros);
    arma::uword entryIndex = 0;
    arma::uword termStart = 0;
    for(const auto & block : cachedBlocks){
      for(auto cell = block.begin(); cell != block.end(); ++cell){
        locations(0, entryIndex) = termStart + cell.row();
        locations(1, entryIndex) = termStart + cell.col();
        values(entryIndex++) = *cell;
      }
      termStart += block.n_rows;
    }
    matrix = arma::sp_mat(locations, values, totalDimension, totalDimension);
    matrixReady = true;
    return true;
  }

  std::vector<LaplacePrecisionDerivative> derivatives(){
    std::vector<LaplacePrecisionDerivative> result;
    arma::uword termStart = 0;
    for(int termIndex = 0; termIndex < descriptors.size(); ++termIndex){
      const Rcpp::List descriptor = descriptors[termIndex];
      const Rcpp::LogicalVector freeFlags = descriptor["free"];
      if(Rcpp::is_true(Rcpp::any(freeFlags)) && !derivativesReady[termIndex]){
        engine.buildDescriptorPrecisionDerivatives(descriptor,
          cachedParameters[termIndex], states[termIndex]);
        derivativesReady[termIndex] = true;
      }
      for(int parameterIndex = 0; parameterIndex < freeFlags.size(); ++parameterIndex){
        if(freeFlags[parameterIndex] != TRUE) continue;
        result.push_back({termStart, arma::kron(
          descriptorMatrix(states[termIndex], parameterIndex), relationships[termIndex]),
          -static_cast<double>(relationships[termIndex].n_rows) *
            states[termIndex].logDetD1[parameterIndex]});
      }
      termStart += cachedBlocks[termIndex].n_rows;
    }
    return result;
  }
};

} // namespace sommer

#endif
