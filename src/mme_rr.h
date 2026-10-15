#ifndef SOMMER_MME_RR_H
#define SOMMER_MME_RR_H

#include <functional>
#include <Eigen/SparseLU>

namespace sommer {

inline arma::mat rrDescriptorDerivative(CovarianceFactorEngine & engine,
                                        const Rcpp::List & descriptor,
                                        const arma::vec & par,
                                        const arma::uword parameter) {
  if(parameter >= par.n_elem) {
    Rcpp::stop("Reduced-rank covariance derivative index is out of range.");
  }
  const arma::mat covariance = engine.evaluateDescriptor(descriptor, par);
  if(parameter == 0) { return covariance; }
  const Rcpp::List factors = descriptor["factors"];
  const double scale = std::exp(par(0));
  for(int axis = 0; axis < factors.size(); ++axis) {
    const Rcpp::List factor = factors[axis];
    const int start = Rcpp::as<int>(factor["par_start"]) - 1;
    const int end = Rcpp::as<int>(factor["par_end"]) - 1;
    if(start < 0 || end < start) { continue; }
    if(parameter < static_cast<arma::uword>(start) ||
       parameter > static_cast<arma::uword>(end)) {
      continue;
    }
    const arma::vec local = par.subvec(static_cast<arma::uword>(start),
                                       static_cast<arma::uword>(end));
    const arma::uword dimension = static_cast<arma::uword>(Rcpp::as<int>(factor["dim"]));
    arma::mat derivative;
    derivative = engine.factorD1(factor, local,
      parameter - static_cast<arma::uword>(start));
    if(derivative.n_rows != dimension || derivative.n_cols != dimension) {
      Rcpp::stop("Covariance factor derivative returned an invalid matrix.");
    }
    arma::mat product(1, 1, arma::fill::ones);
    for(int other = 0; other < factors.size(); ++other) {
      const Rcpp::List otherFactor = factors[other];
      const int otherStart = Rcpp::as<int>(otherFactor["par_start"]) - 1;
      const int otherEnd = Rcpp::as<int>(otherFactor["par_end"]) - 1;
      arma::mat piece;
      if(otherStart < 0 || otherEnd < otherStart) {
        piece = engine.evalFactor(otherFactor, arma::vec());
      } else {
        const arma::vec otherPar = par.subvec(static_cast<arma::uword>(otherStart),
                                              static_cast<arma::uword>(otherEnd));
        if(other == axis) {
          piece = derivative;
        } else {
          piece = engine.evalFactor(otherFactor, otherPar);
        }
      }
      product = arma::kron(product, piece);
    }
    return scale * product;
  }
  Rcpp::stop("Covariance derivative parameter does not belong to a factor.");
  return arma::mat();
}

struct RRAIResult {
  arma::vec score;
  arma::mat information;
  arma::uword randomTraceDirections = 0;
};

inline arma::mat reducedRankCoefficientMap(const Rcpp::List & descriptor,
                                           const arma::vec & parameters) {
  const Rcpp::List factors = descriptor["factors"];
  arma::mat mapping(1, 1, arma::fill::ones);
  for(int axis = 0; axis < factors.size(); ++axis) {
    const Rcpp::List factor = factors[axis];
    const int dimension = Rcpp::as<int>(factor["dim"]);
    arma::mat axisMap = arma::eye<arma::mat>(dimension, dimension);
    if(zeroSpecificFactor(factor)) {
      const int first = Rcpp::as<int>(factor["par_start"]) - 1;
      const int last = Rcpp::as<int>(factor["par_end"]) - 1;
      axisMap = normalizedFactorLoadings(factor, parameters.subvec(first, last));
    }
    mapping = arma::kron(mapping, axisMap);
  }
  return mapping;
}

inline RRAIResult reducedRankInformation(
    const Rcpp::List & fit,
    const arma::sp_mat & X,
    const Rcpp::List & originalZ,
    const arma::vec & Zind,
    const Rcpp::List & Ai,
    const arma::uvec & residualBlock,
    const arma::uvec & residualIndex,
    const arma::sp_mat & H,
    const bool useH,
    const Rcpp::List & descriptors,
    const arma::vec & parameters,
    const std::vector<arma::uword> & offsets,
    const bool reml) {
  const arma::sp_mat W = Rcpp::as<arma::sp_mat>(fit["W"]);
  const arma::sp_mat C = Rcpp::as<arma::sp_mat>(fit["C"]);
  const arma::vec y = Rcpp::as<arma::vec>(fit["yWorking"]);
  const arma::uword n = y.n_elem;
  const arma::uword nParameters = parameters.n_elem;
  const double covarianceScale = Rcpp::as<double>(fit["Cscale"]);
  if(W.n_rows != n || C.n_rows != C.n_cols || residualBlock.n_elem != n ||
     residualIndex.n_elem != n) {
    Rcpp::stop("Reduced-rank AI received inconsistent Henderson dimensions.");
  }
  if(useH && (H.n_rows != n || H.n_cols != n)) {
    Rcpp::stop("Reduced-rank residual transform has incompatible dimensions.");
  }

  typedef Eigen::SparseMatrix<double, Eigen::ColMajor, int> EigenSpMat;
  typedef Eigen::SimplicialLDLT<EigenSpMat, Eigen::Lower, SommerSparseOrdering> EigenLDLT;
  arma::sp_mat residualTransform;
  if(useH) {
    EigenSpMat weightPrecision = armaSparseToEigenGlobal(H);
    Eigen::SimplicialLLT<EigenSpMat, Eigen::Upper, Eigen::NaturalOrdering<int>> weightFactor;
    weightFactor.compute(weightPrecision);
    if(weightFactor.info() != Eigen::Success) {
      Rcpp::stop("Residual weight precision factorization failed in reduced-rank AI evaluation.");
    }
    EigenSpMat upper = weightFactor.matrixU();
    upper.makeCompressed();
    residualTransform = eigenSparseToArmaGlobal(upper);
  }
  EigenSpMat Ce = armaSparseToEigenGlobal(C);
  EigenLDLT Cfactor;
  Cfactor.compute(Ce);
  if(Cfactor.info() != Eigen::Success) {
    Rcpp::stop("Sparse Henderson factorization failed during reduced-rank AI evaluation.");
  }
  auto solveFactor = [](const EigenLDLT & factor, const arma::mat & rhs,
                        const std::string & context) -> arma::mat {
    Eigen::Map<const Eigen::MatrixXd> mapped(rhs.memptr(),
      static_cast<Eigen::Index>(rhs.n_rows), static_cast<Eigen::Index>(rhs.n_cols));
    Eigen::MatrixXd solved = factor.solve(mapped);
    if(factor.info() != Eigen::Success || !solved.allFinite()) {
      Rcpp::stop("Henderson solve failed in " + context + ".");
    }
    return arma::mat(solved.data(), static_cast<arma::uword>(solved.rows()),
                     static_cast<arma::uword>(solved.cols()), true, true);
  };
  EigenSpMat HEigen, HtEigen;
  Eigen::SparseLU<EigenSpMat, Eigen::COLAMDOrdering<int>> Hfactor, HtFactor;
  bool diagonalTransform = true;
  arma::vec inverseTransformDiagonal;
  if(useH) {
    HEigen = armaSparseToEigenGlobal(residualTransform);
    for(arma::sp_mat::const_iterator entry = residualTransform.begin();
        entry != residualTransform.end(); ++entry) {
      if(entry.row() != entry.col() && *entry != 0.0) {
        diagonalTransform = false;
        break;
      }
    }
    if(diagonalTransform) {
      arma::vec diagonal(residualTransform.n_rows);
      for(arma::uword row = 0; row < residualTransform.n_rows; ++row) {
        diagonal(row) = residualTransform(row, row);
      }
      if(arma::any(diagonal == 0.0) || !diagonal.is_finite()) {
        Rcpp::stop("Residual weight transform has a zero or non-finite diagonal.");
      }
      inverseTransformDiagonal = 1.0 / diagonal;
    } else {
      HtEigen = HEigen.transpose();
      Hfactor.analyzePattern(HEigen);
      Hfactor.factorize(HEigen);
      HtFactor.analyzePattern(HtEigen);
      HtFactor.factorize(HtEigen);
      if(Hfactor.info() != Eigen::Success || HtFactor.info() != Eigen::Success) {
        Rcpp::stop("Sparse residual transform factorization failed in reduced-rank AI evaluation.");
      }
    }
  }
  auto solveTransform = [&](const arma::mat & rhs, const bool transpose,
                            const std::string & context) -> arma::mat {
    if(diagonalTransform) {
      arma::mat solved = rhs;
      solved.each_col() %= inverseTransformDiagonal;
      return solved;
    }
    const auto & factor = transpose ? HtFactor : Hfactor;
    Eigen::Map<const Eigen::MatrixXd> mapped(rhs.memptr(),
      static_cast<Eigen::Index>(rhs.n_rows), static_cast<Eigen::Index>(rhs.n_cols));
    Eigen::MatrixXd solved = factor.solve(mapped);
    if(factor.info() != Eigen::Success || !solved.allFinite()) {
      Rcpp::stop("Sparse residual transform solve failed in " + context + ".");
    }
    return arma::mat(solved.data(), static_cast<arma::uword>(solved.rows()),
                     static_cast<arma::uword>(solved.cols()), true, true);
  };

  CovarianceFactorEngine engine;
  const int terms = descriptors.size();
  std::vector<arma::mat> covariance(static_cast<std::size_t>(terms));
  std::vector<arma::vec> standardizedParameters(static_cast<std::size_t>(terms));
  for(int term = 0; term < terms; ++term) {
    standardizedParameters[static_cast<std::size_t>(term)] =
      parameters.subvec(offsets[term], offsets[term + 1] - 1);
    standardizedParameters[static_cast<std::size_t>(term)](0) -= std::log(covarianceScale);
    covariance[static_cast<std::size_t>(term)] = engine.evaluateDescriptor(
      Rcpp::as<Rcpp::List>(descriptors[term]),
      standardizedParameters[static_cast<std::size_t>(term)]);
  }

  struct ResidualBlockState {
    arma::uvec rows;
    arma::uvec index;
    arma::mat covariance;
    arma::mat precision;
  };
  std::vector<ResidualBlockState> residualStates;
  const bool diagonalResidual = arma::accu(arma::abs(covariance.back() -
    arma::diagmat(covariance.back().diag()))) == 0.0;
  arma::vec residualDiagonal(n), residualPrecisionDiagonal(n);
  const arma::uword residualBlocks = residualBlock.max();
  residualStates.reserve(residualBlocks);
  for(arma::uword block = 1; block <= residualBlocks; ++block) {
    ResidualBlockState state;
    state.rows = arma::find(residualBlock == block);
    if(state.rows.is_empty()) { continue; }
    state.index.set_size(state.rows.n_elem);
    for(arma::uword row = 0; row < state.rows.n_elem; ++row) {
      const arma::uword observation = state.rows(row);
      if(residualIndex(observation) < 1 ||
         residualIndex(observation) > covariance.back().n_rows) {
        Rcpp::stop("Residual covariance index is outside its descriptor dimensions.");
      }
      state.index(row) = residualIndex(observation) - 1;
    }
    state.covariance = covariance.back().submat(state.index, state.index);
    if(!arma::inv_sympd(state.precision, state.covariance) ||
       !state.precision.is_finite()) {
      Rcpp::stop("Residual covariance block is not positive definite in reduced-rank AI evaluation.");
    }
    if(diagonalResidual) {
      residualDiagonal.elem(state.rows) = state.covariance.diag();
      residualPrecisionDiagonal.elem(state.rows) = state.precision.diag();
    }
    residualStates.push_back(std::move(state));
  }
  auto applyResidual = [&](const arma::mat & values, const bool precision) -> arma::mat {
    arma::mat input = values;
    if(useH) {
      if(precision) {
        input = arma::mat(residualTransform.t() * input);
      } else {
        input = solveTransform(input, false, "residual covariance derivative");
      }
    }
    arma::mat output(n, values.n_cols, arma::fill::zeros);
    if(diagonalResidual) {
      output = input;
      output.each_col() %= precision ? residualPrecisionDiagonal : residualDiagonal;
    } else {
      for(const ResidualBlockState & state : residualStates) {
        const arma::mat & local = precision ? state.precision : state.covariance;
        output.rows(state.rows) = local * input.rows(state.rows);
      }
    }
    if(useH) {
      output = precision ? arma::mat(residualTransform * output)
                         : solveTransform(output, true, "residual covariance derivative");
    }
    return output;
  };
  auto applyRi = [&](const arma::mat & value) -> arma::mat {
    return applyResidual(value, true);
  };

  // Relationship precisions are factored once per random term and reused for
  // every analytic covariance derivative.
  struct RelationshipFactor {
    std::vector<arma::sp_mat> designs;
    std::shared_ptr<EigenLDLT> factor;
  };
  std::vector<RelationshipFactor> relationship(static_cast<std::size_t>(terms - 1));
  for(int term = 0; term < terms - 1; ++term) {
    const arma::uvec indices = arma::find(Zind == term + 1);
    RelationshipFactor & state = relationship[static_cast<std::size_t>(term)];
    if(indices.is_empty() || indices.n_elem != covariance[static_cast<std::size_t>(term)].n_rows) {
      Rcpp::stop("Reduced-rank covariance dimension does not match its random design blocks.");
    }
    for(const arma::uword index : indices) {
      state.designs.push_back(Rcpp::as<arma::sp_mat>(originalZ[index]));
    }
    arma::sp_mat precision = Rcpp::as<arma::sp_mat>(Ai[term]);
    EigenSpMat precisionEigen = armaSparseToEigenGlobal(precision);
    state.factor = std::make_shared<EigenLDLT>();
    state.factor->compute(precisionEigen);
    if(state.factor->info() != Eigen::Success) {
      Rcpp::stop("Relationship precision factorization failed during reduced-rank AI evaluation.");
    }
  }
  auto applyRandom = [&](const int term, const arma::mat & dSigma,
                         const arma::mat & value) -> arma::mat {
    if(dSigma.is_zero()) { return arma::mat(n, value.n_cols, arma::fill::zeros); }
    const RelationshipFactor & state = relationship[static_cast<std::size_t>(term)];
    const arma::uword q = state.designs.size();
    arma::mat out(n, value.n_cols, arma::fill::zeros);
    std::vector<arma::mat> relationshipApplied(q);
    for(arma::uword trait = 0; trait < q; ++trait) {
      const arma::mat rhs = arma::mat(state.designs[trait].t() * value);
      relationshipApplied[trait] = solveFactor(*state.factor, rhs,
        "relationship covariance action");
    }
    for(arma::uword first = 0; first < q; ++first) {
      arma::mat combined(state.designs[first].n_cols, value.n_cols, arma::fill::zeros);
      for(arma::uword second = 0; second < q; ++second) {
        if(dSigma(first, second) != 0.0) {
          combined += dSigma(first, second) * relationshipApplied[second];
        }
      }
      out += state.designs[first] * combined;
    }
    return out;
  };

  std::vector<std::function<arma::mat(const arma::mat &)>> derivatives;
  std::vector<arma::mat> derivativeCovariances;
  std::vector<bool> zeroDerivatives;
  derivatives.reserve(nParameters);
  zeroDerivatives.reserve(nParameters);
  for(int term = 0; term < terms; ++term) {
    const arma::uword first = offsets[term];
    const arma::uword count = offsets[term + 1] - first;
    const Rcpp::List descriptor = descriptors[term];
    const arma::vec par = standardizedParameters[static_cast<std::size_t>(term)];
    for(arma::uword local = 0; local < count; ++local) {
      const arma::mat dSigma = rrDescriptorDerivative(engine, descriptor, par, local);
      derivativeCovariances.push_back(dSigma);
      zeroDerivatives.push_back(dSigma.is_zero());
      if(term < terms - 1) {
        derivatives.emplace_back([&, term, dSigma](const arma::mat & v) {
          return applyRandom(term, dSigma, v);
        });
      } else {
        const arma::vec diagonal = dSigma.diag();
        const bool diagonalDerivative =
          arma::accu(arma::abs(dSigma - arma::diagmat(diagonal))) == 0.0;
        derivatives.emplace_back([&, dSigma, diagonal, diagonalDerivative](const arma::mat & v) {
          if(dSigma.is_zero()) { return arma::mat(n, v.n_cols, arma::fill::zeros); }
          arma::mat input = v;
          if(useH) { input = solveTransform(input, false, "residual derivative"); }
          arma::mat output(n, v.n_cols, arma::fill::zeros);
          if(diagonalDerivative) {
            output = input;
            output.each_col() %= diagonal.elem(residualIndex - 1);
          } else {
            for(const ResidualBlockState & state : residualStates) {
              output.rows(state.rows) = dSigma.submat(state.index, state.index) *
                input.rows(state.rows);
            }
          }
          if(useH) { output = solveTransform(output, true, "residual derivative"); }
          return output;
        });
      }
    }
  }

  auto applyP = [&](const arma::mat & value) -> arma::mat {
    const arma::mat riValue = applyRi(value);
    const arma::mat rhs = W.t() * riValue;
    const arma::mat solved = solveFactor(Cfactor, rhs, "P-operator coefficient solve");
    return riValue - applyRi(W * solved);
  };
  EigenLDLT randomFactor;
  arma::sp_mat randomW;
  if(!reml && X.n_cols < W.n_cols) {
    randomW = W.cols(X.n_cols, W.n_cols - 1);
    const arma::sp_mat randomC =
      C.submat(X.n_cols, X.n_cols, C.n_rows - 1, C.n_cols - 1);
    EigenSpMat randomEigen = armaSparseToEigenGlobal(randomC);
    randomFactor.compute(randomEigen);
    if(randomFactor.info() != Eigen::Success) {
      Rcpp::stop("Random-effect Henderson factorization failed during ML trace evaluation.");
    }
  }
  auto applyVi = [&](const arma::mat & value) -> arma::mat {
    if(reml) { return applyP(value); }
    if(randomW.n_cols == 0) { return applyRi(value); }
    const arma::mat riValue = applyRi(value);
    const arma::mat solved = solveFactor(randomFactor, randomW.t() * riValue,
      "ML random-effect coefficient solve");
    return riValue - applyRi(randomW * solved);
  };

  const arma::mat yMat(y);
  const arma::mat py = applyP(yMat);
  std::vector<arma::mat> applied(nParameters), projected(nParameters);
  RRAIResult out;
  out.score.zeros(nParameters);
  out.information.zeros(nParameters, nParameters);
  const arma::uword panelWidth = 32;
  arma::vec traces(nParameters, arma::fill::zeros);
  for(int term = 0; term < terms - 1; ++term) {
    const RelationshipFactor & state = relationship[static_cast<std::size_t>(term)];
    const arma::uword q = state.designs.size();
    const arma::uword width = state.designs.front().n_cols;
    const Rcpp::List descriptor = descriptors[term];
    arma::mat span = reducedRankCoefficientMap(descriptor,
      standardizedParameters[static_cast<std::size_t>(term)]);
    if(span.n_cols == q) {
      span.eye(q, q);
    } else {
      arma::mat orthogonal, triangular;
      if(!arma::qr_econ(orthogonal, triangular, span) || !orthogonal.is_finite()) {
        Rcpp::stop("Reduced-rank loading span factorization failed during trace evaluation.");
      }
      span = std::move(orthogonal);
    }
    std::vector<arma::sp_mat> traceDesigns;
    for(arma::uword column = 0; column < span.n_cols; ++column) {
      arma::sp_mat design(n, width);
      for(arma::uword row = 0; row < q; ++row) {
        if(span(row, column) != 0.0) { design += span(row, column) * state.designs[row]; }
      }
      traceDesigns.push_back(std::move(design));
    }
    arma::mat contractions(q, span.n_cols, arma::fill::zeros);
    out.randomTraceDirections += span.n_cols;
    // dG = A Q' + Q A' for an orthonormal loading span Q; only Q needs P solves.
    for(arma::uword first = 0; first < width; first += panelWidth) {
      Rcpp::checkUserInterrupt();
      const arma::uword count = std::min(panelWidth, width - first);
      arma::mat basis(width, count, arma::fill::zeros);
      for(arma::uword column = 0; column < count; ++column) {
        basis(first + column, column) = 1.0;
      }
      const arma::mat relationshipPanel = solveFactor(*state.factor, basis,
        "relationship trace panel");
      for(arma::uword source = 0; source < span.n_cols; ++source) {
        const arma::mat projectedPanel = applyVi(traceDesigns[source] * relationshipPanel);
        for(arma::uword target = 0; target < q; ++target) {
          const arma::sp_mat designPanel = state.designs[target].cols(first, first + count - 1);
          for(arma::sp_mat::const_iterator entry = designPanel.begin();
              entry != designPanel.end(); ++entry) {
            contractions(target, source) += *entry * projectedPanel(entry.row(), entry.col());
          }
        }
      }
    }
    for(arma::uword parameter = offsets[term]; parameter < offsets[term + 1]; ++parameter) {
      const arma::mat derivativeSpan = derivativeCovariances[parameter] * span;
      const arma::mat tangent = derivativeSpan - 0.5 * span * (span.t() * derivativeSpan);
      traces(parameter) = 2.0 * arma::accu(tangent % contractions);
    }
  }
  const arma::uword residualDimension = covariance.back().n_rows;
  arma::mat residualContractions(residualDimension, residualDimension, arma::fill::zeros);
  arma::uvec blockOfRow(n), localRow(n);
  for(std::size_t block = 0; block < residualStates.size(); ++block) {
    const ResidualBlockState & state = residualStates[block];
    for(arma::uword row = 0; row < state.rows.n_elem; ++row) {
      blockOfRow(state.rows(row)) = block;
      localRow(state.rows(row)) = row;
    }
  }
  const arma::uword residualFirst = offsets[terms - 1];
  const bool residualHasShape = nParameters - residualFirst > 1;
  for(arma::uword first = 0; residualHasShape && first < n; first += panelWidth) {
    Rcpp::checkUserInterrupt();
    const arma::uword count = std::min(panelWidth, n - first);
    arma::mat basis(n, count, arma::fill::zeros);
    for(arma::uword column = 0; column < count; ++column) {
      basis(first + column, column) = 1.0;
    }
    if(useH) { basis = solveTransform(basis, true, "residual trace basis"); }
    arma::mat projectedPanel = applyVi(basis);
    if(useH) {
      projectedPanel = solveTransform(projectedPanel, false, "residual trace contraction");
    }
    for(arma::uword column = 0; column < count; ++column) {
      const arma::uword observation = first + column;
      const ResidualBlockState & state = residualStates[blockOfRow(observation)];
      const arma::uword index = state.index(localRow(observation));
      for(arma::uword row = 0; row < state.rows.n_elem; ++row) {
        residualContractions(state.index(row), index) += projectedPanel(state.rows(row), column);
      }
    }
  }
  for(arma::uword parameter = residualFirst + 1; parameter < nParameters; ++parameter) {
    traces(parameter) = arma::accu(derivativeCovariances[parameter] % residualContractions);
  }
  // tr(P V) = n - p (REML), or tr(V^{-1} V) = n (ML).
  traces(residualFirst) = static_cast<double>(n) - (reml ? X.n_cols : 0);
  for(int term = 0; term < terms - 1; ++term) {
    traces(residualFirst) -= traces(offsets[term]);
  }
  for(arma::uword parameter = 0; parameter < nParameters; ++parameter) {
    if(zeroDerivatives[parameter]) {
      applied[parameter].zeros(n, 1);
      projected[parameter].zeros(n, 1);
      continue;
    }
    applied[parameter] = derivatives[parameter](py);
    projected[parameter] = applyP(applied[parameter]);
    const double quadratic = arma::dot(py, applied[parameter]);

    out.score(parameter) = 0.5 * (quadratic - traces(parameter));
  }
  for(arma::uword first = 0; first < nParameters; ++first) {
    for(arma::uword second = first; second < nParameters; ++second) {
      const double value = 0.5 * arma::dot(applied[first], projected[second]);
      out.information(first, second) = value;
      out.information(second, first) = value;
    }
  }
  if(!out.score.is_finite() || !out.information.is_finite()) {
    Rcpp::stop("Reduced-rank analytic score or average information is non-finite.");
  }
  return out;
}

inline bool hasZeroSpecific(const Rcpp::List & descriptors) {
  bool found = false;
  for(int term = 0; term < descriptors.size(); ++term) {
    const Rcpp::List descriptor = descriptors[term];
    const Rcpp::List factors = descriptor["factors"];
    for(int axis = 0; axis < factors.size(); ++axis) {
      if(zeroSpecificFactor(Rcpp::as<Rcpp::List>(factors[axis]))) {
        if(term == descriptors.size() - 1) {
          Rcpp::stop("Zero-specific reduced rank is supported for random effects only; residual covariance must be positive definite.");
        }
        found = true;
      }
    }
  }
  return found;
}

// Conditional fits use ordinary nonsingular Henderson equations for latent
// scores. Only the small loading maps change between likelihood evaluations.
inline Rcpp::List fitReducedRank(
    const arma::sp_mat & X, const Rcpp::List & Z,
    const arma::vec & Zind, const Rcpp::List & descriptors,
    const Rcpp::List & Ai, const arma::sp_mat & y0,
    const arma::sp_mat & H, const bool useH,
    const arma::uvec & residualBlock, const arma::uvec & residualIndex,
    const int iterations, const double tolerance, const double gradientTolerance,
    const int inverseMode, const bool verbose,
    const arma::vec & weightEmInf, const arma::vec & weightInf,
    const bool reml, const bool responsePrepared,
    const double preparedMean, const double preparedSd,
    const bool preparedIntercept,
    const std::function<Rcpp::List(const Rcpp::List &, const arma::vec &,
                                 const Rcpp::List &, int)> & conditionalFit) {
  if(!Rf_isNull(descriptors.attr("vcmap"))) {
    Rcpp::stop("Latent reduced-rank Henderson fitting does not yet support vcc constraints.");
  }
  if(descriptors.hasAttribute("acceleration") &&
     Rcpp::as<std::string>(descriptors.attr("acceleration")) != "none") {
    Rcpp::stop("Latent reduced-rank Henderson average-information fitting does not support covariance-update acceleration.");
  }
  if(iterations < 1) { Rcpp::stop("nIters must be positive."); }
  if(inverseMode < 0 || inverseMode > 2) {
    Rcpp::stop("computeCi must be 0, 1, or 2.");
  }
  const int terms = descriptors.size();
  arma::uword parameterCount = 0;
  std::vector<arma::uword> offsets(terms + 1, 0);
  for(int term = 0; term < terms; ++term) {
    const Rcpp::List descriptor = descriptors[term];
    const arma::vec par = Rcpp::as<arma::vec>(descriptor["par"]);
    if(par.is_empty()) { Rcpp::stop("Covariance descriptor must contain log_sigma2."); }
    parameterCount += par.n_elem;
    offsets[term + 1] = parameterCount;
  }
  const arma::vec response = arma::vec(arma::mat(y0).col(0));
  const double responseVariance = responsePrepared
    ? preparedSd * preparedSd : arma::var(response);
  if(!std::isfinite(responseVariance) || responseVariance <= 0.0) {
    Rcpp::stop("Response variance must be positive and finite in reduced-rank fitting.");
  }
  arma::vec parameters(parameterCount), freeFlags(parameterCount, arma::fill::zeros);
  for(int term = 0; term < terms; ++term) {
    const Rcpp::List descriptor = descriptors[term];
    const arma::vec par = Rcpp::as<arma::vec>(descriptor["par"]);
    const Rcpp::LogicalVector free = descriptor["free"];
    if(free.size() != static_cast<int>(par.n_elem)) {
      Rcpp::stop("covStruct$free and covStruct$par have inconsistent lengths.");
    }
    parameters.subvec(offsets[term], offsets[term + 1] - 1) = par;
    for(arma::uword local = 0; local < par.n_elem; ++local) {
      freeFlags(offsets[term] + local) = free[local] == TRUE;
    }
  }
  const arma::uvec free = arma::find(freeFlags > 0.5);
  CovarianceFactorEngine engine;
  int evaluations = 0;
  std::vector<arma::mat> maps;
  std::vector<arma::uword> widths;
  auto evaluate = [&](const arma::vec & values, const int ci) -> Rcpp::List {
    ++evaluations;
    Rcpp::List covariance = Rcpp::clone(descriptors);
    covariance.attr("rrConditional") = true;
    std::vector<arma::sp_mat> designs;
    std::vector<double> designTerms;
    maps.clear();
    widths.clear();
    for(int term = 0; term < terms; ++term) {
      Rcpp::List descriptor = covariance[term];
      const arma::vec par = values.subvec(offsets[term], offsets[term + 1] - 1);
      descriptor["par"] = par;
      descriptor["free"] = Rcpp::LogicalVector(par.n_elem, false);
      if(term == terms - 1) {
        covariance[term] = descriptor;
        continue;
      }
      Rcpp::List factors = Rcpp::clone(Rcpp::as<Rcpp::List>(descriptor["factors"]));
      arma::mat mapping = reducedRankCoefficientMap(descriptor, par);
      arma::uword reducedDimension = 1;
      for(int axis = 0; axis < factors.size(); ++axis) {
        Rcpp::List factor = factors[axis];
        const int dimension = Rcpp::as<int>(factor["dim"]);
        arma::mat axisMap = arma::eye<arma::mat>(dimension, dimension);
        if(zeroSpecificFactor(factor)) {
          const int first = Rcpp::as<int>(factor["par_start"]);
          const int last = Rcpp::as<int>(factor["par_end"]);
          axisMap = normalizedFactorLoadings(factor, par.subvec(first - 1, last - 1));
          const int rank = axisMap.n_cols;
          factor["dim"] = rank;
          factor["zero_specific"] = false;
          factor["representation"] = Rcpp::List::create(Rcpp::Named("kind") = "diagonal");
          factor["evaluator"] = Rcpp::List::create(
            Rcpp::Named("backend") = "fixed",
            Rcpp::Named("matrix") = arma::eye<arma::mat>(rank, rank));
          factor["derivative"] = Rcpp::List::create(Rcpp::Named("backend") = "numeric");
          factor["structurally_diagonal"] = true;
          factors[axis] = factor;
        }
        reducedDimension *= axisMap.n_cols;
      }
      descriptor["factors"] = factors;
      descriptor["dim"] = static_cast<int>(reducedDimension);
      covariance[term] = descriptor;
      const arma::uvec indices = arma::find(Zind == term + 1);
      if(indices.n_elem != mapping.n_rows || indices.is_empty()) {
        Rcpp::stop("Reduced-rank covariance dimensions do not match the random design blocks.");
      }
      const arma::sp_mat first = Rcpp::as<arma::sp_mat>(Z[indices(0)]);
      const arma::uword width = first.n_cols;
      for(arma::uword column = 0; column < mapping.n_cols; ++column) {
        arma::sp_mat design(first.n_rows, width);
        for(arma::uword row = 0; row < mapping.n_rows; ++row) {
          if(mapping(row, column) != 0.0) {
            const arma::sp_mat original = Rcpp::as<arma::sp_mat>(Z[indices(row)]);
            if(original.n_cols != width) {
              Rcpp::stop("Reduced-rank random design blocks must have the same relationship dimension.");
            }
            design += mapping(row, column) * original;
          }
        }
        designs.push_back(std::move(design));
        designTerms.push_back(term + 1);
      }
      maps.push_back(std::move(mapping));
      widths.push_back(width);
    }
    Rcpp::List latentDesign(designs.size());
    for(std::size_t block = 0; block < designs.size(); ++block) {
      latentDesign[block] = designs[block];
    }
    return conditionalFit(latentDesign, arma::vec(designTerms), covariance, ci);
  };
  auto unpack = [&](const arma::vec & reduced) -> arma::vec {
    arma::vec result = parameters;
    result.elem(free) = reduced;
    return result;
  };
  bool rejectInvalidTrials = true;
  int invalidTrialEvaluations = 0;
  std::string firstInvalidTrialReason;
  std::map<std::string, int> invalidTrialReasons;
  auto rejectTrial = [&](const std::string & reason) -> double {
    if(!rejectInvalidTrials) { Rcpp::stop(reason); }
    if(invalidTrialEvaluations == 0) { firstInvalidTrialReason = reason; }
    ++invalidTrialEvaluations;
    ++invalidTrialReasons[reason];
    return 1.0e100;
  };
  auto objectiveEvaluate = [&](const arma::vec & reduced) -> double {
    Rcpp::checkUserInterrupt();
    try {
      const arma::vec values = unpack(reduced);
      if(!values.is_finite()) {
        return rejectTrial("Reduced-rank trial covariance parameters are not finite.");
      }
      for(int term = 0; term < terms; ++term) {
        if(values(offsets[term]) < -700.0 || values(offsets[term]) > 700.0) {
          return rejectTrial("Reduced-rank trial log variance is outside the finite exponential range.");
        }
      }
      const Rcpp::List result = evaluate(values, 0);
      const arma::vec likelihood = Rcpp::as<arma::vec>(result["llik"]);
      if(likelihood.is_empty() || !likelihood.is_finite()) {
        return rejectTrial("The latent reduced-rank Henderson likelihood is empty or non-finite.");
      }
      return -likelihood(likelihood.n_elem - 1);
    } catch(const Rcpp::internal::InterruptedException &) {
      throw;
    } catch(const Rcpp::LongjumpException &) {
      throw;
    } catch(const Rcpp::exception & error) {
      if(!rejectInvalidTrials) { throw; }
      return rejectTrial(error.what());
    } catch(const std::runtime_error & error) {
      if(!rejectInvalidTrials) { throw; }
      return rejectTrial(error.what());
    }
  };
  arma::vec reduced = parameters.elem(free);
  rejectInvalidTrials = true;
  Rcpp::List result;
  RRAIResult rrInformation;
  arma::vec likelihoodHistory;
  arma::mat parameterHistory(parameterCount, 0);
  arma::mat normHistory(3, 0);
  arma::mat changeRatioHistory(parameterCount, 0);
  arma::mat progressHistory(4, 0);
  bool converged = free.is_empty();
  bool lineSearchStalled = false;
  double previousLikelihood = -std::numeric_limits<double>::infinity();
  double acceptedTrustScale = 1.0, acceptedLineAlpha = 1.0;
  const std::clock_t progressStart = std::clock();
  const std::vector<int> progressWidths = {9, 16, 8, 14, 10, 10, 10, 17};
  auto printProgressRow = [&](const std::vector<std::string> & values) {
    std::ostringstream row;
    for(std::size_t column = 0; column < values.size(); ++column) {
      if(column > 0) { row << "  "; }
      row << std::setw(progressWidths[column]) << values[column];
    }
    Rcpp::Rcout << row.str() << arma::endl;
  };
  auto progressNumber = [](const double value, const int precision,
                           const bool scientific = false) -> std::string {
    std::ostringstream text;
    text << (scientific ? std::scientific : std::fixed)
         << std::setprecision(precision) << value;
    return text.str();
  };
  if(verbose) {
    Rcpp::Rcout << "MME path: exact latent reduced rank + shared analytic AI traces"
                << arma::endl;
  }
  for(int iteration = 0; iteration < iterations; ++iteration) {
    Rcpp::checkUserInterrupt();
    const arma::vec current = unpack(reduced);
    result = evaluate(current, 0);
    const arma::vec likelihood = Rcpp::as<arma::vec>(result["llik"]);
    if(likelihood.is_empty() || !likelihood.is_finite()) {
      Rcpp::stop("Reduced-rank Henderson fit returned a non-finite likelihood.");
    }
    if(verbose && iteration == 0) {
      const Rcpp::List diagnostics = result["engineDiagnostics"];
      if(Rcpp::as<bool>(diagnostics["blockSchurActive"])) {
        const std::string engine = Rcpp::as<bool>(diagnostics["factorSchurActive"])
          ? "Latent-factor Schur" : Rcpp::as<bool>(diagnostics["blockChainActive"])
          ? "Bordered block-chain" : Rcpp::as<bool>(diagnostics["blockBandedActive"])
          ? "Block-banded Schur" : "Dense block-Schur";
        Rcpp::Rcout << engine << " engine active ("
          << Rcpp::as<int>(diagnostics["blockSchurGroups"]) << " latent random-effect blocks, "
          << Rcpp::as<int>(diagnostics["blockSchurBorder"]) << " border effects)."
          << arma::endl;
      } else {
        Rcpp::Rcout << "Latent coefficient solver: "
          << Rcpp::as<std::string>(result["solver"]) << "." << arma::endl;
      }
      printProgressRow({"iteration", "LogLik", "wall", "cpu(sec)", "EMweight",
                        "trustAlpha", "lineAlpha", "scoreNorm"});
    }
    const double currentLikelihood = likelihood(likelihood.n_elem - 1);
    rrInformation = reducedRankInformation(result, X, Z, Zind, Ai,
      residualBlock, residualIndex, H, useH, descriptors, current, offsets,
      reml);
    likelihoodHistory.insert_rows(likelihoodHistory.n_elem,
      arma::vec(1, arma::fill::value(currentLikelihood)));
    parameterHistory.insert_cols(parameterHistory.n_cols, current);
    arma::vec parameterChange(parameterCount, arma::fill::zeros);
    arma::vec previousChange(parameterCount, arma::fill::zeros);
    if(parameterHistory.n_cols > 1) {
      parameterChange = parameterHistory.col(parameterHistory.n_cols - 1) -
        parameterHistory.col(parameterHistory.n_cols - 2);
    }
    if(parameterHistory.n_cols > 2) {
      previousChange = parameterHistory.col(parameterHistory.n_cols - 2) -
        parameterHistory.col(parameterHistory.n_cols - 3);
    }
    arma::vec changeRatio(parameterCount, arma::fill::zeros);
    if(parameterHistory.n_cols > 2) {
      for(arma::uword parameter = 0; parameter < parameterCount; ++parameter) {
        changeRatio(parameter) = std::abs(previousChange(parameter)) >
          std::numeric_limits<double>::epsilon()
          ? parameterChange(parameter) / previousChange(parameter) : NA_REAL;
      }
    }
    changeRatioHistory.insert_cols(changeRatioHistory.n_cols, changeRatio);
    arma::vec norms(3, arma::fill::zeros);
    norms(0) = arma::norm(parameterChange.elem(free), 1);
    norms(1) = arma::norm(rrInformation.score.elem(free), 1);
    norms(2) = arma::norm(parameterChange.elem(free), 2) /
      (arma::norm(current.elem(free), 2) + 1.0e-12);
    normHistory.insert_cols(normHistory.n_cols, norms);
    const double emWeight = weightEmInf.is_empty() ? 0.0
      : weightEmInf(std::min<arma::uword>(iteration, weightEmInf.n_elem - 1));
    const double cpuSeconds = static_cast<double>(std::clock() - progressStart) / CLOCKS_PER_SEC;
    progressHistory.insert_cols(progressHistory.n_cols,
      arma::vec({emWeight, acceptedTrustScale, acceptedLineAlpha, cpuSeconds}));
    if(verbose) {
      const std::time_t now = std::time(nullptr);
      const std::tm * local = std::localtime(&now);
      std::ostringstream wall;
      wall << std::put_time(local, "%H:%M:%S");
      printProgressRow({std::to_string(iteration + 1),
        progressNumber(currentLikelihood, 6), wall.str(), progressNumber(cpuSeconds, 3),
        progressNumber(emWeight, 3), progressNumber(acceptedTrustScale, 3),
        progressNumber(acceptedLineAlpha, 3), progressNumber(norms(1), 6, true)});
    }

    if(free.is_empty()) {
      converged = true;
      break;
    }
    if(iteration == iterations - 1) {
      converged = std::abs(arma::norm(rrInformation.score.elem(free), 1)) <= gradientTolerance ||
        (iteration > 0 && std::abs(currentLikelihood - previousLikelihood) < tolerance);
      break;
    }
    if(arma::norm(rrInformation.score.elem(free), 1) <= gradientTolerance ||
       (iteration > 0 && std::abs(currentLikelihood - previousLikelihood) < tolerance)) {
      converged = true;
      break;
    }

    const double scoreWeight = 1.0;
    if(!std::isfinite(emWeight) || emWeight < 0.0 || emWeight > 1.0 ||
       !std::isfinite(scoreWeight) || scoreWeight <= 0.0) {
      Rcpp::stop("Reduced-rank AI weights must be finite and in their valid ranges.");
    }
    arma::mat freeInformation = rrInformation.information.submat(free, free);
    freeInformation = 0.5 * (freeInformation + freeInformation.t());
    if(emWeight > 0.0) {
      freeInformation = (1.0 - emWeight) * freeInformation +
        emWeight * arma::diagmat(freeInformation.diag());
    }
    arma::vec direction;
    bool solved = arma::solve(direction, freeInformation,
      scoreWeight * rrInformation.score.elem(free), arma::solve_opts::likely_sympd);
    if(!solved) {
      const double ridge = std::max(1.0e-10,
        1.0e-8 * arma::abs(freeInformation.diag()).max());
      freeInformation.diag() += ridge;
      solved = arma::solve(direction, freeInformation,
        scoreWeight * rrInformation.score.elem(free));
    }
    if(!solved || !direction.is_finite()) {
      Rcpp::stop("Unable to solve the reduced-rank average-information system.");
    }
    const double maxStep = arma::abs(direction).max();
    const double trustRadius = 2.5;
    const double trustScale = maxStep > trustRadius ? trustRadius / maxStep : 1.0;
    direction *= trustScale;
    previousLikelihood = currentLikelihood;
    const double initialStep = weightInf.is_empty() ? 1.0
      : weightInf(std::min<arma::uword>(iteration, weightInf.n_elem - 1));
    // weightInf is the stepWeight vector at this native boundary.
    double alpha = initialStep;
    bool accepted = false;
    for(int search = 0; search < 24; ++search) {
      arma::vec candidate = reduced + alpha * direction;
      const double candidateObjective = objectiveEvaluate(candidate);
      if(candidateObjective < 1.0e99 &&
         -candidateObjective >= currentLikelihood - 1.0e-12) {
        reduced = std::move(candidate);
        acceptedTrustScale = trustScale;
        acceptedLineAlpha = alpha;
        accepted = true;
        break;
      }
      alpha *= 0.7;
    }
    if(!accepted) {
      lineSearchStalled = true;
      break;
    }
  }
  const arma::vec fitted = unpack(reduced);
  result = evaluate(fitted, inverseMode == 2 ? 2 : 0);
  const arma::vec latent = Rcpp::as<arma::vec>(result["bu"]);
  arma::uword originalSize = X.n_cols, latentSize = X.n_cols;
  for(std::size_t term = 0; term < maps.size(); ++term) {
    originalSize += maps[term].n_rows * widths[term];
    latentSize += maps[term].n_cols * widths[term];
  }
  arma::sp_mat transform(originalSize, latentSize);
  for(arma::uword fixed = 0; fixed < X.n_cols; ++fixed) { transform(fixed, fixed) = 1.0; }
  arma::uword originalOffset = X.n_cols, latentOffset = X.n_cols;
  Rcpp::List partitions(terms - 1), blups(terms - 1), pevs(terms - 1);
  partitions.attr("dim") = Rcpp::IntegerVector::create(terms - 1, 1);
  blups.attr("dim") = Rcpp::IntegerVector::create(terms - 1, 1);
  pevs.attr("dim") = Rcpp::IntegerVector::create(terms - 1, 1);
  for(std::size_t term = 0; term < maps.size(); ++term) {
    const arma::uword width = widths[term];
    arma::mat partition(maps[term].n_rows, 2);
    for(arma::uword row = 0; row < maps[term].n_rows; ++row) {
      partition(row, 0) = originalOffset + row * width + 1;
      partition(row, 1) = originalOffset + (row + 1) * width;
      for(arma::uword column = 0; column < maps[term].n_cols; ++column) {
        if(maps[term](row, column) == 0.0) { continue; }
        for(arma::uword id = 0; id < width; ++id) {
          transform(originalOffset + row * width + id,
                    latentOffset + column * width + id) = maps[term](row, column);
        }
      }
    }
    partitions[term] = partition;
    originalOffset += maps[term].n_rows * width;
    latentOffset += maps[term].n_cols * width;
  }
  const arma::vec coefficients = transform * latent;
  arma::sp_mat covariance;
  if(inverseMode == 2) {
    const arma::sp_mat latentCovariance = Rcpp::as<arma::sp_mat>(result["Ci"]);
    covariance = transform * latentCovariance * transform.t();
  } else if(inverseMode == 1) {
    const arma::sp_mat C = Rcpp::as<arma::sp_mat>(result["C"]);
    Eigen::SparseMatrix<double> Ce = armaSparseToEigenGlobal(C);
    SommerSpatialLDLT factor;
    factor.compute(Ce);
    if(factor.info() != Eigen::Success) {
      Rcpp::stop("Latent Henderson factorization failed during prediction-error computation.");
    }
    const double scale = Rcpp::as<double>(result["Cscale"]);
    covariance = arma::sp_mat(originalSize, originalSize);
    for(arma::uword first = 0; first < originalSize; first += 32) {
      const arma::uword count = std::min<arma::uword>(32, originalSize - first);
      const arma::mat panel = arma::mat(transform.rows(first, first + count - 1).t());
      const Eigen::Map<const Eigen::MatrixXd> rhs(panel.memptr(), panel.n_rows, panel.n_cols);
      const Eigen::MatrixXd solved = factor.solve(rhs);
      if(factor.info() != Eigen::Success) { Rcpp::stop("Latent prediction-error solve failed."); }
      for(arma::uword column = 0; column < count; ++column) {
        covariance(first + column, first + column) =
          scale * rhs.col(column).dot(solved.col(column));
      }
    }
  }
  originalOffset = X.n_cols;
  for(std::size_t term = 0; term < maps.size(); ++term) {
    arma::mat values(widths[term], maps[term].n_rows);
    arma::mat errors;
    if(inverseMode > 0) { errors.set_size(values.n_rows, values.n_cols); }
    for(arma::uword column = 0; column < values.n_cols; ++column) {
      values.col(column) = coefficients.subvec(originalOffset, originalOffset + values.n_rows - 1);
      if(inverseMode > 0) {
        for(arma::uword row = 0; row < values.n_rows; ++row) {
          errors(row, column) = covariance(originalOffset + row, originalOffset + row);
        }
      }
      originalOffset += values.n_rows;
    }
    blups[term] = values;
    pevs[term] = errors;
  }
  Rcpp::List theta = Rcpp::clone(Rcpp::as<Rcpp::List>(result["theta"]));
  for(int term = 0; term < terms - 1; ++term) {
    const Rcpp::List descriptor = descriptors[term];
    const arma::vec par = fitted.subvec(offsets[term], offsets[term + 1] - 1);
    theta[term] = engine.evaluateDescriptor(descriptor, par);
  }
  result["theta"] = theta;
  result["latentBu"] = latent;
  result["latentW"] = result["W"];
  arma::sp_mat originalDesign(X.n_rows, originalSize);
  if(X.n_cols > 0){ originalDesign.cols(0, X.n_cols - 1) = X; }
  originalOffset = X.n_cols;
  for(int term = 0; term < terms - 1; ++term){
    const arma::uvec indices = arma::find(Zind == term + 1);
    for(arma::uword block = 0; block < indices.n_elem; ++block){
      const arma::sp_mat design = Rcpp::as<arma::sp_mat>(Z[indices(block)]);
      originalDesign.cols(originalOffset, originalOffset + design.n_cols - 1) = design;
      originalOffset += design.n_cols;
    }
  }
  result["W"] = originalDesign;
  result["coefficientMap"] = transform;
  result["bu"] = coefficients;
  result["u"] = arma::vec(coefficients.tail(originalSize - X.n_cols));
  result["uList"] = blups;
  result["uPevList"] = pevs;
  result["partitions"] = partitions;
  result["Ci"] = inverseMode == 0 ? arma::sp_mat() : covariance;
  result["CiMode"] = inverseMode;
  result["CiComputed"] = inverseMode == 2;
  result["convergence"] = converged && !lineSearchStalled;
  result[".ldltCache"] = R_NilValue;
  result[".cholmodCache"] = R_NilValue;
  Rcpp::List diagnostics = result["engineDiagnostics"];
  diagnostics["latentReducedRank"] = true;
  diagnostics["reducedRankActive"] = true;
  diagnostics["latentCoefficients"] = static_cast<double>(latentSize);
  diagnostics["originalCoefficients"] = static_cast<double>(originalSize);
  diagnostics["likelihoodEvaluations"] = evaluations;
  diagnostics["optimizer"] = "AI-Henderson-analytic";
  diagnostics["traceMethod"] = "shared-covariance-contractions";
  diagnostics["randomTraceDirections"] = static_cast<double>(rrInformation.randomTraceDirections);
  result["engineDiagnostics"] = diagnostics;

  arma::mat information = rrInformation.information;
  arma::mat uncertainty(parameterCount, parameterCount, arma::fill::zeros);
  const arma::vec score = rrInformation.score;
  if(!free.is_empty()) {
    arma::mat aiFree = information.submat(free, free);
    aiFree = 0.5 * (aiFree + aiFree.t());
    arma::mat inverse = arma::pinv(aiFree);
    arma::vec jacobian(parameterCount, arma::fill::ones);
    for(int term = 0; term < terms; ++term) {
      const arma::uword first = offsets[term];
      const arma::uword last = offsets[term + 1] - 1;
      jacobian(first) = std::exp(fitted(first));
      const Rcpp::List descriptor = descriptors[term];
      const Rcpp::List factors = descriptor["factors"];
      for(int axis = 0; axis < factors.size(); ++axis) {
        const Rcpp::List factor = factors[axis];
        const int start = Rcpp::as<int>(factor["par_start"]) - 1;
        const int end = Rcpp::as<int>(factor["par_end"]) - 1;
        if(end < start) { continue; }
        const Rcpp::List report = factor["report"];
        if(Rcpp::as<std::string>(report["backend"]) != "builtin") {
          Rcpp::stop("Reduced-rank standard errors require builtin reporting transforms.");
        }
        const Rcpp::CharacterVector transform = report["transform"];
        const Rcpp::NumericVector lower = report["lower"];
        const Rcpp::NumericVector upper = report["upper"];
        for(int local = 0; local <= end - start; ++local) {
          const arma::uword index = first + static_cast<arma::uword>(start + local);
          const double value = fitted(index);
          const std::string code = Rcpp::as<std::string>(transform[local]);
          if(code == "exp") {
            jacobian(index) = std::exp(value);
          } else if(code == "tanh") {
            const double transformed = std::tanh(value);
            jacobian(index) = 1.0 - transformed * transformed;
          } else if(code == "bounded_logit") {
            const double logistic = value >= 0.0 ? 1.0 / (1.0 + std::exp(-value))
              : std::exp(value) / (1.0 + std::exp(value));
            jacobian(index) = (upper[local] - lower[local]) * logistic * (1.0 - logistic);
          } else if(code != "identity") {
            Rcpp::stop("Unknown reduced-rank reporting transform: " + code);
          }
        }
      }
    }
    arma::mat jacobianFree(parameterCount, free.n_elem, arma::fill::zeros);
    for(arma::uword i = 0; i < free.n_elem; ++i) {
      jacobianFree(free(i), i) = jacobian(free(i));
    }
    uncertainty = jacobianFree * inverse * jacobianFree.t();
  }
  result["InfMat"] = 2.0 * information;
  result["InfMatInv"] = uncertainty;
  result["theta_se"] = uncertainty;
  result["dLu"] = 2.0 * score;
  result["monitor"] = parameterHistory;
  result["llik"] = likelihoodHistory;
  result["normMonitor"] = normHistory;
  result["percDelta"] = changeRatioHistory;
  result["toBoundary"] = arma::mat(1, parameterCount, arma::fill::zeros);
  diagnostics["scoreWorking"] = score;
  diagnostics["averageInformationWorking"] = information;
  diagnostics["likelihoodEvaluations"] = evaluations;
  diagnostics["parameterInformation"] = "average";
  diagnostics["lineSearchStalled"] = lineSearchStalled;
  Rcpp::NumericMatrix progress = Rcpp::wrap(arma::mat(progressHistory.t()));
  Rcpp::colnames(progress) = Rcpp::CharacterVector::create(
    "EMweight", "trustAlpha", "lineAlpha", "cpu(sec)");
  diagnostics["iterationProgress"] = progress;
  diagnostics["invalidTrialEvaluations"] = invalidTrialEvaluations;
  Rcpp::IntegerVector invalidCounts(invalidTrialReasons.size());
  Rcpp::CharacterVector invalidReasons(invalidTrialReasons.size());
  int reasonIndex = 0;
  for(const auto & reason : invalidTrialReasons) {
    invalidReasons[reasonIndex] = reason.first;
    invalidCounts[reasonIndex] = reason.second;
    ++reasonIndex;
  }
  invalidCounts.names() = invalidReasons;
  diagnostics["invalidTrialReasons"] = invalidCounts;
  if(invalidTrialEvaluations > 0) {
    Rcpp::warning(
      "Latent reduced-rank Henderson fitting rejected %d invalid trial evaluations. "
      "First reason: %s. All reasons and counts are recorded in "
      "engineDiagnostics$invalidTrialReasons.",
      invalidTrialEvaluations, firstInvalidTrialReason.c_str());
  }
  result["engineDiagnostics"] = diagnostics;
  return result;
}

}

#endif
