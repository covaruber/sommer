#include "mme_glmm_multinomial.h"
#include "mme_glmm_precision.h"
#include <limits>

namespace {

class MultinomialModeContext {
public:
  arma::vec response;
  arma::sp_mat design;
  arma::sp_mat precision;
  arma::uword n;
  arma::uword d;
  arma::uword q;
  double logDetPrecision = 0.0;
  sommer::CholmodState factor;
  sommer::CholmodState precisionFactor;
  Eigen::SparseMatrix<double> hessian;
  std::vector<std::vector<int>> supports;
  double structuralPairs = 0.0;
  double workspaceBytes = 0.0;
  double memoryMB;
  bool precisionReady = false;
  int symbolicAnalyses = 0;
  int precisionSymbolicAnalyses = 0;
  int factorizations = 0;
  int evaluations = 0;
  int scoreEvaluations = 0;
  int derivativeSolves = 0;

  MultinomialModeContext(const arma::vec & suppliedResponse,
      const arma::sp_mat & suppliedDesign, const arma::sp_mat & suppliedPrecision,
      const int dimension, const double memory)
      : response(suppliedResponse), design(suppliedDesign),
        precision(suppliedPrecision), n(response.n_elem),
        d(dimension > 0 ? dimension : 0), q(design.n_cols), memoryMB(memory){
    if(n == 0 || d == 0 || q == 0 ||
       static_cast<double>(n) * d * d > std::numeric_limits<int>::max() ||
       q > static_cast<arma::uword>(std::numeric_limits<int>::max()) ||
       design.n_rows != n * d || precision.n_rows != q || precision.n_cols != q){
      Rcpp::stop("Multinomial mode context dimensions are incompatible or exceed sparse index limits.");
    }
    if(!response.is_finite() || !design.is_finite() || !precision.is_finite() ||
       !precision.is_symmetric(1e-12)){
      Rcpp::stop("Multinomial mode context requires finite inputs and symmetric precision.");
    }
    if(!std::isfinite(memory) || memory <= 0.0){
      Rcpp::stop("Multinomial mode denseMemoryMB must be positive and finite.");
    }
    for(arma::uword row = 0; row < n; ++row){
      if(response(row) < 1 || response(row) > d + 1 ||
         response(row) != std::floor(response(row))){
        Rcpp::stop("Multinomial response codes must be integers from 1 through the declared category count.");
      }
    }
    checkWorkspace(precision);
    supports.resize(n);
    for(auto cell = design.begin(); cell != design.end(); ++cell){
      supports[cell.row() % n].push_back(static_cast<int>(cell.col()));
    }
    for(auto & support : supports){
      std::sort(support.begin(), support.end());
      support.erase(std::unique(support.begin(), support.end()), support.end());
      structuralPairs += static_cast<double>(support.size()) * support.size();
    }
    updatePrecision(precision);
  }

  void checkWorkspace(const arma::sp_mat & suppliedPrecision){
    const double budget = memoryMB * 1048576.0;
    const double baseBytes = 8.0 * (16.0 * n * d * d + 24.0 * n * d +
      16.0 * q + 8.0 * design.n_nonzero + 8.0 * suppliedPrecision.n_nonzero);
    if(!std::isfinite(baseBytes) || baseBytes > budget){
      Rcpp::stop("Multinomial mode workspace exceeds denseMemoryMB.");
    }
    const double pairs = suppliedPrecision.n_nonzero + structuralPairs;
    if(!std::isfinite(pairs) || pairs > std::numeric_limits<int>::max() ||
       baseBytes + 64.0 * pairs > budget){
      Rcpp::stop("Multinomial structural Hessian workspace exceeds denseMemoryMB.");
    }
    workspaceBytes = baseBytes + 64.0 * pairs;
  }

  void updatePrecision(const arma::sp_mat & suppliedPrecision){
    if(suppliedPrecision.n_rows != q || suppliedPrecision.n_cols != q){
      Rcpp::stop("Multinomial random precision dimensions are incompatible.");
    }
    if(!suppliedPrecision.is_finite() || !suppliedPrecision.is_symmetric(1e-12)){
      Rcpp::stop("Multinomial random precision must be finite and symmetric.");
    }
    if(precisionReady){
      const arma::sp_mat difference = suppliedPrecision - precision;
      if(difference.n_nonzero == 0) return;
    }
    checkWorkspace(suppliedPrecision);
    precisionReady = false;
    const auto eigenPrecision = sommer::laplaceEigenSparse(suppliedPrecision);
    double logDet = 0.0;
    if(!precisionFactor.analyzePattern(eigenPrecision, &precisionSymbolicAnalyses) ||
       !precisionFactor.factorize(eigenPrecision) ||
       !sommer::laplaceFactorLogDeterminant(precisionFactor, logDet)){
      Rcpp::stop("Multinomial random precision must be positive definite.");
    }
    std::vector<Eigen::Triplet<double>> entries;
    entries.reserve(static_cast<std::size_t>(suppliedPrecision.n_nonzero + structuralPairs));
    for(auto cell = suppliedPrecision.begin(); cell != suppliedPrecision.end(); ++cell){
      entries.emplace_back(cell.row(), cell.col(), 1.0);
    }
    for(const auto & support : supports){
      for(int column : support) for(int row : support){
        entries.emplace_back(row, column, 1.0);
      }
    }
    hessian.resize(static_cast<int>(q), static_cast<int>(q));
    hessian.setFromTriplets(entries.begin(), entries.end());
    hessian.makeCompressed();
    precision = suppliedPrecision;
    logDetPrecision = logDet;
    precisionReady = true;
  }

  Rcpp::List work() const {
    return Rcpp::List::create(Rcpp::_["evaluations"]=evaluations,
      Rcpp::_["symbolicAnalyses"]=symbolicAnalyses,
      Rcpp::_["precisionSymbolicAnalyses"]=precisionSymbolicAnalyses,
      Rcpp::_["factorizations"]=factorizations,
      Rcpp::_["scoreEvaluations"]=scoreEvaluations,
      Rcpp::_["derivativeSolves"]=derivativeSolves);
  }

  arma::mat solve(sommer::CholmodState & state, const arma::mat & rhs){
    ++derivativeSolves;
    Eigen::Map<const Eigen::MatrixXd> mapped(rhs.memptr(), rhs.n_rows, rhs.n_cols);
    const Eigen::MatrixXd value = state.solve(mapped, "Multinomial marginal score solve failed.");
    if(!value.allFinite()) Rcpp::stop("Multinomial marginal score solve is non-finite.");
    return arma::mat(value.data(), value.rows(), value.cols());
  }

  arma::vec gradient(const Rcpp::List & modeResult, const arma::sp_mat & fixedDesign,
      const std::vector<sommer::LaplacePrecisionDerivative> & precisionDerivatives){
    if(fixedDesign.n_cols == 0 && precisionDerivatives.empty()){
      ++scoreEvaluations;
      return arma::vec();
    }
    const double columnBytes = 8.0 * (4.0 * q + d);
    double derivativeBytes = 8.0 * (12.0 * n * d + 8.0 * d * d +
      12.0 * q + 4.0 * design.n_nonzero + fixedDesign.n_cols +
      precisionDerivatives.size());
    for(const auto & derivative : precisionDerivatives){
      derivativeBytes += 32.0 * derivative.block.n_nonzero;
    }
    const double available = memoryMB * 1048576.0 - workspaceBytes - derivativeBytes;
    if(!std::isfinite(available) || available < columnBytes){
      Rcpp::stop("Multinomial marginal score workspace exceeds denseMemoryMB.");
    }
    const arma::uword panelSize = static_cast<arma::uword>(
      std::min(32.0, std::floor(available / columnBytes)));
    const arma::vec mode = Rcpp::as<arma::vec>(modeResult["mode"]);
    const arma::mat eta = Rcpp::as<arma::mat>(modeResult["linearPredictors"]);
    arma::mat probability;
    arma::vec score;
    arma::sp_mat curvature;
    evaluate(eta, probability, score, &curvature);
    const arma::sp_mat transposed = design.t();
    arma::vec correction(n * d, arma::fill::zeros);
    for(arma::uword row = 0; row < n; ++row){
      Rcpp::checkUserInterrupt();
      arma::uvec indices(d);
      for(arma::uword a = 0; a < d; ++a) indices(a) = row + n * a;
      const arma::sp_mat block = transposed.cols(indices);
      arma::mat leverage(d, d);
      for(arma::uword first = 0; first < d; first += panelSize){
        const arma::uword last = std::min(d - 1, first + panelSize - 1);
        const arma::mat solved = solve(factor, arma::mat(block.cols(first, last)));
        leverage.cols(first, last) = block.t() * solved;
      }
      const arma::vec contracted = sommer::multinomialCurvatureContraction(
        probability.row(row).t(), leverage);
      for(arma::uword a = 0; a < d; ++a) correction(indices(a)) = contracted(a);
    }
    const arma::vec adjoint = solve(factor, arma::mat(design.t() * correction));
    const arma::vec predictorScore = score - correction + curvature * (design * adjoint);
    arma::vec result(fixedDesign.n_cols + precisionDerivatives.size());
    if(fixedDesign.n_cols) result.head(fixedDesign.n_cols) = fixedDesign.t() * predictorScore;
    for(std::size_t parameter = 0; parameter < precisionDerivatives.size(); ++parameter){
      Rcpp::checkUserInterrupt();
      const auto & derivative = precisionDerivatives[parameter];
      double tracePrecision = std::isfinite(derivative.logDet) ? derivative.logDet : 0.0;
      double traceHessian = 0.0;
      const arma::uword blockEnd = derivative.start + derivative.block.n_rows - 1;
      for(arma::uword first = 0; first < derivative.block.n_cols; first += panelSize){
        const arma::uword last = std::min(derivative.block.n_cols - 1, first + panelSize - 1);
        arma::mat rhs(q, last - first + 1, arma::fill::zeros);
        rhs.rows(derivative.start, blockEnd) = arma::mat(derivative.block.cols(first, last));
        if(arma::accu(arma::abs(rhs)) == 0.0) continue;
        if(!std::isfinite(derivative.logDet)){
          const arma::mat solvedPrecision = solve(precisionFactor, rhs);
          for(arma::uword column = first; column <= last; ++column){
            tracePrecision += solvedPrecision(derivative.start + column, column - first);
          }
        }
        const arma::mat solvedHessian = solve(factor, rhs);
        for(arma::uword column = first; column <= last; ++column){
          traceHessian += solvedHessian(derivative.start + column, column - first);
        }
      }
      const arma::vec blockMode = mode.subvec(derivative.start, blockEnd);
      const arma::vec precisionMode = derivative.block * blockMode;
      result(fixedDesign.n_cols + parameter) = -0.5 * arma::dot(blockMode, precisionMode) +
        0.5 * (tracePrecision - traceHessian) +
        arma::dot(adjoint.subvec(derivative.start, blockEnd), precisionMode);
    }
    if(!result.is_finite()) Rcpp::stop("Multinomial marginal score is non-finite.");
    ++scoreEvaluations;
    return result;
  }

  double evaluate(const arma::mat & eta, arma::mat & probability,
      arma::vec & score, arma::sp_mat * weightedCurvature){
    probability.set_size(n, d + 1);
    score.set_size(n * d);
    arma::umat locations;
    arma::vec values;
    if(weightedCurvature != nullptr){
      locations.set_size(2, n * d * d);
      values.set_size(n * d * d);
    }
    arma::uword position = 0;
    double likelihood = 0.0;
    const arma::vec direction(d, arma::fill::zeros);
    for(arma::uword row = 0; row < n; ++row){
      if(row % 128 == 0) Rcpp::checkUserInterrupt();
      const auto value = sommer::multinomialEvaluation(eta.row(row).t(),
        static_cast<arma::uword>(response(row) - 1), direction,
        weightedCurvature == nullptr ? -1 : 0);
      likelihood += value.logDensity;
      probability.row(row) = value.probability.t();
      for(arma::uword a = 0; a < d; ++a){
        score(row + n * a) = value.score(a);
        if(weightedCurvature == nullptr) continue;
        for(arma::uword b = 0; b < d; ++b){
          locations(0, position) = row + n * a;
          locations(1, position) = row + n * b;
          values(position++) = value.curvature(a, b);
        }
      }
    }
    if(weightedCurvature != nullptr){
      *weightedCurvature = arma::sp_mat(locations, values, n * d, n * d);
    }
    return likelihood;
  }

  bool factorize(const arma::sp_mat & curvature){
    const arma::sp_mat numeric = precision + design.t() * curvature * design;
    if(!numeric.is_finite()) return false;
    // Keep structural zeros so saturation cannot invalidate symbolic reuse.
    for(int column = 0; column < hessian.outerSize(); ++column){
      for(Eigen::SparseMatrix<double>::InnerIterator cell(hessian, column); cell; ++cell){
        cell.valueRef() = numeric(cell.row(), cell.col());
      }
    }
    ++factorizations;
    return factor.analyzePattern(hessian, &symbolicAnalyses) && factor.factorize(hessian);
  }
};

MultinomialModeContext & contextFrom(SEXP pointer){
  if(TYPEOF(pointer) != EXTPTRSXP ||
     R_ExternalPtrTag(pointer) != Rf_install("sommer_multinomial_mode_context") ||
     R_ExternalPtrAddr(pointer) == nullptr){
    Rcpp::stop("Invalid multinomial mode context; create a fresh context.");
  }
  return *Rcpp::XPtr<MultinomialModeContext>(pointer);
}

} // namespace

// [[Rcpp::export(rng = false)]]
SEXP sommer_multinomial_mode_create(const arma::vec & response,
    const arma::sp_mat & randomDesign, const arma::sp_mat & randomPrecision,
    const int dimension, const double denseMemoryMB = 1000.0){
  Rcpp::XPtr<MultinomialModeContext> pointer(new MultinomialModeContext(response,
    randomDesign, randomPrecision, dimension, denseMemoryMB), true);
  R_SetExternalPtrTag(pointer, Rf_install("sommer_multinomial_mode_context"));
  return pointer;
}

// [[Rcpp::export(rng = false)]]
Rcpp::List sommer_multinomial_mode_eval(SEXP pointer, const arma::mat & fixedPredictor,
    const arma::vec & modeStart, const int maxModeIterations = 100,
    const double tolerance = 1e-8){
  auto & context = contextFrom(pointer);
  if(fixedPredictor.n_rows != context.n || fixedPredictor.n_cols != context.d ||
     modeStart.n_elem != context.q){
    Rcpp::stop("Multinomial mode evaluation dimensions are incompatible.");
  }
  if(!fixedPredictor.is_finite() || !modeStart.is_finite() ||
     maxModeIterations < 1 || !std::isfinite(tolerance) || tolerance <= 0.0){
    Rcpp::stop("Multinomial mode evaluation requires finite predictors/starts and positive controls.");
  }
  ++context.evaluations;
  int iterations = 0;
  auto invalid = [&](const std::string & reason){
    return Rcpp::List::create(Rcpp::_["valid"]=false, Rcpp::_["reason"]=reason,
      Rcpp::_["iterations"]=iterations,
      Rcpp::_["work"]=context.work());
  };
  arma::vec mode = modeStart;
  const arma::vec fixed = arma::vectorise(fixedPredictor);
  arma::mat eta = arma::reshape(fixed + context.design * mode, context.n, context.d);
  if(!eta.is_finite()) return invalid("non-finite linear predictors");
  arma::mat probability;
  arma::vec score;
  arma::sp_mat curvature;
  double likelihood = context.evaluate(eta, probability, score, &curvature);
  double penalized = likelihood - 0.5 * arma::dot(mode, context.precision * mode);
  if(!std::isfinite(penalized)) return invalid("non-finite conditional likelihood");
  double scoreNorm = NA_REAL;
  double logDetHessian = NA_REAL;
  for(int iteration = 0; iteration <= maxModeIterations; ++iteration){
    Rcpp::checkUserInterrupt();
    const arma::vec modeScore = context.design.t() * score - context.precision * mode;
    if(!modeScore.is_finite()) return invalid("non-finite mode score");
    scoreNorm = arma::norm(modeScore, "inf");
    if(scoreNorm <= tolerance * (1.0 + arma::norm(mode, "inf"))){
      if(!context.factorize(curvature) ||
         !sommer::laplaceFactorLogDeterminant(context.factor, logDetHessian)){
        return invalid("non-positive-definite mode Hessian");
      }
      const double logLik = penalized + 0.5 * context.logDetPrecision -
        0.5 * logDetHessian;
      if(!std::isfinite(logLik)) return invalid("non-finite Laplace likelihood");
      const arma::sp_mat hessian =
        context.precision + context.design.t() * curvature * context.design;
      return Rcpp::List::create(Rcpp::_["valid"]=true, Rcpp::_["reason"]="converged",
        Rcpp::_["logLik"]=logLik, Rcpp::_["conditionalLogLik"]=likelihood,
        Rcpp::_["mode"]=mode, Rcpp::_["linearPredictors"]=eta,
        Rcpp::_["probability"]=probability, Rcpp::_["modeHessian"]=hessian,
        Rcpp::_["observationCurvature"]=curvature,
        Rcpp::_["modeScore"]=scoreNorm, Rcpp::_["iterations"]=iterations,
        Rcpp::_["logDetHessian"]=logDetHessian,
        Rcpp::_["work"]=context.work());
    }
    if(iteration == maxModeIterations) break;
    if(!context.factorize(curvature)) return invalid("non-positive-definite mode Hessian");
    Eigen::MatrixXd rhs(context.q, 1);
    for(arma::uword row = 0; row < context.q; ++row) rhs(row, 0) = modeScore(row);
    const Eigen::MatrixXd solution = context.factor.solve(rhs, "Multinomial mode solve failed.");
    arma::vec step(context.q);
    for(arma::uword row = 0; row < context.q; ++row) step(row) = solution(row, 0);
    if(!step.is_finite()) return invalid("non-finite Newton step");
    bool accepted = false;
    for(double scale = 1.0; scale >= 1.0 / 1024.0; scale *= 0.5){
      const arma::vec candidate = mode + scale * step;
      const arma::mat candidateEta = arma::reshape(fixed + context.design * candidate,
        context.n, context.d);
      if(!candidate.is_finite() || !candidateEta.is_finite()) continue;
      arma::mat candidateProbability;
      arma::vec candidateScore;
      const double candidateLikelihood =
        context.evaluate(candidateEta, candidateProbability, candidateScore, nullptr);
      const double candidatePenalized = candidateLikelihood -
        0.5 * arma::dot(candidate, context.precision * candidate);
      const double roundoff = 64.0 * std::numeric_limits<double>::epsilon() *
        (1.0 + std::abs(penalized));
      if(!std::isfinite(candidatePenalized) || candidatePenalized < penalized - roundoff) continue;
      const arma::vec candidateModeScore =
        context.design.t() * candidateScore - context.precision * candidate;
      if(!candidateModeScore.is_finite()) continue;
      if(candidatePenalized <= penalized + roundoff &&
         arma::norm(candidateModeScore, "inf") >= scoreNorm &&
         arma::norm(candidateModeScore, "inf") >
           tolerance * (1.0 + arma::norm(candidate, "inf"))) continue;
      mode = candidate;
      eta = candidateEta;
      likelihood = context.evaluate(eta, probability, score, &curvature);
      penalized = candidatePenalized;
      accepted = true;
      break;
    }
    iterations = iteration + 1;
    if(!accepted) return invalid("mode line search failed");
  }
  return invalid("mode iteration limit reached");
}

// [[Rcpp::export(rng = false)]]
Rcpp::List sommer_multinomial_laplace_eval(SEXP pointer,
    const arma::mat & fixedPredictor, const arma::sp_mat & randomPrecision,
    const arma::sp_mat & fixedDesign, const Rcpp::List & precisionDerivatives,
    const arma::vec & modeStart, const int maxModeIterations = 100,
    const double tolerance = 1e-8, const bool gradient = true){
  auto & context = contextFrom(pointer);
  if(fixedDesign.n_rows != context.n * context.d || !fixedDesign.is_finite()){
    Rcpp::stop("Multinomial fixed design must be finite with n(K-1) category-major rows.");
  }
  std::vector<sommer::LaplacePrecisionDerivative> derivatives;
  derivatives.reserve(precisionDerivatives.size());
  for(int parameter = 0; parameter < precisionDerivatives.size(); ++parameter){
    const arma::sp_mat derivative = Rcpp::as<arma::sp_mat>(precisionDerivatives[parameter]);
    if(derivative.n_rows != context.q || derivative.n_cols != context.q ||
       !derivative.is_finite() || !derivative.is_symmetric(1e-12)){
      Rcpp::stop("Multinomial precision derivatives must be finite symmetric q by q matrices.");
    }
    derivatives.push_back({0, derivative, NA_REAL});
  }
  context.updatePrecision(randomPrecision);
  Rcpp::List result = sommer_multinomial_mode_eval(pointer, fixedPredictor, modeStart,
    maxModeIterations, tolerance);
  if(gradient && Rcpp::as<bool>(result["valid"])){
    result["gradient"] = context.gradient(result, fixedDesign, derivatives);
    result["work"] = context.work();
  }
  return result;
}

namespace {

class MultinomialObjective {
public:
  arma::sp_mat fixedDesign;
  arma::mat offset;
  sommer::LaplacePrecision precision;
  Rcpp::RObject modePointer;
  arma::uword fixedCount;
  arma::uword covarianceCount;
  int evaluations = 0;
  bool cachedMode = false;
  bool cachedGradient = false;
  arma::vec cachedParameters;
  arma::vec cachedScore;
  Rcpp::List cachedResult;
  int cachedMaxIterations = 0;
  double cachedTolerance = 0.0;

  MultinomialObjective(const arma::vec & response, const arma::sp_mat & fixed,
      const arma::sp_mat & random, const arma::mat & suppliedOffset,
      const Rcpp::List & descriptors, const Rcpp::List & relationships,
      const double memory)
    : fixedDesign(fixed), offset(suppliedOffset), precision(descriptors, relationships),
      fixedCount(fixed.n_cols), covarianceCount(precision.parameterCount()){
    if(offset.n_rows != response.n_elem || offset.n_cols == 0 ||
       fixed.n_rows != offset.n_elem || !fixed.is_finite() || !offset.is_finite()){
      Rcpp::stop("Multinomial objective requires finite category-major fixed design and n by (K-1) offsets.");
    }
    arma::vec initial(covarianceCount);
    arma::uword position = 0;
    for(int term = 0; term < descriptors.size(); ++term){
      const Rcpp::List descriptor = descriptors[term];
      const arma::vec parameters = Rcpp::as<arma::vec>(descriptor["par"]);
      const Rcpp::LogicalVector free = descriptor["free"];
      if(free.size() != static_cast<int>(parameters.n_elem)){
        Rcpp::stop("Multinomial covariance descriptor free flags have incompatible length.");
      }
      for(arma::uword parameter = 0; parameter < parameters.n_elem; ++parameter){
        if(free[parameter] == TRUE) initial(position++) = parameters(parameter);
      }
    }
    std::string reason;
    if(!precision.update(initial, reason)){
      Rcpp::stop("Invalid multinomial initial covariance: " + reason + ".");
    }
    modePointer = sommer_multinomial_mode_create(response, random, precision.matrix,
      static_cast<int>(offset.n_cols), memory);
  }
};

MultinomialObjective & objectiveFrom(SEXP pointer){
  if(TYPEOF(pointer) != EXTPTRSXP ||
     R_ExternalPtrTag(pointer) != Rf_install("sommer_multinomial_objective") ||
     R_ExternalPtrAddr(pointer) == nullptr){
    Rcpp::stop("Invalid multinomial objective; create a fresh context.");
  }
  return *Rcpp::XPtr<MultinomialObjective>(pointer);
}

} // namespace

// [[Rcpp::export(rng = false)]]
SEXP sommer_multinomial_objective_create(const arma::vec & response,
    const arma::sp_mat & fixedDesign, const arma::sp_mat & randomDesign,
    const arma::mat & offset, const Rcpp::List & descriptors,
    const Rcpp::List & relationships, const double denseMemoryMB = 1000.0){
  Rcpp::XPtr<MultinomialObjective> pointer(new MultinomialObjective(response,
    fixedDesign, randomDesign, offset, descriptors, relationships, denseMemoryMB),
    true, Rf_install("sommer_multinomial_objective"));
  return pointer;
}

// [[Rcpp::export(rng = false)]]
Rcpp::List sommer_multinomial_objective_eval(SEXP pointer,
    const arma::vec & parameters, const int maxModeIterations = 100,
    const double tolerance = 1e-8, const bool gradient = true){
  auto & objective = objectiveFrom(pointer);
  if(parameters.n_elem != objective.fixedCount + objective.covarianceCount){
    Rcpp::stop("Multinomial objective parameter vector has incompatible length.");
  }
  if(maxModeIterations < 1 || !std::isfinite(tolerance) || tolerance <= 0.0){
    Rcpp::stop("Multinomial mode controls must be positive.");
  }
  ++objective.evaluations;
  auto & context = contextFrom(objective.modePointer);
  auto work = [&](){
    Rcpp::List result = context.work();
    result["objectiveEvaluations"] = objective.evaluations;
    result["descriptorUpdates"] = objective.precision.descriptorUpdates;
    result["descriptorReuses"] = objective.precision.descriptorReuses;
    result["relationshipFactorizations"] = objective.precision.relationshipFactorizations;
    return result;
  };
  std::string reason;
  if(!parameters.is_finite() ||
     !objective.precision.update(parameters.tail(objective.covarianceCount), reason)){
    objective.cachedMode = false;
    objective.cachedGradient = false;
    return Rcpp::List::create(Rcpp::_["valid"]=false,
      Rcpp::_["reason"]=parameters.is_finite() ? reason : "non-finite parameter proposal",
      Rcpp::_["work"]=work());
  }
  context.updatePrecision(objective.precision.matrix);
  const bool reuse = objective.cachedMode &&
    objective.cachedMaxIterations == maxModeIterations &&
    objective.cachedTolerance == tolerance &&
    arma::all(parameters == objective.cachedParameters);
  Rcpp::List result;
  if(reuse){
    result = Rcpp::clone(objective.cachedResult);
  }else{
    objective.cachedMode = false;
    objective.cachedGradient = false;
    const arma::mat predictor = objective.offset + arma::reshape(
      objective.fixedDesign * parameters.head(objective.fixedCount), context.n, context.d);
    result = sommer_multinomial_mode_eval(objective.modePointer, predictor,
      arma::vec(context.q, arma::fill::zeros), maxModeIterations, tolerance);
    if(Rcpp::as<bool>(result["valid"])){
      objective.cachedResult = Rcpp::clone(result);
      objective.cachedParameters = parameters;
      objective.cachedMaxIterations = maxModeIterations;
      objective.cachedTolerance = tolerance;
      objective.cachedMode = true;
    }
  }
  if(gradient && Rcpp::as<bool>(result["valid"])){
    if(!objective.cachedGradient){
      objective.cachedScore = context.gradient(result, objective.fixedDesign,
        objective.precision.derivatives());
      objective.cachedGradient = true;
    }
    result["gradient"] = objective.cachedScore;
  }
  result["work"] = work();
  return result;
}
