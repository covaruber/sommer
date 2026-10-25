#include "mme_sp2_support.h"
#include "mme_gpu.h"
#include "mme_rr.h"
#include "mme_saddle.h"
#include "mme_exact.h"
#include "mme_em.h"
#include "mme_scale.h"
#include "mme_sparse.h"
#include "mme_cholmod.h"
#include "mme_posterior.h"
#include "mme_pcg.h"
#include "mme_pcg_estimators.h"
#include "mme_sufficient.h"
#include "mme_layout.h"
#include "mme_progress.h"
#include <mutex>
#include <numeric>
#include <unordered_set>

static Rcpp::List ai_mme_sp2_impl(const arma::sp_mat & X, const Rcpp::List & ZI,
                     const arma::vec & Zind, const Rcpp::List & AiI,
                     const arma::sp_mat & y0,
                     const arma::sp_mat & H, const bool & useH,
                     const arma::uvec & residualBlockI,
                     const arma::uvec & residualIndexI,
                     int nIters, double tolParConvLL, double tolParConvNorm,
                     double tolParInv, const Rcpp::List & covStructI,
                     const arma::vec & weightEmInf,
                     const arma::vec & weightInf, const bool & verbose,
                     const int & computeCi = 0,
                     const std::string & solver = "ldlt",
                     const double & pcgTol = 1.0e-8,
                     const int & pcgMaxIters = 0,
                     const int & pcgTraceProbesInput = 8,
                     const int & pcgLanczosStepsInput = 20,
                     const bool & reml = true,
                     const bool & responsePrepared = false,
                     const double & preparedMean = 0.0,
                     const double & preparedSd = 1.0,
                     const bool & preparedIntercept = false,
                     const bool gpu = false,
                     const int gpuMinSize = 256,
                     const sommer::RRConditionalCallback * rrCallback = nullptr
){
  int pcgTraceProbes = pcgTraceProbesInput;
  int pcgLanczosSteps = pcgLanczosStepsInput;
  const auto fitOnce = [&](const Rcpp::List & design, const arma::vec & index,
      const Rcpp::List & covariance, const arma::vec & emWeights,
      const arma::vec & scoreWeights, bool fitVerbose, int inverseMode,
      const sommer::RRConditionalCallback * callback) -> Rcpp::List {
    return ai_mme_sp2_impl(X, design, index, AiI, y0, H, useH,
      residualBlockI, residualIndexI, 1, tolParConvLL, tolParConvNorm,
      tolParInv, covariance, emWeights, scoreWeights, fitVerbose, inverseMode,
      solver, pcgTol, pcgMaxIters, pcgTraceProbes, pcgLanczosSteps, reml,
      responsePrepared, preparedMean, preparedSd, preparedIntercept,
      gpu, gpuMinSize, callback);
  };
  const sommer::control::Transport controls(covStructI);
  sommer::control::Representation representation;
  representation.commonScale = controls.string("commonScale", "off") != "off";
  if(representation.commonScale){
    if((solver != "ldlt" && solver != "cholmod") || rrCallback ||
       sommer::hasZeroSpecific(covStructI) || arma::any(weightEmInf > 0.0) ||
       controls.method().proposal.accelerationRequested ||
       (covStructI.hasAttribute("saddlePoint") &&
        Rcpp::as<std::string>(covStructI.attr("saddlePoint")) == "force")){
      Rcpp::stop("Common-scale profiling requires deterministic Gaussian Henderson AI without positive blending, acceleration, saddle-point or latent reduced rank.");
    }
    return sommer::fitCommonScale(covStructI, nIters,
      tolParConvLL, tolParConvNorm, tolParInv, weightInf, computeCi, verbose, reml,
      [&](const Rcpp::List & covariance, int inverseMode) {
        return fitOnce(ZI, Zind, covariance,
          arma::vec(1, arma::fill::zeros), arma::vec(1, arma::fill::ones),
          false, inverseMode, nullptr);
      });
  }
  sommer::control::Outputs outputs;
  outputs.covarianceInformation = controls.information();
  outputs.retainCoefficientMatrix = controls.logical("retainCoefficientMatrix", true);
  const bool retainCovarianceInformation =
    outputs.covarianceInformation == sommer::control::Information::Full;
  const std::string covarianceInformation = retainCovarianceInformation ? "full" : "none";
  const bool retainCoefficientMatrix = outputs.retainCoefficientMatrix;
  if(!retainCovarianceInformation && sommer::hasZeroSpecific(covStructI)){
    Rcpp::stop("covarianceInformation='none' does not support latent reduced-rank fitting.");
  }

  sommer::control::Proposal proposal;
  proposal.directionBlend = controls.logical("directionBlend");
  if(proposal.directionBlend) {
    if(solver != "ldlt" && solver != "cholmod") {
      Rcpp::stop("EM direction blending requires a deterministic Gaussian Henderson solver.");
    }
    Rcpp::List blended = sommer::fitEMBlend(X, ZI, Zind, AiI, H, useH, residualBlockI,
      residualIndexI, covStructI, nIters, tolParConvLL, tolParConvNorm,
      tolParInv, weightEmInf, weightInf, computeCi, verbose, reml,
      [&](const Rcpp::List & design, const arma::vec & index,
          const Rcpp::List & covariance, int inverseMode,
          const sommer::RRConditionalCallback * callback) {
        return fitOnce(design, index, covariance,
          arma::vec(1, arma::fill::zeros), arma::vec(1, arma::fill::ones),
          false, inverseMode, callback);
      });
    if(!retainCovarianceInformation){ blended["InfMat"] = R_NilValue; }
    return blended;
  }
  const sommer::control::Method method = controls.method();
  sommer::control::CurvaturePolicy curvature = method.curvature;
  proposal.extrapolation = method.proposal.extrapolation;
  proposal.accelerationRequested = method.proposal.accelerationRequested;
  const bool hybridRequested = curvature.method == sommer::control::Curvature::HybridAI;
  if(hybridRequested && (solver != "ldlt" && solver != "cholmod")){
    Rcpp::stop("Hybrid AI requires a deterministic Henderson factorization solver.");
  }
  if(hybridRequested && (rrCallback || sommer::hasZeroSpecific(covStructI))){
    Rcpp::stop("Hybrid AI does not support latent reduced-rank Henderson fitting.");
  }

  if(verbose){
#ifdef _OPENMP
    Rcpp::Rcout
      << "OpenMP configured maximum: "
      << omp_get_max_threads()
      << " threads (not BLAS threads or measured usage)."
      << arma::endl;
#else
    Rcpp::Rcout
      << "OpenMP unavailable: Sommer OpenMP loops are serial; BLAS may still be threaded."
      << arma::endl;
#endif
  }
  representation.conditionalReducedRank = controls.logical("rrConditional");
  const bool rrConditional = representation.conditionalReducedRank;
  curvature.informationOnly = controls.logical("emInformationOnly");
  const bool emInformationOnly = curvature.informationOnly;
  if(emInformationOnly && nIters != 1){
    Rcpp::stop("EM final information evaluation requires exactly one iteration.");
  }
  const bool saddleRequested = covStructI.hasAttribute("saddlePoint") &&
    Rcpp::as<std::string>(covStructI.attr("saddlePoint")) == "force";
  if(!retainCoefficientMatrix && (!saddleRequested || computeCi != 0)){
    Rcpp::stop("retainCoefficientMatrix=FALSE requires saddlePoint='force' and computeCi=0.");
  }
  if(saddleRequested && (rrConditional || sommer::hasZeroSpecific(covStructI) ||
                        hybridRequested)){
    Rcpp::stop("Saddle-point fitting does not support latent reduced rank or hybrid AI.");
  }
  if(!rrConditional && sommer::hasZeroSpecific(covStructI)){
    if(solver != "ldlt" && solver != "cholmod"){
      Rcpp::stop("Exact latent reduced-rank Henderson fitting requires solver='ldlt' or solver='cholmod'.");
    }
    return sommer::fitReducedRank(
      X, ZI, Zind, covStructI, AiI, y0, H, useH, residualBlockI, residualIndexI,
      nIters, tolParConvLL, tolParConvNorm, computeCi, verbose,
      weightEmInf, weightInf, reml, responsePrepared, preparedMean, preparedSd,
      preparedIntercept,
      [&](const Rcpp::List & latentZ, const arma::vec & latentIndex,
          const Rcpp::List & covariance, int inverseMode,
          const sommer::RRConditionalCallback * callback) -> Rcpp::List {
        return fitOnce(latentZ, latentIndex, covariance, weightEmInf, weightInf,
          false, inverseMode, callback);
      });
  }

  sommer::gpu::DenseBackend gpuBackend(gpu, gpuMinSize);

  if(computeCi < 0 || computeCi > 2){
    Rcpp::stop(
      "computeCi must be 0, 1, or 2: "
      "0 = no final C-inverse work; "
      "1 = Takahashi sparse inverse subset; "
      "2 = full C inverse."
    );
  }

  // Linear-system backend used for C x = rhs solves.
  // In solver="pcg" mode the MME matrix C is NEVER factorised by LDLT:
  // log|C| is estimated by deterministic-probe stochastic Lanczos quadrature
  // and C^{-1} trace terms by Hutchinson probes solved with PCG.
  std::string solverName = solver;
  std::transform(solverName.begin(), solverName.end(), solverName.begin(),
                 [](unsigned char c){ return static_cast<char>(std::tolower(c)); });
  if(solverName != "ldlt" && solverName != "pcg" && solverName != "cholmod"){
    Rcpp::stop("solver must be 'ldlt', 'pcg', or 'cholmod'.");
  }
  std::string pcgPreconditioner = "diagonal";
  int pcgNystromRank = static_cast<int>(sommer::WorkspacePolicy::nystromRank());
  bool pcgWarmStart = false;
  int pcgRecycleRank = 0;
  std::string pcgRecycleMethod = "coarse";
  bool pcgExportRecycle = false;
  std::string compressionMethod = "off";
  std::string pcgRhsMethod = "independent";
  int pcgBlockSize = 16;
  int pcgResidualRefresh = 0;
  std::string pcgLogDetScaling = "none";
  std::string pcgTraceEstimator = "hutchinson";
  int pcgTraceRank = 16;
  std::string pcgAccuracy = "fixed";
  int pcgMaxTraceProbes = pcgTraceProbes;
  int pcgMaxLanczosSteps = pcgLanczosSteps;
  double pcgLogDetTol = 1.0e-3;
  double pcgScoreTol = 1.0e-3;
  int pcgSeed = 0;
  if(covStructI.hasAttribute("pcgControl")){
    const Rcpp::List pcgControl = Rcpp::as<Rcpp::List>(covStructI.attr("pcgControl"));
    if(pcgControl.size() > 0){
      auto setString = [&pcgControl](const char * key, std::string & value) {
        if(pcgControl.containsElementNamed(key)){
          value = Rcpp::as<std::string>(pcgControl[key]);
        }
      };
      auto setInteger = [&pcgControl](const char * key, int & value) {
        if(pcgControl.containsElementNamed(key)){
          value = Rcpp::as<int>(pcgControl[key]);
        }
      };
      auto setDouble = [&pcgControl](const char * key, double & value) {
        if(pcgControl.containsElementNamed(key)){
          value = Rcpp::as<double>(pcgControl[key]);
        }
      };
      auto setLogical = [&pcgControl](const char * key, bool & value) {
        if(pcgControl.containsElementNamed(key)){
          value = Rcpp::as<bool>(pcgControl[key]);
        }
      };
      setString("preconditioner", pcgPreconditioner);
      setInteger("nystromRank", pcgNystromRank);
      setLogical("warmStart", pcgWarmStart);
      setInteger("recycleRank", pcgRecycleRank);
      setString("recycleMethod", pcgRecycleMethod);
      setLogical("exportRecycle", pcgExportRecycle);
      setString("compression", compressionMethod);
      setString("rhsMethod", pcgRhsMethod);
      setInteger("blockSize", pcgBlockSize);
      setInteger("residualRefresh", pcgResidualRefresh);
      setString("logDetScaling", pcgLogDetScaling);
      setString("traceEstimator", pcgTraceEstimator);
      setInteger("traceRank", pcgTraceRank);
      setString("accuracy", pcgAccuracy);
      setInteger("maxTraceProbes", pcgMaxTraceProbes);
      setInteger("maxLanczosSteps", pcgMaxLanczosSteps);
      setDouble("logDetTol", pcgLogDetTol);
      setDouble("scoreTol", pcgScoreTol);
      setInteger("seed", pcgSeed);
    }
  }
  if(covStructI.hasAttribute("pcgPreconditioner")){
    pcgPreconditioner = Rcpp::as<std::string>(covStructI.attr("pcgPreconditioner"));
  }
  if(covStructI.hasAttribute("compression")){
    compressionMethod = Rcpp::as<std::string>(covStructI.attr("compression"));
  }
  if(covStructI.hasAttribute("pcgNystromRank")){
    pcgNystromRank = Rcpp::as<int>(covStructI.attr("pcgNystromRank"));
  }
  if(pcgPreconditioner != "diagonal" && pcgPreconditioner != "nystrom" &&
     pcgPreconditioner != "block-jacobi"){
    Rcpp::stop("pcgPreconditioner must be 'diagonal', 'nystrom', or 'block-jacobi'.");
  }
  if(pcgNystromRank < 1){
    Rcpp::stop("pcgNystromRank must be positive.");
  }
  if(pcgBlockSize < 1){
    Rcpp::stop("pcgControl$blockSize must be positive.");
  }
  if(pcgResidualRefresh < 0){
    Rcpp::stop("pcgControl$residualRefresh must be non-negative.");
  }
  if(pcgRecycleRank < 0){
    Rcpp::stop("pcgControl$recycleRank must be non-negative.");
  }
  if(pcgRecycleMethod != "coarse" && pcgRecycleMethod != "krylov"){
    Rcpp::stop("pcgControl$recycleMethod must be 'coarse' or 'krylov'.");
  }
  if(pcgTraceRank < 1 || pcgMaxTraceProbes < 1 || pcgMaxLanczosSteps < 2){
    Rcpp::stop("pcgControl accuracy settings are invalid.");
  }
  if(pcgLogDetTol <= 0.0 || pcgScoreTol <= 0.0){
    Rcpp::stop("pcgControl logDetTol and scoreTol must be positive.");
  }
  if(!std::isfinite(pcgTol) || pcgTol <= 0.0){
    Rcpp::stop("pcgTol must be positive and finite.");
  }
  if(pcgMaxIters < 0){
    Rcpp::stop("pcgMaxIters must be >= 0; 0 selects an automatic limit.");
  }
  if(pcgTraceProbes < 1){
    Rcpp::stop("pcgTraceProbes must be >= 1.");
  }
  if(pcgLanczosSteps < 2){
    Rcpp::stop("pcgLanczosSteps must be >= 2.");
  }
  if(solverName == "pcg" && computeCi == 1){
    Rcpp::stop("computeCi=1 requires LDLT/Takahashi. Use computeCi=0 for a genuinely factorisation-free PCG fit, or computeCi=2 for an explicit PCG full inverse (small systems only).");
  }
  if(!reml && solverName != "ldlt" && solverName != "cholmod"){
    Rcpp::stop("reml=FALSE (maximum likelihood) currently requires solver='ldlt' or solver='cholmod'.");
  }

  const sommer::SolverPlan solverPlan{
    solverName == "ldlt" ? sommer::SolverBackend::LDLT
      : solverName == "cholmod" ? sommer::SolverBackend::Cholmod
      : sommer::SolverBackend::PCG,
    reml, computeCi, pcgPreconditioner == "diagonal"};

  const Rcpp::Function getDenseOption = Rcpp::Environment::base_env()["getOption"];

  auto parseExpMode = [&](const std::string & attrName, const std::string & optName,
      sommer::control::Mode defaultMode) -> sommer::control::Mode {
    return controls.experimentalMode(attrName, optName, defaultMode, getDenseOption);
  };

  sommer::control::Elimination elimination;
  elimination.blockBanded = parseExpMode("blockBanded", "sommer.mme.blockBanded", sommer::control::Mode::Auto);
  elimination.factorSchurMultiFactor = parseExpMode("factorSchurMultiFactor", "sommer.mme.factorSchurMultiFactor", sommer::control::Mode::Auto);
  elimination.dFirstML = parseExpMode("dFirstML", "sommer.mme.dFirstML", sommer::control::Mode::Auto);
  elimination.matrixFreePcg = parseExpMode("matrixFreePcg", "sommer.mme.matrixFreePcg", sommer::control::Mode::Auto);
  elimination.cholmodFinalSelected = parseExpMode("cholmodFinalSelected", "sommer.mme.cholmodFinalSelected", sommer::control::Mode::Auto);
  elimination.pdCertificates = parseExpMode("pdCertificates", "sommer.mme.pdCertificates", sommer::control::Mode::Auto);
  const std::string blockBandedMode = sommer::control::modeName(elimination.blockBanded);
  const std::string factorSchurMultiFactorMode = sommer::control::modeName(elimination.factorSchurMultiFactor);
  const std::string dFirstMLMode = sommer::control::modeName(elimination.dFirstML);
  const std::string matrixFreePcgMode = sommer::control::modeName(elimination.matrixFreePcg);
  const std::string cholmodFinalSelectedMode = sommer::control::modeName(elimination.cholmodFinalSelected);
  const std::string pdCertificatesMode = sommer::control::modeName(elimination.pdCertificates);
  sommer::control::Execution execution;
  execution.profile = parseExpMode("profile", "sommer.mme.profile", sommer::control::Mode::Off);
  const std::string profileMode = sommer::control::modeName(execution.profile);
  sommer::PhaseProfile phaseProfile(execution.profile == sommer::control::Mode::Force);
  using sommer::ProfileTimer;
  double traceDerivativePreparationSeconds = 0.0;
  double traceRandomPrecisionDerivativeSeconds = 0.0;
  double traceFactorCovarianceDerivativeSeconds = 0.0;
  double traceFactorPrecisionDerivativeSeconds = 0.0;
  double traceRandomPrecisionOtherSeconds = 0.0;
  double traceResidualPrecisionDerivativeSeconds = 0.0;
  double traceGridDerivativeBasisSeconds = 0.0;
  double traceResidualDerivativeBasisSeconds = 0.0;
  double traceDerivativeOtherSeconds = 0.0;
  double traceSelectedInverseSeconds = 0.0;
  double sensitivitySolveSeconds = 0.0;
  double selectedInverseLookupRequests = 0.0;
  double selectedInverseLookupHits = 0.0;
  double selectedInverseLookupMisses = 0.0;
  double sensitivityRHSNonzeros = 0.0;
  double sensitivitySolveReachEntries = 0.0;
  double sensitivitySolveReachPossibleEntries = 0.0;
  double scoreTraceFallbackSeconds = 0.0;
  double aiAssemblyWallSeconds = 0.0;
  double aiRandomFirstTermSeconds = 0.0;
  double aiRandomCurvatureSeconds = 0.0;
  double aiResidualBlockSeconds = 0.0;
  execution.memory = controls.memory(getDenseOption);
  const std::string memoryPolicyMode = sommer::control::memoryName(execution.memory);
  const bool memoryStrict = execution.memory == sommer::control::Memory::Speed;
  elimination.blockSchur = parseExpMode("blockSchur", "sommer.mme.blockSchur", sommer::control::Mode::Auto);
  sommer::control::Traces traces;
  traces.batch = parseExpMode("traceBatch", "sommer.mme.traceBatch", sommer::control::Mode::Off);
  const std::string blockSchurMode = sommer::control::modeName(elimination.blockSchur);
  const std::string traceBatchMode = sommer::control::modeName(traces.batch);
  elimination.saddlePoint = parseExpMode("saddlePoint", "sommer.mme.saddlePoint", sommer::control::Mode::Off);
  const bool saddleActive = elimination.saddlePoint == sommer::control::Mode::Force;
  if(saddleActive && (rrConditional || rrCallback || solverPlan.pcg() ||
      proposal.accelerationRequested)){
    Rcpp::stop("Saddle-point fitting requires a deterministic solver without acceleration.");
  }
  sommer::SaddleSystem saddleJoint, saddleRandom;
  Eigen::SparseMatrix<double> saddleResidual, saddleDesign, saddlePrior;
  int saddleCompatibilityAssemblies = 0;
  int normalEquationAssemblies = 0;
  double saddleFactorSeconds = 0.0;

  int panelCi2Batch = static_cast<int>(sommer::WorkspacePolicy::inversePanelColumns());
  if(covStructI.hasAttribute("panelCi2")){
    panelCi2Batch = std::max(1, Rcpp::as<int>(covStructI.attr("panelCi2")));
  } else {
    const Rcpp::RObject opt = getDenseOption("sommer.mme.panelCi2");
    if(!opt.isNULL() && Rf_isNumeric(opt) && Rf_xlength(opt) == 1){
      panelCi2Batch = std::max(1, Rcpp::as<int>(opt));
    }
  }

  sommer::BlockSchurPolicy blockSchurPolicy;
  if(covStructI.hasAttribute("denseMemoryMB")){
    SEXP denseAttr = static_cast<SEXP>(covStructI.attr("denseMemoryMB"));
    if(!Rf_isNumeric(denseAttr) || Rf_xlength(denseAttr) != 1){
      Rcpp::stop("sommer.mme.denseMemoryMB must be a positive finite number.");
    }
    const double denseMemoryMB = Rcpp::as<double>(denseAttr);
    if(!std::isfinite(denseMemoryMB) || denseMemoryMB <= 0.0 ||
       denseMemoryMB > std::numeric_limits<double>::max() / 1.0e6){
      Rcpp::stop("sommer.mme.denseMemoryMB must be a positive finite number.");
    }
    blockSchurPolicy.maxDoubles = denseMemoryMB * 1.0e6 / sizeof(double);
  } else {
    const Rcpp::RObject denseMemoryOption = getDenseOption("sommer.mme.denseMemoryMB");
    if(!denseMemoryOption.isNULL()){
      if(!Rf_isNumeric(denseMemoryOption) || Rf_xlength(denseMemoryOption) != 1){
        Rcpp::stop("sommer.mme.denseMemoryMB must be a positive finite number.");
      }
      const double denseMemoryMB = Rcpp::as<double>(denseMemoryOption);
      if(!std::isfinite(denseMemoryMB) || denseMemoryMB <= 0.0 ||
         denseMemoryMB > std::numeric_limits<double>::max() / 1.0e6){
        Rcpp::stop("sommer.mme.denseMemoryMB must be a positive finite number.");
      }
      blockSchurPolicy.maxDoubles = denseMemoryMB * 1.0e6 / sizeof(double);
    }
  }
  bool enableBlockBandedEngine = (blockBandedMode == "force");


  const int nZs = ZI.size();
  const int nRRe = covStructI.size();
  if(nRRe < 1){
    Rcpp::stop("At least one covariance descriptor (the residual structure) is required.");
  }
  const int nRe = nRRe - 1;
  const int nZsFake = 1;
  const int nReFake = 1;
  const int nX = X.n_cols;
  const int nR = y0.n_rows;

  if(residualBlockI.n_elem != static_cast<arma::uword>(nR) ||
     residualIndexI.n_elem != static_cast<arma::uword>(nR)){
    Rcpp::stop("Residual block/index vectors must have one entry per observation.");
  }

  double vary = responsePrepared ? preparedSd * preparedSd : arma::mean(arma::var(y0));
  double stdy = responsePrepared ? preparedSd : arma::mean(arma::stddev(y0));
  if(!std::isfinite(vary) || vary <= 0.0){
    Rcpp::stop("Response variance must be positive and finite.");
  }

  arma::vec cShift(X.n_cols, arma::fill::zeros);
  bool hasIntercept = false;
  if(responsePrepared){
    hasIntercept = preparedIntercept;
    if(hasIntercept && X.n_cols > 0){
      cShift(0) = 1.0;
    }
  } else {
    hasIntercept = find_intercept_vector(X, cShift);
  }

  double muy = 0.0;
  if(responsePrepared){
    muy = preparedMean;
  } else if(hasIntercept){
    muy = arma::mean(arma::mean(y0));
  }

  arma::mat responseScaled = arma::mat(y0);
  if(!responsePrepared){
    if(hasIntercept){
      const arma::rowvec responseMeans = arma::mean(responseScaled, 0);
      responseScaled.each_row() -= responseMeans;
    }
    const arma::rowvec responseSds = arma::stddev(responseScaled, 0, 0);
    responseScaled.each_row() /= responseSds;
  }
  arma::sp_mat y = arma::sp_mat(responseScaled);
  bool intercept = hasIntercept;

  // ------------------------------------------------------------------
  // Generic arbitrary-Kronecker covariance descriptor interface.
  //
  // Every random/residual covariance structure is represented as
  //
  //   Sigma = exp(log_sigma2) * K1 (x) K2 (x) ... (x) Km
  //
  // with unconstrained working parameters. Each factor is a universal
  // CovarianceFactor descriptor carrying evaluator, derivative, reporting,
  // and trust-cap specifications. Statistical model names are not consumed
  // by the REML/MME optimizer.
  // ------------------------------------------------------------------
  typedef sommer::FactorPrecisionState DescriptorPrecisionState;
  sommer::CovarianceFactorEngine factorEngine;
  arma::field<arma::mat> theta(nRRe);
  arma::field<arma::vec> covPar(nRRe);
  arma::field<arma::vec> covConstraint(nRRe);
  arma::field<arma::vec> covScale(nRRe);
  arma::field<arma::vec> covLower(nRRe);
  arma::field<arma::vec> covUpper(nRRe);

  std::vector<std::string> covType(
    static_cast<std::size_t>(nRRe),
    "kron"
  );

  std::vector<Rcpp::List> covDescriptor;
  covDescriptor.reserve(static_cast<std::size_t>(nRRe));
  std::vector<arma::uword> covarianceDim(static_cast<std::size_t>(nRRe), 0);
  bool compactResidualCovariance = false;
  bool residualSectioned = false;

  int pdCertificates = 0;
  auto descriptorIsPositiveDefinite =
    [&](const Rcpp::List & cs,
        const arma::vec & par,
        const double floor) -> bool {

      if(par.n_elem < 1 || !par.is_finite()){
        return false;
      }

      Rcpp::List factors = cs["factors"];
      double logMinimumEigenvalue = par(0);
      std::vector<bool> certified(static_cast<std::size_t>(factors.size()), false);
      std::vector<double> minima(static_cast<std::size_t>(factors.size()));

      for(int fidx = 0; fidx < factors.size(); ++fidx){
        Rcpp::List f = Rcpp::as<Rcpp::List>(factors[fidx]);
        const int start1 = Rcpp::as<int>(f["par_start"]);
        const int end1 = Rcpp::as<int>(f["par_end"]);

        arma::vec localPar;
        if(end1 >= start1){
          localPar =
            par.subvec(
              static_cast<arma::uword>(start1 - 1),
              static_cast<arma::uword>(end1 - 1)
            );
        }

        arma::mat covariance;
        double factorMinimum = std::numeric_limits<double>::infinity();
        const int nSections = factorEngine.factorSections(f);
        for(int t = 0; t < std::max(nSections, 1); ++t){
          try{
            covariance = nSections > 0
              ? factorEngine.evalSectionFactor(f, localPar, t)
              : factorEngine.evalFactor(f, localPar);
          }catch(...){
            return false;
          }

          double minimum = 0.0;
          arma::mat offDiagonal = covariance;
          offDiagonal.diag().zeros();
          if(arma::accu(arma::abs(offDiagonal)) == 0.0){
            minimum = covariance.diag().min();
          }else if(sommer::precisionFamily(f) == sommer::PrecisionFamily::CompoundSymmetry &&
                   covariance.n_rows >= 2){
            const arma::vec standardDeviation = arma::sqrt(covariance.diag());
            const double rho = covariance(0, 1) / (standardDeviation(0) * standardDeviation(1));
            minimum = covariance.diag().min() * std::min(1.0 - rho,
              1.0 + static_cast<double>(covariance.n_rows - 1) * rho);
          }
          if(std::isfinite(minimum) && minimum > 0.0){
            certified[static_cast<std::size_t>(fidx)] = true;
          }else if(sommer::FactorPolicy::shiftedCholesky(factors.size(), covariance.n_rows)){
            // One factor: lambda_min > floor / sigma^2 is a Cholesky test, avoiding an O(q^3) eigendecomposition.
            const double threshold = std::exp(std::log(floor) - par(0));
            arma::mat shifted = covariance;
            shifted.diag() -= threshold;
            arma::mat shiftedFactor;
            if(!shifted.is_finite() || !gpuBackend.cholesky(shiftedFactor, shifted)){ return false; }
            minimum = 2.0 * threshold;
          }else{
            arma::vec eigenvalues;
            if(!arma::eig_sym(eigenvalues, covariance) || eigenvalues.n_elem == 0 ||
               !eigenvalues.is_finite() || eigenvalues.min() <= 0.0){ return false; }
            minimum = eigenvalues.min();
          }
          factorMinimum = std::min(factorMinimum, minimum);
        }
        minima[static_cast<std::size_t>(fidx)] = factorMinimum;
        logMinimumEigenvalue += std::log(factorMinimum);
      }

      if(std::isfinite(logMinimumEigenvalue) && logMinimumEigenvalue > std::log(floor)){
        if(pdCertificatesMode != "off" &&
           std::any_of(certified.begin(), certified.end(), [](bool value){ return value; })){
          ++pdCertificates;
          return true;
        }
      }
      for(int axis = 0; axis < factors.size(); ++axis){
        if(!certified[static_cast<std::size_t>(axis)]){ continue; }
        const Rcpp::List factor = Rcpp::as<Rcpp::List>(factors[axis]);
        const int start = Rcpp::as<int>(factor["par_start"]);
        const int end = Rcpp::as<int>(factor["par_end"]);
        const arma::vec local = end >= start ? arma::vec(par.subvec(start - 1, end - 1)) : arma::vec();
        const int sections = factorEngine.factorSections(factor);
        double exactMinimum = std::numeric_limits<double>::infinity();
        for(int section = 0; section < std::max(sections, 1); ++section){
          arma::vec eigenvalues;
          const arma::mat covariance = sections > 0 ? factorEngine.evalSectionFactor(factor, local, section) : factorEngine.evalFactor(factor, local);
          if(!arma::eig_sym(eigenvalues, covariance) || !eigenvalues.is_finite() ||
             eigenvalues.n_elem == 0 || eigenvalues.min() <= 0.0){ return false; }
          exactMinimum = std::min(exactMinimum, eigenvalues.min());
        }
        logMinimumEigenvalue += std::log(exactMinimum) - std::log(minima[static_cast<std::size_t>(axis)]);
      }

      return
        std::isfinite(logMinimumEigenvalue)
        &&
        logMinimumEigenvalue > std::log(floor);
    };

  for(int i = 0; i < nRRe; ++i){

    if(Rf_isNull(covStructI[i])){
      Rcpp::stop("NULL covariance descriptor supplied to ai_mme_sp2().");
    }

    Rcpp::List cs =
      Rcpp::as<Rcpp::List>(
        covStructI[i]
      );

    if(!cs.containsElementNamed("type") ||
       Rcpp::as<std::string>(cs["type"]) != "kron"){
      Rcpp::stop("ai_mme_sp2() accepts only the generic type='kron' covariance descriptor.");
    }
    if(!cs.containsElementNamed("descriptor_version") ||
       Rcpp::as<int>(cs["descriptor_version"]) < 2){
      Rcpp::stop("ai_mme_sp2() requires CovarianceFactor descriptor_version >= 2.");
    }

    covDescriptor.push_back(cs);
    covarianceDim[static_cast<std::size_t>(i)] =
      static_cast<arma::uword>(Rcpp::as<int>(cs["dim"]));

    const bool completeResidualBlocks =
      cs.containsElementNamed("complete_residual_blocks") &&
      Rcpp::as<bool>(cs["complete_residual_blocks"]);
    if(i == nRRe - 1){
      Rcpp::List residualFactors = cs["factors"];
      bool nonDiagonalFactor = false;
      for(int fidx = 0; fidx < residualFactors.size(); ++fidx){
        Rcpp::List factor = Rcpp::as<Rcpp::List>(residualFactors[fidx]);
        if(factor.containsElementNamed("section_levels")){
          residualSectioned = true;
        }
        if(!factor.containsElementNamed("structurally_diagonal") ||
           !Rcpp::as<bool>(factor["structurally_diagonal"])){
          nonDiagonalFactor = true;
        }
      }
      if(residualSectioned && useH){
        Rcpp::stop("dsumm() residual structures are not yet available with W weights.");
      }
      compactResidualCovariance = sommer::ResidualPolicy::compactStorage(
        residualSectioned, nonDiagonalFactor, completeResidualBlocks, useH);
    }

    covPar(i) =
      Rcpp::as<arma::vec>(
        cs["par"]
      );

    Rcpp::LogicalVector freeR =
      cs["free"];

    if(freeR.size() != static_cast<int>(covPar(i).n_elem)){
      Rcpp::stop("covStruct$free and covStruct$par have inconsistent lengths.");
    }

    covConstraint(i).set_size(covPar(i).n_elem);
    for(arma::uword k = 0; k < covPar(i).n_elem; ++k){
      covConstraint(i)(k) =
        freeR[static_cast<int>(k)]
        ?
        2.0
        :
        3.0;
    }

    covScale(i) =
      arma::vec(
        covPar(i).n_elem,
        arma::fill::zeros
      );

    covLower(i) = sommer::covarianceVarianceLowerBounds(
      Rcpp::List::create(cs), tolParInv, std::log(vary));

    covUpper(i) =
      arma::vec(
        covPar(i).n_elem,
        arma::fill::value(
          std::numeric_limits<double>::infinity()
        )
      );

    // Standardise only the product-level variance:
    // log(sigma2_internal) = log(sigma2_original) - log(var(y)).
    covPar(i)(0) -= std::log(vary);
    
    // Give the product-level log(sigma2) working parameter (always
    // covPar(i)(0), see vsm()) the same variance floor legacy type-1
    // parameters receive via applyVarianceBounds()/vcFloor, converted to
    // this working/standardised scale. Without this, a boundary-approaching
    // variance component (true value near/at zero) has no lower bound at
    // all, so an unconstrained Newton/AI step can jump arbitrarily far
    // past the PD-check floor in one iteration; the structure then gets
    // flagged "not PD" and, for single-parameter structures, permanently
    // frozen (see the covPar(iStruct).n_elem==1 fallback below).
    covLower(i)(0) =
      std::log(std::max(1.0e-8, tolParInv)) - std::log(vary);

    if(i == nRRe - 1 && compactResidualCovariance){
      theta(i).reset();
    }else{
      theta(i) = factorEngine.evaluateDescriptor(cs, covPar(i));
    }
  }

  auto covarianceD1 =
    [&](const int iStruct,
        const arma::uword k) -> arma::mat {

      return
        factorEngine.descriptorD1(
          covDescriptor[static_cast<std::size_t>(iStruct)],
          covPar(iStruct),
          k
        );
    };

  std::vector< std::vector<arma::mat> > covarianceDerivativeCache(
    static_cast<std::size_t>(nRRe)
  );

  std::vector< std::vector<bool> > covarianceDerivativeReady(
    static_cast<std::size_t>(nRRe)
  );

  auto cachedCovarianceD1 =
    [&](const int iStruct,
        const arma::uword k) -> const arma::mat & {

      const std::size_t structureOffset =
        static_cast<std::size_t>(iStruct);

      if(!covarianceDerivativeReady[structureOffset][k]){
        covarianceDerivativeCache[structureOffset][k] =
          covarianceD1(iStruct, k);
        covarianceDerivativeReady[structureOffset][k] = true;
      }

      return covarianceDerivativeCache[structureOffset][k];
    };

  auto evaluateStructure =
    [&](const int iStruct,
        const arma::vec & par) -> arma::mat {

      return
        factorEngine.evaluateDescriptor(
          covDescriptor[static_cast<std::size_t>(iStruct)],
          par
        );
    };

  auto tryEvaluateStructure =
    [&](const int iStruct,
        const arma::vec & par,
        arma::mat & value) -> bool {

      if(!par.is_finite()){
        return false;
      }

      if(iStruct == nRRe - 1 && compactResidualCovariance){
        value.reset();
        return descriptorIsPositiveDefinite(
          covDescriptor[static_cast<std::size_t>(iStruct)], par, 1.0e-12
        );
      }

      try{
        value =
          evaluateStructure(
            iStruct,
            par
          );
      }catch(const std::exception &){
        return false;
      }

      return
        value.n_rows
        ==
        theta(iStruct).n_rows
        &&
        value.n_cols
        ==
        theta(iStruct).n_cols
        &&
        value.is_finite();
    };

  // move Z to sparse arma objects
  int nZsAl; // integer to define the allocation of Z
  if(nZs > 0){
    nZsAl = nZs; // if there's random effects the nZs to allocate is equal to Z.size
  }else{
    nZsAl = nZsFake; // otherwise at least we allocate 1 element to avoid the program to crash
  }
  arma::field<arma::sp_mat> Z(nZsAl); // allocate size of Z
  if(nZs > 0){ // if there's random effects
    for (int i = 0; i < nZs; ++i) { // for each Z
      Z(i)=convertSparse(ZI(i)); // convert the matrix to sparse and store in the field
    }
  }
  // delete ZI;
  // Residual covariance is now built directly from residual block/local-index
  // metadata plus the same generic Kronecker descriptor used for random G.
  // No residual basis-list expansion is required.
  // move Ai to sparse arma objects
  int nReAl;
  if(nZs > 0){
    nReAl = nRe;
  }else{
    nReAl = nReFake;
  }
  arma::field<arma::sp_mat> Ai(nReAl); // allocate size of Ai field
  if(nZs > 0){
    for (int i = 0; i < nRe; ++i) {
      Ai(i)=convertSparse(AiI(i)); // convert the matrix to sparse and store in the field
    }
  }
  // calculate log determinants of Ai's
  arma::rowvec logDetA(nReAl);
  if(nZs > 0){ // of there's random effects
    for (int i = 0; i < nRe; ++i) { // for each random effect
      // Sparse LDLT log-determinant: avoids densifying the (potentially
      // large) relationship/pedigree precision matrix Ai.
      logDetA(i) = -sparseSpdLogDet(Ai(i));
    }
  }
  
  // define partitions (only used if random effects exist)
  int last = X.n_cols;
  arma::field<arma::mat> partitions(nReAl); // store indices of the random effects
  int Nu = 0;
  if(nZs > 0){ //if there's random effects (Z matrices) check where each starts and ends
    for (int i = 0; i < nRe; ++i) { // for each effect
      arma::uvec indexZind = find(Zind == (i+1) ); // which Z matrices to use , +1 because of the way indeces are used in C++
      int nIndexZind = indexZind.size(); //  number of Z matrices to use
      arma::vec Nus(nIndexZind); // vector to store number of columns in each Z matrix
      // for each matrix in this random effect
      for (int j = 0; j < nIndexZind; ++j) {
        int jj = indexZind(j); // thake the jj matrix
        arma::sp_mat Zprov = Z(jj); // put it in a provisional object
        Nus(j)=Zprov.n_cols; // calculate the number of columns
      }
      arma::vec end = Nus; // define ends and starts
      for (int k = 0; k < nIndexZind; ++k) { // for each effect
        arma::uvec toSum = arma::regspace<arma::uvec>(0,  1,  k); // equivalent to seq()
        end(k)=arma::accu(Nus(toSum));
      }
      arma::vec ones(nIndexZind, arma::fill::ones);
      arma::vec lastM(nIndexZind, arma::fill::value(last));
      arma::vec start = end - Nus + ones;
      start = start + lastM;// adjust start by adding # of fixed effects
      end = end + lastM;//adjust end by adding # of fixed effects
      partitions(i) = arma::join_rows(start,end);
      last = end.max();
      Nu = Nu + accu(Nus);
    }
  }// end of if statement when random effects exist
  
  // define the number of optimizer parameters per covariance structure
  arma::vec nVc(nRRe);
  for (int i = 0; i < nRRe; ++i) {
    nVc(i) =
      static_cast<double>(
        covPar(i).n_elem
      );
  }
  int nVcTotal = accu(nVc); // total number of variance components
  const VarianceParameterMap vcMap =
    readVarianceParameterMap(covStructI, static_cast<arma::uword>(nVcTotal));
  // assign a start and an end index to each covariance structure using the #of VC
  arma::vec nVcEnd = nVc;
  for (int i = 0; i < nRRe; ++i) {
    arma::uvec toSum = arma::regspace<arma::uvec>(0,  1,  i); // equivalent to seq()
    nVcEnd(i)=arma::accu(nVc(toSum));
  }
  arma::vec nVcStart = nVcEnd - nVc + 1;
  // removing complex structures how many effects are really there
  arma::vec nUsTotal(nReAl);
  if(nZs > 0){ //
    for (int i = 0; i < nRe; ++i) {
      arma::mat partitionsProv = partitions(i);
      nUsTotal(i) = partitionsProv(0,1) - partitionsProv(0,0) + 1;
    }
  }
  // define objects to store theta and llik across iterations
  arma::mat monitor(nVcTotal,nIters); // matrix to store variance components
  arma::rowvec llik(nIters); // store log likellihood values
  
  int nEffects = Nu+nX;
  // Stage 2 Eigen migration: W and C are the per-iteration hot-path sparse
  // objects (assembled every iteration and factorised every iteration).
  // They are native Eigen sparse matrices so the LDLT/PCG solvers consume
  // them directly, with no per-iteration Armadillo<->Eigen round trip.
  Eigen::SparseMatrix<double> W(nR,nEffects), C(nEffects,nEffects);
  arma::sp_mat Ci(nEffects,nEffects);
  arma::vec u(Nu), b(nX), bu(nEffects);
  const bool reducedAIAssembly = covStructI.hasAttribute("aiTangent");
  const arma::mat aiTangent = reducedAIAssembly
    ? Rcpp::as<arma::mat>(covStructI.attr("aiTangent")) : arma::mat();
  if(reducedAIAssembly &&
     (!emInformationOnly || retainCovarianceInformation ||
      aiTangent.n_rows != static_cast<arma::uword>(nVcTotal) ||
      aiTangent.n_cols >= aiTangent.n_rows || !aiTangent.is_finite() ||
      (aiTangent.n_cols > 0 && arma::rank(aiTangent) != aiTangent.n_cols))){
    Rcpp::stop("Reduced AI requires a finite independent tangent in an information-only conditional fit without covariance reporting.");
  }
  const arma::uword aiDimension = reducedAIAssembly ? aiTangent.n_cols : nVcTotal;
  const bool profileCommonScale = covStructI.hasAttribute("profileCommonScale") &&
    Rcpp::as<bool>(covStructI.attr("profileCommonScale"));
  if(profileCommonScale && (!reducedAIAssembly || computeCi != 0 ||
      solver == "pcg" || rrCallback || rrConditional)){
    Rcpp::stop("Common-scale conditional evaluation requires reduced exact AI without coefficient inverse work.");
  }
  arma::vec commonScaleQuadratic(profileCommonScale ? nVcTotal : 0, arma::fill::zeros);
  double commonScale = 1.0, commonScaleQuadraticTotal = 0.0;
  const double commonScaleDF = static_cast<double>(nR) - (reml ? nX : 0);
  arma::mat avInf(aiDimension,aiDimension);
  arma::mat emInf(emInformationOnly ? 0 : nVcTotal, emInformationOnly ? 0 : nVcTotal);
  arma::mat InfMat(aiDimension,aiDimension);
  arma::mat InfMatInv(retainCovarianceInformation ? nVcTotal : 0,
                      retainCovarianceInformation ? nVcTotal : 0);
  bool convergence = false;
  arma::vec delta(nVcTotal), delta_minus1(nVcTotal);
  // objects for constraints
  arma::mat percDelta(nVcTotal,nIters,arma::fill::zeros); // store % change of the delta with respect to the previous iteration
  arma::mat normMonitor(3,nIters); // store in each iteration the 3 stopping criteria of Madsen and Jensen
  arma::mat toBoundary(nIters,nVcTotal, arma::fill::zeros ); // store which values have been set to the boundary value
  arma::vec sumToBoundary(nVcTotal, arma::fill::zeros ); // to apply sum across iterations and if a VC goes to the boundary 3 times it is fixed to the boundary
  arma::sp_mat Hs(H.n_cols,H.n_cols); // square of H matrix
  if(useH == true){ // sparse Cholesky decomposition of H if user wants to use weights
    if(!rrConditional && !emInformationOnly){
      Rcpp::Rcout << "Using the weights matrix " << arma::endl;
    }
    // Stage 1 Eigen migration: factorise H directly in sparse form instead
    // of densifying it just to call gpuBackend.cholesky().
    const Eigen::SparseMatrix<double> He = armaSparseToEigenGlobal(H);
    Eigen::SparseMatrix<double> HeU;
    // All optimizer paths must apply weights in the same observation order.
    Eigen::SimplicialLLT<Eigen::SparseMatrix<double>, Eigen::Upper,
      Eigen::NaturalOrdering<int>> HeLLT;
    HeLLT.compute(He);
    if(HeLLT.info() != Eigen::Success){
      Rcpp::stop("Sparse Cholesky factorisation of the weights matrix H failed.");
    }
    HeU = HeLLT.matrixU();
    HeU.makeCompressed();
    Hs = eigenSparseToArmaGlobal(HeU);
  }
  // Hs^{-T}, used to weight residual derivative bases.
  arma::sp_mat HsInvT;
  const std::size_t maxFallbackEntries = sommer::WorkspacePolicy::sparseEntryLimit();
  int weightBlockCount = 0;
  std::vector<std::vector<int>> weightBlockRows;
  std::vector<int> weightBlockOfRow;
  std::vector<Eigen::SparseMatrix<double>> weightUpperBlocks;
  if(useH){
    const Eigen::SparseMatrix<double> upper = armaSparseToEigenGlobal(Hs);
    Rcpp::IntegerVector weightBlocks;
    if(covStructI.hasAttribute("weightBlocks")){
      weightBlocks = Rcpp::as<Rcpp::IntegerVector>(covStructI.attr("weightBlocks"));
      if(weightBlocks.size() != nR){
        Rcpp::stop("The weights block formula has incompatible length for the Cholesky factor.");
      }
    }
    if(weightBlocks.size()){
      std::map<int, std::vector<int>> groupedRows;
      for(int row = 0; row < nR; ++row){
        groupedRows[weightBlocks[row]].push_back(row);
      }
      weightBlockOfRow.assign(static_cast<std::size_t>(nR), -1);
      for(const auto & group : groupedRows){
        const int blockNumber = static_cast<int>(weightBlockRows.size());
        weightBlockRows.push_back(group.second);
        for(int row : group.second){ weightBlockOfRow[static_cast<std::size_t>(row)] = blockNumber; }
      }
      weightBlockCount = static_cast<int>(weightBlockRows.size());
      for(int col = 0; col < upper.outerSize(); ++col){
        for(Eigen::SparseMatrix<double>::InnerIterator entry(upper, col); entry; ++entry){
          if(weightBlockOfRow[static_cast<std::size_t>(entry.row())] !=
             weightBlockOfRow[static_cast<std::size_t>(entry.col())]){
            Rcpp::stop("The Cholesky factor of W couples rows declared independent by weights.");
          }
        }
      }
      std::vector<int> localOfRow(static_cast<std::size_t>(nR), -1);
      for(const std::vector<int> & globalRows : weightBlockRows){
        const int blockSize = static_cast<int>(globalRows.size());
        for(int local = 0; local < blockSize; ++local){
          localOfRow[static_cast<std::size_t>(globalRows[static_cast<std::size_t>(local)])] = local;
        }
        std::vector<Eigen::Triplet<double>> blockTriplets;
        for(int localCol = 0; localCol < blockSize; ++localCol){
          const int globalCol = globalRows[static_cast<std::size_t>(localCol)];
          for(Eigen::SparseMatrix<double>::InnerIterator entry(upper, globalCol); entry; ++entry){
            const int localRow = localOfRow[static_cast<std::size_t>(entry.row())];
            if(localRow >= 0){ blockTriplets.emplace_back(localRow, localCol, entry.value()); }
          }
        }
        Eigen::SparseMatrix<double> blockUpper(blockSize, blockSize);
        blockUpper.setFromTriplets(blockTriplets.begin(), blockTriplets.end());
        blockUpper.makeCompressed();
        weightUpperBlocks.push_back(std::move(blockUpper));
        for(int globalRow : globalRows){ localOfRow[static_cast<std::size_t>(globalRow)] = -1; }
      }
    }else{
      bool diagonalFactor = true;
      for(arma::sp_mat::const_iterator it = Hs.begin(); it != Hs.end(); ++it){
        if(it.row() != it.col() && std::abs(*it) > 0.0){ diagonalFactor = false; break; }
      }
      if(diagonalFactor){
        arma::vec hd(Hs.diag());
        HsInvT = arma::sp_mat(arma::diagmat(1.0 / hd));
      }else{
        const int batchSize = std::max(1, std::min(static_cast<int>(sommer::WorkspacePolicy::panelColumns()),
          static_cast<int>(sommer::WorkspacePolicy::cacheEntries()) / std::max(1, static_cast<int>(nR))));
        std::vector<Eigen::Triplet<double>> inverseTriplets;
        for(int first = 0; first < nR; first += batchSize){
          const int count = std::min(batchSize, static_cast<int>(nR) - first);
          Eigen::MatrixXd rhs = Eigen::MatrixXd::Zero(nR, count);
          for(int local = 0; local < count; ++local){ rhs(first + local, local) = 1.0; }
          const Eigen::MatrixXd solved = upper.triangularView<Eigen::Upper>().solve(rhs);
          for(int local = 0; local < count; ++local){
            for(int row = 0; row <= first + local; ++row){
              const double value = solved(row, local);
              if(value != 0.0){
                if(inverseTriplets.size() >= maxFallbackEntries){
                  Rcpp::stop("Non-diagonal weights produce excessive inverse fill; use weights=~group to declare independent blocks or use a factored residual model.");
                }
                inverseTriplets.emplace_back(first + local, row, value);
              }
            }
          }
        }
        Eigen::SparseMatrix<double> inverseTranspose(nR, nR);
        inverseTranspose.setFromTriplets(inverseTriplets.begin(), inverseTriplets.end());
        HsInvT = eigenSparseToArmaGlobal(inverseTranspose);
      }
    }
  }
  arma::vec dLuOut;//(nVcTotal); // we will join cols

  // ============================================================
  // CHANGE 1: cache objects that never change across iterations
  // ============================================================

  sommer::assembleMMEDesign(W, X, Z, nZs);

  std::vector<Eigen::SparseMatrix<double>> AiEigenCache(
    static_cast<std::size_t>(nRe)
  );
  if(solverName == "pcg" || saddleActive){
    for(int iR = 0; iR < nRe; ++iR){
      AiEigenCache[static_cast<std::size_t>(iR)] =
        armaSparseToEigenGlobal(Ai(iR));
      AiEigenCache[static_cast<std::size_t>(iR)].makeCompressed();
    }
  }

  // Dense response vector is constant.
  const arma::vec yDense = arma::vectorise(arma::mat(y));

  // Fixed/random coefficient indexes are constant.
  const arma::uvec bInd =
    arma::regspace<arma::uvec>(
      static_cast<arma::uword>(0),
      static_cast<arma::uword>(1),
      static_cast<arma::uword>(nX - 1)
    );

  arma::uvec uInd;
  if(nZs > 0){
    uInd =
      arma::regspace<arma::uvec>(
        static_cast<arma::uword>(nX),
        static_cast<arma::uword>(1),
        static_cast<arma::uword>(nX + Nu - 1)
      );
  }

  // Cache the Z indexes belonging to every random covariance structure
  // and the coefficient indexes for every partition.
  arma::field<arma::uvec> useZindCache(nReAl);
  std::vector< std::vector<arma::uvec> > partitionIndexCache(
    static_cast<std::size_t>(nRe)
  );
  std::vector< std::vector<arma::uword> > partitionStartCache(
    static_cast<std::size_t>(nRe)
  );
  std::vector< std::vector<arma::uword> > partitionEndCache(
    static_cast<std::size_t>(nRe)
  );

  if(nZs > 0){
    for(int iR = 0; iR < nRe; ++iR){

      useZindCache(iR) =
        arma::find(
          Zind == (iR + 1)
        );

      arma::mat partitionsP =
        partitions(iR);

      partitionIndexCache[static_cast<std::size_t>(iR)].resize(
        static_cast<std::size_t>(partitionsP.n_rows)
      );
      partitionStartCache[static_cast<std::size_t>(iR)].resize(
        static_cast<std::size_t>(partitionsP.n_rows)
      );
      partitionEndCache[static_cast<std::size_t>(iR)].resize(
        static_cast<std::size_t>(partitionsP.n_rows)
      );

      for(arma::uword j = 0; j < partitionsP.n_rows; ++j){

        const arma::uword first =
          static_cast<arma::uword>(
            partitionsP(j,0) - 1
          );

        const arma::uword lastIndex =
          static_cast<arma::uword>(
            partitionsP(j,1) - 1
          );

        partitionStartCache[static_cast<std::size_t>(iR)][static_cast<std::size_t>(j)] =
          first;

        partitionEndCache[static_cast<std::size_t>(iR)][static_cast<std::size_t>(j)] =
          lastIndex;

        partitionIndexCache[static_cast<std::size_t>(iR)][static_cast<std::size_t>(j)] =
          arma::regspace<arma::uvec>(
            first,
            static_cast<arma::uword>(1),
            lastIndex
          );
      }
    }
  }

  // One-time validation that Ai(iR) is dimensionally consistent with every
  // partition block, so the per-iteration G^-1 accumulation loop does not
  // need to repeat these checks on every REML iteration.
  if(nZs > 0){
    for(int iR = 0; iR < nRe; ++iR){
      const std::size_t iCache = static_cast<std::size_t>(iR);
      for(std::size_t j = 0; j < partitionStartCache[iCache].size(); ++j){
        const arma::uword blockSize =
          partitionEndCache[iCache][j] - partitionStartCache[iCache][j] + 1;
        if(blockSize != Ai(iR).n_rows || blockSize != Ai(iR).n_cols){
          Rcpp::stop(
            "Relationship inverse dimensions are inconsistent with a "
            "random-effect MME partition."
          );
        }
      }
    }
  }

  // ============================================================
  // Generic residual block/layout cache
  // ============================================================
  const int residualStructIndex = nRRe - 1;
  const int residualDim = static_cast<int>(
    covarianceDim[static_cast<std::size_t>(residualStructIndex)]
  );

  if(residualDim < 1){
    Rcpp::stop("Residual covariance dimension must be positive.");
  }

  std::vector<int> residualBlockOfRow(
    static_cast<std::size_t>(nR),
    -1
  );

  std::vector<int> residualLocalOfRow(
    static_cast<std::size_t>(nR),
    -1
  );

  int residualNBlocks = 0;

  for(int rr = 0; rr < nR; ++rr){

    if(residualBlockI(static_cast<arma::uword>(rr)) < 1 ||
       residualIndexI(static_cast<arma::uword>(rr)) < 1){
      Rcpp::stop("Residual block and local indices are 1-based positive integers.");
    }

    const int b =
      static_cast<int>(
        residualBlockI(static_cast<arma::uword>(rr)) - 1
      );

    const int local =
      static_cast<int>(
        residualIndexI(static_cast<arma::uword>(rr)) - 1
      );

    if(local < 0 || local >= residualDim){
      Rcpp::stop("Residual local covariance index exceeds the covariance-product dimension.");
    }

    residualBlockOfRow[static_cast<std::size_t>(rr)] = b;
    residualLocalOfRow[static_cast<std::size_t>(rr)] = local;
    residualNBlocks = std::max(residualNBlocks, b + 1);
  }

  std::vector< std::vector<arma::uword> > residualRowsTmp(
    static_cast<std::size_t>(residualNBlocks)
  );

  for(int rr = 0; rr < nR; ++rr){
    residualRowsTmp[
      static_cast<std::size_t>(
        residualBlockOfRow[static_cast<std::size_t>(rr)]
      )
    ].push_back(
      static_cast<arma::uword>(rr)
    );
  }

  std::vector<arma::uvec> residualKronRows(
    static_cast<std::size_t>(residualNBlocks)
  );

  bool residualKronBlocks = true;

  for(int b = 0; b < residualNBlocks; ++b){

    std::vector< std::pair<int, arma::uword> > ordered;
    ordered.reserve(
      residualRowsTmp[static_cast<std::size_t>(b)].size()
    );

    std::vector<int> seen(
      static_cast<std::size_t>(residualDim),
      0
    );

    for(const arma::uword rr :
        residualRowsTmp[static_cast<std::size_t>(b)]){

      const int local =
        residualLocalOfRow[static_cast<std::size_t>(rr)];

      if(seen[static_cast<std::size_t>(local)] != 0){
        Rcpp::stop("A residual block contains the same covariance-product coordinate more than once.");
      }

      seen[static_cast<std::size_t>(local)] = 1;
      ordered.emplace_back(local, rr);
    }

    std::sort(
      ordered.begin(),
      ordered.end(),
      [](const std::pair<int,arma::uword> & a,
         const std::pair<int,arma::uword> & b){
        return a.first < b.first;
      }
    );

    residualKronRows[static_cast<std::size_t>(b)].set_size(
      ordered.size()
    );

    for(std::size_t j = 0; j < ordered.size(); ++j){
      residualKronRows[static_cast<std::size_t>(b)](
        static_cast<arma::uword>(j)
      ) =
        ordered[j].second;
    }

    if(
        static_cast<int>(ordered.size()) != residualDim
    ){
      residualKronBlocks = false;
    }else{
      for(int local = 0; local < residualDim; ++local){
        if(seen[static_cast<std::size_t>(local)] == 0){
          residualKronBlocks = false;
          break;
        }
      }
    }
  }

  const int residualKronNBlocks =
    residualNBlocks;

  const int residualKronBlockSize =
    residualDim;

  std::vector< std::vector<arma::uword> > residualActiveColumnsTmp(
    static_cast<std::size_t>(residualNBlocks)
  );

  for(int wCol = 0; wCol < W.outerSize(); ++wCol){
    for(Eigen::SparseMatrix<double>::InnerIterator entry(W, wCol); entry; ++entry){
      const int block = residualBlockOfRow[
        static_cast<std::size_t>(entry.row())
      ];
      residualActiveColumnsTmp[static_cast<std::size_t>(block)].push_back(
        static_cast<arma::uword>(entry.col())
      );
    }
  }

  std::vector<arma::uvec> residualActiveColumns(
    static_cast<std::size_t>(residualNBlocks)
  );

  for(int block = 0; block < residualNBlocks; ++block){
    std::vector<arma::uword> & columns =
      residualActiveColumnsTmp[static_cast<std::size_t>(block)];
    std::sort(columns.begin(), columns.end());
    columns.erase(std::unique(columns.begin(), columns.end()), columns.end());
    residualActiveColumns[static_cast<std::size_t>(block)] =
      arma::uvec(columns);
  }

  // Structural diagonality is determined by the covariance factors, not by
  // current parameter values.
  bool residualStructurallyDiagonal = true;

  {
    Rcpp::List rcs =
      covDescriptor[
        static_cast<std::size_t>(
          residualStructIndex
        )
      ];

    Rcpp::List rfactors =
      rcs["factors"];

    for(int fidx = 0; fidx < rfactors.size(); ++fidx){
      Rcpp::List f =
        Rcpp::as<Rcpp::List>(
          rfactors[fidx]
        );

      if(!f.containsElementNamed("structurally_diagonal") ||
         !Rcpp::as<bool>(f["structurally_diagonal"])){
        residualStructurallyDiagonal = false;
        break;
      }
    }
  }

  // Same structural check, per random-effect covariance descriptor. This is
  // a property of the descriptor (which factors it is built from), not of
  // the current parameter values, so it is computed once here rather than
  // every REML iteration. Used to short-circuit the O(q^5) AI-matrix
  // second-derivative loop and the O(q^2) score/trace loop down to their
  // exact closed-form scalar equivalents when a random-effect covariance
  // product is diagonal (e.g. dsm() alone, or a Kronecker product of only
  // diagonal factors).
  std::vector<bool> randomStructurallyDiagonal(
    static_cast<std::size_t>(nRe),
    true
  );

  for(int iR = 0; iR < nRe; ++iR){
    Rcpp::List rcs =
      covDescriptor[static_cast<std::size_t>(iR)];

    Rcpp::List rfactors =
      rcs["factors"];

    bool diagonalHere = true;

    for(int fidx = 0; fidx < rfactors.size(); ++fidx){
      Rcpp::List f =
        Rcpp::as<Rcpp::List>(rfactors[fidx]);

      if(!f.containsElementNamed("structurally_diagonal") ||
         !Rcpp::as<bool>(f["structurally_diagonal"])){
        diagonalHere = false;
        break;
      }
    }

    randomStructurallyDiagonal[static_cast<std::size_t>(iR)] = diagonalHere;
  }

  // Cache whether the optional H Cholesky factor is diagonal.
  bool Hdiag = true;
  arma::vec HdiagSquared(
    nR,
    arma::fill::ones
  );

  if(useH){
    for(arma::sp_mat::const_iterator it = Hs.begin();
        it != Hs.end();
        ++it){

      if(it.row() != it.col() && std::abs(*it) > 0.0){
        Hdiag = false;
      }
    }

    if(Hdiag){
      for(int ii = 0; ii < nR; ++ii){
        const double h = Hs(ii,ii);
        HdiagSquared(ii) = h*h;
      }
    }
  }

  std::vector<int> residualPosInBlock(static_cast<std::size_t>(nR), -1);
  for(int b = 0; b < residualNBlocks; ++b){
    const arma::uvec & rows = residualKronRows[static_cast<std::size_t>(b)];
    for(arma::uword j = 0; j < rows.n_elem; ++j){
      residualPosInBlock[static_cast<std::size_t>(rows(j))] = static_cast<int>(j);
    }
  }

  // Ragged blocks embedded in their factor rectangle; missing cells are absorbed
  // exactly via R_oo^{-1} = Q_oo - Q_om Q_mm^{-1} Q_mo (ASReml missing-value form).
  bool residualGridBlocks = false;
  int gridNFactors = 0;
  std::vector< std::vector<arma::uvec> > gridLevels;
  std::vector<arma::uvec> gridObsPos;
  std::vector<arma::uvec> gridMissPos;
  std::vector<arma::umat> gridMissDigits;
  std::vector<arma::uword> gridSize;
  std::vector< std::vector<int> > gridSectionOf;

  if(sommer::ResidualPolicy::gridCandidate(compactResidualCovariance,
      residualKronBlocks, residualSectioned, residualStructurallyDiagonal, useH)){
    Rcpp::List rfactors =
      covDescriptor[static_cast<std::size_t>(residualStructIndex)]["factors"];
    gridNFactors = rfactors.size();
    std::vector<arma::uword> fdim(static_cast<std::size_t>(gridNFactors));
    std::vector<arma::uword> ftrail(static_cast<std::size_t>(gridNFactors));
    arma::uword trailing = static_cast<arma::uword>(residualDim);
    for(int f = 0; f < gridNFactors; ++f){
      Rcpp::List fl = Rcpp::as<Rcpp::List>(rfactors[f]);
      fdim[static_cast<std::size_t>(f)] = static_cast<arma::uword>(Rcpp::as<int>(fl["dim"]));
      trailing /= fdim[static_cast<std::size_t>(f)];
      ftrail[static_cast<std::size_t>(f)] = trailing;
    }

    residualGridBlocks = gridNFactors > 0;
    gridLevels.resize(static_cast<std::size_t>(residualNBlocks));
    gridObsPos.resize(static_cast<std::size_t>(residualNBlocks));
    gridMissPos.resize(static_cast<std::size_t>(residualNBlocks));
    gridMissDigits.resize(static_cast<std::size_t>(residualNBlocks));
    gridSize.resize(static_cast<std::size_t>(residualNBlocks));
    gridSectionOf.assign(static_cast<std::size_t>(residualNBlocks),
                         std::vector<int>(static_cast<std::size_t>(gridNFactors), -1));

    for(int b = 0; b < residualNBlocks && residualGridBlocks; ++b){
      const std::size_t bb = static_cast<std::size_t>(b);
      const arma::uvec & rows = residualKronRows[bb];
      const arma::uword nb = rows.n_elem;
      std::vector< std::vector<int> > posOfLevel(static_cast<std::size_t>(gridNFactors));
      gridLevels[bb].resize(static_cast<std::size_t>(gridNFactors));
      arma::uword S = 1;

      for(int f = 0; f < gridNFactors; ++f){
        const std::size_t ff = static_cast<std::size_t>(f);
        std::vector<arma::uword> lv;
        lv.reserve(nb);
        for(arma::uword j = 0; j < nb; ++j){
          const arma::uword local = static_cast<arma::uword>(
            residualLocalOfRow[static_cast<std::size_t>(rows(j))]);
          lv.push_back((local / ftrail[ff]) % fdim[ff]);
        }
        std::sort(lv.begin(), lv.end());
        lv.erase(std::unique(lv.begin(), lv.end()), lv.end());
        gridLevels[bb][ff] = arma::uvec(lv);
        posOfLevel[ff].assign(fdim[ff], -1);
        for(std::size_t i = 0; i < lv.size(); ++i){
          posOfLevel[ff][lv[i]] = static_cast<int>(i);
        }
        S *= static_cast<arma::uword>(lv.size());
      }

      if(!sommer::ResidualPolicy::gridFits(static_cast<double>(S), static_cast<double>(nb))){
        residualGridBlocks = false;
        break;
      }
      gridSize[bb] = S;

      std::vector<char> seen(S, 0);
      gridObsPos[bb].set_size(nb);
      for(arma::uword j = 0; j < nb; ++j){
        const arma::uword local = static_cast<arma::uword>(
          residualLocalOfRow[static_cast<std::size_t>(rows(j))]);
        arma::uword g = 0;
        for(int f = 0; f < gridNFactors; ++f){
          const std::size_t ff = static_cast<std::size_t>(f);
          const arma::uword digit = (local / ftrail[ff]) % fdim[ff];
          g = g * gridLevels[bb][ff].n_elem +
            static_cast<arma::uword>(posOfLevel[ff][digit]);
        }
        gridObsPos[bb](j) = g;
        seen[g] = 1;
      }
      std::vector<arma::uword> miss;
      miss.reserve(S - nb);
      for(arma::uword g = 0; g < S; ++g){
        if(!seen[g]){ miss.push_back(g); }
      }
      gridMissPos[bb] = arma::uvec(miss);
      gridMissDigits[bb].set_size(miss.size(), static_cast<arma::uword>(gridNFactors));
      for(std::size_t j = 0; j < miss.size(); ++j){
        arma::uword g = miss[j];
        for(int f = gridNFactors - 1; f >= 0; --f){
          const arma::uword nl = gridLevels[bb][static_cast<std::size_t>(f)].n_elem;
          gridMissDigits[bb](j, static_cast<arma::uword>(f)) = g % nl;
          g /= nl;
        }
      }

      for(int f = 0; f < gridNFactors; ++f){
        Rcpp::List fl = Rcpp::as<Rcpp::List>(rfactors[f]);
        if(!fl.containsElementNamed("section_factor")){
          continue;
        }
        const std::size_t owner =
          static_cast<std::size_t>(Rcpp::as<int>(fl["section_factor"]) - 1);
        if(gridLevels[bb][owner].n_elem != 1){
          Rcpp::stop("A residual block spans several levels of a by= section variable.");
        }
        gridSectionOf[bb][static_cast<std::size_t>(f)] =
          static_cast<int>(gridLevels[bb][owner](0));
      }
    }
  }

  if(residualSectioned && !residualGridBlocks){
    Rcpp::stop("dsumm() residual structures require unweighted, near-rectangular "
               "residual blocks within each section.");
  }

  // ------------------------------------------------------------
  // Cache the block design matrices used by the repeated-block/Kronecker
  // residual C-assembly path (including any diagonal-H row weighting).
  // W, Hs, and the block/column membership are all constant across REML
  // iterations, so this "blockW" is identical every iteration; precomputing
  // it once avoids rescanning every active column of W (previously
  // O(sum of active-column nnz over blocks) per iteration) inside the main
  // loop below.
  // ------------------------------------------------------------
  std::vector<arma::mat> residualKronBlockWCache(
    static_cast<std::size_t>(residualNBlocks)
  );

  // Only complete/repeated blocks (residualKronBlocks) guarantee
  // rows.n_elem == residualDim, which the local-coordinate indexing below
  // requires. When blocks are irregular the Rkron path is never taken at
  // runtime, so the cache is simply left empty.
  if(residualKronBlocks || residualGridBlocks){
  for(int block = 0; block < residualNBlocks; ++block){
    const arma::uvec & rows =
      residualKronRows[static_cast<std::size_t>(block)];
    const arma::uvec & columns =
      residualActiveColumns[static_cast<std::size_t>(block)];

    if(columns.n_elem == 0){ continue; }

    arma::mat blockW(
      rows.n_elem,
      columns.n_elem,
      arma::fill::zeros
    );

    for(arma::uword localCol = 0;
        localCol < columns.n_elem;
        ++localCol){
      for(Eigen::SparseMatrix<double>::InnerIterator designEntry(
            W, static_cast<int>(columns(localCol)));
          designEntry; ++designEntry){
        const arma::uword globalRow =
          static_cast<arma::uword>(designEntry.row());
        if(residualBlockOfRow[static_cast<std::size_t>(globalRow)] == block){
          // For complete blocks the block position equals the local coordinate.
          const int localPosition = residualPosInBlock[
            static_cast<std::size_t>(globalRow)
          ];
          blockW(static_cast<arma::uword>(localPosition), localCol) =
            designEntry.value();
        }
      }
    }

    if(useH && Hdiag){
      arma::vec hDiagonal(rows.n_elem);
      for(arma::uword localRow = 0; localRow < rows.n_elem; ++localRow){
        hDiagonal(localRow) = Hs(rows(localRow), rows(localRow));
      }
      blockW.each_col() %= hDiagonal;
    }

    residualKronBlockWCache[static_cast<std::size_t>(block)] =
      std::move(blockW);
  }
  }

  // ------------------------------------------------------------
  // Persistent residual (gr,gc,lr,lc) index list.
  //
  // residualKronRows/residualLocalOfRow are purely data-driven (fixed for
  // the whole optimisation), so the set of global row/col pairs feeding
  // Rmat and every residualDerivativeBasis(k) never changes across
  // iterations - only the numeric value read at (lr,lc) does.  Building
  // this list once and reusing it every iteration replaces the former
  // per-iteration index bookkeeping plus single-element arma::sp_mat
  // insertion with one batch sparse construction from cached locations.
  // ------------------------------------------------------------
  arma::umat residualPatternLocations;
  std::vector<int> residualPatternLr;
  std::vector<int> residualPatternLc;

  // Weighted fits (useH) assemble R in record space even for complete blocks.
  if(!residualGridBlocks && (!residualKronBlocks || residualStructurallyDiagonal || useH)){
    std::vector<arma::uword> patternRows;
    std::vector<arma::uword> patternCols;

    for(int b = 0; b < residualNBlocks; ++b){
      const arma::uvec & rows =
        residualKronRows[static_cast<std::size_t>(b)];

      for(arma::uword aa = 0; aa < rows.n_elem; ++aa){
        const arma::uword gr = rows(aa);
        const int lr = residualLocalOfRow[static_cast<std::size_t>(gr)];

        for(arma::uword bb = 0; bb < rows.n_elem; ++bb){
          const arma::uword gc = rows(bb);
          const int lc = residualLocalOfRow[static_cast<std::size_t>(gc)];

          patternRows.push_back(gr);
          patternCols.push_back(gc);
          residualPatternLr.push_back(lr);
          residualPatternLc.push_back(lc);
        }
      }
    }

    residualPatternLocations.set_size(2, patternRows.size());
    for(std::size_t k = 0; k < patternRows.size(); ++k){
      residualPatternLocations(0,k) = patternRows[k];
      residualPatternLocations(1,k) = patternCols[k];
    }
  }

  auto buildResidualSparseFromPattern =
    [&](const arma::mat & localMat) -> arma::sp_mat {
      const arma::uword nnz = residualPatternLr.size();
      arma::vec values(nnz);
      for(arma::uword k = 0; k < nnz; ++k){
        values(k) = localMat(
          static_cast<arma::uword>(residualPatternLr[k]),
          static_cast<arma::uword>(residualPatternLc[k])
        );
      }
      return arma::sp_mat(residualPatternLocations, values, nR, nR);
    };

  std::vector< std::vector<arma::uword> > residualDigits;
  if(compactResidualCovariance){
    Rcpp::List rfactors = covDescriptor[static_cast<std::size_t>(residualStructIndex)]["factors"];
    arma::uword trailing = static_cast<arma::uword>(residualDim);
    residualDigits.resize(static_cast<std::size_t>(rfactors.size()));
    for(int fidx = 0; fidx < rfactors.size(); ++fidx){
      Rcpp::List f = Rcpp::as<Rcpp::List>(rfactors[fidx]);
      const arma::uword dim = static_cast<arma::uword>(Rcpp::as<int>(f["dim"]));
      trailing /= dim;
      std::vector<arma::uword> & digits = residualDigits[static_cast<std::size_t>(fidx)];
      digits.resize(static_cast<std::size_t>(residualDim));
      for(int local = 0; local < residualDim; ++local){
        digits[static_cast<std::size_t>(local)] =
          (static_cast<arma::uword>(local) / trailing) % dim;
      }
    }
  }

  // derivativeK < 0 evaluates R itself; otherwise dR/dpar_k.
  auto buildCompactResidualSparse =
    [&](const arma::vec & par, const int derivativeK) -> arma::sp_mat {
      Rcpp::List rfactors = covDescriptor[static_cast<std::size_t>(residualStructIndex)]["factors"];
      std::vector<arma::mat> pieces(static_cast<std::size_t>(rfactors.size()));
      bool found = derivativeK <= 0;
      for(int fidx = 0; fidx < rfactors.size(); ++fidx){
        Rcpp::List f = Rcpp::as<Rcpp::List>(rfactors[fidx]);
        const int start1 = Rcpp::as<int>(f["par_start"]);
        const int end1 = Rcpp::as<int>(f["par_end"]);
        arma::vec localPar;
        if(end1 >= start1){
          localPar = par.subvec(static_cast<arma::uword>(start1 - 1),
                                static_cast<arma::uword>(end1 - 1));
        }
        const int k1 = derivativeK + 1;
        if(derivativeK > 0 && k1 >= start1 && k1 <= end1){
          pieces[static_cast<std::size_t>(fidx)] =
            factorEngine.factorD1(f, localPar, static_cast<arma::uword>(k1 - start1));
          found = true;
        }else{
          pieces[static_cast<std::size_t>(fidx)] = factorEngine.evalFactor(f, localPar);
        }
      }
      if(!found){
        Rcpp::stop("Residual covariance derivative parameter does not belong to any factor.");
      }

      const double scale = std::exp(par(0));
      const arma::uword nnz = residualPatternLr.size();
      arma::vec values(nnz);
      for(arma::uword k = 0; k < nnz; ++k){
        const std::size_t lr = static_cast<std::size_t>(residualPatternLr[k]);
        const std::size_t lc = static_cast<std::size_t>(residualPatternLc[k]);
        double value = scale;
        for(std::size_t fidx = 0; fidx < pieces.size(); ++fidx){
          value *= pieces[fidx](residualDigits[fidx][lr], residualDigits[fidx][lc]);
        }
        values(k) = value;
      }
      return arma::sp_mat(residualPatternLocations, values, nR, nR);
    };


  // CHANGE 2 helper:
  // exact Moore-Penrose inverse of
  // diag(theta^2 / denominator), without forming the diagonal
  // matrix and without an SVD-based pinv().
  // ------------------------------------------------------------

  // ------------------------------------------------------------
  // Information-system helper:
  // solve A x = rhs directly; no inverse is formed.
  // ------------------------------------------------------------
  auto solveInformationSystem =
    [&](arma::vec & solution,
        const arma::mat & A,
        const arma::vec & rhs,
        const std::string & context) -> bool {

      bool ok =
        arma::solve(
          solution,
          A,
          rhs,
          arma::solve_opts::likely_sympd
        );

      if(!ok){
        ok =
          arma::solve(
            solution,
            A,
            rhs
          );
      }

      if(!ok){
        Rcpp::Rcout
          << "Information-system solve failed in "
          << context
          << arma::endl;
      }

      return ok;
    };


  // ============================================================
  // Sparse LDLT / Takahashi selected-inverse infrastructure
  // ============================================================
  typedef Eigen::SparseMatrix<double, Eigen::ColMajor, int> EigenSpMat;
  typedef Eigen::Triplet<double, int> EigenTriplet;
  typedef Eigen::SimplicialLDLT<
    EigenSpMat,
    Eigen::Lower,
    SommerSparseOrdering
  > EigenLDLT;

  struct SelectedInverseSubset {
    std::vector< std::vector<int> > rows;
    std::vector< std::vector<double> > values;
    std::vector< std::unordered_map<int, std::size_t> > rowLookup;
    std::vector<char> contiguous;
    std::vector< std::vector<std::size_t> > dependencyOffsets;
    std::size_t dependencyEntries = 0;
    std::vector<int> originalToPermuted;
  };
  std::size_t selectedDependencyMaxEntries = 0;
  std::unordered_map<const SelectedInverseSubset *, std::vector<std::uint64_t>>
    selectedInverseDemandKeys;
  constexpr std::size_t maxSelectedInverseDemandKeys = 1000000;
  std::size_t selectedInverseDemandKeyCount = 0;
  bool selectedInverseDemandTrackingComplete = true;
  std::mutex selectedInverseDemandMutex;
  std::size_t selectedInverseDemandUniqueEntries = 0;
  std::size_t selectedInverseDemandClosureEntries = 0;
  std::size_t selectedInverseDemandClosureOperations = 0;

  auto finalizeSelectedInverseDemand =
    [&](const SelectedInverseSubset & subset,
        std::vector<std::uint64_t> & keys) -> void {
      std::sort(keys.begin(), keys.end());
      keys.erase(std::unique(keys.begin(), keys.end()), keys.end());
      selectedInverseDemandUniqueEntries += keys.size();
      std::unordered_set<std::uint64_t> closure(keys.begin(), keys.end());
      std::vector<std::uint64_t> pending(keys.begin(), keys.end());
      while(!pending.empty()){
        const std::uint64_t key = pending.back();
        pending.pop_back();
        const int column = static_cast<int>(static_cast<std::uint32_t>(key));
        const int row = static_cast<int>(static_cast<std::uint32_t>(key >> 32));
        const std::vector<int> & columnRows = subset.rows[static_cast<std::size_t>(column)];
        selectedInverseDemandClosureOperations += columnRows.size() - 1;
        for(std::size_t neighbor = 1; neighbor < columnRows.size(); ++neighbor){
          const int first = columnRows[neighbor];
          const int dependencyColumn = std::min(first, row);
          const int dependencyRow = std::max(first, row);
          const std::vector<int> & dependencyRows =
            subset.rows[static_cast<std::size_t>(dependencyColumn)];
          const auto dependency = std::lower_bound(
            dependencyRows.begin(), dependencyRows.end(), dependencyRow);
          if(dependency == dependencyRows.end() || *dependency != dependencyRow){
            Rcpp::stop("Selected-inverse demand closure is outside the filled LDLT pattern.");
          }
          const std::uint64_t dependencyKey =
            (static_cast<std::uint64_t>(static_cast<std::uint32_t>(dependencyRow)) << 32) |
            static_cast<std::uint32_t>(dependencyColumn);
          if(closure.insert(dependencyKey).second){ pending.push_back(dependencyKey); }
        }
      }
      selectedInverseDemandClosureEntries += closure.size();
    };

  // Under memoryPolicy = "speed" a fast path that does not fit the dense budget
  // stops the fit instead of silently using a slower route.
  auto requireMemory = [&](const std::string & what, const double neededDoubles) -> void {
    const double neededMB = neededDoubles * sizeof(double) / 1.0e6;
    const double budgetMB = blockSchurPolicy.maxDoubles * sizeof(double) / 1.0e6;
    Rcpp::stop(what + " needs about " + std::to_string(static_cast<long long>(std::ceil(neededMB))) +
      " MB but denseMemoryMB is " + std::to_string(static_cast<long long>(std::floor(budgetMB))) +
      ". Raise it, e.g. mmes_control(denseMemoryMB = " +
      std::to_string(static_cast<long long>(std::ceil(neededMB * 1.25))) +
      "), to use the faster path, or set memoryPolicy = \"limit\" to accept slower fallbacks.");
  };

  auto buildSelectedInverseSubset =
    [&](const EigenLDLT & factor,
        const std::string & context,
        SelectedInverseSubset & out,
        const bool reuseTopology) -> void {

      const Eigen::VectorXd D = factor.vectorD();
      const int n = static_cast<int>(D.size());
      EigenSpMat Lmat = factor.matrixL();
      Lmat.makeCompressed();

      // Flush demands while their subset address still refers to this pattern.
      if(!reuseTopology && phaseProfile.enabled()){
        std::lock_guard<std::mutex> lock(selectedInverseDemandMutex);
        const auto demand = selectedInverseDemandKeys.find(&out);
        if(demand != selectedInverseDemandKeys.end()){
          const std::size_t recordedKeys = demand->second.size();
          finalizeSelectedInverseDemand(out, demand->second);
          selectedInverseDemandKeyCount -= recordedKeys;
          selectedInverseDemandKeys.erase(demand);
        }
      }

      // When reuseTopology is set, `out` already holds the sparsity pattern
      // (rows/originalToPermuted) from a previous call with the same LDLT
      // fill-in pattern, so only the numeric values need to be reset. This
      // avoids a full deep copy of the (potentially large) pattern vectors
      // on every REML iteration.
      if(reuseTopology){
        if(static_cast<int>(out.rows.size()) != n ||
           static_cast<int>(out.originalToPermuted.size()) != n){
          Rcpp::stop("Cached selected-inverse topology is incompatible with the LDLT factor.");
        }
        if(out.dependencyEntries > selectedDependencyMaxEntries){
          out.dependencyOffsets.clear();
          out.dependencyOffsets.resize(n);
          out.dependencyEntries = 0;
        }
        for(int col = 0; col < n; ++col){
          std::fill(out.values[col].begin(), out.values[col].end(), 0.0);
        }
      }else{
        out.rows.assign(n, std::vector<int>());
        out.values.assign(n, std::vector<double>());
        out.rowLookup.clear();
        out.rowLookup.resize(n);
        out.contiguous.assign(n, 0);
        out.dependencyOffsets.assign(n, std::vector<std::size_t>());
        out.dependencyEntries = 0;
        out.originalToPermuted.assign(n, -1);

        const auto & perm = factor.permutationP();
        // Eigen's convention is P e_i = e_{sigma(i)}, therefore
        // perm.indices()(i) maps ORIGINAL index i -> PERMUTED index.
        for(int i = 0; i < n; ++i){
          const int permutedIndex = perm.indices()(i);
          if(permutedIndex < 0 || permutedIndex >= n){
            Rcpp::stop("Invalid permutation encountered while building selected inverse subset.");
          }
          out.originalToPermuted[i] = permutedIndex;
        }

        for(int col = 0; col < n; ++col){
          out.rows[col].push_back(col);
          for(EigenSpMat::InnerIterator it(Lmat, col); it; ++it){
            const int row = it.row();
            if(row > col){ out.rows[col].push_back(row); }
          }
          std::sort(out.rows[col].begin(), out.rows[col].end());
          out.rows[col].erase(std::unique(out.rows[col].begin(), out.rows[col].end()), out.rows[col].end());
          out.values[col].assign(out.rows[col].size(), 0.0);
          // A dense column maps row -> position by subtraction; no hash needed.
          out.contiguous[col] = static_cast<std::size_t>(out.rows[col].back() - col) + 1 == out.rows[col].size();
          if(!out.contiguous[col] && out.rows[col].size() >= sommer::WorkspacePolicy::indexedColumnRows()){
            std::unordered_map<int, std::size_t> & lookup = out.rowLookup[col];
            lookup.reserve(out.rows[col].size() * 2);
            for(std::size_t position = 0; position < out.rows[col].size(); ++position){
              lookup.emplace(out.rows[col][position], position);
            }
          }
        }

        for(int col = n - 1; col >= 0; --col){
          const std::size_t neighbours = out.rows[col].size() - 1;
          const std::size_t remaining = selectedDependencyMaxEntries - out.dependencyEntries;
          if(neighbours == 0 || neighbours > remaining / neighbours){ continue; }
          std::vector<std::size_t> & offsets = out.dependencyOffsets[col];
          offsets.resize(neighbours * neighbours);
          for(std::size_t target = 0; target < neighbours; ++target){
            for(std::size_t source = 0; source < neighbours; ++source){
              const int first = out.rows[col][target + 1];
              const int second = out.rows[col][source + 1];
              const int dependencyCol = std::min(first, second);
              const int dependencyRow = std::max(first, second);
              const std::vector<int> & rows = out.rows[dependencyCol];
              const auto found = std::lower_bound(rows.begin(), rows.end(), dependencyRow);
              if(found == rows.end() || *found != dependencyRow){
                Rcpp::stop("Selected-inverse dependency is outside the filled LDLT pattern.");
              }
              offsets[target * neighbours + source] = static_cast<std::size_t>(found - rows.begin());
            }
          }
          out.dependencyEntries += offsets.size();
        }
      }

      auto getPermuted = [&](int a, int b, double & value) -> bool {
        const int col = std::min(a,b);
        const int row = std::max(a,b);
        if(out.contiguous[col]){
          if(row > out.rows[col].back()){ return false; }
          value = out.values[col][static_cast<std::size_t>(row - col)];
          return true;
        }
        if(!out.rowLookup.empty() && !out.rowLookup[col].empty()){
          const auto pos = out.rowLookup[col].find(row);
          if(pos == out.rowLookup[col].end()){ return false; }
          value = out.values[col][pos->second];
          return true;
        }
        const std::vector<int> & rr = out.rows[col];
        auto pos = std::lower_bound(rr.begin(), rr.end(), row);
        if(pos == rr.end() || (*pos) != row){ return false; }
        const std::size_t idx = static_cast<std::size_t>(pos - rr.begin());
        value = out.values[col][idx];
        return true;
      };

      // Takahashi recursions for A = P' L D L' P.
      for(int i = n - 1; i >= 0; --i){
        if(!std::isfinite(D(i)) || D(i) <= 0.0){
          Rcpp::stop("Non-positive LDLT pivot while building selected inverse subset in " + context + ".");
        }

        std::vector<int> neighbours;
        std::vector<double> lvalues;
        for(EigenSpMat::InnerIterator it(Lmat, i); it; ++it){
          if(it.row() > i){
            neighbours.push_back(it.row());
            lvalues.push_back(it.value());
          }
        }
        if(out.rows[i].size() != neighbours.size() + 1){
          Rcpp::stop("Cached selected-inverse column is incompatible with the LDLT pattern.");
        }

        // Off-diagonal entries first.
        for(std::size_t jj = 0; jj < neighbours.size(); ++jj){
          const int j = neighbours[jj];
          if(out.rows[i][jj + 1] != j){
            Rcpp::stop("Cached selected-inverse ordering is incompatible with the LDLT pattern.");
          }
          double sum = 0.0;
          for(std::size_t kk = 0; kk < neighbours.size(); ++kk){
            double zkj = 0.0;
            if(!out.dependencyOffsets[i].empty()){
              zkj = out.values[std::min(neighbours[kk], j)][
                out.dependencyOffsets[i][jj * neighbours.size() + kk]];
            }else if(!getPermuted(neighbours[kk], j, zkj)){
              Rcpp::stop("Takahashi recurrence requested an entry outside the filled LDLT pattern in " + context + ".");
            }
            sum += lvalues[kk] * zkj;
          }
          out.values[i][jj + 1] = -sum;
        }

        double diagCorrection = 0.0;
        for(std::size_t kk = 0; kk < neighbours.size(); ++kk){
          diagCorrection += lvalues[kk] * out.values[i][kk + 1];
        }
        out.values[i][0] = (1.0 / D(i)) - diagCorrection;
      }
    };

  auto getSelectedInverseOriginal =
    [&](const SelectedInverseSubset & subset,
        int originalRow,
        int originalCol,
        double & value) -> bool {
      if(originalRow < 0 || originalCol < 0 ||
     originalRow >= static_cast<int>(subset.originalToPermuted.size()) ||
     originalCol >= static_cast<int>(subset.originalToPermuted.size())){
    if(phaseProfile.enabled()){
      std::lock_guard<std::mutex> lock(selectedInverseDemandMutex);
      ++selectedInverseLookupRequests;
      ++selectedInverseLookupMisses;
    }
    return false;
      }
      const int a = subset.originalToPermuted[originalRow];
      const int b = subset.originalToPermuted[originalCol];
      const int col = std::min(a,b);
      const int row = std::max(a,b);
      bool found = false;
      if(!subset.contiguous.empty() && subset.contiguous[col]){
        if(row <= subset.rows[col].back()){
          value = subset.values[col][static_cast<std::size_t>(row - col)];
          found = true;
        }
      }else if(!subset.rowLookup.empty() && !subset.rowLookup[col].empty()){
        const auto pos = subset.rowLookup[col].find(row);
        if(pos != subset.rowLookup[col].end()){
          value = subset.values[col][pos->second];
          found = true;
        }
      }else{
        const std::vector<int> & rr = subset.rows[col];
        auto pos = std::lower_bound(rr.begin(), rr.end(), row);
        if(pos != rr.end() && (*pos) == row){
          const std::size_t idx = static_cast<std::size_t>(pos - rr.begin());
          value = subset.values[col][idx];
          found = true;
        }
      }
      if(phaseProfile.enabled()){
        std::lock_guard<std::mutex> lock(selectedInverseDemandMutex);
        ++selectedInverseLookupRequests;
        if(found){
          ++selectedInverseLookupHits;
          if(selectedInverseDemandTrackingComplete){
            if(selectedInverseDemandKeyCount < maxSelectedInverseDemandKeys){
              const std::uint64_t key =
                (static_cast<std::uint64_t>(static_cast<std::uint32_t>(row)) << 32) |
                static_cast<std::uint32_t>(col);
              selectedInverseDemandKeys[&subset].push_back(key);
              ++selectedInverseDemandKeyCount;
            }else{
              selectedInverseDemandKeys.clear();
              selectedInverseDemandKeyCount = 0;
              selectedInverseDemandTrackingComplete = false;
            }
          }
        }else{
          ++selectedInverseLookupMisses;
        }
      }
      return found;
    };

  // ============================================================
  // CHANGE 5: reusable symbolic factorisation for C
  // ============================================================
  // The numerical values of C change with the variance components,
  // but its sparsity pattern is usually unchanged.  Cache the sparse
  // pattern and reuse Eigen's symbolic analysis whenever possible.
  using sommer::eigenSparsePatternMatches;
  using sommer::cacheEigenSparsePattern;

  // ============================================================
  // MME linear-solver module
  //
  // Cfactor is the exact direct backend and also supplies the LDLT factors
  // needed by log|C| and Takahashi.  PCG is a factorisation-free backend
  // for the repeated C x = rhs solves (BLUE/BLUP, sensitivities and rare
  // selected-block fallbacks).  Keeping all solve dispatch here makes future
  // solver backends local to this module rather than the AI-REML code.
  // ============================================================
  struct LDLTCache {
    EigenLDLT Cfactor;
    EigenLDLT Dfactor;
    bool CsymbolicReady = false;
    bool DsymbolicReady = false;
    bool CselectedTopologyReady = false;
    bool DselectedTopologyReady = false;
    std::vector<int> CouterPattern;
    std::vector<int> CinnerPattern;
    std::vector<int> DouterPattern;
    std::vector<int> DinnerPattern;
    SelectedInverseSubset CselectedTopology;
    SelectedInverseSubset DselectedTopology;
  };
  const bool retainLDLTCache = solverName == "ldlt" && covStructI.hasAttribute("retainLDLTCache") &&
    Rcpp::as<bool>(covStructI.attr("retainLDLTCache"));
  SEXP cacheInput = covStructI.hasAttribute("ldltCache") ? static_cast<SEXP>(covStructI.attr("ldltCache")) : R_NilValue;
  if(cacheInput != R_NilValue &&
     (TYPEOF(cacheInput) != EXTPTRSXP || R_ExternalPtrTag(cacheInput) != Rf_install("sommer_ldlt_cache") ||
      R_ExternalPtrAddr(cacheInput) == nullptr)){
    cacheInput = R_NilValue;
  }
  Rcpp::XPtr<LDLTCache> ldltCache = retainLDLTCache && cacheInput != R_NilValue
    ? Rcpp::XPtr<LDLTCache>(cacheInput)
    : Rcpp::XPtr<LDLTCache>(new LDLTCache(), true, Rf_install("sommer_ldlt_cache"));
  EigenLDLT & Cfactor = ldltCache->Cfactor;
  bool & CsymbolicReady = ldltCache->CsymbolicReady;
  int CsymbolicAnalyses = 0;
  bool CnumericReady = false;
  bool cholmodCFactorCurrent = false;
  bool cholmodFinalSelectedInverse = false;
  bool CpcgReady = false;
  bool matrixFreePCGReady = false;
  std::function<Eigen::VectorXd(const Eigen::VectorXd &)> matrixFreeApplyC;
  std::function<Eigen::MatrixXd(const Eigen::MatrixXd &)> matrixFreeApplyCMatrix;
  Eigen::VectorXd matrixFreePCGDiagonal;
  std::vector<int> & CouterPattern = ldltCache->CouterPattern;
  std::vector<int> & CinnerPattern = ldltCache->CinnerPattern;
  SelectedInverseSubset & CselectedTopology = ldltCache->CselectedTopology;
  bool & CselectedTopologyReady = ldltCache->CselectedTopologyReady;

  // ML (reml=FALSE) support: D is the random-effects-only block of C,
  // D = Z'R^{-1}Z + G^{-1} (rows/cols nX..nEffects-1 of C, re-indexed to
  // 0..Nu-1). log|D| and D's own Takahashi selected inverse replace
  // log|C| and CselectedTopology in the likelihood/score computations,
  // dropping the log|X'V^{-1}X| restriction term that makes REML
  // "restricted". Only used when reml==false (validated solver=="ldlt"
  // above); never touched otherwise.
  EigenLDLT & Dfactor = ldltCache->Dfactor;
  bool & DsymbolicReady = ldltCache->DsymbolicReady;
  std::vector<int> & DldltOuterPattern = ldltCache->DouterPattern;
  std::vector<int> & DldltInnerPattern = ldltCache->DinnerPattern;
  int DsymbolicAnalyses = 0;
  int DselectedTopologyBuilds = 0;
  int DselectedTopologyReuses = 0;
  SelectedInverseSubset & DselectedTopology = ldltCache->DselectedTopology;
  bool & DselectedTopologyReady = ldltCache->DselectedTopologyReady;

  std::string coeffPathDescription = "Unknown";
  std::string residPathDescription = "Unknown";
  std::string assemPathDescription = "Unknown";
  std::string tracePathDescription = "Unknown";
  std::string finalCiPathDescription = computeCi > 0 ? "Unknown" : "None (computeCi=0)";


  // ============================================================
  // CHOLMOD (via R's Matrix package) supernodal backend for solver="cholmod"
  // ============================================================
  using sommer::CholmodState;
  const bool retainCholmodCache = solverName == "cholmod" &&
    covStructI.hasAttribute("retainCholmodCache") &&
    Rcpp::as<bool>(covStructI.attr("retainCholmodCache"));
  SEXP cholmodCacheInput = covStructI.hasAttribute("cholmodCache")
    ? static_cast<SEXP>(covStructI.attr("cholmodCache")) : R_NilValue;
  if(cholmodCacheInput != R_NilValue &&
     (TYPEOF(cholmodCacheInput) != EXTPTRSXP ||
      R_ExternalPtrTag(cholmodCacheInput) != Rf_install("sommer_cholmod_fit_cache") ||
      R_ExternalPtrAddr(cholmodCacheInput) == nullptr)){
    cholmodCacheInput = R_NilValue;
  }
  Rcpp::XPtr<sommer::CholmodFitCache> cholmodCache = retainCholmodCache && cholmodCacheInput != R_NilValue
    ? Rcpp::XPtr<sommer::CholmodFitCache>(cholmodCacheInput)
    : Rcpp::XPtr<sommer::CholmodFitCache>(new sommer::CholmodFitCache(), true,
        Rf_install("sommer_cholmod_fit_cache"));
  CholmodState & cholmodState = cholmodCache->joint;

  const double cholmodInvCacheMaxDoubles = blockSchurPolicy.maxDoubles;
  selectedDependencyMaxEntries = sommer::WorkspacePolicy::cacheCapacity(
    blockSchurPolicy.maxDoubles, sizeof(std::size_t));
  factorEngine.contractionDigitsMaxEntries = sommer::WorkspacePolicy::cacheCapacity(
    blockSchurPolicy.maxDoubles, sizeof(arma::uword));
  factorEngine.factorEvaluationMaxEntries = sommer::WorkspacePolicy::cacheCapacity(
    blockSchurPolicy.maxDoubles, sizeof(double));
  if(factorEngine.factorEvaluationEntries > factorEngine.factorEvaluationMaxEntries){
    factorEngine.factorEvaluationCache.clear();
    factorEngine.factorEvaluationEntries = 0;
  }

  std::vector<sommer::TensorComponentPlan> randomTensorPlans(
    static_cast<std::size_t>(nRe));
  for(int term = 0; term < nRe && nZs > 0; ++term){
    const Rcpp::List descriptor = covDescriptor[static_cast<std::size_t>(term)];
    const Rcpp::List factors = descriptor["factors"];
    std::vector<std::size_t> dimensions;
    std::vector<std::vector<int>> componentLabels;
    for(int axis = 0; axis < factors.size(); ++axis){
      const Rcpp::List factor = Rcpp::as<Rcpp::List>(factors[axis]);
      const int dimension = Rcpp::as<int>(factor["dim"]);
      if(dimension <= 0){ Rcpp::stop("Covariance factor dimensions must be positive."); }
      dimensions.push_back(static_cast<std::size_t>(dimension));
      std::vector<int> labels(static_cast<std::size_t>(dimension), 0);
      const bool structurallyDiagonal =
        factor.containsElementNamed("structurally_diagonal") &&
        Rcpp::as<bool>(factor["structurally_diagonal"]);
      if(structurallyDiagonal){
        for(int vertex = 0; vertex < dimension; ++vertex){ labels[vertex] = vertex; }
      }else if(factor.containsElementNamed("evaluator")){
        const Rcpp::List evaluator = factor["evaluator"];
        if(Rcpp::as<std::string>(evaluator["backend"]) == "fixed"){
          const arma::mat covariance = Rcpp::as<arma::mat>(evaluator["matrix"]);
          if(covariance.n_rows != static_cast<arma::uword>(dimension) ||
             covariance.n_cols != static_cast<arma::uword>(dimension)){
            Rcpp::stop("Fixed covariance factor dimensions are inconsistent.");
          }
          sommer::DisjointComponents graph(static_cast<std::size_t>(dimension));
          for(int column = 0; column < dimension; ++column){
            for(int row = 0; row < dimension; ++row){
              if(covariance(row, column) != 0.0){ graph.unite(row, column); }
            }
          }
          labels = graph.labels();
        }
      }
      componentLabels.push_back(std::move(labels));
    }
    const std::size_t mainSize = Ai(term).n_rows;
    if(mainSize == 0 || partitions(term).n_rows > std::numeric_limits<std::size_t>::max() / mainSize){
      Rcpp::stop("Random-effect tensor dimension overflow.");
    }
    sommer::DisjointComponents relationshipGraph(mainSize);
    for(arma::sp_mat::const_iterator entry = Ai(term).begin(); entry != Ai(term).end(); ++entry){
      if(*entry != 0.0){ relationshipGraph.unite(static_cast<int>(entry.row()), static_cast<int>(entry.col())); }
    }
    dimensions.push_back(mainSize);
    componentLabels.push_back(relationshipGraph.labels());
    const std::size_t coefficientCount = partitions(term).n_rows * mainSize;
    sommer::TensorComponentPlan & plan = randomTensorPlans[static_cast<std::size_t>(term)];
    plan = sommer::planTensorGraphComponents(dimensions, componentLabels, coefficientCount);
    if(!plan.eligible || plan.coordinates != coefficientCount){
      Rcpp::stop("Covariance factor dimensions do not match random-effect partitions.");
    }
  }

  auto buildRandomInverseGroups = [&]() -> std::vector< std::vector<int> > {
    std::vector< std::vector<int> > out;
    for(int iR = 0; iR < nRe && nZs > 0; ++iR){
      const arma::mat partitionsP = partitions(iR);

      const sommer::TensorComponentPlan & plan = randomTensorPlans[static_cast<std::size_t>(iR)];
      const int mainSize = static_cast<int>(Ai(iR).n_rows);
      for(const std::vector<int> & coordinates : plan.groups){
        out.emplace_back();
        for(const int coordinate : coordinates){
          const int partition = coordinate / mainSize;
          const int level = coordinate % mainSize;
          const int first = static_cast<int>(partitionsP(partition, 0)) - 1;
          out.back().push_back(first + level);
        }
      }
    }
    return out;
  };

  // Exact block-Schur engine for C = [Cbb Cbu; Cub blockdiag(D_g)], where the
  // border b holds the fixed effects plus any random groups that couple other
  // groups: log|C| = sum log|D_g| + log|S|, S = Cbb - sum Cbg D_g^{-1} Cgb.
  struct BlockSchurEngine {
    bool active = false;
    bool chain = false;
    bool banded = false;
    std::size_t bandWidth = 0;
    bool latent = false;
    int latentTerm = -1;
    std::vector<int> lowRankDimensions;
    std::size_t originalBorderSize = 0;
    arma::mat latentBasePrecision;
    double latentBaseLogDet = 0.0;
    bool patternChecked = false;
    bool eligible = false;
    std::vector<int> outerPattern;
    std::vector<int> innerPattern;
    std::vector< std::vector<int> > groups;
    std::vector<int> borderCols;
    std::vector<int> borderPos;
    std::vector<int> groupOf;
    std::vector<int> localOf;
    std::uint64_t groupSignature = 0;
    std::vector<arma::mat> Lg;
    std::vector<arma::mat> Fg;
    std::vector<arma::mat> Cgx;
    std::vector<arma::mat> chainTransfers;
    std::vector<std::vector<arma::mat>> bandLower;
    std::vector<int> chainComponents;
    double chainStorageDoubles = 0.0;
    std::size_t sparseFactorNonzeros = 0;
    double sparseFactorFlops = 0.0;
    double schurCostFlops = 0.0;
    double latentWorkFlops = std::numeric_limits<double>::quiet_NaN();
    double marginalWorkFlops = std::numeric_limits<double>::quiet_NaN();
    std::string latentAdmission = "not-applicable";
    arma::mat Ls;
    double logDet = 0.0;
  };
  BlockSchurEngine blockEngine;
  std::vector<DescriptorPrecisionState> randomPrecisionState(static_cast<std::size_t>(nRe));
  bool factorSchurResidualEligible = false;
  int factorSchurFallbacks = 0;
  int factorSchurNumericalFallbacks = 0;
  int traceFallbackSolves = 0;
  int emInverseTraceCalls = 0, emInverseTraceFallbackPanels = 0, emDesignTraceCalls = 0;
  int diagonalResidualTraceCalls = 0, diagonalResidualTraceSolveColumns = 0;
  int traceBatchSolves = 0;
  double traceBatchWorkspaceBytes = 0.0;
  int gridFactorizations = 0;
  int gridFactorReuses = 0;
  int pcgResidualTraceContractions = 0;
  auto crossGroupSignature = [](const std::vector< std::vector<int> > & groups) -> std::uint64_t {
    std::uint64_t hash = 1469598103934665603ULL;
    for(const std::vector<int> & cols : groups){
      hash ^= static_cast<std::uint64_t>(cols.size());
      hash *= 1099511628211ULL;
      for(const int col : cols){
        hash ^= static_cast<std::uint64_t>(static_cast<unsigned int>(col));
        hash *= 1099511628211ULL;
      }
    }
    return hash;
  };

  auto blockEngineCheckPattern = [&](const EigenSpMat & A) -> bool {
    std::vector<int> lowRankDimensions;
    lowRankDimensions.reserve(static_cast<std::size_t>(2 * nRe));
    for(const DescriptorPrecisionState & state : randomPrecisionState){
      lowRankDimensions.push_back(static_cast<int>(state.lowRankSpecific.n_elem));
      lowRankDimensions.push_back(static_cast<int>(state.lowRankLoadings.n_cols));
    }
    bool cachedLatentWellScaled = true;
    if(blockEngine.patternChecked && blockEngine.latent &&
       blockEngine.latentAdmission == "work-estimate"){
      const DescriptorPrecisionState & state =
        randomPrecisionState[static_cast<std::size_t>(blockEngine.latentTerm)];
      cachedLatentWellScaled = sommer::LatentSchurPolicy::specificsWellScaled(
        arma::accu(arma::square(state.lowRankLoadings)), state.lowRankSpecific.min());
    }
    if(blockEngine.patternChecked && cachedLatentWellScaled &&
       eigenSparsePatternMatches(A, blockEngine.outerPattern, blockEngine.innerPattern) &&
       lowRankDimensions == blockEngine.lowRankDimensions){
      return blockEngine.eligible;
    }
    cacheEigenSparsePattern(A, blockEngine.outerPattern, blockEngine.innerPattern);
    blockEngine.patternChecked = true;
    blockEngine.lowRankDimensions = std::move(lowRankDimensions);
    blockEngine.eligible = false;
    blockEngine.chain = false;
    blockEngine.banded = false;
    blockEngine.latent = false;
    blockEngine.latentWorkFlops = std::numeric_limits<double>::quiet_NaN();
    blockEngine.marginalWorkFlops = std::numeric_limits<double>::quiet_NaN();
    blockEngine.latentAdmission = "not-applicable";
    if(blockSchurMode == "off"){ return false; }
    blockEngine.sparseFactorNonzeros = 0;
    blockEngine.sparseFactorFlops = std::numeric_limits<double>::quiet_NaN();
    // Predicted sparse-factor fill and work (AMD order), the basis for pricing sparse versus Schur.
    // Computed only when a dense plan is priced: latent/chain/banded routes do not need it, and
    // AMD plus symbolic counts on a dense marginal C are expensive.
    auto estimateSparseFactor = [&]() {
      if(!A.isCompressed()){ return; }
      Eigen::AMDOrdering<int> amd;
      Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic, int> amdOrder;
      amd(A, amdOrder);
      const std::size_t sizeA = static_cast<std::size_t>(A.rows());
      const std::vector<int> outerA(A.outerIndexPtr(), A.outerIndexPtr() + sizeA + 1);
      const std::vector<int> innerA(A.innerIndexPtr(), A.innerIndexPtr() + A.nonZeros());
      const std::vector<int> orderA(amdOrder.indices().data(), amdOrder.indices().data() + sizeA);
      const sommer::SymbolicFactor symbolic = sommer::symbolicCholeskyCounts(
        sizeA, outerA, innerA, orderA, static_cast<std::size_t>(blockSchurPolicy.maxDoubles));
      blockEngine.sparseFactorNonzeros = symbolic.nonzeros;
      blockEngine.sparseFactorFlops = symbolic.valid ? symbolic.flops : std::numeric_limits<double>::quiet_NaN();
    };
    if(reml && factorSchurResidualEligible){
      for(int term = 0; term < nRe && nZs > 0; ++term){
        const DescriptorPrecisionState & state = randomPrecisionState[static_cast<std::size_t>(term)];
        if(!state.factorWise || state.lowRankSpecific.is_empty() || state.lowRankFactorIndex < 0){ continue; }
        if(state.factorPrecision.size() > 1 && factorSchurMultiFactorMode == "off"){ continue; }
        const arma::mat termPartitions = partitions(term);
        const std::size_t mainSize = Ai(term).n_rows;
        const std::size_t rank = state.lowRankLoadings.n_cols;
        const std::size_t lowRankAxis = static_cast<std::size_t>(state.lowRankFactorIndex);
          if(rank == 0 || lowRankAxis >= state.factorPrecision.size() ||
           state.lowRankSpecific.n_elem != factorEngine.factorPrecisionDimension(state, lowRankAxis) ||
            rank >= state.lowRankSpecific.n_elem ||
            factorEngine.factorPrecisionDimension(state, lowRankAxis) > termPartitions.n_rows ||
           termPartitions.n_rows % state.lowRankSpecific.n_elem != 0 ||
           !blockSchurPolicy.denseEnough(static_cast<double>(Ai(term).n_nonzero),
                                         static_cast<double>(mainSize))){ continue; }
        std::size_t trailing = 1;
        for(std::size_t axis = lowRankAxis + 1; axis < state.factorPrecision.size(); ++axis){
           const arma::uword axisDimension = factorEngine.factorPrecisionDimension(state, axis);
          if(trailing > std::numeric_limits<std::size_t>::max() /
             axisDimension){ trailing = 0; break; }
           trailing *= axisDimension;
        }
        if(trailing == 0){ continue; }
        const std::size_t lowRankDimension = state.lowRankSpecific.n_elem;
        const std::size_t baseFactorCoordinates = termPartitions.n_rows / lowRankDimension;
        if(mainSize > std::numeric_limits<std::size_t>::max() / baseFactorCoordinates){ continue; }
        const std::size_t baseSize = mainSize * baseFactorCoordinates;
        if(baseSize == 0 || rank > std::numeric_limits<std::size_t>::max() / baseSize){ continue; }
        std::vector<std::vector<int>> groups(lowRankDimension);
        std::vector<int> groupOf(static_cast<std::size_t>(nEffects), -1);
        std::vector<int> rowGroup(static_cast<std::size_t>(W.rows()), -1);
        bool disjoint = true;
        for(arma::uword coordinate = 0; coordinate < termPartitions.n_rows && disjoint; ++coordinate){
          const std::size_t group = (coordinate / trailing) % lowRankDimension;
          const int start = static_cast<int>(termPartitions(coordinate, 0)) - 1;
          const int end = static_cast<int>(termPartitions(coordinate, 1)) - 1;
          if(static_cast<std::size_t>(end - start + 1) != mainSize){ disjoint = false; break; }
          for(int column = start; column <= end && disjoint; ++column){
            groupOf[column] = static_cast<int>(group);
            groups[group].push_back(column);
            for(EigenSpMat::InnerIterator entry(W, column); entry; ++entry){
              if(entry.value() == 0.0){ continue; }
              int & owner = rowGroup[entry.row()];
              if(owner >= 0 && owner != static_cast<int>(group)){ disjoint = false; break; }
              owner = static_cast<int>(group);
            }
          }
        }
        if(!disjoint){ continue; }
        const std::size_t originalBorder = static_cast<std::size_t>(nEffects) -
          termPartitions.n_rows * mainSize;
        const std::size_t border = originalBorder + rank * baseSize;
        const double storage = blockSchurPolicy.latentStorageDoubles(groups, border, baseSize);
        if(storage > blockSchurPolicy.maxDoubles){ continue; }
        blockEngine.latentWorkFlops = std::numeric_limits<double>::quiet_NaN();
        blockEngine.marginalWorkFlops = std::numeric_limits<double>::quiet_NaN();
        if(state.factorPrecision.size() > 1 && factorSchurMultiFactorMode == "auto"){
          if(sommer::LatentSchurPolicy::isAdvantageousMultiFactor(lowRankDimension, border)){
            blockEngine.latentAdmission = "small-border";
          }else{
            // Price the marginal alternative's actual connected components, not
            // one dense tensor: diagonal axes may leave much cheaper blocks.
            const auto & priorGroups = randomTensorPlans[static_cast<std::size_t>(term)].groups;
            sommer::DisjointComponents components(priorGroups.size());
            std::vector<int> marginalGroupOf(static_cast<std::size_t>(nEffects), -1);
            for(std::size_t group = 0; group < priorGroups.size(); ++group){
              for(const int coordinate : priorGroups[group]){
                const int partition = coordinate / static_cast<int>(mainSize);
                const int level = coordinate % static_cast<int>(mainSize);
                const int column = static_cast<int>(termPartitions(partition, 0)) - 1 + level;
                marginalGroupOf[column] = static_cast<int>(group);
              }
            }
            for(int column = 0; column < nEffects; ++column){
              const int group = marginalGroupOf[column];
              if(group < 0){ continue; }
              for(EigenSpMat::InnerIterator entry(A, column); entry; ++entry){
                const int other = marginalGroupOf[entry.row()];
                if(other >= 0 && other != group){ components.unite(group, other); }
              }
            }
            const std::vector<int> labels = components.labels();
            std::vector<std::size_t> widths;
            for(std::size_t group = 0; group < priorGroups.size(); ++group){
              const std::size_t component = static_cast<std::size_t>(labels[group]);
              if(widths.size() <= component){ widths.resize(component + 1, 0); }
              widths[component] += priorGroups[group].size();
            }
            blockEngine.latentWorkFlops = sommer::LatentSchurPolicy::latentWork(
              lowRankDimension, baseSize, border);
            blockEngine.marginalWorkFlops = sommer::LatentSchurPolicy::marginalWork(
              widths, originalBorder);
            if(!sommer::LatentSchurPolicy::workFavoursLatent(
                blockEngine.latentWorkFlops, blockEngine.marginalWorkFlops)){
              blockEngine.latentAdmission = "work-rejected";
              continue;
            }
            if(!sommer::LatentSchurPolicy::specificsWellScaled(
                arma::accu(arma::square(state.lowRankLoadings)), state.lowRankSpecific.min())){
              blockEngine.latentAdmission = "work-numerical-rejected";
              ++factorSchurFallbacks;
              ++factorSchurNumericalFallbacks;
              continue;
            }
            blockEngine.latentAdmission = "work-estimate";
          }
        }else{
          blockEngine.latentAdmission = state.factorPrecision.size() > 1 ? "forced" : "single-factor";
        }

        blockEngine.groups = std::move(groups);
        blockEngine.groupOf = std::move(groupOf);
        blockEngine.borderCols.clear();
        blockEngine.borderPos.assign(static_cast<std::size_t>(nEffects), -1);
        blockEngine.localOf.assign(static_cast<std::size_t>(nEffects), -1);
        for(int column = 0; column < nEffects; ++column){
          if(blockEngine.groupOf[column] < 0){
            blockEngine.borderPos[column] = static_cast<int>(blockEngine.borderCols.size());
            blockEngine.borderCols.push_back(column);
          }
        }
        for(std::size_t column = 0; column < rank * baseSize; ++column){
          blockEngine.borderCols.push_back(nEffects + static_cast<int>(column));
        }
        for(const auto & group : blockEngine.groups){
          for(std::size_t local = 0; local < group.size(); ++local){ blockEngine.localOf[group[local]] = static_cast<int>(local); }
        }
        blockEngine.latentTerm = term;
        blockEngine.originalBorderSize = originalBorder;
        blockEngine.latentBasePrecision.reset();
        blockEngine.latentBaseLogDet = 0.0;
        blockEngine.chainStorageDoubles = storage;
        blockEngine.latent = true;
        blockEngine.eligible = true;
        return true;
      }
    }
    if(reml){
      for(int term = 0; term < nRe && nZs > 0; ++term){
        if(randomStructurallyDiagonal[static_cast<std::size_t>(term)]){ continue; }
        const arma::mat termPartitions = partitions(term);
        if(termPartitions.n_rows < 3){ continue; }
        std::vector<std::vector<int>> blocks(termPartitions.n_rows);
        std::vector<int> blockOf(static_cast<std::size_t>(nEffects), -1);
        std::size_t interiorSize = 0;
        for(arma::uword block = 0; block < termPartitions.n_rows; ++block){
          const int start = static_cast<int>(termPartitions(block, 0)) - 1;
          const int end = static_cast<int>(termPartitions(block, 1)) - 1;
          for(int column = start; column <= end; ++column){
            blocks[block].push_back(column);
            blockOf[column] = static_cast<int>(block);
            ++interiorSize;
          }
        }
        const std::size_t borderSize = static_cast<std::size_t>(nEffects) - interiorSize;
        if(borderSize > blockSchurPolicy.maxBorder){ continue; }
        std::vector<std::pair<int, int>> edges;
        std::vector<double> within(blocks.size(), 0.0);
        for(std::size_t block = 0; block < blocks.size(); ++block){
          for(const int column : blocks[block]){
            for(EigenSpMat::InnerIterator entry(A, column); entry; ++entry){
              const int other = blockOf[entry.row()];
              if(other == static_cast<int>(block)){ ++within[block]; }
              else if(other > static_cast<int>(block)){
                edges.emplace_back(static_cast<int>(block), other);
              }
            }
          }
        }
        std::sort(edges.begin(), edges.end());
        edges.erase(std::unique(edges.begin(), edges.end()), edges.end());
        const sommer::BlockChainOrder order = sommer::orderBlockPaths(blocks.size(), edges);
        if(!order.eligible){ continue; }
        std::vector<std::vector<int>> orderedBlocks;
        orderedBlocks.reserve(blocks.size());
        bool dense = true;
        for(std::size_t position = 0; position < order.blocks.size(); ++position){
          const int original = order.blocks[position];
          const double size = static_cast<double>(blocks[original].size());
          dense = dense && blockSchurPolicy.denseEnough(within[original], size);
          orderedBlocks.push_back(std::move(blocks[original]));
        }
        const double storage = blockSchurPolicy.chainStorageDoubles(orderedBlocks, order.components, borderSize);
        if(!dense || storage > blockSchurPolicy.maxDoubles){ continue; }
        blockEngine.groups = std::move(orderedBlocks);
        blockEngine.borderCols.clear();
        blockEngine.borderPos.assign(static_cast<std::size_t>(nEffects), -1);
        blockEngine.groupOf.assign(static_cast<std::size_t>(nEffects), -1);
        blockEngine.localOf.assign(static_cast<std::size_t>(nEffects), -1);
        for(int column = 0; column < nEffects; ++column){
          if(blockOf[column] < 0){
            blockEngine.borderPos[column] = static_cast<int>(blockEngine.borderCols.size());
            blockEngine.borderCols.push_back(column);
          }
        }
        for(std::size_t block = 0; block < blockEngine.groups.size(); ++block){
          for(std::size_t local = 0; local < blockEngine.groups[block].size(); ++local){
            const int column = blockEngine.groups[block][local];
            blockEngine.groupOf[column] = static_cast<int>(block);
            blockEngine.localOf[column] = static_cast<int>(local);
          }
        }
        blockEngine.chainComponents = order.components;
        blockEngine.chainStorageDoubles = storage;
        blockEngine.chain = true;
        blockEngine.eligible = true;
        return true;
      }
    }
    // Never auto-dispatched: measured up to 80x slower than generic CHOLMOD (score traces).
    if(reml && enableBlockBandedEngine){
      const std::size_t maxBlockBandwidth = sommer::BlockBandPolicy().maxBandwidth;
      for(int term = 0; term < nRe && nZs > 0; ++term){
        if(randomStructurallyDiagonal[static_cast<std::size_t>(term)]){ continue; }
        const arma::mat termPartitions = partitions(term);
        if(termPartitions.n_rows < 3){ continue; }
        std::vector<std::vector<int>> blocks(termPartitions.n_rows);
        std::vector<int> blockOf(static_cast<std::size_t>(nEffects), -1);
        std::size_t interiorSize = 0;
        for(arma::uword block = 0; block < termPartitions.n_rows; ++block){
          const int first = static_cast<int>(termPartitions(block, 0)) - 1;
          const int last = static_cast<int>(termPartitions(block, 1)) - 1;
          for(int column = first; column <= last; ++column){
            blocks[block].push_back(column);
            blockOf[column] = static_cast<int>(block);
            ++interiorSize;
          }
        }
        const std::size_t borderSize = static_cast<std::size_t>(nEffects) - interiorSize;
        if(borderSize > blockSchurPolicy.maxBorder){ continue; }
        std::vector<std::pair<int, int>> edges;
        std::vector<double> diagonalEntries(blocks.size(), 0.0);
        for(std::size_t block = 0; block < blocks.size(); ++block){
          for(const int column : blocks[block]){
            for(EigenSpMat::InnerIterator entry(A, column); entry; ++entry){
              const int other = blockOf[entry.row()];
              if(other == static_cast<int>(block)){ ++diagonalEntries[block]; }
              else if(other >= 0 && other != static_cast<int>(block)){
                edges.emplace_back(static_cast<int>(block), other);
              }
            }
          }
        }
        std::sort(edges.begin(), edges.end());
        edges.erase(std::unique(edges.begin(), edges.end()), edges.end());
        const sommer::BlockBandPlan plan =
          sommer::planBlockBand(blocks.size(), edges, maxBlockBandwidth);
        if(!plan.eligible){ continue; }
        bool denseEnough = true;
        for(std::size_t block = 0; block < blocks.size(); ++block){
          const double size = static_cast<double>(blocks[block].size());
          if(blockSchurPolicy.needsDensityCheck(static_cast<double>(nEffects), size)){
            denseEnough = denseEnough && blockSchurPolicy.denseEnough(diagonalEntries[block], size);
          }
        }
        const double storage = blockSchurPolicy.bandStorageDoubles(
          blocks, plan.bandwidth, borderSize);
        if(!denseEnough || storage > blockSchurPolicy.maxDoubles){ continue; }

        blockEngine.groups = std::move(blocks);
        blockEngine.borderCols.clear();
        blockEngine.borderPos.assign(static_cast<std::size_t>(nEffects), -1);
        blockEngine.groupOf.assign(static_cast<std::size_t>(nEffects), -1);
        blockEngine.localOf.assign(static_cast<std::size_t>(nEffects), -1);
        for(int column = 0; column < nEffects; ++column){
          if(blockOf[column] < 0){
            blockEngine.borderPos[column] = static_cast<int>(blockEngine.borderCols.size());
            blockEngine.borderCols.push_back(column);
          }
        }
        for(std::size_t block = 0; block < blockEngine.groups.size(); ++block){
          for(std::size_t local = 0; local < blockEngine.groups[block].size(); ++local){
            const int column = blockEngine.groups[block][local];
            blockEngine.groupOf[column] = static_cast<int>(block);
            blockEngine.localOf[column] = static_cast<int>(local);
          }
        }
        blockEngine.bandWidth = plan.bandwidth;
        blockEngine.chainStorageDoubles = storage;
        blockEngine.banded = true;
        blockEngine.eligible = true;
        return true;
      }
    }
    std::vector< std::vector<int> > candidates = buildRandomInverseGroups();
    if(Nu == 0 || candidates.empty()){ return false; }

    std::vector<int> candidateOf(static_cast<std::size_t>(nEffects), -1);
    for(std::size_t g = 0; g < candidates.size(); ++g){
      for(const int c : candidates[g]){
        if(c < nX || c >= nEffects || candidateOf[static_cast<std::size_t>(c)] != -1){ return false; }
        candidateOf[static_cast<std::size_t>(c)] = static_cast<int>(g);
      }
    }
    for(int c = nX; c < nEffects; ++c){
      if(candidateOf[static_cast<std::size_t>(c)] < 0){ return false; }
    }

    std::vector< std::pair<int,int> > couplings;
    for(int j = nX; j < nEffects; ++j){
      const int g = candidateOf[static_cast<std::size_t>(j)];
      for(EigenSpMat::InnerIterator it(A, j); it; ++it){
        const int i = static_cast<int>(it.row());
        if(i < nX){ continue; }
        const int h = candidateOf[static_cast<std::size_t>(i)];
        if(g < h){ couplings.emplace_back(g, h); }
      }
    }
    std::sort(couplings.begin(), couplings.end());
    couplings.erase(std::unique(couplings.begin(), couplings.end()), couplings.end());

    bool sparseFactorEstimated = false;
    auto densePlanAdmissible = [&](const sommer::BlockSchurPlan & candidate) -> bool {
      if(blockSchurMode == "force"){ return true; }
      if(!sparseFactorEstimated){ estimateSparseFactor(); sparseFactorEstimated = true; }
      return blockSchurPolicy.fillFavoursSchur(
        blockEngine.sparseFactorFlops,
        blockSchurPolicy.cost(candidate.groups, candidate.border.size()));
    };
    sommer::BlockSchurPlan plan = sommer::planBlockSchur(
      candidates, couplings, nX, blockSchurPolicy, densePlanAdmissible);
    blockEngine.schurCostFlops = blockSchurPolicy.cost(plan.groups, plan.border.size());
    if(!plan.eligible){ return false; }
    blockEngine.borderCols.swap(plan.border);
    blockEngine.groups.swap(plan.groups);
    const std::size_t nB = blockEngine.borderCols.size();
    blockEngine.groupSignature = crossGroupSignature(blockEngine.groups);

    blockEngine.borderPos.assign(static_cast<std::size_t>(nEffects), -1);
    blockEngine.groupOf.assign(static_cast<std::size_t>(nEffects), -1);
    blockEngine.localOf.assign(static_cast<std::size_t>(nEffects), -1);
    for(std::size_t b = 0; b < nB; ++b){
      blockEngine.borderPos[static_cast<std::size_t>(blockEngine.borderCols[b])] = static_cast<int>(b);
    }
    for(std::size_t g = 0; g < blockEngine.groups.size(); ++g){
      const std::vector<int> & cols = blockEngine.groups[g];
      for(std::size_t j = 0; j < cols.size(); ++j){
        blockEngine.groupOf[static_cast<std::size_t>(cols[j])] = static_cast<int>(g);
        blockEngine.localOf[static_cast<std::size_t>(cols[j])] = static_cast<int>(j);
      }
    }
    blockEngine.eligible = true;
    return true;
  };

  auto blockEngineSolveD = [&](const std::size_t g, const arma::mat & rhs) -> arma::mat {
    const arma::mat & L = blockEngine.Lg[g];
    arma::mat tmp;
    arma::mat out;
    const bool ok1 = gpuBackend.triangularSolve(tmp, arma::trimatl(L), rhs, arma::solve_opts::fast);
    const bool ok2 = ok1 && gpuBackend.triangularSolve(out, arma::trimatu(L.t()), tmp, arma::solve_opts::fast);
    if(!ok2){ Rcpp::stop("Dense block-Schur triangular solve failed."); }
    return out;
  };

  auto blockChainSolveInterior = [&](std::vector<arma::mat> rhs) -> std::vector<arma::mat> {
    for(std::size_t block = 1; block < rhs.size(); ++block){
      if(blockEngine.chainComponents[block] == blockEngine.chainComponents[block - 1]){
        rhs[block] -= blockEngine.chainTransfers[block].t() * rhs[block - 1];
      }
    }
    for(std::size_t block = 0; block < rhs.size(); ++block){
      rhs[block] = blockEngineSolveD(block, rhs[block]);
    }
    for(std::size_t block = rhs.size(); block-- > 1; ){
      if(blockEngine.chainComponents[block] == blockEngine.chainComponents[block - 1]){
        rhs[block - 1] -= blockEngine.chainTransfers[block] * rhs[block];
      }
    }
    return rhs;
  };

  auto blockBandedSolveInterior = [&](std::vector<arma::mat> rhs) -> std::vector<arma::mat> {
    const std::size_t nBlocks = rhs.size();
    const std::size_t bandwidth = blockEngine.bandWidth;
    for(std::size_t block = 0; block < nBlocks; ++block){
      const std::size_t first = block > bandwidth ? block - bandwidth : 0;
      for(std::size_t previous = first; previous < block; ++previous){
        rhs[block] -= blockEngine.bandLower[block][block - previous] * rhs[previous];
      }
      arma::mat solved;
      if(!gpuBackend.triangularSolve(solved, arma::trimatl(blockEngine.Lg[block]), rhs[block],
                      arma::solve_opts::fast)){
        Rcpp::stop("Block-banded forward triangular solve failed.");
      }
      rhs[block] = std::move(solved);
    }
    for(std::size_t block = nBlocks; block-- > 0;){
      const std::size_t last = std::min(nBlocks - 1, block + bandwidth);
      for(std::size_t following = block + 1; following <= last; ++following){
        rhs[block] -= blockEngine.bandLower[following][following - block].t() * rhs[following];
      }
      arma::mat solved;
      if(!gpuBackend.triangularSolve(solved, arma::trimatu(blockEngine.Lg[block].t()), rhs[block],
                      arma::solve_opts::fast)){
        Rcpp::stop("Block-banded backward triangular solve failed.");
      }
      rhs[block] = std::move(solved);
    }
    return rhs;
  };

  auto blockEngineFactorize = [&](const EigenSpMat & A) -> void {
    const std::size_t nG = blockEngine.groups.size();
    blockEngine.Lg.resize(nG);
    blockEngine.Fg.resize(nG);
    blockEngine.Cgx.resize(nG);
    blockEngine.logDet = 0.0;
    const arma::uword nB = static_cast<arma::uword>(blockEngine.borderCols.size());

    if(blockEngine.latent){
      const int term = blockEngine.latentTerm;
      const DescriptorPrecisionState & state = randomPrecisionState[static_cast<std::size_t>(term)];
      const arma::uword dimension = blockEngine.groups.front().size();
      const arma::uword rank = state.lowRankLoadings.n_cols;
      const arma::uword originalBorder = blockEngine.originalBorderSize;
      const arma::uword lowRankAxis = static_cast<arma::uword>(state.lowRankFactorIndex);
      arma::mat basePrecision(1, 1, arma::fill::ones);
      std::size_t factorMultiplier = 1;
      double baseLogDet = 0.0;
      for(std::size_t axis = 0; axis < state.factorPrecision.size(); ++axis){
        if(axis == lowRankAxis){ continue; }
        const arma::uword axisDimension = factorEngine.factorPrecisionDimension(state, axis);
        if(factorMultiplier > std::numeric_limits<std::size_t>::max() / axisDimension){
          Rcpp::stop("Latent base precision dimensions overflow.");
        }
        factorMultiplier *= axisDimension;
        basePrecision = arma::kron(basePrecision, factorEngine.factorPrecisionDense(state, axis));
        baseLogDet -= static_cast<double>(dimension / axisDimension) *
          state.factorLogDet[axis];
      }
      baseLogDet -= static_cast<double>(factorMultiplier) * logDetA(term);
      basePrecision = arma::kron(basePrecision, arma::mat(Ai(term)));
      if(basePrecision.n_rows != dimension || basePrecision.n_cols != dimension ||
         !basePrecision.is_finite() || !std::isfinite(baseLogDet)){
        blockEngine.patternChecked = false;
        ++factorSchurFallbacks;
        return;
      }
      blockEngine.latentBasePrecision = std::move(basePrecision);
      blockEngine.latentBaseLogDet = baseLogDet;
      const arma::mat & basePrecisionCurrent = blockEngine.latentBasePrecision;
      arma::mat weightedLoadings = state.lowRankLoadings;
      weightedLoadings.each_col() /= state.lowRankSpecific;
      const arma::mat core = arma::eye<arma::mat>(rank, rank) + state.lowRankLoadings.t() * weightedLoadings;
      arma::mat schur(nB, nB, arma::fill::zeros);
      for(arma::uword column = 0; column < originalBorder; ++column){
        for(EigenSpMat::InnerIterator entry(A, blockEngine.borderCols[column]); entry; ++entry){
          const int row = blockEngine.borderPos[entry.row()];
          if(row >= 0){ schur(row, column) = entry.value(); }
        }
      }
      for(arma::uword first = 0; first < rank; ++first){
        for(arma::uword second = 0; second < rank; ++second){
          schur.submat(originalBorder + first * dimension, originalBorder + second * dimension,
            originalBorder + (first + 1) * dimension - 1, originalBorder + (second + 1) * dimension - 1) =
              state.inverseScale * core(first, second) * basePrecisionCurrent;
        }
      }
      for(std::size_t group = 0; group < nG; ++group){
        arma::mat diagonal(dimension, dimension, arma::fill::zeros);
        arma::mat & cross = blockEngine.Cgx[group];
        cross.zeros(dimension, nB);
        for(arma::uword column = 0; column < dimension; ++column){
          for(EigenSpMat::InnerIterator entry(A, blockEngine.groups[group][column]); entry; ++entry){
            const int border = blockEngine.borderPos[entry.row()];
            if(border >= 0){ cross(column, border) = entry.value(); }
            else if(blockEngine.groupOf[entry.row()] == static_cast<int>(group)){
              diagonal(blockEngine.localOf[entry.row()], column) = entry.value();
            }
          }
        }
        const double marginalDiagonal = factorEngine.factorPrecisionValue(state, lowRankAxis, group, group);
        diagonal += state.inverseScale *
          (1.0 / state.lowRankSpecific(group) - marginalDiagonal) * basePrecisionCurrent;
        for(arma::uword component = 0; component < rank; ++component){
          cross.cols(originalBorder + component * dimension, originalBorder + (component + 1) * dimension - 1) =
            -state.inverseScale * weightedLoadings(group, component) * basePrecisionCurrent;
        }
        diagonal = 0.5 * (diagonal + diagonal.t());
        if(!gpuBackend.cholesky(blockEngine.Lg[group], diagonal, "lower")){
          blockEngine.patternChecked = false;
          ++factorSchurFallbacks;
          return;
        }
        blockEngine.logDet += 2.0 * arma::accu(arma::log(blockEngine.Lg[group].diag()));
        arma::mat forward;
        if(!gpuBackend.triangularSolve(forward, arma::trimatl(blockEngine.Lg[group]), cross, arma::solve_opts::fast) ||
           !gpuBackend.triangularSolve(blockEngine.Fg[group], arma::trimatu(blockEngine.Lg[group].t()),
             forward, arma::solve_opts::fast)){
          Rcpp::stop("Latent block-Schur triangular solve failed.");
        }
        schur -= gpuBackend.product(forward, forward, true);
      }
      schur = 0.5 * (schur + schur.t());
      if(!gpuBackend.cholesky(blockEngine.Ls, schur, "lower")){
        blockEngine.patternChecked = false;
        ++factorSchurFallbacks;
        return;
      }
      blockEngine.logDet += 2.0 * arma::accu(arma::log(blockEngine.Ls.diag()));
      blockEngine.logDet -= static_cast<double>(dimension) * state.lowRankCoreLogDet
        + static_cast<double>(dimension * rank) * std::log(state.inverseScale)
        + static_cast<double>(rank) * blockEngine.latentBaseLogDet;
      if(!std::isfinite(blockEngine.logDet)){
        blockEngine.patternChecked = false;
        ++factorSchurFallbacks;
        return;
      }
      blockEngine.active = true;
      return;
    }

    if(blockEngine.chain){
      std::vector<arma::mat> lower(nG);
      blockEngine.chainTransfers.resize(nG);
      arma::mat schur(nB, nB, arma::fill::zeros);
      for(arma::uword column = 0; column < nB; ++column){
        for(EigenSpMat::InnerIterator entry(A, blockEngine.borderCols[column]); entry; ++entry){
          const int row = blockEngine.borderPos[entry.row()];
          if(row >= 0){ schur(row, column) = entry.value(); }
        }
      }
      for(std::size_t block = 0; block < nG; ++block){
        const arma::uword size = blockEngine.groups[block].size();
        blockEngine.Lg[block].zeros(size, size);
        blockEngine.Cgx[block].zeros(size, nB);
        blockEngine.chainTransfers[block].reset();
        if(block > 0 && blockEngine.chainComponents[block] == blockEngine.chainComponents[block - 1]){
          lower[block].zeros(size, blockEngine.groups[block - 1].size());
        }
      }
      for(std::size_t block = 0; block < nG; ++block){
        for(arma::uword column = 0; column < blockEngine.groups[block].size(); ++column){
          for(EigenSpMat::InnerIterator entry(A, blockEngine.groups[block][column]); entry; ++entry){
            const int border = blockEngine.borderPos[entry.row()];
            const int other = blockEngine.groupOf[entry.row()];
            if(border >= 0){ blockEngine.Cgx[block](column, border) = entry.value(); }
            else if(other == static_cast<int>(block)){
              blockEngine.Lg[block](blockEngine.localOf[entry.row()], column) = entry.value();
            }else if(other == static_cast<int>(block + 1)){
              lower[block + 1](blockEngine.localOf[entry.row()], column) = entry.value();
            }
          }
        }
      }
      for(std::size_t block = 0; block < nG; ++block){
        arma::mat pivot = std::move(blockEngine.Lg[block]);
        if(!lower[block].is_empty()){
          blockEngine.chainTransfers[block] = blockEngineSolveD(block - 1, lower[block].t());
          pivot -= gpuBackend.product(lower[block], blockEngine.chainTransfers[block]);
          lower[block].reset();
        }
        pivot = 0.5 * (pivot + pivot.t());
        if(!gpuBackend.cholesky(blockEngine.Lg[block], pivot, "lower")){
          throw sommer::CovarianceFactorizationError("Block-chain pivot Cholesky failed (matrix not positive definite).");
        }
        blockEngine.logDet += 2.0 * arma::accu(arma::log(blockEngine.Lg[block].diag()));
      }
      blockEngine.Fg = blockChainSolveInterior(blockEngine.Cgx);
      for(std::size_t block = 0; block < nG; ++block){
        schur -= gpuBackend.product(blockEngine.Cgx[block], blockEngine.Fg[block], true);
      }
      if(nB > 0){
        schur = 0.5 * (schur + schur.t());
        if(!gpuBackend.cholesky(blockEngine.Ls, schur, "lower")){
          throw sommer::CovarianceFactorizationError("Block-chain border Cholesky failed (matrix not positive definite).");
        }
        blockEngine.logDet += 2.0 * arma::accu(arma::log(blockEngine.Ls.diag()));
      }else{ blockEngine.Ls.reset(); }
      blockEngine.active = true;
      return;
    }

    if(blockEngine.banded){
      const std::size_t bandwidth = blockEngine.bandWidth;
      std::vector<arma::mat> diagonal(nG);
      blockEngine.bandLower.assign(nG, std::vector<arma::mat>(bandwidth + 1));
      arma::mat schur(nB, nB, arma::fill::zeros);
      for(arma::uword column = 0; column < nB; ++column){
        for(EigenSpMat::InnerIterator entry(A, blockEngine.borderCols[column]); entry; ++entry){
          const int row = blockEngine.borderPos[entry.row()];
          if(row >= 0){ schur(row, column) = entry.value(); }
        }
      }
      for(std::size_t block = 0; block < nG; ++block){
        const arma::uword size = blockEngine.groups[block].size();
        diagonal[block].zeros(size, size);
        blockEngine.Lg[block].reset();
        blockEngine.Cgx[block].zeros(size, nB);
        for(std::size_t offset = 1; offset <= bandwidth && offset <= block; ++offset){
          blockEngine.bandLower[block][offset].zeros(
            size, blockEngine.groups[block - offset].size());
        }
      }
      for(std::size_t columnBlock = 0; columnBlock < nG; ++columnBlock){
        const std::vector<int> & columns = blockEngine.groups[columnBlock];
        for(std::size_t localColumn = 0; localColumn < columns.size(); ++localColumn){
          const int globalColumn = columns[localColumn];
          for(EigenSpMat::InnerIterator entry(A, globalColumn); entry; ++entry){
            const int border = blockEngine.borderPos[entry.row()];
            if(border >= 0){
              blockEngine.Cgx[columnBlock](localColumn, border) = entry.value();
              continue;
            }
            const int rowBlock = blockEngine.groupOf[entry.row()];
            if(rowBlock < static_cast<int>(columnBlock)){ continue; }
            const std::size_t offset = static_cast<std::size_t>(rowBlock) - columnBlock;
            if(offset > bandwidth){
              Rcpp::stop("Block-banded admission disagrees with the assembled coefficient graph.");
            }
            if(offset == 0){
              diagonal[columnBlock](blockEngine.localOf[entry.row()], localColumn) = entry.value();
            }else{
              blockEngine.bandLower[static_cast<std::size_t>(rowBlock)][offset](
                blockEngine.localOf[entry.row()], localColumn) = entry.value();
            }
          }
        }
      }

      for(std::size_t block = 0; block < nG; ++block){
        const std::size_t first = block > bandwidth ? block - bandwidth : 0;
        for(std::size_t previous = first; previous < block; ++previous){
          arma::mat & lower = blockEngine.bandLower[block][block - previous];
          const std::size_t firstCommon = previous > bandwidth ? previous - bandwidth : 0;
          for(std::size_t common = std::max(first, firstCommon); common < previous; ++common){
            lower -= gpuBackend.product(blockEngine.bandLower[block][block - common],
              blockEngine.bandLower[previous][previous - common], false, true);
          }
          arma::mat solvedTranspose;
          if(!gpuBackend.triangularSolve(solvedTranspose, arma::trimatl(blockEngine.Lg[previous]),
                          lower.t(), arma::solve_opts::fast)){
            Rcpp::stop("Block-banded factor triangular solve failed.");
          }
          lower = solvedTranspose.t();
        }
        arma::mat pivot = std::move(diagonal[block]);
        for(std::size_t previous = first; previous < block; ++previous){
          const arma::mat & lower = blockEngine.bandLower[block][block - previous];
          pivot -= gpuBackend.product(lower, lower, false, true);
        }
        pivot = 0.5 * (pivot + pivot.t());
        if(!gpuBackend.cholesky(blockEngine.Lg[block], pivot, "lower")){
          throw sommer::CovarianceFactorizationError("Block-banded pivot Cholesky failed (matrix not positive definite).");
        }
        blockEngine.logDet += 2.0 * arma::accu(arma::log(blockEngine.Lg[block].diag()));
      }

      blockEngine.Fg = blockBandedSolveInterior(blockEngine.Cgx);
      for(std::size_t block = 0; block < nG; ++block){
        schur -= gpuBackend.product(blockEngine.Cgx[block], blockEngine.Fg[block], true);
      }
      if(nB > 0){
        schur = 0.5 * (schur + schur.t());
        if(!gpuBackend.cholesky(blockEngine.Ls, schur, "lower")){
          throw sommer::CovarianceFactorizationError("Block-banded border Schur complement is not positive definite.");
        }
        blockEngine.logDet += 2.0 * arma::accu(arma::log(blockEngine.Ls.diag()));
      }else{
        blockEngine.Ls.reset();
      }
      blockEngine.active = true;
      return;
    }

    arma::mat S(nB, nB, arma::fill::zeros);
    for(arma::uword jb = 0; jb < nB; ++jb){
      for(EigenSpMat::InnerIterator it(A, blockEngine.borderCols[jb]); it; ++it){
        const int ib = blockEngine.borderPos[static_cast<std::size_t>(it.row())];
        if(ib >= 0){ S(static_cast<arma::uword>(ib), jb) = it.value(); }
      }
    }

    for(std::size_t g = 0; g < nG; ++g){
      const std::vector<int> & cols = blockEngine.groups[g];
      const arma::uword m = static_cast<arma::uword>(cols.size());
      arma::mat D(m, m, arma::fill::zeros);
      arma::mat & Cx = blockEngine.Cgx[g];
      Cx.zeros(m, nB);
      for(arma::uword jj = 0; jj < m; ++jj){
        for(EigenSpMat::InnerIterator it(A, cols[jj]); it; ++it){
          const std::size_t i = static_cast<std::size_t>(it.row());
          const int ib = blockEngine.borderPos[i];
          if(ib >= 0){
            Cx(jj, static_cast<arma::uword>(ib)) = it.value();
          }else{
            D(static_cast<arma::uword>(blockEngine.localOf[i]), jj) = it.value();
          }
        }
      }
      if(!gpuBackend.cholesky(blockEngine.Lg[g], D, "lower")){
        throw sommer::CovarianceFactorizationError("Dense block Cholesky of the MME coefficient matrix C failed (matrix not positive definite).");
      }
      blockEngine.logDet += 2.0 * arma::accu(arma::log(blockEngine.Lg[g].diag()));
      if(nB > 0){
        std::vector<arma::uword> activeColumns;
        for(arma::uword b = 0; b < nB; ++b){
          if(arma::any(Cx.col(b) != 0.0)){ activeColumns.push_back(b); }
        }
        blockEngine.Fg[g].zeros(m, nB);
        if(!activeColumns.empty()){
          arma::uvec active(activeColumns.size());
          for(std::size_t j = 0; j < activeColumns.size(); ++j){
            active(static_cast<arma::uword>(j)) = activeColumns[j];
          }
          const arma::mat CxActive = Cx.cols(active);
          arma::mat forward;
          arma::mat FActive;
          if(!gpuBackend.triangularSolve(forward, arma::trimatl(blockEngine.Lg[g]), CxActive, arma::solve_opts::fast) ||
             !gpuBackend.triangularSolve(FActive, arma::trimatu(blockEngine.Lg[g].t()), forward, arma::solve_opts::fast)){
            Rcpp::stop("Dense block-Schur triangular solve failed.");
          }
          blockEngine.Fg[g].cols(active) = FActive;
          S.submat(active, active) -= gpuBackend.product(forward, forward, true);
        }
      }else{
        blockEngine.Fg[g].reset();
      }
    }

    if(nB > 0){
      S = 0.5 * (S + S.t());
      if(!gpuBackend.cholesky(blockEngine.Ls, S, "lower")){
        throw sommer::CovarianceFactorizationError("Dense block-Schur complement of the fixed effects is not positive definite.");
      }
      blockEngine.logDet += 2.0 * arma::accu(arma::log(blockEngine.Ls.diag()));
    }else{
      blockEngine.Ls.reset();
    }
    blockEngine.active = true;
  };

  auto blockEngineSolve = [&](const Eigen::Ref<const Eigen::MatrixXd> & rhs) -> Eigen::MatrixXd {
    const std::size_t nG = blockEngine.groups.size();
    const arma::uword k = static_cast<arma::uword>(rhs.cols());
    Eigen::MatrixXd out(rhs.rows(), rhs.cols());
    const std::vector<int> & border = blockEngine.borderCols;
    const arma::uword nB = static_cast<arma::uword>(border.size());

    arma::mat Rx(nB, k, arma::fill::zeros);
    for(arma::uword c = 0; c < k; ++c){
      for(arma::uword i = 0; i < nB; ++i){
        if(border[i] < rhs.rows()){ Rx(i, c) = rhs(border[i], static_cast<Eigen::Index>(c)); }
      }
    }
    std::vector<arma::mat> Tg(nG);
    for(std::size_t g = 0; g < nG; ++g){
      const std::vector<int> & cols = blockEngine.groups[g];
      arma::mat Rg(static_cast<arma::uword>(cols.size()), k);
      for(arma::uword c = 0; c < k; ++c){
        for(std::size_t jj = 0; jj < cols.size(); ++jj){
          Rg(static_cast<arma::uword>(jj), c) = rhs(cols[jj], static_cast<Eigen::Index>(c));
        }
      }
      Tg[g] = (blockEngine.chain || blockEngine.banded)
        ? std::move(Rg) : blockEngineSolveD(g, Rg);
      if(nB > 0 && !blockEngine.chain && !blockEngine.banded){
        Rx -= blockEngine.Cgx[g].t() * Tg[g];
      }
    }
    if(blockEngine.chain || blockEngine.banded){
      Tg = blockEngine.chain
        ? blockChainSolveInterior(std::move(Tg))
        : blockBandedSolveInterior(std::move(Tg));
      if(nB > 0){
        for(std::size_t block = 0; block < nG; ++block){ Rx -= blockEngine.Cgx[block].t() * Tg[block]; }
      }
    }

    arma::mat B;
    if(nB > 0){
      arma::mat tmp;
      const bool ok =
        gpuBackend.triangularSolve(tmp, arma::trimatl(blockEngine.Ls), Rx, arma::solve_opts::fast) &&
        gpuBackend.triangularSolve(B, arma::trimatu(blockEngine.Ls.t()), tmp, arma::solve_opts::fast);
      if(!ok){ Rcpp::stop("Dense block-Schur border solve failed."); }
      for(arma::uword c = 0; c < k; ++c){
        for(arma::uword i = 0; i < nB; ++i){
          if(border[i] < rhs.rows()){ out(border[i], static_cast<Eigen::Index>(c)) = B(i, c); }
        }
      }
    }
    for(std::size_t g = 0; g < nG; ++g){
      if(nB > 0){ Tg[g] -= blockEngine.Fg[g] * B; }
      const std::vector<int> & cols = blockEngine.groups[g];
      for(arma::uword c = 0; c < k; ++c){
        for(std::size_t jj = 0; jj < cols.size(); ++jj){
          out(cols[jj], static_cast<Eigen::Index>(c)) = Tg[g](static_cast<arma::uword>(jj), c);
        }
      }
    }
    return out;
  };

  auto solveCVectorCholmod =
    [&](const Eigen::Ref<const Eigen::VectorXd> & rhs,
        const std::string & context) -> Eigen::VectorXd {
      if(blockEngine.active){
        return blockEngineSolve(rhs).col(0);
      }
      return cholmodState.solve(rhs, "CHOLMOD solve failed in " + context + ".").col(0);
    };

  auto solveCMatrixCholmod =
    [&](const Eigen::Ref<const Eigen::MatrixXd> & rhs,
        const std::string & context) -> Eigen::MatrixXd {
      if(blockEngine.active){
        return blockEngineSolve(rhs);
      }
      return cholmodState.solve(rhs, "CHOLMOD multi-RHS solve failed in " + context + ".");
    };

  // CHOLMOD backend for D (the random-effects-only block used by ML,
  // reml==false), kept in a separate state/factor from C's - D is a
  // different matrix (Nu x Nu, no fixed-effect rows/cols) with its own
  // sparsity pattern, mirroring how cholmodRState is kept separate from
  // cholmodState for R.
  CholmodState & cholmodDState = cholmodCache->random;

  auto solveDMatrixCholmod =
    [&](const Eigen::Ref<const Eigen::MatrixXd> & rhs,
        const std::string & context) -> Eigen::MatrixXd {
      return cholmodDState.solve(rhs,
        "CHOLMOD multi-RHS solve for the random-effects-only matrix D failed in " +
        context + " (reml=FALSE).");
    };

  bool dFirstMLActive = false;
  int dFirstMLEvaluations = 0;
  EigenSpMat dFirstCross;
  Eigen::MatrixXd dFirstTransfer;
  Eigen::LLT<Eigen::MatrixXd> dFirstFixedFactor;
  auto solveRandomBlock = [&](const Eigen::Ref<const Eigen::MatrixXd> & rhs,
                               const std::string & context) -> Eigen::MatrixXd {
    if(saddleActive){ return saddleRandom.solveCoefficients(rhs, context); }
    if(solverName == "cholmod"){ return solveDMatrixCholmod(rhs, context); }
    Eigen::MatrixXd result = Dfactor.solve(rhs);
    if(Dfactor.info() != Eigen::Success || !result.allFinite()){
      Rcpp::stop("Random-block LDLT solve failed in " + context + ".");
    }
    return result;
  };
  auto solveDfirstJoint = [&](const Eigen::Ref<const Eigen::MatrixXd> & rhs,
                               const std::string & context) -> Eigen::MatrixXd {
    Eigen::MatrixXd result(nEffects, rhs.cols());
    if(rhs.cols() == 0){ return result; }
    Eigen::MatrixXd random = solveRandomBlock(Eigen::MatrixXd(rhs.bottomRows(Nu)), context);
    if(nX > 0){
      const Eigen::MatrixXd fixedRhs = rhs.topRows(nX) - dFirstCross * random;
      result.topRows(nX) = dFirstFixedFactor.solve(fixedRhs);
      random.noalias() -= dFirstTransfer * result.topRows(nX);
    }
    result.bottomRows(Nu) = random;
    if(!result.allFinite()){ Rcpp::stop("D-first joint solve produced non-finite values in " + context + "."); }
    return result;
  };

  // Exact dense blocks of C^{-1} (border-border, border-u_g, u_g-u_g), shared
  // by the random and residual score traces. The border is X for the CHOLMOD
  // factor path and blockEngine.borderCols for the block-Schur engine.
  struct CholmodSelectedInverseState {
    bool selReady = false;
    const cholmod_factor * selL = nullptr;
    std::vector<double> selZ;
    std::vector<int> selColSuper;
    std::vector<int> selIperm;
  };
  struct CholmodInverseBlockCache : CholmodSelectedInverseState {
    bool ready = false;
    bool patternOnly = false;
    arma::mat Cxx;
    std::vector<arma::mat> Cuu;
    std::vector<arma::mat> Cxu;
    std::vector<int> borderPos;
    std::vector<int> groupOf;
    std::vector<int> localOf;
    std::vector<arma::mat> FS;
    std::map< std::pair<int,int>, arma::mat > crossBlocks;
    double crossBlockDoubles = 0.0;
    bool crossAvailable = false;
    // Supernodal selected inverse: C^{-1} on the pattern of L, laid out like L->x.
  };
  CholmodInverseBlockCache cholmodInvCache;
  int patternInverseCacheBuilds = 0;
  int patternInverseCacheFallbacks = 0;

  struct CrossGroupTraceDemand {
    std::vector<std::size_t> groupSizes;
    std::vector<std::size_t> entries;
    std::uint64_t groupSignature = 0;
  };
  CrossGroupTraceDemand previousCrossGroupDemand;
  CrossGroupTraceDemand currentCrossGroupDemand;
  bool reportedCrossGroupCache = false;
  arma::uword currentResidualParameterCount = 0;
  const std::size_t crossGroupCacheMaxGroups = sommer::WorkspacePolicy::crossGroups();
  const double crossGroupCacheMaxDoubles = cholmodInvCacheMaxDoubles / 16.0;
  // Y = L^{-1} P E_cols; supernodes whose RHS rows are still zero are skipped,
  // so only the elimination-tree reach of `cols` is processed.
  auto cholmodForwardSolveUnit =
    [&](const cholmod_factor * L,
        const std::vector<int> & iperm,
        const std::vector<int> & cols,
        arma::mat & Y,
        arma::uvec & touchedRows) -> void {
      const int n = static_cast<int>(L->n);
      const int m = static_cast<int>(cols.size());
      Y.zeros(n, m);
      std::vector<char> touched(static_cast<std::size_t>(n), 0);
      for(int j = 0; j < m; ++j){
        Y(static_cast<arma::uword>(iperm[static_cast<std::size_t>(cols[static_cast<std::size_t>(j)])]),
          static_cast<arma::uword>(j)) = 1.0;
      }
      const int * super = static_cast<const int *>(L->super);
      const int * pi = static_cast<const int *>(L->pi);
      const int * px = static_cast<const int *>(L->px);
      const int * srows = static_cast<const int *>(L->s);
      const double * Lx = static_cast<const double *>(L->x);

      for(std::size_t sn = 0; sn < L->nsuper; ++sn){
        const int k1 = super[sn];
        const int k2 = super[sn + 1];
        const int nscol = k2 - k1;
        const int psi = pi[sn];
        const int nsrow = pi[sn + 1] - psi;

        bool active = false;
        for(int j = 0; j < m && !active; ++j){
          const double * column = Y.colptr(static_cast<arma::uword>(j)) + k1;
          for(int r = 0; r < nscol; ++r){
            if(column[r] != 0.0){ active = true; break; }
          }
        }
        if(!active){ continue; }
        std::fill(touched.begin() + k1, touched.begin() + k2, 1);

        const arma::mat Ls(
          const_cast<double *>(Lx + static_cast<std::ptrdiff_t>(px[sn])),
          static_cast<arma::uword>(nsrow),
          static_cast<arma::uword>(nscol),
          false,
          true
        );
        const arma::mat Yin = Y.rows(static_cast<arma::uword>(k1), static_cast<arma::uword>(k2 - 1));
        const arma::mat L11 = arma::trimatl(Ls.head_rows(static_cast<arma::uword>(nscol)));
        arma::mat Yk;
        if(!gpuBackend.triangularSolve(Yk, arma::trimatl(L11), Yin, arma::solve_opts::fast)){
          Rcpp::stop("Triangular solve failed while forming CHOLMOD inverse blocks.");
        }
        Y.rows(static_cast<arma::uword>(k1), static_cast<arma::uword>(k2 - 1)) = Yk;

        const int nrest = nsrow - nscol;
        if(nrest > 0){
          const arma::mat U = Ls.tail_rows(static_cast<arma::uword>(nrest)) * Yk;
          const int * rest = srows + psi + nscol;
          for(int j = 0; j < m; ++j){
            double * yc = Y.colptr(static_cast<arma::uword>(j));
            const double * uc = U.colptr(static_cast<arma::uword>(j));
            for(int ii = 0; ii < nrest; ++ii){
              yc[rest[ii]] -= uc[ii];
            }
          }
        }
      }

      std::vector<arma::uword> rowsOut;
      for(int r = 0; r < n; ++r){
        if(touched[static_cast<std::size_t>(r)]){ rowsOut.push_back(static_cast<arma::uword>(r)); }
      }
      touchedRows = arma::uvec(rowsOut);
    };

  // Takahashi recursion over supernodes (last to first), with unit factor
  // [I; Y], Y = L_RJ L_JJ^{-1}: Z_RJ = -Z_RR Y, Z_JJ = (L_JJ L_JJ')^{-1} - Y' Z_RJ.
  auto buildCholmodSelectedInverse = [&](const cholmod_factor * L, const std::vector<int> & iperm,
                                         CholmodSelectedInverseState & selected) -> bool {
    selected.selReady = false;
    if(L == nullptr || !L->is_super || !L->is_ll || L->xtype != CHOLMOD_REAL ||
       L->dtype != CHOLMOD_DOUBLE || L->itype != CHOLMOD_INT){ return false; }
    const int n = static_cast<int>(L->n);
    const int nsuper = static_cast<int>(L->nsuper);
    const int * super = static_cast<const int *>(L->super);
    const int * pi = static_cast<const int *>(L->pi);
    const int * px = static_cast<const int *>(L->px);
    const int * srows = static_cast<const int *>(L->s);
    const double * Lx = static_cast<const double *>(L->x);
    std::vector<double> & Z = selected.selZ;
    Z.assign(static_cast<std::size_t>(L->xsize), 0.0);
    std::vector<int> & colSuper = selected.selColSuper;
    colSuper.assign(static_cast<std::size_t>(n), -1);
    for(int sn = 0; sn < nsuper; ++sn){
      for(int k = super[sn]; k < super[sn + 1]; ++k){ colSuper[static_cast<std::size_t>(k)] = sn; }
    }
    std::vector<int> rowPos(static_cast<std::size_t>(n), -1);

    for(int sn = nsuper - 1; sn >= 0; --sn){
      const int k1 = super[sn];
      const int nscol = super[sn + 1] - k1;
      const int psi = pi[sn];
      const int nsrow = pi[sn + 1] - psi;
      const int nrest = nsrow - nscol;
      const arma::mat Ls(const_cast<double *>(Lx + static_cast<std::ptrdiff_t>(px[sn])),
                         static_cast<arma::uword>(nsrow), static_cast<arma::uword>(nscol), false, true);
      arma::mat L11inv;
      if(!arma::inv(L11inv, arma::trimatl(Ls.head_rows(static_cast<arma::uword>(nscol))))){ return false; }
      arma::mat Zjj = L11inv.t() * L11inv;
      arma::mat Zout(Z.data() + static_cast<std::ptrdiff_t>(px[sn]),
                     static_cast<arma::uword>(nsrow), static_cast<arma::uword>(nscol), false, true);
      if(nrest > 0){
        const int * R = srows + psi + nscol;
        const arma::mat Y = Ls.tail_rows(static_cast<arma::uword>(nrest)) * L11inv;
        arma::mat Zrr(static_cast<arma::uword>(nrest), static_cast<arma::uword>(nrest));
        int a = 0;
        while(a < nrest){
          const int t = colSuper[static_cast<std::size_t>(R[a])];
          const int tpsi = pi[t];
          const int tnsrow = pi[t + 1] - tpsi;
          for(int q = 0; q < tnsrow; ++q){ rowPos[static_cast<std::size_t>(srows[tpsi + q])] = q; }
          int b = a;
          bool ok = true;
          for(; b < nrest && colSuper[static_cast<std::size_t>(R[b])] == t; ++b){
            const double * base = Z.data() + static_cast<std::ptrdiff_t>(px[t]) +
              static_cast<std::ptrdiff_t>(R[b] - super[t]) * tnsrow;
            for(int r = b; r < nrest; ++r){
              const int p = rowPos[static_cast<std::size_t>(R[r])];
              if(p < 0){ ok = false; break; }
              Zrr(static_cast<arma::uword>(r), static_cast<arma::uword>(b)) = base[p];
            }
            if(!ok){ break; }
          }
          for(int q = 0; q < tnsrow; ++q){ rowPos[static_cast<std::size_t>(srows[tpsi + q])] = -1; }
          if(!ok){ return false; }
          a = b;
        }
        Zrr = arma::symmatl(Zrr);
        const arma::mat Zrj = -Zrr * Y;
        Zjj -= Y.t() * Zrj;
        Zout.tail_rows(static_cast<arma::uword>(nrest)) = Zrj;
      }
      Zout.head_rows(static_cast<arma::uword>(nscol)) = Zjj;
    }
    selected.selL = L;
    selected.selIperm = iperm;
    selected.selReady = true;
    return true;
  };

  int cholmodDSelectedBuilds = 0;
  int cholmodRSelectedBuilds = 0;
  auto buildCholmodTraceSubset = [&](const cholmod_factor * factor,
                                      SelectedInverseSubset & subset) -> bool {
    subset = SelectedInverseSubset();
    if(factor == nullptr || !factor->is_super || !factor->is_ll ||
       factor->itype != CHOLMOD_INT){ return false; }
    const int size = static_cast<int>(factor->n);
    const int * rowStarts = static_cast<const int *>(factor->pi);
    double largestFront = 0.0;
    for(std::size_t node = 0; node < factor->nsuper; ++node){
      largestFront = std::max(largestFront,
        static_cast<double>(rowStarts[node + 1] - rowStarts[node]));
    }
    const double selectedNeededDoubles =
      2.0 * static_cast<double>(factor->xsize) + 4.0 * largestFront * largestFront;
    if(selectedNeededDoubles > cholmodInvCacheMaxDoubles / 3.0){
      if(memoryStrict){
        requireMemory("The CHOLMOD selected inverse", 3.0 * selectedNeededDoubles);
      }
      return false;
    }
    const int * permutation = static_cast<const int *>(factor->Perm);
    std::vector<int> inversePermutation(static_cast<std::size_t>(size));
    for(int index = 0; index < size; ++index){
      inversePermutation[permutation[index]] = index;
    }
    CholmodSelectedInverseState selected;
    if(!buildCholmodSelectedInverse(factor, inversePermutation, selected)){ return false; }
    subset.originalToPermuted = std::move(inversePermutation);
    subset.rows.resize(size);
    subset.values.resize(size);
    const int * super = static_cast<const int *>(factor->super);
    const int * starts = static_cast<const int *>(factor->px);
    const int * rows = static_cast<const int *>(factor->s);
    for(std::size_t node = 0; node < factor->nsuper; ++node){
      const int width = rowStarts[node + 1] - rowStarts[node];
      for(int column = super[node]; column < super[node + 1]; ++column){
        const std::size_t offset = static_cast<std::size_t>(starts[node]) +
          static_cast<std::size_t>(column - super[node]) * width;
        for(int local = 0; local < width; ++local){
          const int row = rows[rowStarts[node] + local];
          if(row >= column){
            subset.rows[column].push_back(row);
            subset.values[column].push_back(selected.selZ[offset + local]);
          }
        }
      }
    }
    return true;
  };

  auto cholmodSelectedLookup = [&](int i, int j, double & value) -> bool {
    const cholmod_factor * L = cholmodInvCache.selL;
    int a = cholmodInvCache.selIperm[static_cast<std::size_t>(i)];
    int b = cholmodInvCache.selIperm[static_cast<std::size_t>(j)];
    if(a > b){ std::swap(a, b); }
    const int t = cholmodInvCache.selColSuper[static_cast<std::size_t>(a)];
    const int * super = static_cast<const int *>(L->super);
    const int * pi = static_cast<const int *>(L->pi);
    const int * px = static_cast<const int *>(L->px);
    const int * srows = static_cast<const int *>(L->s);
    const int nscol = super[t + 1] - super[t];
    const int nsrow = pi[t + 1] - pi[t];
    int p;
    if(b < super[t + 1]){
      p = b - super[t];
    }else{
      const int * first = srows + pi[t] + nscol;
      const int * last = srows + pi[t + 1];
      const int * hit = std::lower_bound(first, last, b);
      if(hit == last || *hit != b){ return false; }
      p = static_cast<int>(hit - (srows + pi[t]));
    }
    value = cholmodInvCache.selZ[static_cast<std::size_t>(px[t]) +
      static_cast<std::size_t>(a - super[t]) * static_cast<std::size_t>(nsrow) + static_cast<std::size_t>(p)];
    return true;
  };

  auto buildCholmodInverseBlockCache = [&](bool patternLocal = false) -> void {
    cholmodInvCache.ready = false;
    cholmodInvCache.patternOnly = false;
    cholmodInvCache.selReady = false;
    cholmodInvCache.crossBlocks.clear();
    cholmodInvCache.crossBlockDoubles = 0.0;
    cholmodInvCache.Cuu.clear();
    cholmodInvCache.Cxu.clear();
    cholmodInvCache.FS.clear();
    cholmodInvCache.Cxx.reset();
    cholmodInvCache.crossAvailable = false;

    if(blockEngine.active){
      if(blockEngine.banded){ return; }
      // C^{-1}_bb = S^{-1}, C^{-1}_bg = -S^{-1} F_g', C^{-1}_gh = [g==h] D_g^{-1} + F_g S^{-1} F_h'.
      const bool haveBorder = !blockEngine.borderCols.empty();
      arma::mat Sinv;
      if(haveBorder){
        arma::mat LsInv;
        if(!arma::inv(LsInv, arma::trimatl(blockEngine.Ls))){
          Rcpp::stop("Dense block-Schur inverse of the border Schur complement failed.");
        }
        Sinv = LsInv.t() * LsInv;
        cholmodInvCache.Cxx = Sinv;
      }
      const std::size_t nG = blockEngine.groups.size();
      cholmodInvCache.Cuu.resize(nG);
      cholmodInvCache.Cxu.resize(nG);
      cholmodInvCache.FS.resize(nG);
      if(blockEngine.chain){
        for(std::size_t block = nG; block-- > 0; ){
          arma::mat inverseLower;
          if(!arma::inv(inverseLower, arma::trimatl(blockEngine.Lg[block]))){
            Rcpp::stop("Block-chain pivot inverse failed.");
          }
          cholmodInvCache.Cuu[block] = inverseLower.t() * inverseLower;
          if(block + 1 < nG && blockEngine.chainComponents[block] == blockEngine.chainComponents[block + 1]){
            arma::mat adjacent = -blockEngine.chainTransfers[block + 1] * cholmodInvCache.Cuu[block + 1];
            cholmodInvCache.Cuu[block] -= adjacent * blockEngine.chainTransfers[block + 1].t();
            cholmodInvCache.crossBlocks.emplace(std::make_pair(static_cast<int>(block), static_cast<int>(block + 1)), std::move(adjacent));
          }
        }
        for(std::size_t block = 0; block < nG; ++block){
          if(haveBorder){
            cholmodInvCache.FS[block] = blockEngine.Fg[block] * Sinv;
            cholmodInvCache.Cuu[block] += cholmodInvCache.FS[block] * blockEngine.Fg[block].t();
            cholmodInvCache.Cxu[block] = -cholmodInvCache.FS[block].t();
          }
        }
        if(haveBorder){
          for(auto & adjacent : cholmodInvCache.crossBlocks){
            adjacent.second += cholmodInvCache.FS[adjacent.first.first] * blockEngine.Fg[adjacent.first.second].t();
          }
        }
        cholmodInvCache.borderPos = blockEngine.borderPos;
        cholmodInvCache.groupOf = blockEngine.groupOf;
        cholmodInvCache.localOf = blockEngine.localOf;
        cholmodInvCache.crossAvailable = true;
        cholmodInvCache.ready = true;
        return;
      }
      for(std::size_t g = 0; g < nG; ++g){
        arma::mat LgInv;
        if(!arma::inv(LgInv, arma::trimatl(blockEngine.Lg[g]))){
          Rcpp::stop("Dense block-Schur inverse of a random-effect block failed.");
        }
        cholmodInvCache.Cuu[g] = LgInv.t() * LgInv;
        if(haveBorder){
          cholmodInvCache.FS[g] = blockEngine.Fg[g] * Sinv;
          cholmodInvCache.Cuu[g] += cholmodInvCache.FS[g] * blockEngine.Fg[g].t();
          cholmodInvCache.Cxu[g] = -cholmodInvCache.FS[g].t();
        }
      }
      if(blockEngine.latent){
        cholmodInvCache.borderPos = blockEngine.borderPos;
        cholmodInvCache.groupOf = blockEngine.groupOf;
        cholmodInvCache.localOf = blockEngine.localOf;
        cholmodInvCache.crossAvailable = haveBorder;
        cholmodInvCache.ready = true;
        return;
      }
      currentCrossGroupDemand.groupSizes.clear();
      currentCrossGroupDemand.groupSizes.reserve(nG);
      for(const std::vector<int> & cols : blockEngine.groups){
        currentCrossGroupDemand.groupSizes.push_back(cols.size());
      }
      currentCrossGroupDemand.groupSignature = crossGroupSignature(blockEngine.groups);
      currentCrossGroupDemand.entries.clear();
      if(currentResidualParameterCount > 1 && nG <= crossGroupCacheMaxGroups){
        currentCrossGroupDemand.entries.assign(nG * nG, 0);
      }
      if(haveBorder && nG <= crossGroupCacheMaxGroups &&
         previousCrossGroupDemand.groupSignature == currentCrossGroupDemand.groupSignature &&
         previousCrossGroupDemand.groupSizes == currentCrossGroupDemand.groupSizes &&
         previousCrossGroupDemand.entries.size() == nG * nG){
        for(std::size_t g = 0; g < nG; ++g){
          for(std::size_t h = g + 1; h < nG; ++h){
            const double rows = static_cast<double>(blockEngine.groups[g].size());
            const double cols = static_cast<double>(blockEngine.groups[h].size());
            const double cells = rows * cols;
            const std::size_t requests =
              previousCrossGroupDemand.entries[g * nG + h] +
              previousCrossGroupDemand.entries[h * nG + g];
            const double materializeFlops =
              2.0 * rows * cols * static_cast<double>(blockEngine.borderCols.size());
            const double dotFlops =
              2.0 * static_cast<double>(requests) *
              static_cast<double>(blockEngine.borderCols.size());
            if(!sommer::WorkspacePolicy::crossBlockWorthCaching(dotFlops, materializeFlops) ||
               cholmodInvCache.crossBlockDoubles + cells > crossGroupCacheMaxDoubles){
              continue;
            }
            cholmodInvCache.crossBlocks.emplace(
              std::make_pair(static_cast<int>(g), static_cast<int>(h)),
              cholmodInvCache.FS[g] * blockEngine.Fg[h].t()
            );
            cholmodInvCache.crossBlockDoubles += cells;
            if(verbose && !reportedCrossGroupCache){
              Rcpp::Rcout
                << "Materialized cross-group inverse block ("
                << blockEngine.groups[g].size() << " x "
                << blockEngine.groups[h].size()
                << ") from prior-iteration trace demand."
                << arma::endl;
              reportedCrossGroupCache = true;
            }
          }
        }
      }
      cholmodInvCache.borderPos = blockEngine.borderPos;
      cholmodInvCache.groupOf = blockEngine.groupOf;
      cholmodInvCache.localOf = blockEngine.localOf;
      cholmodInvCache.crossAvailable = haveBorder;
      cholmodInvCache.ready = true;
      return;
    }

    const cholmod_factor * L = cholmodState.factor;
    if(L == nullptr || !L->is_super || !L->is_ll ||
       L->xtype != CHOLMOD_REAL || L->dtype != CHOLMOD_DOUBLE ||
       L->itype != CHOLMOD_INT || static_cast<int>(L->n) != nEffects){
      return;
    }
    const int n = nEffects;
    const int * perm = static_cast<const int *>(L->Perm);
    std::vector<int> iperm(static_cast<std::size_t>(n));
    for(int k = 0; k < n; ++k){
      iperm[static_cast<std::size_t>(perm[k])] = k;
    }

    double largestFront = 0.0;
    const int * supernodeRows = static_cast<const int *>(L->pi);
    for(std::size_t node = 0; node < L->nsuper; ++node){
      largestFront = std::max(largestFront,
        static_cast<double>(supernodeRows[node + 1] - supernodeRows[node]));
    }
    const double patternSelectedDoubles = static_cast<double>(L->xsize) +
      4.0 * largestFront * largestFront;
    if(sommer::InverseCachePolicy::preferSelected(
        patternLocal, patternSelectedDoubles, cholmodInvCacheMaxDoubles)){
      if(buildCholmodSelectedInverse(L, iperm, cholmodInvCache)){
        cholmodInvCache.patternOnly = true;
        cholmodInvCache.ready = true;
        ++patternInverseCacheBuilds;
        return;
      }
      ++patternInverseCacheFallbacks;
    }

    cholmodInvCache.borderPos.assign(static_cast<std::size_t>(nEffects), -1);
    cholmodInvCache.groupOf.assign(static_cast<std::size_t>(nEffects), -1);
    cholmodInvCache.localOf.assign(static_cast<std::size_t>(nEffects), -1);
    arma::mat YX;
    double retainedDoubles = 0.0;
    const double fixedPeakDoubles = 2.0 * static_cast<double>(n) * nX +
      static_cast<double>(nX) * nX;
    if(nX > 0 && fixedPeakDoubles <= cholmodInvCacheMaxDoubles){
      std::vector<int> xcols(static_cast<std::size_t>(nX));
      for(int j = 0; j < nX; ++j){ xcols[static_cast<std::size_t>(j)] = j; }
      arma::uvec rowsX;
      cholmodForwardSolveUnit(L, iperm, xcols, YX, rowsX);
      const arma::mat YXT = YX.rows(rowsX);
      cholmodInvCache.Cxx = YXT.t() * YXT;
      retainedDoubles = static_cast<double>(nX) * nX;
      for(int j = 0; j < nX; ++j){ cholmodInvCache.borderPos[static_cast<std::size_t>(j)] = j; }
    }

    bool skippedGroup = false;
    // A sparse group is cheaper through the sparse selected inverse than a dense solve cache.
    std::vector<char> inCacheGroup(static_cast<std::size_t>(n), 0);
    auto groupIsSparse = [&](const std::vector<int> & cols) -> bool {
      if(cols.size() <= static_cast<std::size_t>(blockSchurPolicy.smallGroup)){ return false; }
      for(const int column : cols){ inCacheGroup[static_cast<std::size_t>(column)] = 1; }
      double entries = 0.0;
      for(const int column : cols){
        for(EigenSpMat::InnerIterator entry(C, column); entry; ++entry){
          if(inCacheGroup[static_cast<std::size_t>(entry.row())]){ ++entries; }
        }
      }
      for(const int column : cols){ inCacheGroup[static_cast<std::size_t>(column)] = 0; }
      return !blockSchurPolicy.denseEnough(entries, static_cast<double>(cols.size()));
    };
    for(const std::vector<int> & cols : buildRandomInverseGroups()){
      const int m = static_cast<int>(cols.size());
      const double groupSize = static_cast<double>(m);
      const double fixedSize = static_cast<double>(YX.n_cols);
      const double groupPeakDoubles = retainedDoubles +
        2.0 * static_cast<double>(n) * (groupSize + fixedSize) +
        2.0 * groupSize * groupSize + 2.0 * fixedSize * groupSize;
      if(m == 0 ||
        !blockSchurPolicy.inverseCacheFits(static_cast<double>(n), groupSize) ||
        groupPeakDoubles > cholmodInvCacheMaxDoubles ||
        groupIsSparse(cols)){
        skippedGroup = skippedGroup || m > 0;
        continue;
      }
      arma::mat Y;
      arma::uvec rowsT;
      cholmodForwardSolveUnit(L, iperm, cols, Y, rowsT);
      const arma::mat YT = Y.rows(rowsT);
      Y.reset();
      const int g = static_cast<int>(cholmodInvCache.Cuu.size());
      cholmodInvCache.Cuu.push_back(YT.t() * YT);
      if(!YX.is_empty()){
        cholmodInvCache.Cxu.push_back(arma::mat(YX.rows(rowsT).t() * YT));
      }else{
        cholmodInvCache.Cxu.emplace_back();
      }
      retainedDoubles += groupSize * groupSize + fixedSize * groupSize;
      for(int j = 0; j < m; ++j){
        const std::size_t c = static_cast<std::size_t>(cols[static_cast<std::size_t>(j)]);
        cholmodInvCache.groupOf[c] = g;
        cholmodInvCache.localOf[c] = j;
      }
    }
    YX.reset();
    const double selectedPeakDoubles = retainedDoubles + static_cast<double>(L->xsize) +
      4.0 * largestFront * largestFront;
    if(skippedGroup && selectedPeakDoubles <= cholmodInvCacheMaxDoubles){
      buildCholmodSelectedInverse(L, iperm, cholmodInvCache);
    }else if(skippedGroup && memoryStrict){
      requireMemory("The CHOLMOD inverse-block cache", selectedPeakDoubles);
    }
    cholmodInvCache.ready = true;
  };

  // Original-index lookup of C^{-1}(i,j); false if the entry is not cached.
  auto cholmodInverseLookup = [&](int i, int j, double & value) -> bool {
    if(!cholmodInvCache.ready || i < 0 || j < 0 || i >= nEffects || j >= nEffects){
      return false;
    }
    if(cholmodInvCache.selReady && cholmodSelectedLookup(i, j, value)){ return true; }
    if(cholmodInvCache.patternOnly){ return false; }
    const std::size_t si = static_cast<std::size_t>(i);
    const std::size_t sj = static_cast<std::size_t>(j);
    const int bi = cholmodInvCache.borderPos[si];
    const int bj = cholmodInvCache.borderPos[sj];
    if(bi >= 0 && bj >= 0){
      value = cholmodInvCache.Cxx(static_cast<arma::uword>(bi), static_cast<arma::uword>(bj));
      return true;
    }
    if(bi >= 0 || bj >= 0){
      const int b = bi >= 0 ? bi : bj;
      const std::size_t u = bi >= 0 ? sj : si;
      const int g = cholmodInvCache.groupOf[u];
      if(g < 0){ return false; }
      value = cholmodInvCache.Cxu[static_cast<std::size_t>(g)](
        static_cast<arma::uword>(b),
        static_cast<arma::uword>(cholmodInvCache.localOf[u]));
      return true;
    }
    const int g = cholmodInvCache.groupOf[si];
    const int h = cholmodInvCache.groupOf[sj];
    if(g < 0 || h < 0){ return false; }
    if(g != h){
      if(!cholmodInvCache.crossAvailable){ return false; }
      if(!cholmodInvCache.crossBlocks.empty()){
        const int lo = std::min(g, h);
        const int hi = std::max(g, h);
        const auto cached = cholmodInvCache.crossBlocks.find(std::make_pair(lo, hi));
        if(cached != cholmodInvCache.crossBlocks.end()){
          const arma::uword ri = static_cast<arma::uword>(cholmodInvCache.localOf[si]);
          const arma::uword rj = static_cast<arma::uword>(cholmodInvCache.localOf[sj]);
          value = g < h ? cached->second(ri, rj) : cached->second(rj, ri);
          return true;
        }
      }
      if(blockEngine.chain){
        if(blockEngine.chainComponents[g] == blockEngine.chainComponents[h]){ return false; }
        if(blockEngine.borderCols.empty()){ value = 0.0; return true; }
      }
      value = arma::dot(
        cholmodInvCache.FS[static_cast<std::size_t>(g)].row(static_cast<arma::uword>(cholmodInvCache.localOf[si])),
        blockEngine.Fg[static_cast<std::size_t>(h)].row(static_cast<arma::uword>(cholmodInvCache.localOf[sj])));
      return true;
    }
    value = cholmodInvCache.Cuu[static_cast<std::size_t>(g)](
      static_cast<arma::uword>(cholmodInvCache.localOf[si]),
      static_cast<arma::uword>(cholmodInvCache.localOf[sj]));
    return true;
  };

  bool pcgNystromReady = false;
  int pcgNystromRankUsed = 0;
  Eigen::VectorXd pcgPreconditionerDiagonal;
  Eigen::MatrixXd pcgNystromBasis;
  Eigen::MatrixXd pcgNystromWeightedBasis;
  Eigen::LLT<Eigen::MatrixXd> pcgNystromCoreFactor;
  std::vector<std::vector<int>> pcgBlockIndices;
  std::vector<Eigen::LLT<Eigen::MatrixXd>> pcgBlockFactors;
  std::atomic<std::size_t> pcgOperatorCalls(0), pcgVectorEquivalent(0);
  auto applyCurrentPCG = [&](const Eigen::MatrixXd & values) -> Eigen::MatrixXd {
    ++pcgOperatorCalls;
    pcgVectorEquivalent += static_cast<std::size_t>(values.cols());
    return matrixFreePCGReady ? matrixFreeApplyCMatrix(values) : Eigen::MatrixXd(C * values);
  };
  auto applyPCGPreconditioner = [&](const Eigen::Ref<const Eigen::MatrixXd> & rhs)
      -> Eigen::MatrixXd {
    if(!pcgBlockFactors.empty()){
      Eigen::MatrixXd result = Eigen::MatrixXd::Zero(rhs.rows(), rhs.cols());
      for(std::size_t block = 0; block < pcgBlockIndices.size(); ++block){
        const std::vector<int> & indices = pcgBlockIndices[block];
        Eigen::MatrixXd local(indices.size(), rhs.cols());
        for(std::size_t row = 0; row < indices.size(); ++row){
          local.row(row) = rhs.row(indices[row]);
        }
        const Eigen::MatrixXd solved = pcgBlockFactors[block].solve(local);
        if(!solved.allFinite()){ Rcpp::stop("PCG block-Jacobi solve is non-finite."); }
        for(std::size_t row = 0; row < indices.size(); ++row){
          result.row(indices[row]) = solved.row(row);
        }
      }
      return result;
    }
    if(!pcgNystromReady){
      Eigen::MatrixXd result = rhs;
      for(Eigen::Index col = 0; col < result.cols(); ++col){
        result.col(col).array() /= pcgPreconditionerDiagonal.array();
      }
      return result;
    }
    Eigen::MatrixXd scaled = rhs;
    for(Eigen::Index col = 0; col < scaled.cols(); ++col){
      scaled.col(col).array() /= pcgPreconditionerDiagonal.array();
    }
    const Eigen::MatrixXd projected = pcgNystromBasis.transpose() * scaled;
    const Eigen::MatrixXd correction = pcgNystromCoreFactor.solve(projected);
    if(pcgNystromCoreFactor.info() != Eigen::Success || !correction.allFinite()){
      Rcpp::stop("Nyström PCG preconditioner core solve failed.");
    }
    scaled.noalias() -= pcgNystromWeightedBasis * correction;
    return scaled;
  };

  auto preparePCG = [&](const arma::field<arma::sp_mat> & lambda){
    const Eigen::Index n = C.rows();
    pcgBlockIndices.clear();
    pcgBlockFactors.clear();
    pcgPreconditionerDiagonal = matrixFreePCGReady ? matrixFreePCGDiagonal : Eigen::VectorXd(C.diagonal());
    if(n == 0 || !pcgPreconditionerDiagonal.allFinite() ||
       pcgPreconditionerDiagonal.minCoeff() <= 0.0){
      Rcpp::stop("PCG requires a positive finite coefficient diagonal.");
    }
    CpcgReady = true;
    pcgNystromReady = false;
    pcgNystromRankUsed = 0;
    if(pcgPreconditioner == "block-jacobi"){
      double blockStorage = 0.0;
      auto addBlock = [&](const std::vector<int> & indices, int random, int level,
                          int coordinateStart){
        const int size = static_cast<int>(indices.size());
        blockStorage += static_cast<double>(size) * size;
        const double needed = blockStorage + 3.0 * size * size;
        if(needed > blockSchurPolicy.maxDoubles){
          requireMemory("PCG block-Jacobi factors", needed);
        }
        Eigen::MatrixXd block(size, size);
        for(int col = 0; col < size; ++col){
          for(int row = 0; row < size; ++row){
            block(row,col) = C.coeff(indices[row], indices[col]);
            if(matrixFreePCGReady && random >= 0){
              block(row,col) += lambda(random)(coordinateStart + row, coordinateStart + col) *
                Ai(random)(level, level);
            }
          }
        }
        Eigen::LLT<Eigen::MatrixXd> factor(0.5 * (block + block.transpose()));
        if(factor.info() != Eigen::Success){
          Rcpp::stop("PCG block-Jacobi principal block is not positive definite.");
        }
        pcgBlockIndices.push_back(indices);
        pcgBlockFactors.push_back(factor);
      };
      for(int first = 0; first < nX; first += pcgBlockSize){
        std::vector<int> indices;
        for(int j = first; j < std::min(nX, first + pcgBlockSize); ++j){ indices.push_back(j); }
        addBlock(indices, -1, 0, 0);
      }
      for(int random = 0; random < nRe; ++random){
        const std::vector<arma::uvec> & partitions = partitionIndexCache[random];
        for(arma::uword level = 0; level < Ai(random).n_rows; ++level){
          for(int first = 0; first < static_cast<int>(partitions.size()); first += pcgBlockSize){
            std::vector<int> indices;
            for(int j = first; j < std::min<int>(partitions.size(), first + pcgBlockSize); ++j){
              indices.push_back(partitions[j](level));
            }
            addBlock(indices, random, static_cast<int>(level), first);
          }
        }
      }
      return;
    }
    if(pcgPreconditioner == "nystrom"){
      const Eigen::VectorXd diagonal = pcgPreconditionerDiagonal;
      const int rank = std::min<int>(pcgNystromRank, static_cast<int>(n));
      const double needed = 4.0 * static_cast<double>(n) * rank + 5.0 * rank * rank;
      if(needed > blockSchurPolicy.maxDoubles){ requireMemory("PCG Nystrom workspace", needed); }
      std::vector<int> landmarks(static_cast<std::size_t>(rank));
      for(int landmark = 0; landmark < rank; ++landmark){
        landmarks[static_cast<std::size_t>(landmark)] = static_cast<int>(
          (static_cast<long long>(2 * landmark + 1) * static_cast<long long>(n)) /
          static_cast<long long>(2 * rank));
      }
      Eigen::MatrixXd unit = Eigen::MatrixXd::Zero(n, rank);
      for(int local = 0; local < rank; ++local){
        const int source = landmarks[static_cast<std::size_t>(local)];
        unit(source,local) = 1.0;
      }
      const Eigen::MatrixXd columns = applyCurrentPCG(unit);
      Eigen::MatrixXd core(rank, rank);
      for(int col = 0; col < rank; ++col){
        for(int row = 0; row < rank; ++row){
          core(row, col) = columns(landmarks[static_cast<std::size_t>(row)], col);
        }
      }
      core = 0.5 * (core + core.transpose());
      Eigen::LLT<Eigen::MatrixXd> landmarkFactor(core);
      if(landmarkFactor.info() != Eigen::Success){
        Rcpp::stop("Nyström landmark submatrix is not positive definite; reduce pcgNystromRank or use diagonal preconditioning.");
      }
      pcgNystromBasis = landmarkFactor.matrixL().solve(columns.transpose()).transpose();
      Eigen::VectorXd lowRankDiagonal = pcgNystromBasis.rowwise().squaredNorm();
      const double diagonalFloor = 1.0e-4 * diagonal.maxCoeff();
      pcgPreconditionerDiagonal = (diagonal - lowRankDiagonal).cwiseMax(diagonalFloor);
      if(!pcgPreconditionerDiagonal.allFinite() || pcgPreconditionerDiagonal.minCoeff() <= 0.0){
        Rcpp::stop("Nyström diagonal remainder is not positive; reduce pcgNystromRank or use diagonal preconditioning.");
      }
      pcgNystromWeightedBasis = pcgNystromBasis;
      for(Eigen::Index col = 0; col < pcgNystromWeightedBasis.cols(); ++col){
        pcgNystromWeightedBasis.col(col).array() /= pcgPreconditionerDiagonal.array();
      }
      Eigen::MatrixXd preconditionerCore =
        Eigen::MatrixXd::Identity(rank, rank) + pcgNystromBasis.transpose() * pcgNystromWeightedBasis;
      pcgNystromCoreFactor.compute(0.5 * (preconditionerCore + preconditionerCore.transpose()));
      if(pcgNystromCoreFactor.info() != Eigen::Success){
        Rcpp::stop("Nyström Woodbury core is not positive definite.");
      }
      pcgNystromReady = true;
      pcgNystromRankUsed = rank;
      CpcgReady = true;
      return;
    }
  };

  double pcgSolveMaxRelativeResidual = 0.0;
  std::function<Eigen::MatrixXd(const Eigen::Ref<const Eigen::MatrixXd> &,
    const Eigen::MatrixXd *, const std::string &)> solvePCGBatchFn;
  auto solveCVector = [&](const Eigen::Ref<const Eigen::VectorXd> & rhs,
                          const std::string & context) -> Eigen::VectorXd {
    if(saddleActive){ return saddleJoint.solveCoefficients(Eigen::MatrixXd(rhs), context).col(0); }
    if(dFirstMLActive){
      const Eigen::MatrixXd matrixRhs = rhs;
      return solveDfirstJoint(matrixRhs, context).col(0);
    }
    if(solverName == "ldlt"){
      Eigen::VectorXd ans = Cfactor.solve(rhs);
      if(Cfactor.info() != Eigen::Success){
        Rcpp::stop("Sparse LDLT solve failed in " + context + ".");
      }
      return ans;
    }
    if(solverName == "cholmod"){
      return solveCVectorCholmod(rhs, context);
    }
    if(!CpcgReady && !matrixFreePCGReady){
      Rcpp::stop("PCG solve requested before the PCG backend was prepared.");
    }
    Eigen::MatrixXd rhsMatrix(rhs.size(), 1);
    rhsMatrix.col(0) = rhs;
    return solvePCGBatchFn(rhsMatrix, nullptr, context).col(0);
  };

  int pcgBatchSolves = 0;
  int pcgBatchIterations = 0;
  sommer::pcg::Metrics pcgMetrics;
  const bool compressedData = compressionMethod != "off";
  sommer::SufficientStatistics sufficient;
  Eigen::VectorXd compressedVariance, compressedFitted, compressedSquares;
  Eigen::MatrixXd compressedDerivative;
  if(compressedData){
    if(compressionMethod != "groups" && compressionMethod != "classes"){
      Rcpp::stop("compression must be 'off', 'groups', or 'classes'.");
    }
    if(!residualStructurallyDiagonal || (useH && !Hdiag) ||
       compactResidualCovariance || saddleRequested || rrConditional || rrCallback){
      Rcpp::stop("Compression requires independent residuals and fixed diagonal weights without compact residual, saddle-point or latent reduced-rank modes.");
    }
    const auto compressionStart = std::chrono::steady_clock::now();
    std::vector<int> classes(static_cast<std::size_t>(nR));
    Eigen::VectorXd precision = Eigen::VectorXd::Ones(nR);
    for(int row = 0; row < nR; ++row){
      classes[row] = residualLocalOfRow[row];
      if(useH){ precision(row) = HdiagSquared(row); }
    }
    sufficient.build(W, Eigen::Map<const Eigen::VectorXd>(yDense.memptr(), nR),
      classes, precision, residualDim, compressionMethod == "classes",
      static_cast<std::size_t>(std::min<long double>(
        blockSchurPolicy.maxDoubles * sizeof(double),
        std::numeric_limits<std::size_t>::max())));
    sufficient.setupSeconds = std::chrono::duration<double>(
      std::chrono::steady_clock::now() - compressionStart).count();
  }

  std::map<std::string, Eigen::MatrixXd> pcgAcceptedSolutions, pcgCandidateSolutions;
  const bool krylovPcg = solverName == "pcg" && pcgRecycleMethod == "krylov" && pcgRecycleRank > 0;
  const int krylovRank = std::max(1, std::min(nEffects, pcgRecycleRank));
  const int krylovCandidates = static_cast<int>(std::min<long long>(nEffects, 4LL * krylovRank));
  sommer::pcg::RecycleSpace pcgKrylov(std::max(1, krylovCandidates), krylovRank);
  sommer::pcg::RecycleSpace::Selection pcgPendingKrylov;
  bool pcgPendingKrylovReady = false;
  std::size_t pcgBalancedApplications = 0;
  if(krylovPcg){
    const double workspace = 32.0 * nEffects * krylovRank + 200.0 * krylovRank * krylovRank;
    if(workspace > blockSchurPolicy.maxDoubles){ requireMemory("PCG Krylov recycle space", workspace); }
  }
  Eigen::MatrixXd pcgImportedBasis;
  int pcgImportedRank = 0;
  if(covStructI.hasAttribute("pcgControl")){
    const Rcpp::List control = Rcpp::as<Rcpp::List>(covStructI.attr("pcgControl"));
    if(control.containsElementNamed("recycleState") &&
       !Rf_isNull(control["recycleState"])){
      if(!krylovPcg){ Rcpp::stop("PCG state import requires positive-rank Krylov recycling."); }
      const Rcpp::List state = Rcpp::as<Rcpp::List>(control["recycleState"]);
      const Rcpp::NumericMatrix basis = Rcpp::as<Rcpp::NumericMatrix>(state["basis"]);
      if(basis.nrow() != nEffects || basis.ncol() < 1 || basis.ncol() > krylovRank){
        Rcpp::stop("PCG imported recycle basis has incompatible dimensions.");
      }
      pcgImportedBasis = Eigen::Map<const Eigen::MatrixXd>(
        basis.begin(), basis.nrow(), basis.ncol());
      if(!pcgImportedBasis.allFinite()){
        Rcpp::stop("PCG imported recycle basis must be finite.");
      }
      pcgImportedRank = basis.ncol();
    }
  }
  if(pcgExportRecycle && !krylovPcg){
    Rcpp::stop("PCG state export requires positive-rank Krylov recycling.");
  }
  bool pcgEvaluationAccepted = false;
  Eigen::MatrixXd pcgRecycleBasis, pcgRecycleApplied;
  Eigen::LLT<Eigen::MatrixXd> pcgRecycleFactor;
  int pcgWarmStarts = 0, pcgCoarseCorrections = 0;
  double pcgSetupSeconds = 0.0, pcgSolveSeconds = 0.0;
  double pcgLogDetSeconds = 0.0, pcgTraceSetupSeconds = 0.0;
  auto solvePCGBatch = [&](const Eigen::Ref<const Eigen::MatrixXd> & rhs,
                           const Eigen::MatrixXd * initialGuess,
                           const std::string & context) -> Eigen::MatrixXd {
    const auto started = std::chrono::steady_clock::now();
    if(!CpcgReady && !matrixFreePCGReady){
      Rcpp::stop("Batched PCG requested without a prepared coefficient operator.");
    }
    if(rhs.rows() != C.rows()){
      Rcpp::stop("Batched PCG right-hand side has incompatible dimensions in " + context + ".");
    }
    ++pcgBatchSolves;
    const Eigen::Index n = rhs.rows();
    const Eigen::Index nRhs = rhs.cols();
    const bool cacheable = context == "BLUE/BLUP calculation" || context == "Hutchinson trace probes";
    Eigen::MatrixXd guess = initialGuess != nullptr &&
      initialGuess->rows() == n && initialGuess->cols() == nRhs
      ? *initialGuess : Eigen::MatrixXd::Zero(n, nRhs);
    const auto cached = pcgAcceptedSolutions.find(context);
    if(pcgWarmStart && cacheable && cached != pcgAcceptedSolutions.end() &&
       cached->second.rows() == n && cached->second.cols() == nRhs){
      guess = cached->second;
      ++pcgWarmStarts;
    }
    if(pcgRecycleBasis.cols() > 0){
      const Eigen::MatrixXd residual = rhs - applyCurrentPCG(guess);
      guess.noalias() += pcgRecycleBasis *
        pcgRecycleFactor.solve(pcgRecycleBasis.transpose() * residual);
      ++pcgCoarseCorrections;
    }
    if(krylovPcg && pcgKrylov.prepared()){
      const Eigen::MatrixXd residual = rhs - applyCurrentPCG(guess);
      guess += pcgKrylov.coarseCorrectionPanel(residual);
      ++pcgCoarseCorrections;
    }
    sommer::pcg::Options options;
    options.tolerance = pcgTol;
    options.maxIterations = pcgMaxIters;
    options.residualRefresh = pcgResidualRefresh;
    options.blockSize = std::min<std::size_t>(pcgBlockSize, sommer::WorkspacePolicy::panelColumns());
    options.maxWorkspaceBytes = static_cast<std::size_t>(std::min<double>(
      blockSchurPolicy.maxDoubles * sizeof(double),
      static_cast<double>(std::numeric_limits<std::size_t>::max() / 2)));
    const double workspace = static_cast<double>(n) *
      (4.0 * nRhs + 12.0 * std::min<Eigen::Index>(nRhs, options.blockSize));
    if(workspace > blockSchurPolicy.maxDoubles){ requireMemory("PCG solve panels", workspace); }
    auto applyVector = [&](const Eigen::VectorXd & value) -> Eigen::VectorXd {
      Rcpp::checkUserInterrupt();
      return applyCurrentPCG(Eigen::MatrixXd(value)).col(0);
    };
    auto preconditionVector = [&](const Eigen::VectorXd & value) -> Eigen::VectorXd {
      return applyPCGPreconditioner(Eigen::MatrixXd(value)).col(0);
    };
    auto applyPanel = [&](const Eigen::MatrixXd & values) -> Eigen::MatrixXd {
      Rcpp::checkUserInterrupt();
      return applyCurrentPCG(values);
    };
    if(krylovPcg){ options.recycleSpace = &pcgKrylov; }
    const sommer::pcg::Result result = pcgRhsMethod == "block"
      ? sommer::pcg::solveBlock(rhs, applyPanel, applyPCGPreconditioner,
          applyVector, preconditionVector, options, &guess)
      : sommer::pcg::solvePackedIndependent(rhs, applyPanel, applyPCGPreconditioner, options, &guess);
    if(krylovPcg && pcgKrylov.prepared()){
      pcgBalancedApplications += result.metrics.preconditionerApplications;
    }
    pcgMetrics.add(result.metrics);
    pcgBatchIterations += result.metrics.iterations;
    const Eigen::MatrixXd & solution = result.solution;
    const Eigen::MatrixXd trueResidual = rhs - applyCurrentPCG(solution);
    for(Eigen::Index col = 0; col < nRhs; ++col){
      if(rhs.col(col).norm() > 0.0){
        const double relativeResidual = trueResidual.col(col).norm() / rhs.col(col).norm();
        pcgSolveMaxRelativeResidual = std::max(pcgSolveMaxRelativeResidual, relativeResidual);
        if(relativeResidual > pcgTol){
          Rcpp::stop("Batched PCG failed the true-residual check for RHS column " +
            std::to_string(static_cast<long long>(col + 1)) + " in " + context + ".");
        }
      }
    }
    if(!solution.allFinite()){
      Rcpp::stop("Batched PCG produced a non-finite solution in " + context + ".");
    }
    if(cacheable && (pcgWarmStart || (pcgRecycleRank > 0 && !krylovPcg))){
      pcgCandidateSolutions[context] = solution;
    }
    pcgSolveSeconds += std::chrono::duration<double>(
      std::chrono::steady_clock::now() - started).count();
    return solution;
  };
  solvePCGBatchFn = solvePCGBatch;

  auto solveCMatrix = [&](const Eigen::Ref<const Eigen::MatrixXd> & rhs,
                          const std::string & context) -> Eigen::MatrixXd {
    if(saddleActive){ return saddleJoint.solveCoefficients(rhs, context); }
    if(dFirstMLActive){ return solveDfirstJoint(rhs, context); }
    Eigen::MatrixXd ans(rhs.rows(), rhs.cols());
    if(solverName == "ldlt"){
      ans = Cfactor.solve(rhs);
      if(Cfactor.info() != Eigen::Success){
        Rcpp::stop("Sparse LDLT multi-RHS solve failed in " + context + ".");
      }
      return ans;
    }
    if(solverName == "cholmod"){
      return solveCMatrixCholmod(rhs, context);
    }
    if(matrixFreePCGReady){
      return solvePCGBatch(rhs, nullptr, context);
    }
    if(!CpcgReady){
      Rcpp::stop("PCG solve requested before the PCG backend was prepared.");
    }
    return solvePCGBatch(rhs, nullptr, context);
  };

  // Computes tr(B * factor^{-1}) using the selected-inverse fast path when
  // available (LDLT only). genericSolve, when non-empty, routes the
  // fallback through it instead of factor.solve() directly - required for
  // a solver="cholmod" factor (its supernodal representation has no
  // Takahashi subset, so `factor`/`subset` are never populated for it).
  auto sparseTraceInverseTimes =
    [&](const arma::sp_mat & B,
        const EigenLDLT & factor,
        const SelectedInverseSubset & subset,
        const std::string & context,
        bool & usedFallback,
        const std::function<Eigen::MatrixXd(const Eigen::Ref<const Eigen::MatrixXd> &, const std::string &)> & genericSolve = nullptr,
        const bool useCholmodInverseCache = false) -> double {
      usedFallback = false;
      double traceValue = 0.0;
      bool allAvailable = true;
      const std::size_t nGroups = currentCrossGroupDemand.groupSizes.size();
      const bool trackCrossDemand =
        useCholmodInverseCache && blockEngine.active && !blockEngine.chain &&
        !blockEngine.banded && !blockEngine.latent && nGroups > 0 &&
        nGroups <= crossGroupCacheMaxGroups &&
        currentCrossGroupDemand.groupSignature == blockEngine.groupSignature &&
        currentCrossGroupDemand.entries.size() == nGroups * nGroups;
      for(arma::sp_mat::const_iterator it = B.begin(); it != B.end(); ++it){
        double zij = 0.0;
        const bool have = useCholmodInverseCache
          ? cholmodInverseLookup(static_cast<int>(it.row()), static_cast<int>(it.col()), zij)
          : getSelectedInverseOriginal(subset, static_cast<int>(it.col()), static_cast<int>(it.row()), zij);
        if(trackCrossDemand){
          const int g = cholmodInvCache.groupOf[static_cast<std::size_t>(it.row())];
          const int h = cholmodInvCache.groupOf[static_cast<std::size_t>(it.col())];
          if(g >= 0 && h >= 0 && g != h){
            const std::size_t pairIndex =
              static_cast<std::size_t>(g) * nGroups + static_cast<std::size_t>(h);
            ++currentCrossGroupDemand.entries[pairIndex];
          }
        }
        if(!have){
          allAvailable = false;
          break;
        }
        traceValue += (*it) * zij;
      }
      if(allAvailable){ return traceValue; }

      usedFallback = true;
      traceValue = 0.0;
      std::vector<arma::uword> activeColumns;
      activeColumns.reserve(static_cast<std::size_t>(B.n_cols));
      for(arma::uword j = 0; j < B.n_cols; ++j){
        if(B.begin_col(j) != B.end_col(j)){ activeColumns.push_back(j); }
      }
      const std::size_t batchSize = sommer::WorkspacePolicy::panelColumns();
      for(std::size_t batchStart = 0; batchStart < activeColumns.size(); batchStart += batchSize){
        const std::size_t batchEnd = std::min(batchStart + batchSize, activeColumns.size());
        const std::size_t currentBatchSize = batchEnd - batchStart;
        Eigen::MatrixXd rhs = Eigen::MatrixXd::Zero(static_cast<Eigen::Index>(B.n_rows), static_cast<Eigen::Index>(currentBatchSize));
        for(std::size_t k = 0; k < currentBatchSize; ++k){
          const arma::uword sourceCol = activeColumns[batchStart + k];
          for(arma::sp_mat::const_col_iterator it = B.begin_col(sourceCol); it != B.end_col(sourceCol); ++it){
            rhs(static_cast<Eigen::Index>(it.row()), static_cast<Eigen::Index>(k)) = (*it);
          }
        }
        Eigen::MatrixXd solution;
        if(genericSolve){
          solution = genericSolve(rhs, context);
        }else{
          solution = factor.solve(rhs);
          if(factor.info() != Eigen::Success){ Rcpp::stop("Sparse LDLT trace fallback solve failed in " + context + "."); }
        }
        for(std::size_t k = 0; k < currentBatchSize; ++k){
          const arma::uword sourceCol = activeColumns[batchStart + k];
          traceValue += solution(static_cast<Eigen::Index>(sourceCol), static_cast<Eigen::Index>(k));
        }
      }
      return traceValue;
    };

  // ============================================================
  // Factorisation-free PCG trace/log-determinant helpers for C
  // ============================================================
  // Deterministic Rademacher signs give common random numbers across REML
  // iterations.  This is important for the likelihood line search: changes
  // in the approximate log determinant then reflect C, not fresh Monte Carlo
  // noise at every iteration.
  auto pcgProbeSign = [pcgSeed](const unsigned long long row,
                         const unsigned long long probe) -> double {
    unsigned long long x = row + 0x9e3779b97f4a7c15ULL * (probe + 1ULL) +
      0xd1b54a32d192ed03ULL * static_cast<unsigned long long>(pcgSeed);
    x ^= x >> 30;
    x *= 0xbf58476d1ce4e5b9ULL;
    x ^= x >> 27;
    x *= 0x94d049bb133111ebULL;
    x ^= x >> 31;
    return (x & 1ULL) ? 1.0 : -1.0;
  };

  Eigen::MatrixXd pcgLanczosInverseGuess;
  double pcgLogDetProbeStandardError = std::numeric_limits<double>::quiet_NaN();
  std::vector<double> pcgLogDetProbeValues;
  double pcgLogDetDepthDifference = std::numeric_limits<double>::quiet_NaN();
  Eigen::MatrixXd pcgScoreProbeValues;
  double pcgScoreMaxSE = std::numeric_limits<double>::quiet_NaN();
  int pcgAccuracyRefinements = 0;
  Eigen::VectorXd pcgLastSparseTraceSamples;
  std::vector<double> pcgAcceptedLogDetValues;
  double pcgPairedLogLikSE = std::numeric_limits<double>::quiet_NaN();

  auto pcgApproxLogDet =
    [&](const Eigen::Index n,
        const std::function<Eigen::VectorXd(const Eigen::VectorXd &)> & rawOperator) -> double {
    if(n <= 0){ return 0.0; }
    const auto started = std::chrono::steady_clock::now();
    int workers = 1;
#ifdef _OPENMP
    workers = std::min(pcgTraceProbes, omp_get_max_threads());
#endif
    const double needed = 3.0 * static_cast<double>(n) *
      (std::min<Eigen::Index>(pcgLanczosSteps, n) * static_cast<double>(workers) +
       pcgTraceProbes + std::min<Eigen::Index>(n, pcgRecycleRank) +
       (pcgTraceEstimator == "control-variate" ? std::min<Eigen::Index>(n, pcgTraceRank) : 0));
    if(needed > blockSchurPolicy.maxDoubles){ requireMemory("PCG Lanczos and probe workspace", needed); }
    const Eigen::VectorXd diagonal = matrixFreePCGReady ? matrixFreePCGDiagonal : Eigen::VectorXd(C.diagonal());
    const bool scaled = pcgLogDetScaling == "diagonal";
    Eigen::VectorXd root = Eigen::VectorXd::Ones(n);
    if(scaled){ root = diagonal.array().sqrt().matrix(); }
    const double correction = scaled ? diagonal.array().log().sum() : 0.0;
    auto applyC = [&](const Eigen::VectorXd & value) -> Eigen::VectorXd {
      ++pcgOperatorCalls;
      ++pcgVectorEquivalent;
      return scaled ? rawOperator(value.cwiseQuotient(root)).cwiseQuotient(root).eval() : rawOperator(value);
    };
    const int mMax = std::min<int>(pcgLanczosSteps, static_cast<int>(n));
    std::vector<double> probeEstimates(
      static_cast<std::size_t>(pcgTraceProbes),
      0.0
    );
    pcgLanczosInverseGuess.resize(n, pcgTraceProbes);
    std::vector<double> shallowEstimates(static_cast<std::size_t>(pcgTraceProbes), 0.0);
    std::atomic<bool> invalidProbe(false);

    // Each deterministic Rademacher probe is independent.  Keeping one
    // accumulator per probe makes the final reduction reproducible across
    // OpenMP schedules.
    #pragma omp parallel for if(pcgTraceProbes > 1)
    for(int probe = 0; probe < pcgTraceProbes; ++probe){
      Eigen::VectorXd q(n);
      for(Eigen::Index i = 0; i < n; ++i){
        q(i) = pcgProbeSign(static_cast<unsigned long long>(i),
                           static_cast<unsigned long long>(probe));
      }
      const double normz = q.norm();
      q /= normz;
      Eigen::MatrixXd lanczosBasis(n, mMax);

      std::vector<double> alpha;
      alpha.reserve(static_cast<std::size_t>(mMax));
      Eigen::MatrixXd projected = Eigen::MatrixXd::Zero(mMax, mMax);

      for(int j = 0; j < mMax; ++j){
        lanczosBasis.col(j) = q;
        Eigen::VectorXd w = applyC(q);
        const Eigen::MatrixXd basis = lanczosBasis.leftCols(j + 1);
        Eigen::VectorXd coefficients = basis.transpose() * w;
        w.noalias() -= basis * coefficients;
        const Eigen::VectorXd correction = basis.transpose() * w;
        w.noalias() -= basis * correction;
        coefficients += correction;
        for(int row = 0; row <= j; ++row){
          projected(row, j) = coefficients(row);
          projected(j, row) = coefficients(row);
        }
        const double a = coefficients(j);
        const double b = w.norm();
        alpha.push_back(a);
        if(j + 1 >= mMax || b <= 1.0e-14 * std::max(1.0, std::abs(a))){ break; }
        q = w / b;
      }

      const int m = static_cast<int>(alpha.size());
      const Eigen::MatrixXd T = 0.5 *
        (projected.topLeftCorner(m, m) + projected.topLeftCorner(m, m).transpose());
      Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(T);
      if(es.info() != Eigen::Success){
        invalidProbe.store(true, std::memory_order_relaxed);
        continue;
      }
      const Eigen::VectorXd eval = es.eigenvalues();
      const Eigen::MatrixXd evec = es.eigenvectors();
      double quad = 0.0;
      Eigen::VectorXd inverseCoefficients = Eigen::VectorXd::Zero(m);
      for(int j = 0; j < m; ++j){
        if(!std::isfinite(eval(j)) || eval(j) <= 0.0){
          invalidProbe.store(true, std::memory_order_relaxed);
          quad = 0.0;
          break;
        }
        const double w0 = evec(0,j);
        quad += w0*w0*std::log(eval(j));
        inverseCoefficients.noalias() +=
          (normz * w0 / eval(j)) * evec.col(j);
      }
      probeEstimates[static_cast<std::size_t>(probe)] = normz*normz*quad;
      if(pcgAccuracy == "adaptive"){
        if(m == n){
          shallowEstimates[probe] = probeEstimates[probe];
        }else{
          const int shallow = std::min(m, std::max(2, mMax / 2));
          Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> nested(T.topLeftCorner(shallow, shallow));
          if(nested.info() != Eigen::Success || nested.eigenvalues().minCoeff() <= 0.0){
            invalidProbe.store(true, std::memory_order_relaxed);
            continue;
          }
          shallowEstimates[probe] = normz * normz *
            (nested.eigenvectors().row(0).array().square() *
             nested.eigenvalues().transpose().array().log()).sum();
        }
      }
      pcgLanczosInverseGuess.col(probe).noalias() =
        lanczosBasis.leftCols(m) * inverseCoefficients;
    }

    if(invalidProbe.load(std::memory_order_relaxed)){
      Rcpp::stop("PCG/SLQ failed while estimating log|C|; C may not be positive definite or Lanczos accuracy is insufficient.");
    }

    double total = 0.0;
    for(int probe = 0; probe < pcgTraceProbes; ++probe){
      total += probeEstimates[static_cast<std::size_t>(probe)];
    }
    pcgLogDetProbeValues = probeEstimates;
    for(double & estimate : pcgLogDetProbeValues){ estimate += correction; }
    if(pcgAccuracy == "adaptive"){
      const double shallowMean = std::accumulate(shallowEstimates.begin(), shallowEstimates.end(), 0.0) /
        static_cast<double>(pcgTraceProbes);
      pcgLogDetDepthDifference = std::abs(total / pcgTraceProbes - shallowMean);
    }
    if(pcgTraceProbes > 1){
      const double mean = total / static_cast<double>(pcgTraceProbes);
      double squaredDeviation = 0.0;
      for(const double estimate : probeEstimates){
        const double difference = estimate - mean;
        squaredDeviation += difference * difference;
      }
      pcgLogDetProbeStandardError = std::sqrt(
        squaredDeviation /
        (static_cast<double>(pcgTraceProbes - 1) * pcgTraceProbes));
    }else{
      pcgLogDetProbeStandardError = std::numeric_limits<double>::quiet_NaN();
    }
    pcgLogDetSeconds += std::chrono::duration<double>(
      std::chrono::steady_clock::now() - started).count();
    return total / static_cast<double>(pcgTraceProbes) + correction;
  };

  // Z and X=C^{-1}Z are allocated only in PCG mode and reused by every
  // random/residual trace in the current REML iteration.
  Eigen::MatrixXd pcgTraceZ;
  Eigen::MatrixXd pcgTraceX;
  bool reportedOpenMpSlq = false;
  bool reportedBlockEngine = false;
  bool reportedOpenMpResidualBlocks = false;

  auto preparePCGTraceProbes = [&](const Eigen::Index n){
    if(solverName != "pcg"){ return; }
    const auto started = std::chrono::steady_clock::now();
    pcgTraceZ.resize(n, pcgTraceProbes);
    pcgTraceX.resize(n, pcgTraceProbes);
    for(int p = 0; p < pcgTraceProbes; ++p){
      for(Eigen::Index i = 0; i < n; ++i){
        pcgTraceZ(i,p) = pcgProbeSign(static_cast<unsigned long long>(i),
                                     static_cast<unsigned long long>(p));
      }
    }
    const bool haveLanczosGuess =
      pcgLanczosInverseGuess.rows() == n &&
      pcgLanczosInverseGuess.cols() == pcgTraceProbes && pcgLogDetScaling == "none";
    pcgTraceX = solvePCGBatch(
      pcgTraceZ,
      haveLanczosGuess ? &pcgLanczosInverseGuess : nullptr,
      "Hutchinson trace probes"
    );
    for(int p = 0; p < pcgTraceProbes; ++p){
      if(!pcgTraceX.col(p).allFinite()){
        Rcpp::stop("PCG failed while preparing Hutchinson trace probes for C^{-1}.");
      }
    }
    if(pcgTraceEstimator == "control-variate"){
      const int rank = std::min<int>(pcgTraceRank, n);
      const double needed = 8.0 * static_cast<double>(n) * (rank + pcgTraceProbes) +
        5.0 * rank * rank;
      if(needed > blockSchurPolicy.maxDoubles){ requireMemory("PCG trace control variate", needed); }
      auto applyVector = [&](const Eigen::VectorXd & value) -> Eigen::VectorXd {
        Rcpp::checkUserInterrupt();
        return applyCurrentPCG(Eigen::MatrixXd(value)).col(0);
      };
      const sommer::pcg_estimators::LandmarkControlVariate control =
        sommer::pcg_estimators::makeLandmarkControlVariate(n, applyVector, rank,
          static_cast<unsigned>(pcgSeed),
          static_cast<std::size_t>(blockSchurPolicy.maxDoubles * sizeof(double)));
      pcgTraceX.noalias() -= control.basis *
        (control.projectedInverse * (control.basis.transpose() * pcgTraceZ));
      pcgTraceZ.conservativeResize(n, pcgTraceProbes + rank);
      pcgTraceX.conservativeResize(n, pcgTraceProbes + rank);
      const double rootCount = std::sqrt(static_cast<double>(pcgTraceProbes));
      pcgTraceZ.rightCols(rank) = rootCount * control.basis;
      pcgTraceX.rightCols(rank).noalias() = rootCount * control.basis * control.projectedInverse;
    }
    pcgTraceSetupSeconds += std::chrono::duration<double>(
      std::chrono::steady_clock::now() - started).count();
  };


    if(verbose && !reportedOpenMpSlq && solverName == "pcg"){
#ifdef _OPENMP
      Rcpp::Rcout
        << "OpenMP active: parallel SLQ log-determinant probes ("
        << pcgTraceProbes
        << " probes)."
        << arma::endl;
#else
      Rcpp::Rcout
        << "OpenMP unavailable: SLQ log-determinant probes run serially."
        << arma::endl;
#endif
      reportedOpenMpSlq = true;
    }
  auto pcgTraceCInverseTimesSparse = [&](const arma::sp_mat & B) -> double {
    if(solverName != "pcg" || pcgTraceX.cols() < pcgTraceProbes){
      Rcpp::stop("PCG trace probes are not available.");
    }
    double out = 0.0;
    pcgLastSparseTraceSamples = Eigen::VectorXd::Zero(pcgTraceProbes);
    for(Eigen::Index p = 0; p < pcgTraceX.cols(); ++p){
      double one = 0.0;
      for(arma::sp_mat::const_iterator it = B.begin(); it != B.end(); ++it){
        one += pcgTraceZ(static_cast<Eigen::Index>(it.row()), p)
             * (*it)
             * pcgTraceX(static_cast<Eigen::Index>(it.col()), p);
      }
      out += one;
      if(p < pcgTraceProbes){ pcgLastSparseTraceSamples(p) = one; }
    }
    return out / static_cast<double>(pcgTraceProbes);
  };


  // ============================================================
  // CHANGE 7: reusable factorisation of the small Kronecker block R0
  // ============================================================
  EigenLDLT RkronFactor;
  bool RkronSymbolicReady = false;
  std::vector<int> RkronOuterPattern;
  std::vector<int> RkronInnerPattern;

  auto factorizeRkronWithCachedPattern =
    [&](const EigenSpMat & Re) -> bool {
      return sommer::factorizeLDLTWithCachedPattern(
        Re, RkronFactor, RkronSymbolicReady, RkronOuterPattern, RkronInnerPattern);
    };

  // ============================================================
  // CHANGE 6: reusable symbolic factorisation for generic non-diagonal R
  // ============================================================
  EigenLDLT Rfactor;
  bool RsymbolicReady = false;
  std::vector<int> RouterPattern;
  std::vector<int> RinnerPattern;

  // CHOLMOD backend for R, used when solverName=="cholmod" (Phase 2). R is
  // usually far smaller than C (this path is only reached for irregular
  // residual designs that are neither diagonal nor repeated-block), but the
  // supernodal factor is reused with the same symbolic-pattern-caching
  // discipline as C for consistency.
  CholmodState cholmodRState;

  auto solveRMatrixCholmod =
    [&](const Eigen::Ref<const Eigen::MatrixXd> & rhs,
        const std::string & context) -> Eigen::MatrixXd {
      return cholmodRState.solve(rhs,
        "CHOLMOD solve for the residual covariance matrix R failed in " + context + ".");
    };

  auto factorizeRWithCachedPattern =
    [&](const EigenSpMat & Re) -> bool {

      if(solverName == "cholmod"){
        return cholmodRState.analyzePattern(Re) && cholmodRState.factorize(Re);
      }
      return sommer::factorizeLDLTWithCachedPattern(
        Re, Rfactor, RsymbolicReady, RouterPattern, RinnerPattern);
    };


  // ============================================================
  // GLOBAL OPTIMIZER SAFEGUARDS FOR TRANSFORMED COVARIANCE COORDINATES
  // ============================================================
  //
  // The generic covariance descriptors use log / tanh / Cholesky working
  // coordinates.  Those mappings guarantee valid covariance structures over
  // very large regions of parameter space, so "is Sigma positive definite?"
  // is no longer an adequate proxy for "is this AI step reasonable?".
  //
  // We therefore use three complementary safeguards:
  //
  //   1) likelihood convergence uses |Delta logLik|, never a signed decrease;
  //   2) every proposed global covariance step is checked at the next pass
  //      through the exact REML likelihood and backtracked geometrically when
  //      it decreases the likelihood;
  //   3) transformed-coordinate trust caps limit one-iteration movement before
  //      the likelihood line search is even attempted.
  //
  // Backtracking is implemented without duplicating the expensive likelihood
  // code: a rejected trial is replaced by a half-step, iIter is decremented,
  // and the same iteration slot is recomputed.  Consequently retries do not
  // consume the user's requested nIters.
  // ============================================================

  bool haveAcceptedLikelihood = false;
  double acceptedLikelihood =
    -std::numeric_limits<double>::infinity();

  bool pendingLineSearch = false;
  arma::vec lineSearchBase;
  arma::vec lineSearchTarget;
  double lineSearchAlpha = 1.0;
  int lineSearchHalvings = 0;
  int coefficientEvaluations = 0;
  int totalLineSearchHalvings = 0;
  int acceleratedProposals = 0;
  struct FactorizationTiming {
    double elapsedSeconds = 0.0;
    double cpuSeconds = 0.0;
    int calls = 0;
  };
  struct NumericFactorTimer {
    FactorizationTiming & timing;
    std::chrono::steady_clock::time_point started;
    double cpuStarted;
    static double processCPUSeconds(){
#ifdef _WIN32
      FILETIME created, exited, kernelTime, userTime;
      if(!GetProcessTimes(GetCurrentProcess(), &created, &exited, &kernelTime, &userTime)){
        return std::numeric_limits<double>::quiet_NaN();
      }
      ULARGE_INTEGER kernelTicks, userTicks;
      kernelTicks.LowPart = kernelTime.dwLowDateTime;
      kernelTicks.HighPart = kernelTime.dwHighDateTime;
      userTicks.LowPart = userTime.dwLowDateTime;
      userTicks.HighPart = userTime.dwHighDateTime;
      return (static_cast<double>(kernelTicks.QuadPart) +
        static_cast<double>(userTicks.QuadPart)) / 1.0e7;
#else
      const std::clock_t ticks = std::clock();
      return ticks == static_cast<std::clock_t>(-1)
        ? std::numeric_limits<double>::quiet_NaN()
        : static_cast<double>(ticks) / CLOCKS_PER_SEC;
#endif
    }
    explicit NumericFactorTimer(FactorizationTiming & value)
      : timing(value), started(std::chrono::steady_clock::now()), cpuStarted(processCPUSeconds()) {}
    ~NumericFactorTimer(){
      const double cpuFinished = processCPUSeconds();
      timing.elapsedSeconds += std::chrono::duration<double>(
        std::chrono::steady_clock::now() - started).count();
      timing.cpuSeconds += cpuFinished - cpuStarted;
      ++timing.calls;
    }
  };
  FactorizationTiming cholmodTiming;
  FactorizationTiming ldltTiming;
  FactorizationTiming blockSchurTiming;
  std::size_t maxIrregularDesignEntries = 0;
  std::size_t maxIrregularMMEEntries = 0;
  arma::vec previousAcceptedParameters;
  const bool useAitken = proposal.extrapolation == sommer::control::Extrapolation::Aitken;
  const bool useHybridAI = hybridRequested && !vcMap.active;
  HybridAIState hybridState;
  bool hybridForceExact = false;
  int lastExactAIIteration = -1;
  int exactAIRefreshes = 0;
  int hybridAIUpdates = 0;
  int hybridAIResets = 0;
  int hybridAICertifications = 0;
  const int hybridRefreshInterval = 5;
  std::size_t sensitivityRHSColumns = 0;
  std::size_t sensitivityRHSColumnsSkipped = 0;
  bool finalAIExact = false;
  bool pendingHybridProposal = false;
  double finalProjectedScoreNorm = std::numeric_limits<double>::quiet_NaN();
  arma::uvec hybridFree;
  {
    arma::vec constraints;
    for(int term = 0; term < nRRe; ++term){
      constraints = arma::join_cols(constraints, covConstraint(term));
    }
    hybridFree = arma::find(constraints != 3);
  }
  const arma::uword modelFreeCount = hybridFree.n_elem;

  const int maxLineSearchHalvings = 24;

  // Permit only a negligible numerical drop when deciding whether a trial
  // likelihood is acceptable.
  const double likelihoodAcceptTolerance =
    std::max(
      1.0e-10,
      0.1 * tolParConvLL
    );

  bool lineSearchStalled = false;

  // ------------------------------------------------------------
  // Per-parameter trust caps in WORKING coordinates.
  // ------------------------------------------------------------
  // The covariance model owns these caps through its CovarianceFactor
  // descriptor.  The optimizer no longer contains model-specific branches.
  const arma::vec workingStepCap = sommer::covarianceWorkingStepCaps(covStructI);

  std::vector<arma::mat> lastEvaluatedRandomPrecision(
    static_cast<std::size_t>(nRe)
  );

  const sommer::ResidualPlan residualPlan = sommer::ResidualPlan::select(
    residualStructurallyDiagonal, residualKronBlocks, residualGridBlocks,
    useH, Hdiag);
  factorSchurResidualEligible = residualPlan.elementwisePrecision();
  if(saddleActive && residualPlan.grid()){
    Rcpp::stop("Saddle-point fitting does not yet support compact grid/ragged residual layouts.");
  }

  ////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////
  // START ITERATIVE ALGORITHM
  ////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////

sommer::MMEProgress progress(nRRe);
for (int iIter = 0; iIter < nIters; ++iIter) {
  phaseProfile.startIteration();
  ++coefficientEvaluations;
  std::swap(previousCrossGroupDemand, currentCrossGroupDemand);
  currentCrossGroupDemand.groupSizes.clear();
  currentCrossGroupDemand.entries.clear();
  currentCrossGroupDemand.groupSignature = 0;

    // ###########################
    // # 1) absorption of M into y to obtain y'Py and logDetC
    // ###########################

    const int residualStruct =
      nRRe - 1;

    const arma::uword nResidualPar =
      covPar(residualStruct).n_elem;
    currentResidualParameterCount = nResidualPar;

    const arma::mat residualSigma = compactResidualCovariance
      ? arma::mat()
      : theta(residualStruct);

    randomPrecisionState.assign(static_cast<std::size_t>(nRe), DescriptorPrecisionState());

    for(int iStruct = 0; iStruct < nRRe; ++iStruct){
      const std::size_t structureOffset =
        static_cast<std::size_t>(iStruct);
      const std::size_t parameterCount =
        static_cast<std::size_t>(covPar(iStruct).n_elem);
      covarianceDerivativeCache[structureOffset].resize(parameterCount);
      covarianceDerivativeReady[structureOffset].assign(parameterCount, false);
    }

    if(!compactResidualCovariance &&
       (residualSigma.n_rows != static_cast<arma::uword>(residualDim) ||
      residualSigma.n_cols != static_cast<arma::uword>(residualDim))){
      Rcpp::stop("Residual covariance dimension changed unexpectedly.");
    }

    arma::field<arma::mat> residualLocalD1(
      nResidualPar
    );

    arma::field<arma::sp_mat> residualDerivativeBasis(
      nResidualPar
    );

    // ============================================================
    // STRUCTURED R HANDLING
    // ============================================================
    // R is assembled sparsely from the residual covariance bases.
    // If R (and the optional weight transform) are diagonal, all
    // R^{-1} applications are elementwise.  Otherwise R is kept
    // sparse and factorised by Eigen::SimplicialLDLT; no dense
    // inv_sympd(R) is formed.
    // ============================================================

    arma::sp_mat Rmat;

    const bool Rdiag = residualPlan.diagonal();
    const bool Rkron = residualPlan.kronecker();

    arma::vec RdiagInv;
    arma::vec effectiveRiDiag;
    double logDetR = 0.0;
    SelectedInverseSubset Rselected;
    arma::mat RkronInv;
    DescriptorPrecisionState residualPrecisionState;

    const bool Rgrid = residualPlan.grid();
    std::vector< std::vector<arma::mat> > gridF;
    std::vector< std::vector<arma::mat> > gridP;
    std::vector<arma::mat> gridQcolM;
    std::vector<arma::mat> gridQmmChol;
    std::vector<arma::mat> gridQmmInv;
    std::vector<arma::mat> gridDFull;
    std::vector<int> gridDFactor;
    std::vector<int> gridDSection;
    double gridScale = 1.0;

    auto gridKron =
      [&](const arma::mat & X, const std::vector<arma::mat> & factors,
          const double scale) -> arma::mat {
        return factorEngine.applyKroneckerRight(X.t(), factors, scale).t();
      };

    auto gridEmbed =
      [&](const int b, const arma::mat & V) -> arma::mat {
        arma::mat X(gridSize[static_cast<std::size_t>(b)], V.n_cols, arma::fill::zeros);
        X.rows(gridObsPos[static_cast<std::size_t>(b)]) = V;
        return X;
      };

    auto gridApplyRi =
      [&](const int b, const arma::mat & V) -> arma::mat {
        const std::size_t bb = static_cast<std::size_t>(b);
        arma::mat out(V.n_rows, V.n_cols, arma::fill::zeros);
        // Columns that vanish on this block (sectioned derivatives) stay zero.
        std::vector<arma::uword> active;
        active.reserve(V.n_cols);
        for(arma::uword c = 0; c < V.n_cols; ++c){
          if(arma::any(V.col(c) != 0.0)){ active.push_back(c); }
        }
        if(active.empty()){
          return out;
        }
        const bool allActive = active.size() == V.n_cols;
        const arma::uvec activeCols(active);
        arma::mat Y = gridKron(
          gridEmbed(b, allActive ? V : arma::mat(V.cols(activeCols))),
          gridP[bb], 1.0 / gridScale);
        if(gridMissPos[bb].n_elem > 0){
          const arma::mat Ym = Y.rows(gridMissPos[bb]);
          const arma::mat & U = gridQmmChol[bb];
          Y -= gridQcolM[bb] *
            arma::solve(arma::trimatu(U), arma::solve(arma::trimatl(U.t()), Ym));
        }
        if(allActive){
          return Y.rows(gridObsPos[bb]);
        }
        out.cols(activeCols) = Y.rows(gridObsPos[bb]);
        return out;
      };

    // Returns false when dR_S/dpar_k is structurally zero on this block.
    auto gridDerivFactors =
      [&](const int b, const arma::uword k, std::vector<arma::mat> & factors) -> bool {
        const std::size_t bb = static_cast<std::size_t>(b);
        const int f = gridDFactor[static_cast<std::size_t>(k)];
        if(f < 0){
          factors = gridF[bb];
          return true;
        }
        const int section = gridDSection[static_cast<std::size_t>(k)];
        if(section >= 0 && gridSectionOf[bb][static_cast<std::size_t>(f)] != section){
          return false;
        }
        const arma::uvec & lev = gridLevels[bb][static_cast<std::size_t>(f)];
        arma::mat d = gridDFull[static_cast<std::size_t>(k)].submat(lev, lev);
        if(!arma::any(arma::vectorise(d) != 0.0)){
          return false;
        }
        factors = gridF[bb];
        factors[static_cast<std::size_t>(f)] = std::move(d);
        return true;
      };

    auto gridApplyDR =
      [&](const int b, const arma::uword k, const arma::mat & V) -> arma::mat {
        std::vector<arma::mat> factors;
        if(!gridDerivFactors(b, k, factors)){
          return arma::zeros<arma::mat>(V.n_rows, V.n_cols);
        }
        return gridKron(gridEmbed(b, V), factors, gridScale)
          .rows(gridObsPos[static_cast<std::size_t>(b)]);
      };

    // tr(R_oo^{-1} dR_oo) = dlog|R_S| - tr(Q_mm^{-1} [Q dR_S Q]_mm); the
    // Kronecker form of Q dR_S Q gives its (m x m) block entrywise.
    auto gridTraceDR =
      [&](const arma::uword k) -> double {
        double total = 0.0;
        for(int b = 0; b < residualNBlocks; ++b){
          const std::size_t bb = static_cast<std::size_t>(b);
          std::vector<arma::mat> factors;
          if(!gridDerivFactors(b, k, factors)){
            continue;
          }
          const double S = static_cast<double>(gridSize[bb]);
          const int f = gridDFactor[static_cast<std::size_t>(k)];
          const arma::uword m = gridMissPos[bb].n_elem;
          double trace = S;
          if(f < 0){
            trace -= static_cast<double>(m);
            total += trace;
            continue;
          }
          const std::size_t ff = static_cast<std::size_t>(f);
          trace = (S / static_cast<double>(factors[ff].n_rows)) *
            arma::accu(gridP[bb][ff] % factors[ff]);
          if(m > 0){
            if(gridQmmInv[bb].is_empty()){
              const arma::mat Uinv = arma::inv(arma::trimatu(gridQmmChol[bb]));
              gridQmmInv[bb] = Uinv * Uinv.t();
            }
            const arma::mat PdP = gridP[bb][ff] * factors[ff] * gridP[bb][ff];
            const arma::umat & dig = gridMissDigits[bb];
            const arma::mat & Qinv = gridQmmInv[bb];
            double corr = 0.0;
            for(arma::uword j = 0; j < m; ++j){
              for(arma::uword i = 0; i < m; ++i){
                double v = PdP(dig(i, ff), dig(j, ff));
                for(int g = 0; g < gridNFactors && v != 0.0; ++g){
                  if(g == f){ continue; }
                  const std::size_t gg = static_cast<std::size_t>(g);
                  v *= gridP[bb][gg](dig(i, gg), dig(j, gg));
                }
                corr += Qinv(i, j) * v;
              }
            }
            trace -= corr / gridScale;
          }
          total += trace;
        }
        return total;
      };

    if(Rdiag && compressedData){
      logDetR = 0.0;
      for(int category = 0; category < residualDim; ++category){
        const double variance = residualSigma(category, category);
        if(!std::isfinite(variance) || variance <= 0.0){
          Rcpp::stop("Compressed residual covariance produced a non-positive/non-finite value.");
        }
        logDetR += sufficient.classCounts[category] * std::log(variance);
      }
    }else if(Rdiag){

      arma::vec rd(
        nR,
        arma::fill::zeros
      );

      for(int rr = 0; rr < nR; ++rr){
        const int local =
          residualLocalOfRow[
            static_cast<std::size_t>(rr)
          ];

        rd(static_cast<arma::uword>(rr)) =
          residualSigma(
            static_cast<arma::uword>(local),
            static_cast<arma::uword>(local)
          );
      }

      const double minRd = rd.min();

      if(!rd.is_finite() || minRd <= 0.0){
        Rcpp::stop("Diagonal residual covariance produced a non-positive/non-finite value.");
      }

      RdiagInv = 1.0 / rd;
      if(saddleActive){
        const Eigen::Map<const Eigen::VectorXd> diagonal(rd.memptr(), nR);
        saddleResidual = diagonal.asDiagonal();
      }
      logDetR = arma::accu(arma::log(rd));

      if(!useH || Hdiag){
        effectiveRiDiag = RdiagInv;
        if(useH){
          effectiveRiDiag %= HdiagSquared;
        }
      }

    }else if(Rgrid){

      Rcpp::List rfactors =
        covDescriptor[static_cast<std::size_t>(residualStruct)]["factors"];
      const arma::vec & par = covPar(residualStruct);
      std::vector< std::vector<arma::mat> > Ffull(static_cast<std::size_t>(gridNFactors));
      for(int f = 0; f < gridNFactors; ++f){
        Rcpp::List fl = Rcpp::as<Rcpp::List>(rfactors[f]);
        const int start1 = Rcpp::as<int>(fl["par_start"]);
        const int end1 = Rcpp::as<int>(fl["par_end"]);
        arma::vec localPar;
        if(end1 >= start1){
          localPar = par.subvec(static_cast<arma::uword>(start1 - 1),
                                static_cast<arma::uword>(end1 - 1));
        }
        const int nSections = factorEngine.factorSections(fl);
        std::vector<arma::mat> & evaluated = Ffull[static_cast<std::size_t>(f)];
        if(nSections > 0){
          evaluated.resize(static_cast<std::size_t>(nSections));
          for(int t = 0; t < nSections; ++t){
            evaluated[static_cast<std::size_t>(t)] = factorEngine.evalSectionFactor(fl, localPar, t);
          }
        }else{
          evaluated.assign(1, factorEngine.evalFactor(fl, localPar));
        }
      }

      gridScale = std::exp(par(0));
      gridF.assign(static_cast<std::size_t>(residualNBlocks), std::vector<arma::mat>());
      gridP.assign(static_cast<std::size_t>(residualNBlocks), std::vector<arma::mat>());
      gridQcolM.assign(static_cast<std::size_t>(residualNBlocks), arma::mat());
      gridQmmChol.assign(static_cast<std::size_t>(residualNBlocks), arma::mat());
      gridQmmInv.assign(static_cast<std::size_t>(residualNBlocks), arma::mat());
      std::vector<std::map<std::vector<arma::uword>, std::pair<int, double>>> gridFactorCache(
        static_cast<std::size_t>(gridNFactors));

      for(int b = 0; b < residualNBlocks; ++b){
        const std::size_t bb = static_cast<std::size_t>(b);
        const double S = static_cast<double>(gridSize[bb]);
        double blockLogDet = S * par(0);
        gridF[bb].resize(static_cast<std::size_t>(gridNFactors));
        gridP[bb].resize(static_cast<std::size_t>(gridNFactors));

        for(int f = 0; f < gridNFactors; ++f){
          const std::size_t ff = static_cast<std::size_t>(f);
          const arma::uvec & lev = gridLevels[bb][ff];
          const int section = gridSectionOf[bb][ff];
          std::vector<arma::uword> key;
          key.reserve(lev.n_elem + 1);
          key.push_back(static_cast<arma::uword>(section + 1));
          key.insert(key.end(), lev.begin(), lev.end());
          const auto cached = gridFactorCache[ff].find(key);
          if(cached != gridFactorCache[ff].end()){
            const std::size_t source = static_cast<std::size_t>(cached->second.first);
            gridF[bb][ff] = gridF[source][ff];
            gridP[bb][ff] = gridP[source][ff];
            blockLogDet += (S / static_cast<double>(lev.n_elem)) * cached->second.second;
            ++gridFactorReuses;
            continue;
          }
          gridF[bb][ff] = Ffull[ff][section >= 0 ? static_cast<std::size_t>(section) : 0].submat(lev, lev);
          arma::mat U;
          if(!gpuBackend.cholesky(U, gridF[bb][ff])){
            Rcpp::stop("Residual covariance factor is not positive definite on a residual block.");
          }
          const double factorLogDet = 2.0 * arma::accu(arma::log(U.diag()));
          blockLogDet += (S / static_cast<double>(lev.n_elem)) *
            factorLogDet;
          const arma::mat Uinv = arma::inv(arma::trimatu(U));
          gridP[bb][ff] = Uinv * Uinv.t();
          gridFactorCache[ff].emplace(std::move(key), std::make_pair(b, factorLogDet));
          ++gridFactorizations;
        }

        const arma::uvec & miss = gridMissPos[bb];
        if(miss.n_elem > 0){
          const arma::umat & dig = gridMissDigits[bb];
          gridQcolM[bb].set_size(gridSize[bb], miss.n_elem);
          for(arma::uword j = 0; j < miss.n_elem; ++j){
            arma::vec column(1);
            column(0) = 1.0 / gridScale;
            for(int f = 0; f < gridNFactors; ++f){
              const std::size_t ff = static_cast<std::size_t>(f);
              column = arma::kron(column, gridP[bb][ff].col(dig(j, ff)));
            }
            gridQcolM[bb].col(j) = column;
          }
          arma::mat Qmm = gridQcolM[bb].rows(miss);
          Qmm = 0.5 * (Qmm + Qmm.t());
          if(!gpuBackend.cholesky(gridQmmChol[bb], Qmm)){
            Rcpp::stop("Missing-cell residual precision block is not positive definite.");
          }
          blockLogDet += 2.0 * arma::accu(arma::log(gridQmmChol[bb].diag()));
        }

        if(!std::isfinite(blockLogDet)){
          Rcpp::stop("Non-finite residual log-determinant on a residual block.");
        }
        logDetR += blockLogDet;
      }

      if(verbose && iIter == 0){
        Rcpp::Rcout
          << "Using exact grid residual engine (missing cells absorbed): "
          << residualNBlocks
          << " blocks"
          << arma::endl;
      }

    }else if(Rkron){

      const bool residualPrecisionOK =
        factorEngine.buildDescriptorPrecision(
          covDescriptor[static_cast<std::size_t>(residualStruct)],
          covPar(residualStruct),
          residualPrecisionState
        )
        &&
        residualPrecisionState.totalDim
          == static_cast<arma::uword>(residualKronBlockSize);

      if(residualPrecisionOK){
        logDetR =
          static_cast<double>(residualKronNBlocks)
          *
          residualPrecisionState.logDet;
        // Weighted paths read the dense block inverse directly.
        if(useH){
          arma::mat blockPrecision(1, 1, arma::fill::ones);
          for(std::size_t factor = 0; factor < residualPrecisionState.factorPrecision.size(); ++factor){
            blockPrecision = arma::kron(blockPrecision,
              factorEngine.factorPrecisionDense(residualPrecisionState, factor));
          }
          RkronInv = residualPrecisionState.inverseScale * blockPrecision;
          RkronInv = 0.5 * (RkronInv + RkronInv.t());
        }

      }else{

      arma::sp_mat R0 =
        arma::sp_mat(
          residualSigma
        );

      EigenSpMat R0e(
        residualKronBlockSize,
        residualKronBlockSize
      );

      std::vector<EigenTriplet> tripsR0;
      tripsR0.reserve(
        static_cast<std::size_t>(
          R0.n_nonzero
        )
      );

      for(arma::sp_mat::const_iterator it = R0.begin();
          it != R0.end();
          ++it){

        tripsR0.emplace_back(
          static_cast<int>(it.row()),
          static_cast<int>(it.col()),
          (*it)
        );
      }

      R0e.setFromTriplets(
        tripsR0.begin(),
        tripsR0.end()
      );

      R0e.makeCompressed();

      bool factorOK =
        factorizeRkronWithCachedPattern(
          R0e
        );

      if(!factorOK){

        bool recovered = false;

        for(int k = 0; k < 3 && !recovered; ++k){

          const double jitter =
            tolParInv
            *
            std::pow(
              10.0,
              static_cast<double>(k)
            );

          arma::sp_mat R0try =
            R0
            +
            arma::speye<arma::sp_mat>(
              residualKronBlockSize,
              residualKronBlockSize
            )
            *
            jitter;

          std::vector<EigenTriplet> retryTrips;
          retryTrips.reserve(
            static_cast<std::size_t>(
              R0try.n_nonzero
            )
          );

          EigenSpMat R0tryE(
            residualKronBlockSize,
            residualKronBlockSize
          );

          for(arma::sp_mat::const_iterator it = R0try.begin();
              it != R0try.end();
              ++it){

            retryTrips.emplace_back(
              static_cast<int>(it.row()),
              static_cast<int>(it.col()),
              (*it)
            );
          }

          R0tryE.setFromTriplets(
            retryTrips.begin(),
            retryTrips.end()
          );

          R0tryE.makeCompressed();

          if(
              factorizeRkronWithCachedPattern(
                R0tryE
              )
          ){
            R0 = R0try;
            R0e = R0tryE;
            recovered = true;
          }
        }

        if(!recovered){
          Rcpp::stop(
            "Kronecker residual block factorisation failed."
          );
        }
      }

      const Eigen::VectorXd Dr0 =
        RkronFactor.vectorD();

      double logDetR0 = 0.0;

      for(Eigen::Index j = 0;
          j < Dr0.size();
          ++j){

        const double dj = Dr0(j);

        if(!std::isfinite(dj) || dj <= 0.0){
          Rcpp::stop(
            "Kronecker residual block LDLT produced a non-positive/non-finite pivot."
          );
        }

        logDetR0 += std::log(dj);
      }

      logDetR =
        static_cast<double>(
          residualKronNBlocks
        )
        *
        logDetR0;

      Eigen::MatrixXd Iq =
        Eigen::MatrixXd::Identity(
          residualKronBlockSize,
          residualKronBlockSize
        );

      Eigen::MatrixXd R0invEig =
        RkronFactor.solve(
          Iq
        );

      if(RkronFactor.info() != Eigen::Success){
        Rcpp::stop(
          "Kronecker residual block inverse solve failed."
        );
      }

      RkronInv.set_size(
        static_cast<arma::uword>(
          residualKronBlockSize
        ),
        static_cast<arma::uword>(
          residualKronBlockSize
        )
      );

      std::copy(
        R0invEig.data(),
        R0invEig.data() + R0invEig.size(),
        RkronInv.memptr()
      );

      RkronInv =
        0.5
        *
        (
          RkronInv
          +
          RkronInv.t()
        );

      }

      if(verbose && iIter == 0){
        Rcpp::Rcout
          << "Using exact descriptor-driven repeated-block/Kronecker residual engine: "
          << residualKronNBlocks
          << " blocks x "
          << residualKronBlockSize
          << " covariance coordinates"
          << arma::endl;
      }

    }else{

      // Generic sparse residual path for incomplete/unbalanced blocks.
      Rmat = compactResidualCovariance
        ? buildCompactResidualSparse(covPar(residualStruct), -1)
        : buildResidualSparseFromPattern(residualSigma);

      if(Rmat.n_rows > static_cast<arma::uword>(std::numeric_limits<int>::max())){
        Rcpp::stop("R is too large for the 32-bit Eigen sparse index type used in ai_mme_sp2().");
      }

      auto armaSparseToEigen =
        [&](const arma::sp_mat & A) -> EigenSpMat {

          EigenSpMat out(
            static_cast<int>(A.n_rows),
            static_cast<int>(A.n_cols)
          );

          std::vector<EigenTriplet> trips;
          trips.reserve(
            static_cast<std::size_t>(
              A.n_nonzero
            )
          );

          for(arma::sp_mat::const_iterator it = A.begin();
              it != A.end();
              ++it){

            trips.emplace_back(
              static_cast<int>(it.row()),
              static_cast<int>(it.col()),
              *it
            );
          }

          out.setFromTriplets(
            trips.begin(),
            trips.end()
          );

          out.makeCompressed();
          return out;
        };

      EigenSpMat Re =
        armaSparseToEigen(
          Rmat
        );

      bool RfactorOK =
        factorizeRWithCachedPattern(
          Re
        );

      if(!RfactorOK){

        bool recovered = false;

        for(int k = 0; k < 3 && !recovered; ++k){

          const double jitter =
            tolParInv
            *
            std::pow(
              10.0,
              static_cast<double>(k)
            );

          arma::sp_mat Rtry =
            Rmat
            +
            arma::speye<arma::sp_mat>(
              nR,
              nR
            )
            *
            jitter;

          Re =
            armaSparseToEigen(
              Rtry
            );

          if(
              factorizeRWithCachedPattern(
                Re
              )
          ){
            Rmat = Rtry;
            recovered = true;
          }
        }

        if(!recovered){
          Rcpp::stop(
            "Sparse LDLT factorisation of the residual covariance matrix R failed."
          );
        }
      }

      if(solverName == "cholmod"){
        logDetR += M_cholmod_factor_ldetA(cholmodRState.factor);
        if(!std::isfinite(logDetR)){
          Rcpp::stop("CHOLMOD produced a non-finite log-determinant for R.");
        }
        // No Takahashi subset for a supernodal factor: Rselected stays
        // empty, and sparseTraceInverseTimes()'s fallback is routed through
        // solveRMatrixCholmod (see the R-trace call site below).
      }else{
        const Eigen::VectorXd Dr =
          Rfactor.vectorD();

        for(Eigen::Index j = 0; j < Dr.size(); ++j){

          if(!std::isfinite(Dr(j)) || Dr(j) <= 0.0){
            Rcpp::stop(
              "Residual sparse LDLT produced a non-positive/non-finite pivot."
            );
          }

          logDetR += std::log(Dr(j));
        }

        if(!saddleActive && !rrConditional){
          buildSelectedInverseSubset(
            Rfactor,
            "R",
            Rselected,
            false
          );
        }
      }
    }

    // Apply the effective residual precision Hs * R^{-1} * Hs'.
    auto applyRiDense = [&](const arma::mat & B) -> arma::mat {
      if(B.n_rows != static_cast<arma::uword>(nR)){
        Rcpp::stop("R^{-1} application received an incompatible RHS.");
      }

      if(effectiveRiDiag.n_elem == static_cast<arma::uword>(nR)){
        arma::mat out = B;
        out.each_col() %= effectiveRiDiag;
        return out;
      }

      arma::mat rhs = B;
      if(useH){ rhs = arma::mat(Hs.t() * rhs); }

      arma::mat solved;
      if(Rdiag){
        solved = rhs;
        solved.each_col() %= RdiagInv;
      }else if(Rgrid){

        solved.zeros(rhs.n_rows, rhs.n_cols);
        #pragma omp parallel for if(residualNBlocks > 1)
        for(int b = 0; b < residualNBlocks; ++b){
          const arma::uvec & rows = residualKronRows[static_cast<std::size_t>(b)];
          const arma::mat solvedBlock = gridApplyRi(b, rhs.rows(rows));
          for(arma::uword localRow = 0; localRow < rows.n_elem; ++localRow){
            solved.row(rows(localRow)) = solvedBlock.row(localRow);
          }
        }

      }else if(Rkron){

        solved.zeros(
          rhs.n_rows,
          rhs.n_cols
        );

        if(verbose && !reportedOpenMpResidualBlocks){
#ifdef _OPENMP
          Rcpp::Rcout
            << "OpenMP active: parallel repeated residual-block precision applications ("
            << residualKronNBlocks
            << " blocks)."
            << arma::endl;
#else
          Rcpp::Rcout
            << "OpenMP unavailable: repeated residual-block precision applications run serially."
            << arma::endl;
#endif
          reportedOpenMpResidualBlocks = true;
        }

        // Apply the same small R0^{-1} to every independent block.
        #pragma omp parallel for if(residualKronNBlocks > 1)
        for(int b = 0;
            b < residualKronNBlocks;
            ++b){

          const arma::uvec & rows =
            residualKronRows[
              static_cast<std::size_t>(b)
            ];

          arma::mat solvedBlock;
          if(residualPrecisionState.factorWise){
            solvedBlock = factorEngine.applyKroneckerRight(
                rhs.rows(rows).t(),
                residualPrecisionState.factorPrecision,
                residualPrecisionState.inverseScale,
                &residualPrecisionState.sparseFactorPrecision
              ).t();
          }else{
            solvedBlock = RkronInv * rhs.rows(rows);
          }

          for(arma::uword localRow = 0;
              localRow < rows.n_elem;
              ++localRow){
            solved.row(rows(localRow)) = solvedBlock.row(localRow);
          }
        }

      }else{
        Eigen::Map<const Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor> > rhsEig(
          rhs.memptr(), static_cast<Eigen::Index>(rhs.n_rows), static_cast<Eigen::Index>(rhs.n_cols)
        );
        Eigen::MatrixXd solEig;
        if(solverName == "cholmod"){
          solEig = solveRMatrixCholmod(rhsEig, "R^{-1} application");
        }else{
          solEig = Rfactor.solve(rhsEig);
          if(Rfactor.info() != Eigen::Success){
            Rcpp::stop("Sparse residual LDLT solve failed.");
          }
        }
        solved.set_size(rhs.n_rows, rhs.n_cols);
        std::copy(solEig.data(), solEig.data()+solEig.size(), solved.memptr());
      }

      if(useH){ solved = arma::mat(Hs * solved); }
      return solved;
    };

    auto applyRiVec = [&](const arma::vec & v) -> arma::vec {
      arma::mat tmp(v);
      arma::mat ans = applyRiDense(tmp);
      return ans.col(0);
    };

    phaseProfile.enter(sommer::FitPhase::Assembly);
    // ------------------------------------------------------------
    // Build the data part of C and the MME RHS without an explicit Ri.
    // ------------------------------------------------------------
    arma::vec Riy;
    double yRiy = 0.0;
    if(compressedData){
      compressedVariance.resize(residualDim);
      for(int category = 0; category < residualDim; ++category){
        compressedVariance(category) = residualSigma(category, category);
      }
    }else{
      Riy = applyRiVec(yDense);
      yRiy = arma::dot(yDense, Riy);
    }
    arma::vec rhsMME;

    Eigen::SparseMatrix<double> RiWsp;
    const sommer::AssemblyPath assemblyPath = residualPlan.assembly(
      useH && !Hdiag && Hs.n_nonzero <= 4.0 * nR && HsInvT.n_nonzero <= 4.0 * nR);
    bool RiWisSparse = assemblyPath == sommer::AssemblyPath::Sparse;
    std::vector<arma::mat> residualKronBlockRiW(
      static_cast<std::size_t>(residualKronNBlocks)
    );

    if(!saddleActive){ ++normalEquationAssemblies; }
    if(compressedData){
      Eigen::VectorXd compressedRhs;
      sufficient.assemble(compressedVariance, C, compressedRhs, yRiy);
      rhsMME = arma::vec(compressedRhs.data(), compressedRhs.size());
      RiWisSparse = true;
    }else if(saddleActive){
      C.resize(nEffects, nEffects);
      C.setZero();
      RiWisSparse = false;
      const Eigen::Map<const Eigen::VectorXd> RiyEig(Riy.memptr(), nR);
      const Eigen::VectorXd rhs = W.transpose() * RiyEig;
      rhsMME = arma::vec(rhs.data(), rhs.size());
      if(!Rdiag && !Rkron){
        saddleResidual = armaSparseToEigenGlobal(Rmat);
      }else if(Rkron){
        arma::mat covariance = residualSigma;
        if(residualPrecisionState.factorWise){
          covariance = factorEngine.evaluateDescriptor(
            covDescriptor[static_cast<std::size_t>(residualStruct)], covPar(residualStruct));
        }else if(!eigenSpdInverse(RkronInv, covariance)){
          Rcpp::stop("Cannot recover the evaluated residual block for saddle-point fitting.");
        }
        std::vector<EigenTriplet> entries;
        for(const arma::uvec & rows : residualKronRows){
          for(arma::uword column = 0; column < rows.n_elem; ++column){
            for(arma::uword row = 0; row < rows.n_elem; ++row){
              const double value = covariance(row, column);
              if(value != 0.0){
                if(entries.size() >= maxFallbackEntries){
                  Rcpp::stop("Saddle-point residual covariance exceeds the sparse entry limit.");
                }
                entries.emplace_back(rows(row), rows(column), value);
              }
            }
          }
        }
        saddleResidual.resize(nR, nR);
        saddleResidual.setFromTriplets(entries.begin(), entries.end());
      }
      saddlePrior.resize(nEffects, nEffects);
      saddlePrior.setZero();
    }else if(RiWisSparse){
      if(effectiveRiDiag.n_elem == static_cast<arma::uword>(nR)){
        const Eigen::Map<const Eigen::VectorXd> effectiveRiDiagEig(
          effectiveRiDiag.memptr(), static_cast<Eigen::Index>(effectiveRiDiag.n_elem)
        );
        RiWsp = effectiveRiDiagEig.asDiagonal() * W;
      }else{
        const Eigen::SparseMatrix<double> He = armaSparseToEigenGlobal(Hs);
        const Eigen::Map<const Eigen::VectorXd> diagonal(
          RdiagInv.memptr(), static_cast<Eigen::Index>(RdiagInv.n_elem)
        );
        const Eigen::SparseMatrix<double> weightedDesign = He.transpose() * W;
        RiWsp = He * (diagonal.asDiagonal() * weightedDesign);
      }
      C = Eigen::SparseMatrix<double>(W.transpose()) * RiWsp;
      C.makeCompressed();
      const Eigen::Map<const Eigen::VectorXd> RiyEig(
        Riy.memptr(), static_cast<Eigen::Index>(Riy.n_elem)
      );
      const Eigen::VectorXd rhsMMEEig = W.transpose() * RiyEig;
      rhsMME.set_size(rhsMMEEig.size());
      std::copy(rhsMMEEig.data(), rhsMMEEig.data() + rhsMMEEig.size(), rhsMME.memptr());
    }else if(assemblyPath == sommer::AssemblyPath::Blocks){
      // A repeated residual block does not couple observations from other
      // blocks.  Restricting each operation to the MME columns active in a
      // block avoids the global nR x nEffects dense RiW temporary.
      C.resize(nEffects, nEffects);
      RiWsp.resize(nR, nEffects);
      std::vector<Eigen::Triplet<double>> Ctriplets;
      std::vector<Eigen::Triplet<double>> RiWspTriplets;

      for(int block = 0; block < residualKronNBlocks; ++block){
        const arma::uvec & rows =
          residualKronRows[static_cast<std::size_t>(block)];
        const arma::uvec & columns =
          residualActiveColumns[static_cast<std::size_t>(block)];

        if(columns.n_elem == 0){
          continue;
        }

        // blockW (including any diagonal-H weighting) is constant across
        // REML iterations and was precomputed once before the main loop.
        const arma::mat & blockW =
          residualKronBlockWCache[static_cast<std::size_t>(block)];

        arma::mat blockRiW;
        if(Rgrid){
          blockRiW = gridApplyRi(block, blockW);
        }else if(residualPrecisionState.factorWise){
          blockRiW = factorEngine.applyKroneckerRight(
              blockW.t(),
              residualPrecisionState.factorPrecision,
              residualPrecisionState.inverseScale,
              &residualPrecisionState.sparseFactorPrecision
            ).t();
        }else{
          blockRiW = RkronInv * blockW;
        }

        // blockW already carries the H diagonal, so the cross product is
        // (H W)' R^{-1} (H W); the record-space R^{-1}W needs one more H.
        const arma::mat blockCross = blockW.t() * blockRiW;
        if(useH){
          arma::vec hDiagonal(rows.n_elem);
          for(arma::uword localRow = 0; localRow < rows.n_elem; ++localRow){
            hDiagonal(localRow) = Hs(rows(localRow), rows(localRow));
          }
          blockRiW.each_col() %= hDiagonal;
        }
        residualKronBlockRiW[static_cast<std::size_t>(block)] = blockRiW;

        for(arma::uword localCol = 0; localCol < columns.n_elem; ++localCol){
          const arma::uword globalCol = columns(localCol);
          for(arma::uword localRow = 0; localRow < columns.n_elem; ++localRow){
            const double crossValue = blockCross(localRow, localCol);
            if(crossValue != 0.0){
              Ctriplets.emplace_back(
                static_cast<int>(columns(localRow)),
                static_cast<int>(globalCol),
                crossValue
              );
            }
          }
        }

        for(arma::uword localRow = 0; localRow < rows.n_elem; ++localRow){
          for(arma::uword localCol = 0; localCol < columns.n_elem; ++localCol){
            const double precisionDesignValue = blockRiW(localRow, localCol);
            if(precisionDesignValue != 0.0){
              RiWspTriplets.emplace_back(
                static_cast<int>(rows(localRow)),
                static_cast<int>(columns(localCol)),
                precisionDesignValue
              );
            }
          }
        }
      }

      C.setFromTriplets(Ctriplets.begin(), Ctriplets.end());
      C.makeCompressed();
      RiWsp.setFromTriplets(RiWspTriplets.begin(), RiWspTriplets.end());
      RiWsp.makeCompressed();

      {
        const Eigen::Map<const Eigen::VectorXd> RiyEig(
          Riy.memptr(), static_cast<Eigen::Index>(Riy.n_elem)
        );
        const Eigen::VectorXd rhsMMEEig = W.transpose() * RiyEig;
        rhsMME.set_size(rhsMMEEig.size());
        std::copy(rhsMMEEig.data(), rhsMMEEig.data() + rhsMMEEig.size(), rhsMME.memptr());
      }
      RiWisSparse = true;
    }else{
      const int batchSize = std::max(1, std::min(static_cast<int>(sommer::WorkspacePolicy::panelColumns()),
        static_cast<int>(sommer::WorkspacePolicy::cacheEntries()) /
          std::max(static_cast<int>(nR), static_cast<int>(nEffects))));
      std::vector<Eigen::Triplet<double>> Ctriplets;
      std::vector<Eigen::Triplet<double>> precisionTriplets;
      for(int first = 0; first < nEffects; first += batchSize){
        const int count = std::min(batchSize, static_cast<int>(nEffects) - first);
        arma::mat design(nR, count, arma::fill::zeros);
        for(int local = 0; local < count; ++local){
          for(Eigen::SparseMatrix<double>::InnerIterator entry(W, first + local); entry; ++entry){
            design(entry.row(), local) = entry.value();
          }
        }
        const arma::mat precisionDesign = applyRiDense(design);
        const Eigen::Map<const Eigen::MatrixXd> precisionMap(precisionDesign.memptr(), nR, count);
        const Eigen::MatrixXd cross = W.transpose() * precisionMap;
        for(int local = 0; local < count; ++local){
          for(int row = 0; row < nR; ++row){
            const double value = precisionDesign(row, local);
            if(value != 0.0){
              if(precisionTriplets.size() >= maxFallbackEntries){
                Rcpp::stop("Irregular residual precision produces excessive design fill; use independent residual blocks or a structured residual model.");
              }
              precisionTriplets.emplace_back(row, first + local, value);
            }
          }
          for(int row = 0; row < nEffects; ++row){
            const double value = cross(row, local);
            if(value != 0.0){
              if(Ctriplets.size() >= maxFallbackEntries){
                Rcpp::stop("Irregular residual precision produces excessive MME fill; use independent residual blocks or a structured residual model.");
              }
              Ctriplets.emplace_back(row, first + local, value);
            }
          }
        }
      }
      C.resize(nEffects, nEffects);
      C.setFromTriplets(Ctriplets.begin(), Ctriplets.end());
      C.makeCompressed();
      RiWsp.resize(nR, nEffects);
      RiWsp.setFromTriplets(precisionTriplets.begin(), precisionTriplets.end());
      RiWsp.makeCompressed();
      maxIrregularDesignEntries = std::max(maxIrregularDesignEntries, precisionTriplets.size());
      maxIrregularMMEEntries = std::max(maxIrregularMMEEntries, Ctriplets.size());
      RiWisSparse = true;
      const Eigen::Map<const Eigen::VectorXd> RiyEig(Riy.memptr(), nR);
      const Eigen::VectorXd rhsMMEEig = W.transpose() * RiyEig;
      rhsMME.set_size(rhsMMEEig.size());
      std::copy(rhsMMEEig.data(), rhsMMEEig.data() + rhsMMEEig.size(), rhsMME.memptr());
    }

    // ------------------------------------------------------------
    // Construct inverses of random-effect covariance matrices
    // and add G^-1 blocks to M
    //
    // IMPORTANT:
    // lambda must remain in THIS scope because it is subsequently
    // used to construct Wu and the score equations.
    // ------------------------------------------------------------

    arma::field<arma::sp_mat> lambda(nReAl);

    if(nZs > 0){

        for(int i = 0; i < nRe; ++i){

            arma::mat thetaSym =
                arma::symmatu(
                    theta(i)
                );

            DescriptorPrecisionState & precisionState =
              randomPrecisionState[static_cast<std::size_t>(i)];

            bool lambdaOK =
              factorEngine.buildDescriptorPrecision(
                covDescriptor[static_cast<std::size_t>(i)],
                covPar(i),
                precisionState
              );

            arma::mat lambdaDense(1, 1, arma::fill::ones);
            for(std::size_t factor = 0; factor < precisionState.factorPrecision.size(); ++factor){
              lambdaDense = arma::kron(lambdaDense, factorEngine.factorPrecisionDense(precisionState, factor));
            }
            lambdaDense *= precisionState.inverseScale;
            lambdaDense = 0.5 * (lambdaDense + lambdaDense.t());

            if(!lambdaOK){

              precisionState.factorWise = false;

                // The covariance update machinery normally guarantees
                // positive definiteness.  Repair only when the fast
                // SPD inverse actually fails.
                arma::mat bend =
                    nearPDcpp(
                        thetaSym,
                        100,
                        1e-06,
                        1e-07
                    );

                lambdaOK =
                    eigenSpdInverse(
                        bend,
                        lambdaDense
                    );

                if(!lambdaOK){
                    Rcpp::stop(
                        "Unable to invert a random-effect covariance matrix "
                        "even after the nearPD fallback."
                    );
                }
            }

            lambda(i) =
                arma::sp_mat(
                    lambdaDense
                );

            lastEvaluatedRandomPrecision[
              static_cast<std::size_t>(i)
            ] = lambdaDense;

            const std::size_t iCache =
                static_cast<std::size_t>(i);

            if(
                lambdaDense.n_rows
                !=
                partitionStartCache[iCache].size()
                ||
                lambdaDense.n_cols
                !=
                partitionStartCache[iCache].size()
            ){
                Rcpp::stop(
                    "Random-effect covariance dimensions are inconsistent "
                    "with the cached MME partitions."
                );
            }

        }

        bool useMatrixFreePCG = (matrixFreePcgMode != "off");
        if(useMatrixFreePCG){
          useMatrixFreePCG = solverName == "pcg" && reml && computeCi == 0 && nRe > 0;
        }

        for(int i = 0; i < nRe && useMatrixFreePCG; ++i){
          useMatrixFreePCG =
            randomPrecisionState[static_cast<std::size_t>(i)].factorWise;
        }

        if(useMatrixFreePCG){
          matrixFreePCGReady = true;
          CpcgReady = false;

          matrixFreeApplyCMatrix =
            [&](const Eigen::MatrixXd & values) -> Eigen::MatrixXd {
              if(values.rows() != nEffects){
                Rcpp::stop("Matrix-free C application received an incompatible panel.");
              }

              Eigen::MatrixXd result = C * values;
              for(int iR = 0; iR < nRe; ++iR){
                const std::size_t iCache = static_cast<std::size_t>(iR);
                const arma::uword q =
                  randomPrecisionState[iCache].totalDim;
                const arma::uword nLevels = Ai(iR).n_rows;
                if(values.cols() == 0){ continue; }
                if(q > std::numeric_limits<arma::uword>::max() / static_cast<arma::uword>(values.cols())){
                  Rcpp::stop("Matrix-free tensor panel dimension overflow.");
                }
                arma::mat coefficients(nLevels, q * values.cols(), arma::fill::zeros);
                for(Eigen::Index rhs = 0; rhs < values.cols(); ++rhs){
                  for(arma::uword coordinate = 0; coordinate < q; ++coordinate){
                    const arma::uvec & indices = partitionIndexCache[iCache][static_cast<std::size_t>(coordinate)];
                    for(arma::uword level = 0; level < nLevels; ++level){
                      coefficients(level, coordinate + q * rhs) = values(indices(level), rhs);
                    }
                  }
                }
                const arma::mat related = arma::mat(Ai(iR) * coefficients);
                for(Eigen::Index rhs = 0; rhs < values.cols(); ++rhs){
                  const arma::mat priorApplied = factorEngine.applyKroneckerRight(
                    related.cols(q * rhs, q * (rhs + 1) - 1),
                    randomPrecisionState[iCache].factorPrecision, randomPrecisionState[iCache].inverseScale,
                    &randomPrecisionState[iCache].sparseFactorPrecision);
                  for(arma::uword coordinate = 0; coordinate < q; ++coordinate){
                    const arma::uvec & indices = partitionIndexCache[iCache][static_cast<std::size_t>(coordinate)];
                    for(arma::uword level = 0; level < nLevels; ++level){
                      result(indices(level), rhs) += priorApplied(level, coordinate);
                    }
                  }
                }
              }
              return result;
            };
          matrixFreeApplyC = [&](const Eigen::VectorXd & value) -> Eigen::VectorXd {
            return matrixFreeApplyCMatrix(Eigen::MatrixXd(value)).col(0);
          };

          matrixFreePCGDiagonal.resize(nEffects);
          for(int effect = 0; effect < nEffects; ++effect){
            matrixFreePCGDiagonal(effect) = C.coeff(effect,effect);
          }
          for(int iR = 0; iR < nRe; ++iR){
            const std::size_t iCache = static_cast<std::size_t>(iR);
            const arma::vec precisionDiagonal =
              factorEngine.kroneckerDiagonal(
                factorEngine.factorPrecisionFactorsDense(randomPrecisionState[iCache]),
                randomPrecisionState[iCache].inverseScale);
            const arma::vec relationshipDiagonal =
              arma::vec(Ai(iR).diag());
            for(arma::uword coordinate = 0;
                coordinate < precisionDiagonal.n_elem;
                ++coordinate){
              const arma::uvec & indices =
                partitionIndexCache[iCache][static_cast<std::size_t>(coordinate)];
              for(arma::uword level = 0; level < relationshipDiagonal.n_elem; ++level){
                matrixFreePCGDiagonal(
                  static_cast<Eigen::Index>(indices(level))
                ) += precisionDiagonal(coordinate) * relationshipDiagonal(level);
              }
            }
          }

          if(!matrixFreePCGDiagonal.allFinite() ||
             matrixFreePCGDiagonal.minCoeff() <= 0.0){
            Rcpp::stop("Matrix-free PCG produced a non-positive/non-finite Jacobi diagonal.");
          }
          pcgPreconditionerDiagonal = matrixFreePCGDiagonal;
          pcgNystromReady = false;

        }else{
          matrixFreePCGReady = false;
          matrixFreeApplyC = nullptr;
          matrixFreeApplyCMatrix = nullptr;

          EigenSpMat Gcontribution(nEffects, nEffects);
          Eigen::VectorXi columnCapacity = Eigen::VectorXi::Zero(nEffects);
          for(int i = 0; i < nRe; ++i){
            const arma::mat lambdaDense = arma::mat(lambda(i));
            const std::size_t iCache = static_cast<std::size_t>(i);
            for(arma::uword iCol = 0; iCol < lambdaDense.n_cols; ++iCol){
              const int nonzeroFactors = static_cast<int>(arma::accu(lambdaDense.col(iCol) != 0.0));
              const arma::uword colStart = partitionStartCache[iCache][static_cast<std::size_t>(iCol)];
              for(arma::uword relCol = 0; relCol < Ai(i).n_cols; ++relCol){
                columnCapacity(static_cast<int>(colStart + relCol)) =
                  nonzeroFactors * static_cast<int>(Ai(i).col(relCol).n_nonzero);
              }
            }
          }
          Gcontribution.reserve(columnCapacity);
          for(int i = 0; i < nRe; ++i){
            const arma::mat lambdaDense = arma::mat(lambda(i));
            const std::size_t iCache = static_cast<std::size_t>(i);
            for(arma::uword iCol = 0; iCol < lambdaDense.n_cols; ++iCol){
              const arma::uword colStart = partitionStartCache[iCache][static_cast<std::size_t>(iCol)];
              for(arma::uword relCol = 0; relCol < Ai(i).n_cols; ++relCol){
                const int globalCol = static_cast<int>(colStart + relCol);
                for(arma::uword iRow = 0; iRow < lambdaDense.n_rows; ++iRow){
                  const double coefficient = lambdaDense(iRow,iCol);
                  if(coefficient == 0.0){ continue; }
                  const arma::uword rowStart = partitionStartCache[iCache][static_cast<std::size_t>(iRow)];
                  for(arma::sp_mat::const_col_iterator relEntry = Ai(i).begin_col(relCol);
                      relEntry != Ai(i).end_col(relCol); ++relEntry){
                    Gcontribution.insert(static_cast<int>(rowStart + relEntry.row()), globalCol) =
                      coefficient * (*relEntry);
                  }
                }
              }
            }
          }
          Gcontribution.makeCompressed();
          if(saddleActive){ saddlePrior = Gcontribution; }
          else { C = C + Gcontribution; }
          C.makeCompressed();
        }
    }


    phaseProfile.enter(sommer::FitPhase::Factorization);
    // ------------------------------------------------------------
    // CHANGE #2: sparse symmetric factorisation of C
    //
    // We use Eigen::SimplicialLDLT because C is the symmetric
    // mixed-model coefficient matrix.  Unlike Armadillo's
    // spsolve_factoriser interface, this factorisation exposes the
    // diagonal D, so log|C| is available without making C dense.
    // ------------------------------------------------------------

    if(
        C.rows() > static_cast<Eigen::Index>(std::numeric_limits<int>::max())
        ||
        C.cols() > static_cast<Eigen::Index>(std::numeric_limits<int>::max())
    ){
      Rcpp::stop(
        "C is too large for the 32-bit Eigen sparse index type used in ai_mme_sp()."
      );
    }

    // Stage 2 Eigen migration: C is already native Eigen sparse, so there is
    // no per-iteration Armadillo<->Eigen copy left before factorisation.
    C.makeCompressed();

    double logDetC = 0.0;
    cholmodCFactorCurrent = false;
    double logDetD = 0.0; // only computed/used when reml==false (Nu==0 => 0)
    // Diagnostic only: minimum LDLT pivot in direct mode; minimum diagonal
    // entry of C in PCG mode. This must not be used by the PCG algorithm.
    double minD = std::numeric_limits<double>::quiet_NaN();
    bool reuseCselectedTopology = false;
    bool reuseDselectedTopology = false;

    // D = Z'R^{-1}Z + G^{-1}, the bottom-right (nX..nEffects-1) block of C,
    // re-indexed to 0..Nu-1. This block never involves X, so it is exactly
    // the same matrix regardless of the fixed-effects design - a single
    // pass over C's stored entries suffices. Shared by both the ldlt and
    // cholmod reml=FALSE (ML) branches below.
    auto buildDFromC = [&]() -> EigenSpMat {
      std::vector<Eigen::Triplet<double>> Dtriplets;
      for(int col = 0; col < static_cast<int>(C.outerSize()); ++col){
        if(col < nX){ continue; }
        for(EigenSpMat::InnerIterator it(C, col); it; ++it){
          if(it.row() < nX){ continue; }
          Dtriplets.emplace_back(
            static_cast<int>(it.row() - nX),
            col - nX,
            it.value()
          );
        }
      }
      EigenSpMat Dmat(Nu, Nu);
      Dmat.setFromTriplets(Dtriplets.begin(), Dtriplets.end());
      Dmat.makeCompressed();
      return Dmat;
    };

    bool randomBlockReady = false;
    auto factorRandomBlock = [&]() -> void {
      if(randomBlockReady){ return; }
      EigenSpMat Dmat = buildDFromC();
      if(solverPlan.ldlt()){
        const auto factorization = sommer::factorizeLDLTWithCachedPattern(
          Dmat, Dfactor, DsymbolicReady, DldltOuterPattern, DldltInnerPattern,
          [&](){ ++DsymbolicAnalyses; DselectedTopologyReady = false; },
          [&](){
            NumericFactorTimer timer(ldltTiming);
            Dfactor.factorize(Dmat);
          });
        if(factorization.status == sommer::SparseFactorizationStatus::SymbolicFailure){
          Rcpp::stop("Sparse symbolic analysis of D failed.");
        }
        if(factorization.status != sommer::SparseFactorizationStatus::Success){
          Rcpp::stop("Sparse LDLT factorisation of the random-effects-only matrix D failed (reml=FALSE).");
        }
        reuseDselectedTopology = factorization.samePattern && DselectedTopologyReady;
        const Eigen::VectorXd pivots = Dfactor.vectorD();
        if(!pivots.allFinite() || pivots.minCoeff() <= 0.0){
          Rcpp::stop("Random-effects-only LDLT produced non-positive/non-finite pivots.");
        }
        logDetD = pivots.array().log().sum();
        minD = pivots.minCoeff();
      }else{
        if(!cholmodDState.analyzePattern(Dmat, &DsymbolicAnalyses)){
          Rcpp::stop("CHOLMOD symbolic analysis of D failed (reml=FALSE).");
        }
        {
          NumericFactorTimer timer(cholmodTiming);
          if(!cholmodDState.factorize(Dmat)){
            Rcpp::stop("CHOLMOD factorisation of D failed (reml=FALSE).");
          }
        }
        logDetD = M_cholmod_factor_ldetA(cholmodDState.factor);
        if(!std::isfinite(logDetD)){ Rcpp::stop("CHOLMOD produced a non-finite D log-determinant."); }
      }
      randomBlockReady = true;
    };

    dFirstMLActive = false;
    bool eligibleDFirst = false;
    const bool needsJointCache = (retainLDLTCache || retainCholmodCache) &&
      !rrConditional && !emInformationOnly;
    if(saddleActive || dFirstMLMode == "off"){
      eligibleDFirst = false;
    } else if(dFirstMLMode == "force"){
      eligibleDFirst = !solverPlan.pcg() && !reml && solverPlan.inverseMode != 1 && Nu > 0 &&
                       !needsJointCache;
    } else {
      eligibleDFirst = solverPlan.dFirstML(static_cast<std::size_t>(nX), static_cast<std::size_t>(Nu),
         blockSchurPolicy.maxDoubles, needsJointCache);
    }
    if(eligibleDFirst){
      factorRandomBlock();
      dFirstCross = C.block(0, nX, nX, Nu);
      double fixedLogDet = 0.0;
      if(nX > 0){
        const Eigen::MatrixXd randomCross = Eigen::MatrixXd(dFirstCross.transpose());
        dFirstTransfer = solveRandomBlock(randomCross, "ML fixed Schur complement");
        Eigen::MatrixXd fixed = Eigen::MatrixXd(C.topLeftCorner(nX, nX));
        fixed.noalias() -= dFirstCross * dFirstTransfer;
        dFirstFixedFactor.compute(0.5 * (fixed + fixed.transpose()));
        if(dFirstFixedFactor.info() == Eigen::Success){
          const Eigen::MatrixXd lower = dFirstFixedFactor.matrixL();
          fixedLogDet = 2.0 * lower.diagonal().array().log().sum();
          dFirstMLActive = std::isfinite(fixedLogDet);
        }
      }else{
        dFirstTransfer.resize(Nu, 0);
        dFirstMLActive = true;
      }
      if(dFirstMLActive){
        ++dFirstMLEvaluations;
        blockEngine.active = false;
        cholmodInvCache.ready = false;
        CselectedTopologyReady = false;
        CnumericReady = true;
        logDetC = logDetD + fixedLogDet;
      }
    }

    if(saddleActive){
      const auto started = std::chrono::steady_clock::now();
      saddleDesign = useH
        ? EigenSpMat(armaSparseToEigenGlobal(Hs).transpose() * W) : W;
      const double augmentedEntries = static_cast<double>(saddleResidual.nonZeros()) +
        2.0 * saddleDesign.nonZeros() + saddlePrior.nonZeros();
      if(augmentedEntries > static_cast<double>(maxFallbackEntries)){
        Rcpp::stop("Saddle-point augmented matrix exceeds the sparse entry limit.");
      }
      saddleJoint.factorize(saddleResidual, saddleDesign, saddlePrior, "joint saddle-point MME");
      logDetC = saddleJoint.logAbsDeterminant() - logDetR;
      if(!reml && Nu > 0){
        const EigenSpMat randomDesign = saddleDesign.rightCols(Nu);
        const EigenSpMat randomPrior = saddlePrior.bottomRightCorner(Nu, Nu);
        saddleRandom.factorize(saddleResidual, randomDesign, randomPrior, "ML random saddle-point MME");
        logDetD = saddleRandom.logAbsDeterminant() - logDetR;
      }
      saddleFactorSeconds += std::chrono::duration<double>(
        std::chrono::steady_clock::now() - started).count();
      CnumericReady = true;
      blockEngine.active = false;
    }else if(!dFirstMLActive && solverPlan.ldlt()){
      const auto factorization = sommer::factorizeLDLTWithCachedPattern(
        C, Cfactor, CsymbolicReady, CouterPattern, CinnerPattern,
        [&](){ ++CsymbolicAnalyses; CselectedTopologyReady = false; },
        [&](){
          NumericFactorTimer timer(ldltTiming);
          Cfactor.factorize(C);
        });
      if(factorization.status == sommer::SparseFactorizationStatus::SymbolicFailure){
        Rcpp::stop("Sparse symbolic analysis of the MME coefficient matrix C failed.");
      }
      if(factorization.status != sommer::SparseFactorizationStatus::Success){
        throw sommer::CovarianceFactorizationError("Sparse LDLT factorisation of the MME coefficient matrix C failed.");
      }
      CnumericReady = true;

      const Eigen::VectorXd Dldlt = Cfactor.vectorD();
      if(Dldlt.size() > 0){ minD = Dldlt.minCoeff(); }
      for(Eigen::Index j = 0; j < Dldlt.size(); ++j){
        const double dj = Dldlt(j);
        if(!std::isfinite(dj) || dj <= 0.0){
          throw sommer::CovarianceFactorizationError("Sparse LDLT produced a non-positive/non-finite pivot in C.");
        }
        logDetC += std::log(dj);
      }
      reuseCselectedTopology = factorization.samePattern && CselectedTopologyReady;

      if(solverPlan.factorRandomOnly(Nu > 0) && !randomBlockReady){
        EigenSpMat Dmat = buildDFromC();

        const auto randomFactorization = sommer::factorizeLDLTWithCachedPattern(
          Dmat, Dfactor, DsymbolicReady, DldltOuterPattern, DldltInnerPattern,
          [&](){ ++DsymbolicAnalyses; DselectedTopologyReady = false; },
          [&](){
            NumericFactorTimer timer(ldltTiming);
            Dfactor.factorize(Dmat);
          });
        if(randomFactorization.status == sommer::SparseFactorizationStatus::SymbolicFailure){
          Rcpp::stop("Sparse symbolic analysis of the random-effects-only matrix D failed (reml=FALSE).");
        }
        if(randomFactorization.status != sommer::SparseFactorizationStatus::Success){
          Rcpp::stop("Sparse LDLT factorisation of the random-effects-only matrix D failed (reml=FALSE).");
        }
        reuseDselectedTopology = randomFactorization.samePattern && DselectedTopologyReady;

        const Eigen::VectorXd DdiagLDLT = Dfactor.vectorD();
        for(Eigen::Index j = 0; j < DdiagLDLT.size(); ++j){
          const double dj = DdiagLDLT(j);
          if(!std::isfinite(dj) || dj <= 0.0){
            Rcpp::stop("Sparse LDLT produced a non-positive/non-finite pivot in D (reml=FALSE).");
          }
          logDetD += std::log(dj);
        }
      }
    }else if(!dFirstMLActive && solverPlan.cholmod()){
      // Supernodal (BLAS-3) direct factorisation via R's Matrix package, or
      // the exact dense block-Schur engine when C's random block is
      // block-diagonal with dense blocks (e.g. diag(env) x dense Gu).
      cholmodInvCache.ready = false;
      blockEngine.active = false;
      if(blockEngineCheckPattern(C)){
        NumericFactorTimer timer(blockSchurTiming);
        blockEngineFactorize(C);
      }
      if(blockEngine.active){
        logDetC = blockEngine.logDet;
        if(!std::isfinite(logDetC)){
          Rcpp::stop("Dense block-Schur engine produced a non-finite log-determinant for C.");
        }
        CnumericReady = true;
        if(verbose && !reportedBlockEngine){
          Rcpp::Rcout
            << (blockEngine.latent ? "Latent-factor Schur engine active (" :
              blockEngine.chain ? "Bordered block-chain engine active (" :
              blockEngine.banded ? "Block-banded Schur engine active (" :
              "Dense block-Schur engine active (")
            << blockEngine.groups.size()
            << " random-effect blocks, "
            << blockEngine.borderCols.size()
            << " border effects)."
            << arma::endl;
          reportedBlockEngine = true;
        }
      }else{
        if(!cholmodState.analyzePattern(C, &CsymbolicAnalyses)){
          Rcpp::stop("CHOLMOD symbolic analysis of the MME coefficient matrix C failed.");
        }

        bool factorizeOk;
        {
          NumericFactorTimer timer(cholmodTiming);
          factorizeOk = cholmodState.factorize(C);
        }
        if(!factorizeOk){
          if(cholmodState.common.status == CHOLMOD_NOT_POSDEF){
            throw sommer::CovarianceFactorizationError(
              "CHOLMOD supernodal factorisation of the MME coefficient matrix C failed.");
          }
          Rcpp::stop("CHOLMOD supernodal factorisation of the MME coefficient matrix C failed.");
        }
        CnumericReady = true;
        cholmodCFactorCurrent = true;

        logDetC = M_cholmod_factor_ldetA(cholmodState.factor);
        if(!std::isfinite(logDetC)){
          throw sommer::CovarianceFactorizationError("CHOLMOD produced a non-finite log-determinant for C.");
        }
      }

      if(solverPlan.factorRandomOnly(Nu > 0) && !randomBlockReady){
        EigenSpMat Dmat = buildDFromC();

        if(!cholmodDState.analyzePattern(Dmat, &DsymbolicAnalyses)){
          Rcpp::stop("CHOLMOD symbolic analysis of the random-effects-only matrix D failed (reml=FALSE).");
        }

        bool factorizeDOk;
        {
          NumericFactorTimer timer(cholmodTiming);
          factorizeDOk = cholmodDState.factorize(Dmat);
        }
        if(!factorizeDOk){
          Rcpp::stop("CHOLMOD supernodal factorisation of the random-effects-only matrix D failed (reml=FALSE).");
        }

        logDetD = M_cholmod_factor_ldetA(cholmodDState.factor);
        if(!std::isfinite(logDetD)){
          Rcpp::stop("CHOLMOD produced a non-finite log-determinant for D (reml=FALSE).");
        }
      }
    }else if(!dFirstMLActive){
      // Genuine factorisation-free MME path: no analyzePattern(), factorize(),
      // vectorD(), matrixL(), or Takahashi call is made for C.
      CnumericReady = false;
      if(pcgEvaluationAccepted){
        pcgAcceptedSolutions = std::move(pcgCandidateSolutions);
        if(krylovPcg && pcgPendingKrylovReady){ pcgKrylov.accept(pcgPendingKrylov); }
      }
      pcgPendingKrylovReady = false;
      if(krylovPcg){ pcgKrylov.clearCandidates(); }
      pcgCandidateSolutions.clear();
      pcgEvaluationAccepted = false;
      const auto pcgSetupStarted = std::chrono::steady_clock::now();
      preparePCG(lambda);
      pcgRecycleBasis.resize(nEffects, 0);
      if(krylovPcg){
        if(pcgImportedBasis.cols() > 0){
          pcgKrylov.capturePanel(pcgImportedBasis);
          if(pcgKrylov.candidateCount() == 0){
            Rcpp::stop("PCG imported recycle basis contains no independent directions.");
          }
          const sommer::pcg::RecycleSpace::Selection imported =
            pcgKrylov.selectCandidates(
              [&](const Eigen::VectorXd & value) -> Eigen::VectorXd {
                Rcpp::checkUserInterrupt();
                return applyCurrentPCG(Eigen::MatrixXd(value)).col(0);
              }, pcgPreconditionerDiagonal, krylovRank);
          if(imported.basis.cols() == 0){
            Rcpp::stop("PCG imported recycle basis contains no independent directions.");
          }
          pcgKrylov.accept(imported);
          pcgImportedBasis.resize(0, 0);
        }
        pcgKrylov.prepare(
          [&](const Eigen::VectorXd & value) -> Eigen::VectorXd {
            Rcpp::checkUserInterrupt();
            return applyCurrentPCG(Eigen::MatrixXd(value)).col(0);
          },
          [&](const Eigen::VectorXd & value) -> Eigen::VectorXd {
            return applyPCGPreconditioner(Eigen::MatrixXd(value)).col(0);
          });
        pcgSetupSeconds += std::chrono::duration<double>(
          std::chrono::steady_clock::now() - pcgSetupStarted).count();
      }else if(pcgRecycleRank > 0 && !pcgAcceptedSolutions.empty()){
        const int requestedRank = std::min(nEffects, pcgRecycleRank);
        const double needed = 4.0 * nEffects * requestedRank + 3.0 * requestedRank * requestedRank;
        if(needed > blockSchurPolicy.maxDoubles){ requireMemory("PCG recycled coarse space", needed); }
        Eigen::MatrixXd candidates(nEffects, requestedRank);
        int count = 0;
        for(const auto & cached : pcgAcceptedSolutions){
          for(Eigen::Index col = 0; col < cached.second.cols() && count < requestedRank; ++col){
            candidates.col(count++) = cached.second.col(col);
          }
        }
        Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(candidates.leftCols(count));
        const int rank = std::min<int>(pcgRecycleRank, qr.rank());
        if(rank > 0){
          pcgRecycleBasis = qr.householderQ() * Eigen::MatrixXd::Identity(nEffects, rank);
          pcgRecycleApplied = applyCurrentPCG(pcgRecycleBasis);
          const Eigen::MatrixXd gram = pcgRecycleBasis.transpose() * pcgRecycleApplied;
          pcgRecycleFactor.compute(0.5 * (gram + gram.transpose()));
          if(pcgRecycleFactor.info() != Eigen::Success){
            Rcpp::stop("PCG recycled coarse operator is not positive definite.");
          }
        }
        pcgSetupSeconds += std::chrono::duration<double>(
          std::chrono::steady_clock::now() - pcgSetupStarted).count();
      }
      // Keep the historical verbose diagnostic column defined without
      // introducing a factorisation in PCG mode.
      if(matrixFreePCGReady){
        minD = matrixFreePCGDiagonal.minCoeff();
      }else if(C.rows() > 0){
        minD = std::numeric_limits<double>::infinity();
        for(Eigen::Index jj = 0; jj < C.rows(); ++jj){
          minD = std::min(minD, C.coeff(jj,jj));
        }
      }
      for(;;){
        logDetC = pcgApproxLogDet(
          C.rows(), [&](const Eigen::VectorXd & value) -> Eigen::VectorXd {
              return matrixFreePCGReady ? matrixFreeApplyC(value) : C * value;
            });
        if(pcgAccuracy != "adaptive" || iIter != 0 || modelFreeCount != 0 ||
           (std::isfinite(pcgLogDetProbeStandardError) && pcgLogDetProbeStandardError <= pcgLogDetTol &&
            pcgLogDetDepthDifference <= pcgLogDetTol) ||
           (pcgTraceProbes == pcgMaxTraceProbes && pcgLanczosSteps == pcgMaxLanczosSteps)){
          break;
        }
        pcgTraceProbes = std::min<long long>(pcgMaxTraceProbes, 2LL * pcgTraceProbes);
        pcgLanczosSteps = std::min<long long>(pcgMaxLanczosSteps, 2LL * pcgLanczosSteps);
        ++pcgAccuracyRefinements;
      }
    }

    // ------------------------------------------------------------
    // Solve C * bu = W' R^-1 y
    // ------------------------------------------------------------

    Eigen::Map<const Eigen::VectorXd> rhsBu(
      rhsMME.memptr(),
      static_cast<Eigen::Index>(rhsMME.n_elem)
    );

    Eigen::VectorXd buEig;
    if(saddleActive){
      const arma::vec transformedResponse = useH ? arma::vec(Hs.t() * yDense) : yDense;
      const Eigen::Map<const Eigen::VectorXd> transformedResponseEig(
        transformedResponse.memptr(), nR);
      buEig = saddleJoint.solveResponse(transformedResponseEig, "BLUE/BLUP calculation");
    }else{
      buEig = solveCVector(rhsBu, "BLUE/BLUP calculation");
    }

    bu.set_size(nEffects);

    for(arma::uword j = 0; j < bu.n_elem; ++j){
      bu(j) = buEig(static_cast<Eigen::Index>(j));
    }

    // ------------------------------------------------------------
    // y' P y from the Schur complement
    //
    // M = [ C  r ]
    //     [ r' d ]
    //
    // y'Py = d - r' C^{-1} r
    // ------------------------------------------------------------

    const double yPy =
      yRiy
      -
      arma::dot(
        rhsMME,
        bu
      );

    if(!std::isfinite(yPy)){
      Rcpp::stop(
        "Non-finite y'Py obtained from the sparse Schur-complement calculation."
      );
    }
    if(profileCommonScale){
      if(commonScaleDF <= 0.0 || !std::isfinite(commonScaleDF) || yPy <= 0.0){
        Rcpp::stop("Common-scale profiling requires positive residual degrees of freedom and positive y'Py.");
      }
      commonScaleQuadraticTotal = yPy;
      commonScale = yPy / commonScaleDF;
      if(!std::isfinite(commonScale) || commonScale <= 0.0){
        Rcpp::stop("Common-scale profiling produced a non-positive or non-finite scale.");
      }
    }

    // ###########################
    // # 1.1) calculate the log-likelihood
    // ###########################

    double llikp = 0.0;

    if(nZs > 0){
      for(int i = 0; i < nRe; ++i){

        double covarianceLogDet =
          randomPrecisionState[static_cast<std::size_t>(i)].logDet;

        if(!randomPrecisionState[static_cast<std::size_t>(i)].factorWise){
          double val;
          double sign;

          bool ok1 =
            arma::log_det(
              val,
              sign,
              theta(i)
            );

          if(ok1 == false){
            Rcpp::Rcout
              << "log determinant failed "
              << arma::endl;
          }
          covarianceLogDet = val * sign;
        }

        llikp =
          llikp
          +
          (
            nUsTotal(i)
            *
            covarianceLogDet
          )
          +
          (
            logDetA(i)
            *
            theta(i).n_rows
          );
      }
    }

    llik(iIter) =
      (-0.5)
      *
      (
        llikp
        +
        (reml ? logDetC : logDetD)
        +
        logDetR
        +
        yPy
      );
    if(profileCommonScale){
      llik(iIter) -= 0.5 * (commonScaleDF * std::log(commonScale) +
        commonScaleDF - yPy);
    }

    phaseProfile.enter(sommer::FitPhase::Acceptance);
    // ==========================================================
    // GLOBAL EXACT-LIKELIHOOD STEP ACCEPTANCE / BACKTRACKING
    // ==========================================================
    //
    // The current likelihood belongs to the CURRENT covPar/theta point.
    // If that point is a proposal from the previous accepted iteration,
    // compare it with the previous accepted likelihood before computing any
    // new score/information quantities.
    //
    // A likelihood decrease is NOT convergence.  It triggers global
    // step-halving in the complete free covariance-parameter vector.
    // ==========================================================

    if(!haveAcceptedLikelihood){
      if(solverName == "pcg"){ pcgAcceptedLogDetValues = pcgLogDetProbeValues; }

      acceptedLikelihood =
        llik(iIter);

      haveAcceptedLikelihood =
        true;

    }else if(pendingLineSearch){
      if(solverName == "pcg" && pcgTraceProbes > 1 &&
         pcgAcceptedLogDetValues.size() == pcgLogDetProbeValues.size()){
        Eigen::VectorXd differences(pcgTraceProbes);
        for(int probe = 0; probe < pcgTraceProbes; ++probe){
          differences(probe) = pcgLogDetProbeValues[probe] - pcgAcceptedLogDetValues[probe];
        }
        pcgPairedLogLikSE = 0.5 * std::sqrt(
          (differences.array() - differences.mean()).square().sum() /
          (static_cast<double>(pcgTraceProbes - 1) * pcgTraceProbes));
      }

      const double trialLikelihood =
        llik(iIter);

      const bool likelihoodAcceptable =
        std::isfinite(trialLikelihood)
        &&
        (
          trialLikelihood
          +
          likelihoodAcceptTolerance
          >=
          acceptedLikelihood
        );

      if(!likelihoodAcceptable){

        if(lineSearchHalvings < maxLineSearchHalvings){

          ++lineSearchHalvings;
          ++totalLineSearchHalvings;
          lineSearchAlpha *= 0.5;

          arma::vec trialParameters =
            lineSearchBase
            +
            lineSearchAlpha
            *
            (
              lineSearchTarget
              -
              lineSearchBase
            );

          // Respect model-defined fixed parameters exactly.
          for(int iStruct = 0; iStruct < nRRe; ++iStruct){

            const arma::uword g0 =
              static_cast<arma::uword>(
                nVcStart(iStruct) - 1
              );

            const arma::uword g1 =
              static_cast<arma::uword>(
                nVcEnd(iStruct) - 1
              );

            arma::uvec local =
              arma::regspace<arma::uvec>(
                g0,
                static_cast<arma::uword>(1),
                g1
              );

            covPar(iStruct) =
              trialParameters(local);

            if(!(iStruct == residualStruct && compactResidualCovariance)){
              theta(iStruct) = evaluateStructure(iStruct, covPar(iStruct));
            }
          }

          // Recompute this SAME optimizer iteration at the smaller trial.
          --iIter;
          continue;

        }else{

          // No improving step was found.  Restore the previous accepted point
          // and evaluate it once more so all final C/b/u objects correspond to
          // an accepted covariance parameter vector.
          for(int iStruct = 0; iStruct < nRRe; ++iStruct){

            const arma::uword g0 =
              static_cast<arma::uword>(
                nVcStart(iStruct) - 1
              );

            const arma::uword g1 =
              static_cast<arma::uword>(
                nVcEnd(iStruct) - 1
              );

            arma::uvec local =
              arma::regspace<arma::uvec>(
                g0,
                static_cast<arma::uword>(1),
                g1
              );

            covPar(iStruct) =
              lineSearchBase(local);

            if(!(iStruct == residualStruct && compactResidualCovariance)){
              theta(iStruct) = evaluateStructure(iStruct, covPar(iStruct));
            }
          }

          pendingLineSearch = false;
          lineSearchStalled = true;
          lineSearchAlpha = 0.0;
          if(useHybridAI && pendingHybridProposal){
            hybridForceExact = true;
            pendingHybridProposal = false;
            lineSearchStalled = false;
            ++hybridAIResets;
          }

          if(verbose){
            Rcpp::Rcout
              << "  No likelihood-improving step found after "
              << maxLineSearchHalvings
              << " halvings; restoring the previous accepted covariance parameters."
              << arma::endl;
          }

          --iIter;
          continue;
        }
      }else{

        acceptedLikelihood =
          trialLikelihood;
        if(solverName == "pcg"){ pcgAcceptedLogDetValues = pcgLogDetProbeValues; }

        pendingLineSearch =
          false;
      }
    }

    if(solverName == "pcg"){ pcgEvaluationAccepted = true; }
    if(!rrConditional){ phaseProfile.enter(sommer::FitPhase::TracePrep); }
    const double traceDerivativePhaseStart = traceDerivativePreparationSeconds;
    const double traceRandomPrecisionStart = traceRandomPrecisionDerivativeSeconds;
    const double traceFactorCovarianceStart = traceFactorCovarianceDerivativeSeconds;
    const double traceFactorPrecisionStart = traceFactorPrecisionDerivativeSeconds;
    const double traceResidualPrecisionStart = traceResidualPrecisionDerivativeSeconds;
    const double traceGridDerivativeStart = traceGridDerivativeBasisSeconds;
    const double traceResidualBasisStart = traceResidualDerivativeBasisSeconds;
    const bool traceDerivativeProfile = phaseProfile.enabled() && !rrConditional;
    ProfileTimer traceDerivativeTimer(traceDerivativePreparationSeconds,
      traceDerivativeProfile);
    // Derivatives and trace/AI preparation are unnecessary for likelihood
    // trials rejected above. Reuse the factor covariance/precision state
    // already evaluated for the accepted likelihood rather than evaluating
    // descriptor factors again.
    {
      ProfileTimer timer(traceRandomPrecisionDerivativeSeconds, traceDerivativeProfile);
      for(int iR = 0; !rrConditional && iR < nRe; ++iR){
        DescriptorPrecisionState & precisionState =
          randomPrecisionState[static_cast<std::size_t>(iR)];
        if(precisionState.factorWise){
          factorEngine.buildDescriptorPrecisionDerivatives(
            covDescriptor[static_cast<std::size_t>(iR)],
            covPar(iR),
            precisionState,
            traceDerivativeProfile ? &traceFactorCovarianceDerivativeSeconds : nullptr,
            traceDerivativeProfile ? &traceFactorPrecisionDerivativeSeconds : nullptr
          );
        }
      }
    }
    traceRandomPrecisionOtherSeconds += std::max(
      0.0,
      (traceRandomPrecisionDerivativeSeconds - traceRandomPrecisionStart) -
        (traceFactorCovarianceDerivativeSeconds - traceFactorCovarianceStart) -
        (traceFactorPrecisionDerivativeSeconds - traceFactorPrecisionStart)
    );

    {
      ProfileTimer timer(traceResidualPrecisionDerivativeSeconds, traceDerivativeProfile);
      if(!rrConditional && Rkron && residualPrecisionState.factorWise){
        factorEngine.buildDescriptorPrecisionDerivatives(
          covDescriptor[static_cast<std::size_t>(residualStruct)],
          covPar(residualStruct),
          residualPrecisionState
        );
      }
    }

    const bool residualDerivativesBlockLocal = (Rkron || Rgrid) && !useH;

    const bool emMomentsOnly = covStructI.hasAttribute("emMomentsOnly") &&
      Rcpp::as<bool>(covStructI.attr("emMomentsOnly"));
    ProfileTimer traceGridDerivativeTimer(traceGridDerivativeBasisSeconds,
      traceDerivativeProfile);
    if(Rgrid && (!rrConditional || (rrCallback && !emMomentsOnly))){
      Rcpp::List rfactors =
        covDescriptor[static_cast<std::size_t>(residualStruct)]["factors"];
      const arma::vec & par = covPar(residualStruct);
      gridDFactor.assign(static_cast<std::size_t>(nResidualPar), -1);
      gridDSection.assign(static_cast<std::size_t>(nResidualPar), -1);
      gridDFull.assign(static_cast<std::size_t>(nResidualPar), arma::mat());
      for(int f = 0; f < gridNFactors; ++f){
        Rcpp::List fl = Rcpp::as<Rcpp::List>(rfactors[f]);
        const int start1 = Rcpp::as<int>(fl["par_start"]);
        const int end1 = Rcpp::as<int>(fl["par_end"]);
        if(end1 < start1){
          continue;
        }
        const arma::vec localPar =
          par.subvec(static_cast<arma::uword>(start1 - 1),
                     static_cast<arma::uword>(end1 - 1));
        const int nSections = factorEngine.factorSections(fl);
        const int nPerSection = nSections > 0 ? Rcpp::as<int>(fl["section_npar"]) : 0;
        std::vector<int> listedToSection;
        for(int t = 0; t < nSections; ++t){
          const int listed = factorEngine.sectionListed(fl, t);
          if(listed >= 0){
            if(static_cast<int>(listedToSection.size()) <= listed){
              listedToSection.resize(static_cast<std::size_t>(listed + 1), -1);
            }
            listedToSection[static_cast<std::size_t>(listed)] = t;
          }
        }
        for(int k1 = start1; k1 <= end1; ++k1){
          const std::size_t k = static_cast<std::size_t>(k1 - 1);
          const int local = k1 - start1;
          gridDFactor[k] = f;
          if(nSections > 0){
            const int t = listedToSection[static_cast<std::size_t>(local / nPerSection)];
            gridDSection[k] = t;
            gridDFull[k] = factorEngine.factorD1(fl, factorEngine.sectionSlice(fl, localPar, t),
                                    static_cast<arma::uword>(local % nPerSection));
          }else{
            gridDFull[k] = factorEngine.factorD1(fl, localPar, static_cast<arma::uword>(local));
          }
        }
      }
    }
    traceGridDerivativeTimer.stop();

    auto invokeConditionalCallback = [&]() {
      if(rrCallback){
        sommer::RRConditionalOperators operators;
        operators.W = eigenSparseToArmaGlobal(W);
        // Callbacks only use C's dimensions; solves go through operators.solve.
        operators.C = arma::sp_mat(static_cast<arma::uword>(C.rows()),
                                   static_cast<arma::uword>(C.cols()));
        operators.y = yDense;
        operators.mean = bu;
        if(useH) { operators.residualTransform = Hs.t(); }
        operators.scale = vary;
        if(emInformationOnly) {
          operators.score = -0.5 * dLuOut;
          operators.information = 0.5 * InfMat;
          if(reducedAIAssembly){ operators.informationTangent = aiTangent; }
          operators.sensitivityColumns = sensitivityRHSColumns;
        }
        operators.residualEngine = Rgrid ? "grid" : Rdiag ? "diagonal"
          : Rkron ? "kronecker" : "sparse";
        operators.precision = applyRiDense;
        operators.precisionDiagonal = effectiveRiDiag;
        operators.solve = [&](const arma::mat & rhs) -> arma::mat {
          Eigen::Map<const Eigen::MatrixXd> mapped(rhs.memptr(), rhs.n_rows, rhs.n_cols);
          const Eigen::MatrixXd solved = solveCMatrix(mapped, "reduced-rank AI");
          if(!solved.allFinite()){ Rcpp::stop("Non-finite reduced-rank AI coefficient solve."); }
          return arma::mat(solved.data(), solved.rows(), solved.cols());
        };
        operators.solveRandom = [&](const arma::mat & rhs) -> arma::mat {
          Eigen::Map<const Eigen::MatrixXd> mapped(rhs.memptr(), rhs.n_rows, rhs.n_cols);
          const Eigen::MatrixXd solved = solveRandomBlock(mapped, "reduced-rank ML AI");
          if(!solved.allFinite()){ Rcpp::stop("Non-finite reduced-rank ML AI coefficient solve."); }
          return arma::mat(solved.data(), solved.rows(), solved.cols());
        };
        bool rrInversePrepared = false, rrInverseAvailable = false;
        if(reml && !saddleActive && !dFirstMLActive && CnumericReady &&
           (solverName == "cholmod" || solverName == "ldlt")){
          operators.inverseEntry = [&](arma::uword row, arma::uword column, double & value) -> bool {
            if(!rrInversePrepared){
              rrInversePrepared = true;
              try {
                if(solverName == "ldlt"){
                  buildSelectedInverseSubset(Cfactor, "C", CselectedTopology, reuseCselectedTopology);
                  CselectedTopologyReady = true;
                  rrInverseAvailable = true;
                } else {
                  if(!cholmodInvCache.ready){ buildCholmodInverseBlockCache(true); }
                  rrInverseAvailable = cholmodInvCache.ready;
                }
              } catch(const std::exception &) {
                rrInverseAvailable = false;
              }
            }
            if(!rrInverseAvailable){ return false; }
            const int i = static_cast<int>(row), j = static_cast<int>(column);
            if(solverName == "ldlt"){
              return getSelectedInverseOriginal(CselectedTopology, i, j, value);
            }
            if(cholmodInverseLookup(i, j, value)){ return true; }
            if(cholmodInvCache.patternOnly){
              // Zero loading coordinates can request entries outside C's
              // pattern. Restore the general cache before exact solve fallback.
              ++patternInverseCacheFallbacks;
              buildCholmodInverseBlockCache();
              return cholmodInverseLookup(i, j, value);
            }
            return false;
          };
        }
        operators.posteriorSolve = [&](const arma::mat & rhs) -> arma::mat {
          Eigen::Map<const Eigen::MatrixXd> mapped(rhs.memptr(), rhs.n_rows, rhs.n_cols);
          const Eigen::MatrixXd solved = sommer::posteriorCoefficientSolve(mapped, nX, reml,
            [&](const Eigen::Ref<const Eigen::MatrixXd> & values) {
              return solveCMatrix(values, "REML posterior moments");
            },
            [&](const Eigen::Ref<const Eigen::MatrixXd> & values) {
              return solveRandomBlock(values, "ML posterior moments");
            });
          if(!solved.allFinite()){ Rcpp::stop("Non-finite posterior moment solve."); }
          return arma::mat(solved.data(), solved.rows(), solved.cols());
        };
        const bool genericEMMoments = covStructI.hasAttribute("emGenericMoments") &&
          Rcpp::as<bool>(covStructI.attr("emGenericMoments"));
        bool emInverseCachePrepared = false;
        auto prepareEMInverseCache = [&]() {
          if(emInverseCachePrepared) { return; }
          if(emInformationOnly) {
            emInverseCachePrepared = true;
            return;
          } else if(solverName == "ldlt") {
            if(reml) {
              buildSelectedInverseSubset(Cfactor, "C", CselectedTopology, reuseCselectedTopology);
              CselectedTopologyReady = true;
            } else if(Nu > 0) {
              buildSelectedInverseSubset(Dfactor, "D", DselectedTopology, reuseDselectedTopology);
              if(reuseDselectedTopology) { ++DselectedTopologyReuses; }
              else { ++DselectedTopologyBuilds; }
              DselectedTopologyReady = true;
            }
          } else if(solverName == "cholmod") {
            if(reml) { buildCholmodInverseBlockCache(); }
            else if(Nu > 0 && buildCholmodTraceSubset(cholmodDState.factor, DselectedTopology)) {
              ++cholmodDSelectedBuilds;
            }
          }
          emInverseCachePrepared = true;
        };
        if(!genericEMMoments) {
          operators.posteriorRandomTraces = [&](int term, bool diagonal) -> arma::mat {
            prepareEMInverseCache();
            ++emInverseTraceCalls;
            EigenSpMat & relationship = AiEigenCache[static_cast<std::size_t>(term)];
            if(relationship.rows() != static_cast<Eigen::Index>(Ai(term).n_rows)) {
              relationship = armaSparseToEigenGlobal(Ai(term));
              relationship.makeCompressed();
            }
            const arma::mat & termPartitions = partitions(term);
            std::vector<Eigen::Index> starts;
            for(arma::uword column = 0; column < termPartitions.n_rows; ++column) {
              starts.push_back(static_cast<Eigen::Index>(termPartitions(column, 0)) - 1 -
                (reml ? 0 : nX));
            }
            const Eigen::MatrixXd traces = sommer::exactInverseBlockTraces(
              reml ? nEffects : Nu, starts, relationship,
              sommer::WorkspacePolicy::panelColumns(),
              [&](const Eigen::MatrixXd & rhs) -> Eigen::MatrixXd {
                Rcpp::checkUserInterrupt();
                ++traceFallbackSolves;
                ++emInverseTraceFallbackPanels;
                return reml ? solveCMatrix(rhs, "REML inverse-block trace fallback")
                  : solveRandomBlock(rhs, "ML inverse-block trace fallback");
              },
              [&](std::size_t row, std::size_t column) { return !diagonal || row == column; },
              [&](Eigen::Index row, Eigen::Index column, double & value) {
                if(!reml) {
                  return getSelectedInverseOriginal(DselectedTopology, row, column, value);
                }
                return solverName == "cholmod"
                  ? cholmodInverseLookup(row, column, value)
                  : getSelectedInverseOriginal(CselectedTopology, row, column, value);
              });
            return arma::mat(traces.data(), traces.rows(), traces.cols());
          };
          operators.posteriorTrace = [&](const arma::sp_mat & matrix) {
            prepareEMInverseCache();
            ++emDesignTraceCalls;
            bool fallback = false;
            if(!reml) {
              if(Nu == 0) { return 0.0; }
              const arma::sp_mat randomMatrix = matrix.submat(nX, nX, nEffects - 1, nEffects - 1);
              return sparseTraceInverseTimes(randomMatrix, Dfactor, DselectedTopology,
                "ML posterior design trace fallback", fallback,
                [&](const Eigen::Ref<const Eigen::MatrixXd> & rhs,
                    const std::string & context) -> Eigen::MatrixXd {
                  Rcpp::checkUserInterrupt();
                  ++emInverseTraceFallbackPanels;
                  return solveRandomBlock(rhs, context);
                });
            }
            return sparseTraceInverseTimes(matrix, Cfactor, CselectedTopology,
              "EM design trace fallback", fallback,
              [&](const Eigen::Ref<const Eigen::MatrixXd> & rhs,
                  const std::string & context) -> Eigen::MatrixXd {
                Rcpp::checkUserInterrupt();
                ++emInverseTraceFallbackPanels;
                return solveCMatrix(rhs, context);
              }, solverName == "cholmod");
          };
        }
        Eigen::SparseMatrix<double, Eigen::RowMajor> rrRowDesign;
        if(residualPlan.elementwisePrecision() && !emMomentsOnly) {
          operators.residualDiagonalDerivative = [&](arma::uword parameter) -> arma::vec {
            const arma::mat & derivative = cachedCovarianceD1(residualStruct, parameter);
            arma::vec diagonal(nR);
            for(int row = 0; row < nR; ++row) {
              const int coordinate = residualLocalOfRow[static_cast<std::size_t>(row)];
              diagonal(row) = derivative(coordinate, coordinate);
            }
            return diagonal;
          };
          operators.residualDiagonalContraction = [&]() -> arma::vec {
            prepareEMInverseCache();
            ++diagonalResidualTraceCalls;
            EigenSpMat design = reml ? W : EigenSpMat(W.rightCols(Nu));
            const Eigen::VectorXd quadratic = sommer::exactInverseRowQuadratics(
              design, sommer::WorkspacePolicy::panelColumns(),
              [&](const Eigen::MatrixXd & rhs) -> Eigen::MatrixXd {
                Rcpp::checkUserInterrupt();
                diagonalResidualTraceSolveColumns += rhs.cols();
                return reml ? solveCMatrix(rhs, "diagonal residual contraction")
                  : solveRandomBlock(rhs, "ML diagonal residual contraction");
              },
              [&](Eigen::Index row, Eigen::Index column, double & value) {
                if(!reml) {
                  return getSelectedInverseOriginal(DselectedTopology, row, column, value);
                }
                return solverName == "cholmod"
                  ? cholmodInverseLookup(row, column, value)
                  : getSelectedInverseOriginal(CselectedTopology, row, column, value);
              });
            const arma::vec correction(quadratic.data(), quadratic.size());
            return RdiagInv - (RdiagInv % effectiveRiDiag % correction);
          };
        }
        if(Rgrid && !emMomentsOnly){
          rrRowDesign = W;
          operators.derivative = [&](arma::uword parameter, const arma::mat & values) -> arma::mat {
            arma::mat output(values.n_rows, values.n_cols, arma::fill::zeros);
            for(int block = 0; block < residualNBlocks; ++block){
              const arma::uvec & rows = residualKronRows[static_cast<std::size_t>(block)];
              output.rows(rows) = gridApplyDR(block, parameter, values.rows(rows));
            }
            return output;
          };
          operators.residualContraction = [&](arma::uword block) -> arma::mat {
            const arma::uvec & rows = residualKronRows[block];
            arma::mat design(rows.n_elem, reml ? nEffects : Nu, arma::fill::zeros);
            for(arma::uword local = 0; local < rows.n_elem; ++local){
              for(Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator entry(rrRowDesign, rows(local)); entry; ++entry){
                if(reml || entry.col() >= nX){
                  design(local, entry.col() - (reml ? 0 : nX)) = entry.value();
                }
              }
            }
            const arma::mat precisionDesign = gridApplyRi(block, design);
            arma::mat contraction(rows.n_elem, rows.n_elem);
            const arma::uword panel = sommer::WorkspacePolicy::panelColumns();
            for(arma::uword first = 0; first < rows.n_elem; first += panel){
              Rcpp::checkUserInterrupt();
              const arma::uword count = std::min(panel, rows.n_elem - first);
              arma::mat basis(rows.n_elem, count, arma::fill::zeros);
              for(arma::uword column = 0; column < count; ++column){ basis(first + column, column) = 1.0; }
              const arma::mat rhs = precisionDesign.rows(first, first + count - 1).t();
              const arma::mat solved = reml ? operators.solve(rhs) : operators.solveRandom(rhs);
              contraction.cols(first, first + count - 1) =
                gridApplyRi(block, basis) - precisionDesign * solved;
            }
            return contraction;
          };
        }
        (*rrCallback)(operators);
      }
    };
    if(rrConditional){
      invokeConditionalCallback();
      b = bu(bInd);
      if(nZs > 0){ u = bu(uInd); }
      arma::uword offset = 0;
      for(int term = 0; term < nRRe; ++term){
        monitor.submat(offset, 0, offset + covPar(term).n_elem - 1, 0) = covPar(term);
        offset += covPar(term).n_elem;
      }
      avInf.zeros();
      InfMat.zeros();
      InfMatInv.zeros();
      normMonitor.zeros();
      dLuOut.zeros(nVcTotal);
      convergence = true;
      break;
    }

    auto applyGroupedWeightDerivative = [&](const arma::sp_mat & derivative) -> arma::sp_mat {
      std::map<std::pair<int,int>, Eigen::MatrixXd> sourceBlocks;
      std::vector<int> localOfRow(static_cast<std::size_t>(nR), -1);
      for(int block = 0; block < weightBlockCount; ++block){
        const std::vector<int> & rows = weightBlockRows[static_cast<std::size_t>(block)];
        for(int local = 0; local < static_cast<int>(rows.size()); ++local){
          localOfRow[static_cast<std::size_t>(rows[static_cast<std::size_t>(local)])] = local;
        }
      }
      for(arma::sp_mat::const_iterator entry = derivative.begin(); entry != derivative.end(); ++entry){
        const int globalRow = static_cast<int>(entry.row());
        const int globalCol = static_cast<int>(entry.col());
        if(globalRow > globalCol){ continue; }
        int blockRow = weightBlockOfRow[static_cast<std::size_t>(globalRow)];
        int blockCol = weightBlockOfRow[static_cast<std::size_t>(globalCol)];
        int localRow = localOfRow[static_cast<std::size_t>(globalRow)];
        int localCol = localOfRow[static_cast<std::size_t>(globalCol)];
        if(blockRow > blockCol){
          std::swap(blockRow, blockCol);
          std::swap(localRow, localCol);
        }
        const auto key = std::make_pair(blockRow, blockCol);
        auto found = sourceBlocks.find(key);
        if(found == sourceBlocks.end()){
          found = sourceBlocks.emplace(key, Eigen::MatrixXd::Zero(
            weightUpperBlocks[static_cast<std::size_t>(blockRow)].rows(),
            weightUpperBlocks[static_cast<std::size_t>(blockCol)].rows())).first;
        }
        found->second(localRow, localCol) += *entry;
        if(blockRow == blockCol && localRow != localCol){
          found->second(localCol, localRow) += *entry;
        }
      }
      std::vector<Eigen::Triplet<double>> transformedTriplets;
      for(const auto & source : sourceBlocks){
        const int blockRow = source.first.first;
        const int blockCol = source.first.second;
        const Eigen::MatrixXd left = weightUpperBlocks[static_cast<std::size_t>(blockRow)]
          .transpose().triangularView<Eigen::Lower>().solve(source.second);
        const Eigen::MatrixXd transformed = weightUpperBlocks[static_cast<std::size_t>(blockCol)]
          .transpose().triangularView<Eigen::Lower>().solve(left.transpose()).transpose();
        const std::vector<int> & rows = weightBlockRows[static_cast<std::size_t>(blockRow)];
        const std::vector<int> & cols = weightBlockRows[static_cast<std::size_t>(blockCol)];
        for(Eigen::Index col = 0; col < transformed.cols(); ++col){
          for(Eigen::Index row = 0; row < transformed.rows(); ++row){
            const double value = transformed(row, col);
            if(value == 0.0){ continue; }
            if(transformedTriplets.size() >= maxFallbackEntries){
              Rcpp::stop("Weighted residual derivatives produce excessive block fill; reduce weight-block sizes or use a structured residual model.");
            }
            transformedTriplets.emplace_back(rows[static_cast<std::size_t>(row)],
                                              cols[static_cast<std::size_t>(col)], value);
            if(blockRow != blockCol){
              transformedTriplets.emplace_back(cols[static_cast<std::size_t>(col)],
                                                rows[static_cast<std::size_t>(row)], value);
            }
          }
        }
      }
      Eigen::SparseMatrix<double> transformed(derivative.n_rows, derivative.n_cols);
      transformed.setFromTriplets(transformedTriplets.begin(), transformedTriplets.end());
      transformed.makeCompressed();
      return eigenSparseToArmaGlobal(transformed);
    };

    {
      ProfileTimer timer(traceResidualDerivativeBasisSeconds, traceDerivativeProfile);
      for(arma::uword k = 0; k < nResidualPar && !Rgrid; ++k){
        if(compressedData){
          residualLocalD1(k) = cachedCovarianceD1(residualStruct, k);
          continue;
        }
        if(compactResidualCovariance && !residualDerivativesBlockLocal){
          residualDerivativeBasis(k) =
            buildCompactResidualSparse(covPar(residualStruct), static_cast<int>(k));
        }else{
          if(!(residualDerivativesBlockLocal && residualPrecisionState.factorWise)){
            residualLocalD1(k) =
              cachedCovarianceD1(
                residualStruct,
                k
              );
          }
          if(!residualDerivativesBlockLocal){
            residualDerivativeBasis(k) =
              buildResidualSparseFromPattern(residualLocalD1(k));
          }
        }
        // The effective residual is R = Hs^{-T} R0 Hs^{-1}; its derivative
        // carries the same weighting (the trace code relies on this).
        if(useH){
          residualDerivativeBasis(k) = weightBlockCount > 0
            ? applyGroupedWeightDerivative(residualDerivativeBasis(k))
            : HsInvT * residualDerivativeBasis(k) * HsInvT.t();
        }
      }
    }

    traceDerivativeTimer.stop();
    traceDerivativeOtherSeconds += std::max(
      0.0,
      (traceDerivativePreparationSeconds - traceDerivativePhaseStart) -
        (traceRandomPrecisionDerivativeSeconds - traceRandomPrecisionStart) -
        (traceResidualPrecisionDerivativeSeconds - traceResidualPrecisionStart) -
        (traceGridDerivativeBasisSeconds - traceGridDerivativeStart) -
        (traceResidualDerivativeBasisSeconds - traceResidualBasisStart)
    );
    ProfileTimer tracePreparationTimer(traceSelectedInverseSeconds, phaseProfile.enabled());
    if(saddleActive){
      // Exact augmented inverse actions supply traces without an SPD C factor.
    }else if(solverPlan.traces() == sommer::TracePreparation::SelectedInverse){
      if(!dFirstMLActive){
        buildSelectedInverseSubset(
          Cfactor,
          "C",
          CselectedTopology,
          reuseCselectedTopology
        );
        CselectedTopologyReady = true;
      }

      if(!reml && Nu > 0){
        buildSelectedInverseSubset(
          Dfactor,
          "D",
          DselectedTopology,
          reuseDselectedTopology
        );
        if(reuseDselectedTopology){ ++DselectedTopologyReuses; }
        else{ ++DselectedTopologyBuilds; }
        DselectedTopologyReady = true;
      }
    }else if(solverPlan.traces() == sommer::TracePreparation::InverseBlocks){
      buildCholmodInverseBlockCache();
    }else if(solverPlan.traces() == sommer::TracePreparation::Probes){
      preparePCGTraceProbes(static_cast<Eigen::Index>(nEffects));
    }
    if(!saddleActive && solverName == "cholmod"){
      if(!reml && Nu > 0 && buildCholmodTraceSubset(cholmodDState.factor, DselectedTopology)){
        ++cholmodDSelectedBuilds;
      }
      if(!Rdiag && !Rgrid && !Rkron && buildCholmodTraceSubset(cholmodRState.factor, Rselected)){
        ++cholmodRSelectedBuilds;
      }
    }
    tracePreparationTimer.stop();

    b = bu(bInd); // move BLUEs to a different vector
    if(nZs > 0){ // move BLUPs to a different vector
      u = bu(uInd);
    }
    phaseProfile.enter(sommer::FitPhase::Sensitivities);
    // ============================================================
    // CHANGE 8: analytic sparse-factor differentiation for the AI matrix
    // ============================================================
    //
    // The previous implementation constructed one observation-space
    // working variate per variance component,
    //
    //        w_i = V_i P y,
    //
    // and then evaluated
    //
    //        AI_ij = w_i' P w_j.
    //
    // That formulation is exact, but for many RANDOM covariance
    // parameters it creates an nR x nVC working-variate matrix and sends
    // all of those columns through R^{-1}.
    //
    // Here we forward-differentiate the already-factorised mixed-model
    // equations instead.  If
    //
    //        C b = r,
    //
    // then for variance parameter theta_i
    //
    //        C b_i = r_i - C_i b.
    //
    // The SAME sparse LDLT factorisation of C is reused for all sensitivity
    // right-hand sides.  For random-effect covariance parameters r_i = 0
    // and C_i is confined to the corresponding G^{-1} block, so no
    // observation-space random working variate is needed.
    //
    // The AI matrix is recovered from one half of the exact second
    // derivative of y' P y.  Because the marginal covariance V is linear
    // in the variance parameters,
    //
    //        0.5 d^2(y'Py)/(d theta_i d theta_j)
    //          = y' P V_i P V_j P y
    //          = AI_ij.
    //
    // For residual covariance parameters we retain the compact residual
    // working variates only.  This is intentional: the number of residual
    // covariance parameters is normally small, and this avoids explicitly
    // constructing second derivatives of R^{-1}.  Cross random/residual
    // terms are obtained from the differentiated coefficient system.
    //
    // This is an exact analytic forward differentiation of the sparse
    // factorised MME solve; no finite differences or numerical derivatives
    // are used.
    // ============================================================

    // Number/index at which residual variance components begin.
    const arma::uword residualVcStart =
      static_cast<arma::uword>(
        nVcStart(nRRe - 1) - 1
      );

    const arma::uword nRandomVc =
      residualVcStart;
    const arma::mat residualAITangent = reducedAIAssembly && aiDimension > 0
      ? aiTangent.rows(residualVcStart, nVcTotal - 1) : arma::mat();
    const arma::uvec residualAICoordinates = reducedAIAssembly && aiDimension > 0
      ? arma::find(arma::sum(arma::abs(residualAITangent), 0).t() > 0.0) : arma::uvec();
    const arma::mat packedResidualTangent = reducedAIAssembly && residualAICoordinates.n_elem > 0
      ? residualAITangent.cols(residualAICoordinates) : arma::mat(nResidualPar, 0);
    const arma::uword residualAIDimension = reducedAIAssembly
      ? residualAICoordinates.n_elem : nResidualPar;

    arma::vec e;
    if(!compressedData){
      const Eigen::Map<const Eigen::VectorXd> buEig(
        bu.memptr(), static_cast<Eigen::Index>(bu.n_elem)
      );
      const Eigen::VectorXd WbuEig = W * buEig;
      e = yDense - arma::vec(WbuEig.data(), WbuEig.size());
    }
    arma::vec Rie;
    if(compressedData){
      const Eigen::Map<const Eigen::VectorXd> coefficients(bu.memptr(), nEffects);
      compressedFitted = sufficient.design * coefficients;
      compressedSquares = sufficient.residualSquares(compressedFitted);
      compressedDerivative.resize(residualDim, nResidualPar);
      for(arma::uword parameter = 0; parameter < nResidualPar; ++parameter){
        for(int category = 0; category < residualDim; ++category){
          compressedDerivative(category, parameter) = residualLocalD1(parameter)(category, category);
        }
      }
    }else{
      Rie = applyRiVec(e);
    }
    std::vector<arma::mat> randomQuadraticBase(static_cast<std::size_t>(nRe));

    bool exactAIThisIteration = !useHybridAI || !hybridState.ready ||
      hybridForceExact || lineSearchStalled || iIter == nIters - 1 ||
      iIter - lastExactAIIteration >= hybridRefreshInterval;
    const bool fullAIOutput = !useHybridAI || hybridForceExact ||
      lineSearchStalled || iIter == nIters - 1;
    bool completeAIThisIteration = fullAIOutput ||
      hybridFree.n_elem == static_cast<arma::uword>(nVcTotal);
    hybridForceExact = false;

    auto computeAverageInformation = [&](bool required, bool reduceFixed){
    phaseProfile.enter(sommer::FitPhase::Sensitivities);

    // Forward-sensitivity right-hand sides:
    //
    //   random component i:    -(dC_i) b
    //   residual component i:  -W'Ri S_i Ri e
    //
    // Solving C * db = sensitivityRHS gives db/dtheta.
    arma::mat sensitivityRHS(
      nEffects,
      required ? aiDimension : 0,
      arma::fill::zeros
    );

    // Cache the small covariance-space objects needed for the random/random
    // curvature correction:
    //
    // d Lambda_i = -Lambda B_i Lambda
    //
    // d2 Lambda_ij =
    //     Lambda B_j Lambda B_i Lambda
    //   + Lambda B_i Lambda B_j Lambda.
    std::vector< std::vector<arma::mat> > randomDLambda(
      static_cast<std::size_t>(nRe)
    );

    if(nZs > 0){

      for(int iR = 0; iR < nRe; ++iR){

        const std::size_t iCache =
          static_cast<std::size_t>(iR);

        arma::mat partitionsP =
          partitions(iR);

        const arma::uword nLevels =
          static_cast<arma::uword>(
            partitionsP(0,1)
            -
            partitionsP(0,0)
            +
            1
          );

        const arma::uword nTraits =
          static_cast<arma::uword>(
            partitionsP.n_rows
          );

        arma::mat U(
          nLevels,
          nTraits,
          arma::fill::zeros
        );

        for(arma::uword iCol = 0; iCol < nTraits; ++iCol){

          const arma::uvec & idx =
            partitionIndexCache[
              iCache
            ][
              static_cast<std::size_t>(iCol)
            ];

          if(idx.n_elem != nLevels){
            Rcpp::stop(
              "Random-effect partition length mismatch in the "
              "factor-differentiation AI engine."
            );
          }

          U.col(iCol) =
            bu(idx);
        }

        arma::mat lambdaDense =
          arma::mat(
            lambda(iR)
          );

        // A * U is reused for every local covariance derivative.
        arma::mat AU =
          arma::mat(
            Ai(iR)
            *
            U
          );

        randomQuadraticBase[iCache] =
          U.t()
          *
          AU;

        const arma::uword localNvc =
          static_cast<arma::uword>(
            nVc(iR)
          );

        randomDLambda[iCache].resize(
          static_cast<std::size_t>(localNvc)
        );

        DescriptorPrecisionState & precisionState =
          randomPrecisionState[iCache];

        if(!precisionState.factorWise){
          precisionState.precisionD1.resize(
            static_cast<std::size_t>(localNvc)
          );
          precisionState.logDetD1.resize(
            static_cast<std::size_t>(localNvc)
          );
        }

        for(arma::uword localK = 0; localK < localNvc; ++localK){

          arma::mat dLambda;
          if(!precisionState.factorWise){
            const arma::mat & Bk =
              cachedCovarianceD1(
                iR,
                localK
              );
            dLambda =
              -lambdaDense
              *
              Bk
              *
              lambdaDense;
            precisionState.precisionD1[
              static_cast<std::size_t>(localK)
            ] = dLambda;
            precisionState.logDetD1[
              static_cast<std::size_t>(localK)
            ] =
              arma::accu(
                lambdaDense
                %
                Bk.t()
              );
            dLambda = 0.5 * (dLambda + dLambda.t());
            randomDLambda[iCache][
              static_cast<std::size_t>(localK)
            ] = dLambda;
          }

          if(!required || aiDimension == 0){ continue; }

          // (dC_k) b for this random covariance structure is
          //
          //      (dLambda_k kron A) b
          //
          // and can be evaluated blockwise as
          //
          //      (A U) dLambda_k'.
          arma::mat dCbuLocal =
            precisionState.factorWise
            ? factorEngine.applyDescriptorPrecisionD1Right(
                AU,
                precisionState,
                localK
              )
            : AU * dLambda.t();

          const arma::uword globalK =
            static_cast<arma::uword>(
              nVcStart(iR) - 1
            )
            +
            localK;

          for(arma::uword trait = 0; trait < nTraits; ++trait){

            const arma::uvec & idx =
              partitionIndexCache[
                iCache
              ][
                static_cast<std::size_t>(trait)
              ];

            for(arma::uword rr = 0; rr < idx.n_elem; ++rr){
              if(reducedAIAssembly){
                sensitivityRHS.row(idx(rr)) -= dCbuLocal(rr, trait) * aiTangent.row(globalK);
              }else{
                sensitivityRHS(idx(rr), globalK) = -dCbuLocal(rr, trait);
              }
            }
          }
        }
      }
    }

    if(!required){ return; }
    if(aiDimension == 0){
      avInf.zeros(0, 0);
      sensitivityRHSColumnsSkipped += nVcTotal;
      return;
    }

    // ------------------------------------------------------------
    // Residual sensitivities.
    // ------------------------------------------------------------
    // Only residual-parameter working variates are retained.  AR1 has
    // two columns (sigma2, rho) no matter how many lag bases define R.
    arma::mat residualWorking(
      compressedData ? 0 : nR,
      residualAIDimension,
      arma::fill::zeros
    );

    for(arma::uword iP = 0; !compressedData && residualAIDimension > 0 && iP < nResidualPar; ++iP){
      if(residualDerivativesBlockLocal){
        for(int block = 0; block < residualKronNBlocks; ++block){
          const arma::uvec & rows =
            residualKronRows[static_cast<std::size_t>(block)];
          const arma::vec localRie = Rie.elem(rows);
          arma::vec localWorking;
          if(Rgrid){
            localWorking = gridApplyDR(block, iP, localRie);
          }else if(residualPrecisionState.factorWise){
            localWorking = factorEngine.applyDescriptorCovarianceD1Right(
              localRie.t(), residualPrecisionState, iP
            ).t();
          }else{
            localWorking = residualLocalD1(iP) * localRie;
          }
          for(arma::uword localRow = 0; localRow < rows.n_elem; ++localRow){
            if(reducedAIAssembly){
              residualWorking.row(rows(localRow)) +=
                localWorking(localRow) * packedResidualTangent.row(iP);
            }else{
              residualWorking(rows(localRow), iP) = localWorking(localRow);
            }
          }
        }
      }else{
        const arma::vec working =
          arma::vec(
            residualDerivativeBasis(iP)
            *
            Rie
          );
        if(reducedAIAssembly){
          residualWorking += working * packedResidualTangent.row(iP);
        }else{
          residualWorking.col(iP) = working;
        }
      }
    }

    arma::mat RiResidualWorking;
    if(!compressedData){ RiResidualWorking = applyRiDense(residualWorking); }

    arma::mat residualCrossRHS;
    if(compressedData){
      Eigen::MatrixXd groupedRhs(sufficient.groups.size(),
        residualAIDimension);
      for(std::size_t group = 0; residualAIDimension > 0 && group < sufficient.groups.size(); ++group){
        const auto & statistics = sufficient.groups[group];
        const int category = statistics.residualClass;
        const double variance = compressedVariance(category);
        const Eigen::RowVectorXd derivative = compressedDerivative.row(category);
        if(reducedAIAssembly){
          const Eigen::Map<const Eigen::MatrixXd> mapping(
            packedResidualTangent.memptr(), packedResidualTangent.n_rows,
            packedResidualTangent.n_cols);
          groupedRhs.row(group) = (derivative * mapping) *
            (statistics.weight * (statistics.mean - compressedFitted(group)) / (variance * variance));
        }else{
          groupedRhs.row(group) = derivative *
            (statistics.weight * (statistics.mean - compressedFitted(group)) / (variance * variance));
        }
      }
      const Eigen::MatrixXd cross = sufficient.design.transpose() * groupedRhs;
      residualCrossRHS = arma::mat(cross.data(), cross.rows(), cross.cols());
    }else{
      const Eigen::Map<const Eigen::MatrixXd> RiResidualWorkingEig(
        RiResidualWorking.memptr(),
        static_cast<Eigen::Index>(RiResidualWorking.n_rows),
        static_cast<Eigen::Index>(RiResidualWorking.n_cols)
      );
      const Eigen::MatrixXd residualCrossRHSEig = W.transpose() * RiResidualWorkingEig;
      residualCrossRHS.set_size(residualCrossRHSEig.rows(), residualCrossRHSEig.cols());
      std::copy(residualCrossRHSEig.data(), residualCrossRHSEig.data() + residualCrossRHSEig.size(), residualCrossRHS.memptr());
    }

    if(reducedAIAssembly){
      sensitivityRHS.cols(residualAICoordinates) -= residualCrossRHS;
    }else for(arma::uword iP = 0; iP < nResidualPar; ++iP){

      const arma::uword globalK =
        residualVcStart
        +
        iP;

      sensitivityRHS.col(globalK) =
        -residualCrossRHS.col(iP);
    }

    // ------------------------------------------------------------
    // Differentiate the sparse MME solve.
    //
    // C is already numerically factorised.  Solve all derivative RHSs
    // simultaneously using that same sparse LDLT factor.
    // ------------------------------------------------------------
    Eigen::Map<
      const Eigen::Matrix<
        double,
        Eigen::Dynamic,
        Eigen::Dynamic,
        Eigen::ColMajor
      >
    > sensitivityRHSEig(
      sensitivityRHS.memptr(),
      static_cast<Eigen::Index>(
        sensitivityRHS.n_rows
      ),
      static_cast<Eigen::Index>(
        sensitivityRHS.n_cols
      )
    );

    if(phaseProfile.enabled()){
      for(arma::uword column = 0; column < sensitivityRHS.n_cols; ++column){
        for(arma::uword row = 0; row < sensitivityRHS.n_rows; ++row){
          sensitivityRHSNonzeros += sensitivityRHS(row, column) != 0.0;
        }
      }
      if(solverName == "ldlt" && !saddleActive && !dFirstMLActive &&
         Cfactor.info() == Eigen::Success &&
         Cfactor.rows() == static_cast<Eigen::Index>(sensitivityRHS.n_rows)){
        const std::size_t dimension = static_cast<std::size_t>(sensitivityRHS.n_rows);
        EigenSpMat lower = Cfactor.matrixL();
        lower.makeCompressed();
        std::vector<std::vector<int>> forwardRows(dimension);
        std::vector<std::vector<int>> backwardColumns(dimension);
        for(int column = 0; column < lower.outerSize(); ++column){
          for(EigenSpMat::InnerIterator entry(lower, column); entry; ++entry){
            if(entry.row() > column){
              forwardRows[static_cast<std::size_t>(column)].push_back(entry.row());
              backwardColumns[static_cast<std::size_t>(entry.row())].push_back(column);
            }
          }
        }
        std::vector<int> forwardSeen(dimension, 0);
        std::vector<int> backwardSeen(dimension, 0);
        const auto & permutation = Cfactor.permutationP();
        const bool reducedSensitivityColumns =
          reduceFixed && hybridFree.n_elem < static_cast<arma::uword>(nVcTotal);
        const arma::uword reachColumnCount = reducedSensitivityColumns
          ? hybridFree.n_elem
          : sensitivityRHS.n_cols;
        for(arma::uword reachColumn = 0; reachColumn < reachColumnCount; ++reachColumn){
          const arma::uword column = reducedSensitivityColumns
            ? hybridFree(reachColumn)
            : reachColumn;
          const int marker = static_cast<int>(reachColumn + 1);
          std::vector<int> forward;
          for(arma::uword row = 0; row < sensitivityRHS.n_rows; ++row){
            if(sensitivityRHS(row, column) == 0.0){ continue; }
            const int permuted = permutation.indices()(static_cast<int>(row));
            if(forwardSeen[static_cast<std::size_t>(permuted)] != marker){
              forwardSeen[static_cast<std::size_t>(permuted)] = marker;
              forward.push_back(permuted);
            }
          }
          for(std::size_t cursor = 0; cursor < forward.size(); ++cursor){
            for(const int row : forwardRows[static_cast<std::size_t>(forward[cursor])]){
              if(forwardSeen[static_cast<std::size_t>(row)] != marker){
                forwardSeen[static_cast<std::size_t>(row)] = marker;
                forward.push_back(row);
              }
            }
          }
          std::vector<int> solveReach = forward;
          for(const int row : solveReach){
            backwardSeen[static_cast<std::size_t>(row)] = marker;
          }
          for(std::size_t cursor = 0; cursor < solveReach.size(); ++cursor){
            for(const int parent : backwardColumns[
                  static_cast<std::size_t>(solveReach[cursor])]){
              if(backwardSeen[static_cast<std::size_t>(parent)] != marker){
                backwardSeen[static_cast<std::size_t>(parent)] = marker;
                solveReach.push_back(parent);
              }
            }
          }
          sensitivitySolveReachEntries += solveReach.size();
          sensitivitySolveReachPossibleEntries += dimension;
        }
      }
    }
    ProfileTimer sensitivitySolveTimer(sensitivitySolveSeconds, phaseProfile.enabled());
    Eigen::MatrixXd dBuEig;
    if(reduceFixed && !reducedAIAssembly && hybridFree.n_elem < static_cast<arma::uword>(nVcTotal)){
      // Only the free/free block is needed here; final refreshes solve every column.
      dBuEig = Eigen::MatrixXd::Zero(nEffects, nVcTotal);
      if(hybridFree.n_elem > 0){
        const arma::mat packed = sensitivityRHS.cols(hybridFree);
        const Eigen::Map<const Eigen::MatrixXd> mapped(packed.memptr(), packed.n_rows, packed.n_cols);
        const Eigen::MatrixXd solved = solveCMatrix(mapped, "free AI sensitivity equations");
        for(arma::uword column = 0; column < hybridFree.n_elem; ++column){
          dBuEig.col(hybridFree(column)) = solved.col(column);
        }
      }
      sensitivityRHSColumns += hybridFree.n_elem;
      sensitivityRHSColumnsSkipped += nVcTotal - hybridFree.n_elem;
    }else{
      dBuEig = solveCMatrix(
        sensitivityRHSEig,
        "factor-differentiation AI sensitivity equations"
      );
      sensitivityRHSColumns += aiDimension;
      sensitivityRHSColumnsSkipped += nVcTotal - aiDimension;
    }
    sensitivitySolveTimer.stop();

    // Zero-copy Armadillo view of Eigen's column-major result.  dBuEig owns
    // the memory and remains alive for the complete AI assembly below.
    arma::mat dBu(
      dBuEig.data(),
      static_cast<arma::uword>(dBuEig.rows()),
      static_cast<arma::uword>(dBuEig.cols()),
      false,
      true
    );

    // ------------------------------------------------------------
    // Assemble the exact Average Information matrix.
    // ------------------------------------------------------------
    phaseProfile.enter(sommer::FitPhase::AIAssembly);
    ProfileTimer aiAssemblyTimer(aiAssemblyWallSeconds, phaseProfile.enabled());
    avInf.zeros(aiDimension, aiDimension);
    if(reducedAIAssembly){
      // Pull back -B' C^{-1} B before adding local curvature and residual Gram terms.
      avInf = -sensitivityRHS.t() * dBu;
    }
    std::vector<arma::uvec> aiCoordinates;
    if(reducedAIAssembly){
      aiCoordinates.resize(nVcTotal);
      for(int parameter = 0; parameter < nVcTotal; ++parameter){
        aiCoordinates[parameter] = arma::find(aiTangent.row(parameter).t() != 0.0);
      }
    }
    const auto addCurvature = [&](arma::uword i, arma::uword j, double value){
      if(!reducedAIAssembly){ avInf(i, j) += value; return; }
      for(const arma::uword row : aiCoordinates[i]){
        for(const arma::uword column : aiCoordinates[j]){
          avInf(row, column) += value * aiTangent(i, row) * aiTangent(j, column);
        }
      }
    };

    // Random rows:
    //
    // AI_ij =
    //      b' C_i b_j
    //    + 0.5 b' C_ij b
    //
    // The first term for every j is obtained in one dense BLAS product,
    // because sensitivityRHS_i = -C_i b.
    if(nRandomVc > 0){

      ProfileTimer randomFirstTermTimer(aiRandomFirstTermSeconds, phaseProfile.enabled());
      if(!reducedAIAssembly){
      arma::mat randomFirstTerm =
        (
          -sensitivityRHS.cols(
            static_cast<arma::uword>(0),
            nRandomVc - 1
          )
        ).t()
        *
        dBu;

      avInf.submat(
        static_cast<arma::uword>(0),
        static_cast<arma::uword>(0),
        nRandomVc - 1,
        static_cast<arma::uword>(nVcTotal - 1)
      ) =
        randomFirstTerm;
      }
      randomFirstTermTimer.stop();

      // Add 0.5 b' C_ij b for pairs of covariance parameters belonging
      // to the SAME random covariance structure.  Cross-structure second
      // derivatives are exactly zero.
      ProfileTimer randomCurvatureTimer(aiRandomCurvatureSeconds, phaseProfile.enabled());
    #ifdef _OPENMP
      #pragma omp parallel for if(nRe > 1 && !reducedAIAssembly) schedule(static)
    #endif
      for(int iR = 0; iR < nRe; ++iR){

        const std::size_t iCache =
          static_cast<std::size_t>(iR);

        const arma::mat & quadBase =
          randomQuadraticBase[iCache];

        const arma::uword localNvc =
          static_cast<arma::uword>(
            nVc(iR)
          );

        const arma::uword globalStart =
          static_cast<arma::uword>(
            nVcStart(iR) - 1
          );

        if(randomStructurallyDiagonal[iCache]){
          // Closed-form fast path: for a structurally diagonal covariance
          // product, Sigma and every dLambda/dphi_k are diagonal, so
          //   0.5 d2Lambda_ij = dLambda_i Sigma dLambda_j
          // is an elementwise product.
          // This replaces localNvc^2 dense q x q matrix-chain products
          // (O(q^3) each) with localNvc^2 O(q) dot products.
          const arma::vec covarianceDiag =
            theta(iR).diag();
          const arma::vec quadDiag =
            quadBase.diag();

          std::vector<arma::vec> dLambdaDiag(localNvc);
          const DescriptorPrecisionState & precisionState =
            randomPrecisionState[iCache];
          for(arma::uword localK = 0; localK < localNvc; ++localK){
            if(precisionState.factorWise){
              std::vector<arma::mat> derivativeFactors;
              double derivativeScale = 1.0;
              factorEngine.descriptorPrecisionD1Factors(
                precisionState,
                localK,
                derivativeFactors,
                derivativeScale
              );
              dLambdaDiag[localK] =
                factorEngine.kroneckerDiagonal(derivativeFactors, derivativeScale);
            }else{
              dLambdaDiag[localK] =
                randomDLambda[iCache][
                  static_cast<std::size_t>(localK)
                ].diag();
            }
          }

#ifdef _OPENMP
          #pragma omp parallel for if(nRe <= 1 && localNvc > 1 && !reducedAIAssembly) schedule(static)
#endif
          for(arma::uword localI = 0; localI < localNvc; ++localI){
            for(arma::uword localJ = 0; localJ < localNvc; ++localJ){
              if(reducedAIAssembly &&
                 (aiCoordinates[globalStart + localI].is_empty() ||
                  aiCoordinates[globalStart + localJ].is_empty())){ continue; }

              const double secondCorrection =
                arma::accu(
                  dLambdaDiag[localI]
                  %
                  dLambdaDiag[localJ]
                  %
                  covarianceDiag
                  %
                  quadDiag
                );

              addCurvature(globalStart + localI, globalStart + localJ, secondCorrection);
            }
          }

          continue;
        }

        const DescriptorPrecisionState & precisionState = randomPrecisionState[iCache];
      #ifdef _OPENMP
        #pragma omp parallel for if(nRe <= 1 && localNvc > 1 && !reducedAIAssembly) schedule(static)
      #endif
        for(arma::uword localI = 0; localI < localNvc; ++localI){
          const int factorI = precisionState.derivativeFactor[static_cast<std::size_t>(localI)];

            // The covariance curvature correction is symmetric in each parameter pair.
            for(arma::uword localJ = localI; localJ < localNvc; ++localJ){
            if(reducedAIAssembly &&
               (aiCoordinates[globalStart + localI].is_empty() ||
                aiCoordinates[globalStart + localJ].is_empty())){ continue; }

            // Average-information metric in generic coordinates:
            // J' AI(Sigma) J.  Build dLambda_j Sigma dLambda_i from
            // factor-level products when the descriptor precision path is
            // available, avoiding a full q x q dense matrix chain.
            double secondCorrection = 0.0;

            if(precisionState.factorWise){
              const int factorJ =
                precisionState.derivativeFactor[
                  static_cast<std::size_t>(localJ)
                ];

              std::vector<arma::mat> crossFactors;
              double crossScale = precisionState.inverseScale;

              if(factorI < 0 && factorJ < 0){
                crossFactors = factorEngine.factorPrecisionFactorsDense(precisionState);
              }else if(factorI < 0){
                factorEngine.descriptorPrecisionD1Factors(
                  precisionState,
                  localJ,
                  crossFactors,
                  crossScale
                );
                crossScale = -crossScale;
              }else if(factorJ < 0){
                factorEngine.descriptorPrecisionD1Factors(
                  precisionState,
                  localI,
                  crossFactors,
                  crossScale
                );
                crossScale = -crossScale;
              }else{
                crossFactors = factorEngine.factorPrecisionFactorsDense(precisionState);
                for(std::size_t factor = 0;
                    factor < precisionState.factorPrecision.size();
                    ++factor){
                  arma::mat piece;
                  if(static_cast<int>(factor) == factorI &&
                     static_cast<int>(factor) == factorJ){
                    const std::size_t parameterI = static_cast<std::size_t>(localI);
                    const std::size_t parameterJ = static_cast<std::size_t>(localJ);
                    if(parameterI < precisionState.sparseLocalPrecisionD1.size() &&
                       parameterJ < precisionState.sparseLocalPrecisionD1.size() &&
                       !precisionState.sparseLocalPrecisionD1[parameterI].is_empty() &&
                       !precisionState.sparseLocalPrecisionD1[parameterJ].is_empty()){
                      piece = precisionState.sparseLocalPrecisionD1[parameterJ] *
                        precisionState.factorCovariance[factor] *
                        precisionState.sparseLocalPrecisionD1[parameterI];
                    }else{
                      piece = factorEngine.localPrecisionD1Dense(precisionState, parameterJ) *
                        precisionState.factorCovariance[factor] *
                        factorEngine.localPrecisionD1Dense(precisionState, parameterI);
                    }
                  }else if(static_cast<int>(factor) == factorI){
                    piece = factorEngine.localPrecisionD1Dense(precisionState, localI);
                  }else if(static_cast<int>(factor) == factorJ){
                    piece = factorEngine.localPrecisionD1Dense(precisionState, localJ);
                  }else{
                    piece = factorEngine.factorPrecisionDense(precisionState, factor);
                  }
                  crossFactors[factor] = std::move(piece);
                }
              }

              secondCorrection =
                factorEngine.contractKronecker(
                  quadBase,
                  crossFactors,
                  crossScale
                );
            }else{
              const arma::mat & dLambdaI =
                randomDLambda[iCache][
                  static_cast<std::size_t>(localI)
                ];
              const arma::mat & dLambdaJ =
                randomDLambda[iCache][
                  static_cast<std::size_t>(localJ)
                ];
              const arma::mat crossTerm =
                dLambdaJ
                *
                theta(iR)
                *
                dLambdaI;
              secondCorrection =
                0.5
                *
                arma::accu(
                  (crossTerm + crossTerm.t())
                  %
                  quadBase
                );
            }

            addCurvature(globalStart + localI, globalStart + localJ, secondCorrection);
            if(localI != localJ){
              addCurvature(globalStart + localJ, globalStart + localI, secondCorrection);
            }
          }
        }
      }

      // Numerical roundoff can make the random/random block very slightly
      // asymmetric even though the analytical expression is symmetric.
      if(!reducedAIAssembly){
      arma::mat randomBlock =
        avInf.submat(
          static_cast<arma::uword>(0),
          static_cast<arma::uword>(0),
          nRandomVc - 1,
          nRandomVc - 1
        );

      randomBlock =
        0.5
        *
        (
          randomBlock
          +
          randomBlock.t()
        );

      avInf.submat(
        static_cast<arma::uword>(0),
        static_cast<arma::uword>(0),
        nRandomVc - 1,
        nRandomVc - 1
      ) =
        randomBlock;
      }
      randomCurvatureTimer.stop();
    }

    // Residual/residual block:
    //
    // w_i' R^{-1} w_j - q_i' C^{-1} q_j
    //
    // and because db_j = -C^{-1} q_j,
    //
    // = w_i' R^{-1} w_j + q_i' db_j.
    ProfileTimer residualBlockTimer(aiResidualBlockSeconds, phaseProfile.enabled());
    arma::mat residualGram;
    if(compressedData){
      Eigen::MatrixXd gram = Eigen::MatrixXd::Zero(residualAIDimension, residualAIDimension);
      for(std::size_t group = 0; group < sufficient.groups.size(); ++group){
        const int category = sufficient.groups[group].residualClass;
        const double variance = compressedVariance(category);
        Eigen::RowVectorXd derivative = compressedDerivative.row(category);
        if(reducedAIAssembly){
          const Eigen::Map<const Eigen::MatrixXd> mapping(
            packedResidualTangent.memptr(), packedResidualTangent.n_rows,
            packedResidualTangent.n_cols);
          derivative = (derivative * mapping).eval();
        }
        gram.noalias() += (compressedSquares(group) / (variance * variance * variance)) *
          derivative.transpose() * derivative;
      }
      residualGram = arma::mat(gram.data(), gram.rows(), gram.cols());
    }else{
      residualGram = residualWorking.t() * RiResidualWorking;
    }

    if(reducedAIAssembly){
      avInf.submat(residualAICoordinates, residualAICoordinates) += residualGram;
      avInf = 0.5 * (avInf + avInf.t());
    }else{
    arma::mat residualSensitivity =
      dBu.cols(
        residualVcStart,
        static_cast<arma::uword>(
          nVcTotal - 1
        )
      );

    arma::mat residualAI =
      residualGram
      +
      residualCrossRHS.t()
      *
      residualSensitivity;

    residualAI =
      0.5
      *
      (
        residualAI
        +
        residualAI.t()
      );

    avInf.submat(
      residualVcStart,
      residualVcStart,
      static_cast<arma::uword>(
        nVcTotal - 1
      ),
      static_cast<arma::uword>(
        nVcTotal - 1
      )
    ) =
      residualAI;

    // Random/residual cross block was already calculated in the random
    // rows via b' C_i b_j.  Copy it explicitly to the residual/random
    // block instead of halving it through a blanket symmetrisation.
    if(nRandomVc > 0 && nResidualPar > 0){

      avInf.submat(
        residualVcStart,
        static_cast<arma::uword>(0),
        static_cast<arma::uword>(
          nVcTotal - 1
        ),
        nRandomVc - 1
      ) =
        avInf.submat(
          static_cast<arma::uword>(0),
          residualVcStart,
          nRandomVc - 1,
          static_cast<arma::uword>(
            nVcTotal - 1
          )
        ).t();
    }
    }
    residualBlockTimer.stop();

    if(!avInf.is_finite()){
      Rcpp::stop(
        "Non-finite Average Information matrix produced by the "
        "factor-differentiation engine."
      );
    }
    ++exactAIRefreshes;
    aiAssemblyTimer.stop();
    };

    if(!exactAIThisIteration){
      avInf = hybridState.metric;
      sensitivityRHSColumnsSkipped += nVcTotal;
    }
    computeAverageInformation(exactAIThisIteration, !fullAIOutput);

    // ##########################
    // # 5) get 1st derivatives (dL/ds2i) from MME-version
    // # PAPER FORMULA (Lee and Van der Werf, 2006)
    // # dL/ds2u = -0.5 [(Nu/s2u) - (tr(AiCuu)/s4u) -  (e/s2e)'(Zu/s2u)]
    // # dL/ds2e = -0.5 [((Nr-Nb)/s2e) - [(Nu - (tr(AiCuu)/s2u))*(1/s2e)] - ... - (e/s2e)'(e/s2e)]
    // #
    // # PAPER FORMULA (Jensen and Madsen, 1997)
    // # dL/ds2u = (q.i * lambda) - (lambda * (T + S) * lambda)  Eq. 18
    // # dL/ds2e = tr(Rij*Ri) - tr(Ci*W'*Ri*Rij*Ri*W) - (e'*Ri*Rij*Ri*e)
    // ###########################
    
    phaseProfile.enter(sommer::FitPhase::ScoreTraces);
    // Score traces use the sparse inverse subset; no full C^{-1}
    // is formed during iterations.
    arma::field<arma::mat> emInfList(nRRe);
    arma::vec dLu;
    double saddleInversePriorTrace = 0.0;
    if(solverName == "pcg"){
      pcgScoreProbeValues = Eigen::MatrixXd::Zero(nVcTotal, pcgTraceProbes);
    }

    // ============================================================
    // CHANGE #3A: selected inverse blocks for random-effect traces
    // ============================================================
    if(nZs > 0){

      for(int iR = 0; iR < nRe; ++iR){

        arma::mat traces(
          lambda(iR).n_rows,
          lambda(iR).n_cols,
          arma::fill::zeros
        );

        arma::mat partitionsP = partitions(iR);
        std::vector<arma::mat> pcgRandomProbeTraces;
        if(solverName == "pcg"){
          const double needed = static_cast<double>(pcgTraceProbes) * traces.n_elem;
          if(needed > blockSchurPolicy.maxDoubles){ requireMemory("PCG random trace diagnostics", needed); }
          pcgRandomProbeTraces.assign(pcgTraceProbes,
            arma::mat(traces.n_rows, traces.n_cols, arma::fill::zeros));
        }

        arma::umat traceSupport;
        {
          const DescriptorPrecisionState & state =
            randomPrecisionState[static_cast<std::size_t>(iR)];
          traceSupport.zeros(lambda(iR).n_rows, lambda(iR).n_cols);
          for(arma::uword parameter = 0; parameter < covPar(iR).n_elem; ++parameter){
            if(state.factorWise){
              std::vector<arma::mat> derivativeFactors;
              double derivativeScale = 1.0;
              factorEngine.descriptorPrecisionD1Factors(state, parameter, derivativeFactors, derivativeScale);
              arma::umat support(1, 1, arma::fill::ones);
              for(const arma::mat & factor : derivativeFactors){
                support = arma::kron(support, arma::umat(factor != 0.0));
              }
              traceSupport += support;
            }else{
              traceSupport += arma::umat(state.precisionD1[static_cast<std::size_t>(parameter)] != 0.0);
            }
          }
        }

        const SelectedInverseSubset & traceSubset = reml ? CselectedTopology : DselectedTopology;
        const EigenLDLT & traceFactor = reml ? Cfactor : Dfactor;
        const int traceDimension = reml ? nEffects : Nu;
        const int traceRowOffset = reml ? 0 : nX;
        bool batchTraces = !saddleActive && solverName == "ldlt" && traceBatchMode == "force";
        int traceWorkers = 1;
      #ifdef _OPENMP
        traceWorkers = omp_get_max_threads();
      #endif
        if(batchTraces){
          const int maxWidth = std::max(static_cast<int>(sommer::WorkspacePolicy::panelColumns()),
            static_cast<int>(partitionsP(0, 1) - partitionsP(0, 0) + 1));
          const double neededBytes =
            3.0 * traceDimension * static_cast<double>(maxWidth) * traceWorkers * sizeof(double);
          const double retainedBytes = static_cast<double>(CselectedTopology.dependencyEntries +
            DselectedTopology.dependencyEntries) * sizeof(std::size_t);
          if(neededBytes + retainedBytes > blockSchurPolicy.maxDoubles * sizeof(double)){
            if(memoryStrict){ requireMemory("Batched exact trace solves", (neededBytes + retainedBytes) / sizeof(double)); }
            batchTraces = false;
          }else{
            traceBatchWorkspaceBytes = std::max(traceBatchWorkspaceBytes, neededBytes);
          }
        }
        const int traceColumnCount = static_cast<int>(lambda(iR).n_cols);
        // Dense block-Schur cross blocks are C^{-1}_gh = F_g S^{-1} F_h', so
        // tr(Ai C^{-1}_gh) = accu(F_g S^{-1} % (Ai F_h)): one BLAS-3 product per
        // block instead of a border-length dot product per stored Ai entry.
        arma::umat crossTraceReady(traceColumnCount, traceColumnCount, arma::fill::zeros);
        arma::mat crossTraceValues(traceColumnCount, traceColumnCount, arma::fill::zeros);
        if(reml && !saddleActive && solverName == "cholmod" && blockEngine.active &&
           !blockEngine.chain && !blockEngine.banded && cholmodInvCache.ready &&
           !cholmodInvCache.selReady && cholmodInvCache.crossAvailable &&
           !randomStructurallyDiagonal[static_cast<std::size_t>(iR)] && traceColumnCount > 1){
          std::vector<int> columnGroup(static_cast<std::size_t>(traceColumnCount), -1);
          std::vector<arma::uvec> columnLocal(static_cast<std::size_t>(traceColumnCount));
          for(int column = 0; column < traceColumnCount; ++column){
            const arma::uword start = static_cast<arma::uword>(partitionsP(column, 0) - 1);
            const arma::uword width = static_cast<arma::uword>(partitionsP(column, 1) - partitionsP(column, 0) + 1);
            if(width != Ai(iR).n_rows || Ai(iR).n_rows != Ai(iR).n_cols){ columnGroup.assign(columnGroup.size(), -1); break; }
            arma::uvec local(width);
            int group = cholmodInvCache.groupOf[start];
            for(arma::uword j = 0; j < width && group >= 0; ++j){
              const std::size_t index = static_cast<std::size_t>(start + j);
              if(cholmodInvCache.borderPos[index] >= 0 || cholmodInvCache.groupOf[index] != group){
                group = -1;
              }else{
                local(j) = static_cast<arma::uword>(cholmodInvCache.localOf[index]);
              }
            }
            columnGroup[static_cast<std::size_t>(column)] = group;
            columnLocal[static_cast<std::size_t>(column)] = std::move(local);
          }
          std::vector<char> needColumn(static_cast<std::size_t>(traceColumnCount), 0);
          for(int column = 0; column < traceColumnCount; ++column){
            for(int row = 0; row < column; ++row){
              const int g = columnGroup[static_cast<std::size_t>(row)];
              const int h = columnGroup[static_cast<std::size_t>(column)];
              if(traceSupport(row, column) != 0 && g >= 0 && h >= 0 && g != h){
                needColumn[static_cast<std::size_t>(row)] = needColumn[static_cast<std::size_t>(column)] = 1;
              }
            }
          }
          const double border = static_cast<double>(blockEngine.borderCols.size());
          const double neededDoubles = 2.0 * std::count(needColumn.begin(), needColumn.end(), 1) *
            static_cast<double>(Ai(iR).n_rows) * border;
          if(border > 0.0 && neededDoubles <= blockSchurPolicy.maxDoubles){
            const bool denseAi = static_cast<double>(Ai(iR).n_nonzero) >
              0.1 * static_cast<double>(Ai(iR).n_rows) * static_cast<double>(Ai(iR).n_cols);
            const arma::mat AiDense = denseAi ? arma::mat(Ai(iR)) : arma::mat();
            std::vector<arma::mat> leftPanels(static_cast<std::size_t>(traceColumnCount));
            std::vector<arma::mat> rightPanels(static_cast<std::size_t>(traceColumnCount));
            for(int column = 0; column < traceColumnCount; ++column){
              if(!needColumn[static_cast<std::size_t>(column)]){ continue; }
              const std::size_t g = static_cast<std::size_t>(columnGroup[static_cast<std::size_t>(column)]);
              const arma::uvec & local = columnLocal[static_cast<std::size_t>(column)];
              leftPanels[static_cast<std::size_t>(column)] = cholmodInvCache.FS[g].rows(local);
              const arma::mat factorRows = blockEngine.Fg[g].rows(local);
              rightPanels[static_cast<std::size_t>(column)] = denseAi
                ? arma::mat(AiDense * factorRows) : arma::mat(Ai(iR) * factorRows);
            }
            for(int column = 0; column < traceColumnCount; ++column){
              for(int row = 0; row < column; ++row){
                const int g = columnGroup[static_cast<std::size_t>(row)];
                const int h = columnGroup[static_cast<std::size_t>(column)];
                if(traceSupport(row, column) == 0 || g < 0 || h < 0 || g == h){ continue; }
                // sum_{a,b} Ai(a,b) C^{-1}(row a, column b).
                crossTraceValues(row, column) = arma::accu(
                  leftPanels[static_cast<std::size_t>(row)] %
                  rightPanels[static_cast<std::size_t>(column)]);
                crossTraceReady(row, column) = 1;
              }
            }
          }
        }
        if(batchTraces){
          const int blockWidth = static_cast<int>(partitionsP(0, 1) - partitionsP(0, 0) + 1);
          for(int column = 0; column < traceColumnCount; ++column){
            if(partitionsP(column, 1) - partitionsP(column, 0) + 1 != blockWidth ||
               Ai(iR).n_rows != static_cast<arma::uword>(blockWidth) ||
               Ai(iR).n_cols != static_cast<arma::uword>(blockWidth)){
              Rcpp::stop("Random-effect inverse block dimensions are inconsistent with Ai.");
            }
          }
          if(traceFactor.info() != Eigen::Success){ Rcpp::stop("LDLT factor is not ready for batched exact traces."); }
        }
        const int traceColumnBatch = batchTraces
            ? std::max(1, static_cast<int>(sommer::WorkspacePolicy::panelColumns()) /
              static_cast<int>(partitionsP(0, 1) - partitionsP(0, 0) + 1))
          : std::max(1, traceColumnCount);
        if(saddleActive){
          std::vector<Eigen::Index> starts;
          for(int column = 0; column < traceColumnCount; ++column){
            const Eigen::Index width = partitionsP(column, 1) - partitionsP(column, 0) + 1;
            if(width != static_cast<Eigen::Index>(Ai(iR).n_rows)){
              Rcpp::stop("Random-effect saddle trace block dimensions are inconsistent with Ai.");
            }
            starts.push_back(static_cast<Eigen::Index>(partitionsP(column, 0)) - 1 - traceRowOffset);
          }
          const Eigen::MatrixXd contracted = sommer::exactInverseBlockTraces(
            traceDimension, starts, AiEigenCache[static_cast<std::size_t>(iR)],
            sommer::WorkspacePolicy::panelColumns(),
            [&](const Eigen::MatrixXd & rhs) -> Eigen::MatrixXd {
              Rcpp::checkUserInterrupt();
              ++traceFallbackSolves;
              return reml ? solveCMatrix(rhs, "REML saddle random traces")
                : solveRandomBlock(rhs, "ML saddle random traces");
            },
            [&](std::size_t row, std::size_t column) -> bool {
              return traceSupport(row, column) != 0 &&
                (!randomStructurallyDiagonal[static_cast<std::size_t>(iR)] || row == column);
            });
          traces = arma::mat(contracted.data(), contracted.rows(), contracted.cols());
        }else{
        const int traceBatches = (traceColumnCount + traceColumnBatch - 1) / traceColumnBatch;
        std::vector<int> traceBatchFailed(static_cast<std::size_t>(traceBatches), 0);
#ifdef _OPENMP
        #pragma omp parallel for if(batchTraces && traceBatches > 1) schedule(static) reduction(+:traceFallbackSolves,traceBatchSolves)
#endif
        for(int traceBatch = 0; traceBatch < traceBatches; ++traceBatch){
          const int traceBatchStart = traceBatch * traceColumnBatch;
          const int traceBatchEnd = std::min(traceBatchStart + traceColumnBatch, traceColumnCount);
          std::vector<int> traceOffsets(static_cast<std::size_t>(traceBatchEnd - traceBatchStart), -1);
          Eigen::MatrixXd tracePanel;
          int panelWidth = 0;
          if(batchTraces){
            for(int column = traceBatchStart; column < traceBatchEnd; ++column){
              bool missing = false;
              const int columnStart = static_cast<int>(partitionsP(column, 0)) - 1 - traceRowOffset;
              for(int row = 0; row <= column && !missing; ++row){
                if(traceSupport(row, column) == 0 ||
                      (randomStructurallyDiagonal[static_cast<std::size_t>(iR)] && row != column)){
                  continue;
                }
                const int rowStart = static_cast<int>(partitionsP(row, 0)) - 1 - traceRowOffset;
                for(arma::sp_mat::const_iterator entry = Ai(iR).begin(); entry != Ai(iR).end(); ++entry){
                  double value = 0.0;
                  if(!getSelectedInverseOriginal(traceSubset, columnStart + entry.col(), rowStart + entry.row(), value)){
                    missing = true;
                    break;
                  }
                }
              }
              if(missing){
                traceOffsets[static_cast<std::size_t>(column - traceBatchStart)] = panelWidth;
                panelWidth += static_cast<int>(partitionsP(column, 1) - partitionsP(column, 0) + 1);
                ++traceFallbackSolves;
              }
            }
            if(panelWidth > 0){
              Eigen::MatrixXd traceRHS = Eigen::MatrixXd::Zero(traceDimension, panelWidth);
              for(int column = traceBatchStart; column < traceBatchEnd; ++column){
                const int offset = traceOffsets[static_cast<std::size_t>(column - traceBatchStart)];
                if(offset < 0){ continue; }
                const int columnStart = static_cast<int>(partitionsP(column, 0)) - 1 - traceRowOffset;
                const int width = static_cast<int>(partitionsP(column, 1) - partitionsP(column, 0) + 1);
                for(int local = 0; local < width; ++local){
                  traceRHS(columnStart + local, offset + local) = 1.0;
                }
              }
              tracePanel = traceFactor.solve(traceRHS);
              if(!tracePanel.allFinite()){
                traceBatchFailed[static_cast<std::size_t>(traceBatch)] = 1;
                continue;
              }
              ++traceBatchSolves;
            }
          }

      #ifdef _OPENMP
        #pragma omp parallel for if(!saddleActive && !batchTraces && traceBatchEnd - traceBatchStart > 1) schedule(static)
      #endif
        for(int iCol = traceBatchStart; iCol < traceBatchEnd; ++iCol){

          const arma::uword colStart =
            static_cast<arma::uword>(partitionsP(iCol, 0) - 1);
          const arma::uword colEnd =
            static_cast<arma::uword>(partitionsP(iCol, 1) - 1);
          const arma::uword blockWidth = colEnd - colStart + 1;

          bool fallbackBlockAvailable = false;
          Eigen::MatrixXd fallbackBlockSolution;
          Eigen::MatrixXd pcgAiZ;

          if(solverName == "pcg"){
            const Eigen::SparseMatrix<double> & AiEigen =
              AiEigenCache[static_cast<std::size_t>(iR)];
            if(AiEigen.cols() != static_cast<Eigen::Index>(blockWidth)){
              Rcpp::stop("Random-effect inverse block dimensions are inconsistent with Ai.");
            }
            pcgAiZ.noalias() =
              AiEigen
              *
              pcgTraceZ.middleRows(
                static_cast<Eigen::Index>(colStart),
                static_cast<Eigen::Index>(blockWidth)
              );
          }

          // traces is symmetric and symmatu() below copies this upper
          // triangle into the lower triangle.
          for(int iRow = 0; iRow <= iCol; ++iRow){

            if(traceSupport(iRow, iCol) == 0){
              continue;
            }
            if(crossTraceReady(iRow, iCol) != 0){
              traces(iRow, iCol) = crossTraceValues(iRow, iCol);
              continue;
            }

            // For a structurally diagonal covariance product, lambda and
            // every dSigma/dphi_k are diagonal matrices, so the score's
            // lambda*traces*lambda term only ever reads traces' diagonal
            // entries. Off-diagonal (iRow != iCol) blocks are exactly
            // unused and skipping them avoids the expensive Ai/selected-
            // inverse work below for q^2-q of the q^2 (row,col) pairs.
            if(
                randomStructurallyDiagonal[static_cast<std::size_t>(iR)]
                &&
                iRow != iCol
            ){
              continue;
            }

            const arma::uword rowStart =
              static_cast<arma::uword>(partitionsP(iRow, 0) - 1);
            const arma::uword rowEnd =
              static_cast<arma::uword>(partitionsP(iRow, 1) - 1);
            const arma::uword blockHeight = rowEnd - rowStart + 1;

            double trAiCuu = 0.0;

            if(solverName == "pcg"){
              if(pcgAiZ.rows() != static_cast<Eigen::Index>(blockHeight)){
                Rcpp::stop("Random-effect inverse block dimensions are inconsistent with Ai.");
              }
              trAiCuu =
                (
                  pcgAiZ.array()
                  *
                  pcgTraceX.middleRows(
                    static_cast<Eigen::Index>(rowStart),
                    static_cast<Eigen::Index>(blockHeight)
                  ).array()
                ).sum()
                /
                static_cast<double>(pcgTraceProbes);
              for(int probe = 0; probe < pcgTraceProbes; ++probe){
                const double sample = pcgAiZ.col(probe).dot(
                  pcgTraceX.col(probe).segment(rowStart, blockHeight));
                pcgRandomProbeTraces[probe](iRow,iCol) = sample;
                pcgRandomProbeTraces[probe](iCol,iRow) = sample;
              }
            }else{
              bool subsetComplete = true;
              for(arma::sp_mat::const_iterator ait = Ai(iR).begin();
                  ait != Ai(iR).end(); ++ait){
                const arma::uword ar = ait.row();
                const arma::uword ac = ait.col();
                if(ar >= blockHeight || ac >= blockWidth){
                  Rcpp::stop("Random-effect inverse block dimensions are inconsistent with Ai.");
                }
                double cij = 0.0;
                const bool haveEntry = solverName == "cholmod"
                  ? (reml ? cholmodInverseLookup(
                      static_cast<int>(rowStart + ar),
                      static_cast<int>(colStart + ac), cij)
                    : getSelectedInverseOriginal(DselectedTopology,
                      static_cast<int>(rowStart + ar - nX),
                      static_cast<int>(colStart + ac - nX), cij))
                  : getSelectedInverseOriginal(
                      traceSubset,
                      static_cast<int>(colStart + ac - (reml ? 0 : nX)),
                      static_cast<int>(rowStart + ar - (reml ? 0 : nX)),
                      cij);
                if(!haveEntry){
                  subsetComplete = false;
                  break;
                }
                trAiCuu += (*ait) * cij;
              }

              if(!subsetComplete){
                const int panelOffset = traceOffsets[static_cast<std::size_t>(iCol - traceBatchStart)];
                if(panelOffset >= 0){
                  trAiCuu = 0.0;
                  for(arma::sp_mat::const_iterator entry = Ai(iR).begin(); entry != Ai(iR).end(); ++entry){
                    trAiCuu += (*entry) * tracePanel(
                      rowStart + entry.row() - traceRowOffset,
                      panelOffset + entry.col());
                  }
                  traces(iRow, iCol) = trAiCuu;
                  continue;
                }
                if(!fallbackBlockAvailable){
#ifdef _OPENMP
                  #pragma omp critical(sommer_factor_solve)
#endif
                  {
                    ProfileTimer fallbackTimer(scoreTraceFallbackSeconds, phaseProfile.enabled());
                    ++traceFallbackSolves;
                    if(reml){
                      Eigen::MatrixXd selectedRHS = Eigen::MatrixXd::Zero(
                        static_cast<Eigen::Index>(nEffects),
                        static_cast<Eigen::Index>(blockWidth));
                      for(arma::uword j = 0; j < blockWidth; ++j){
                        selectedRHS(static_cast<Eigen::Index>(colStart + j),
                                    static_cast<Eigen::Index>(j)) = 1.0;
                      }
                      fallbackBlockSolution = solveCMatrix(
                        selectedRHS,
                        "selected-block random-effect score-trace fallback");
                    }else{
                      Eigen::MatrixXd selectedRHS = Eigen::MatrixXd::Zero(
                        static_cast<Eigen::Index>(Nu),
                        static_cast<Eigen::Index>(blockWidth));
                      for(arma::uword j = 0; j < blockWidth; ++j){
                        selectedRHS(static_cast<Eigen::Index>(colStart + j - nX),
                                    static_cast<Eigen::Index>(j)) = 1.0;
                      }
                      if(solverName == "cholmod"){
                        fallbackBlockSolution = solveDMatrixCholmod(
                          selectedRHS,
                          "selected-block random-effect score-trace fallback");
                      }else{
                        fallbackBlockSolution = Dfactor.solve(selectedRHS);
                        if(Dfactor.info() != Eigen::Success){
                          Rcpp::stop("Sparse LDLT trace fallback solve failed in selected-block random-effect score-trace fallback (reml=FALSE).");
                        }
                      }
                    }
                    fallbackTimer.stop();
                  }
                  fallbackBlockAvailable = true;
                }
                arma::mat inverseBlock(blockHeight, blockWidth);
                const arma::uword rowOffset = reml ? 0 : nX;
                for(arma::uword rr = 0; rr < blockHeight; ++rr){
                  for(arma::uword cc = 0; cc < blockWidth; ++cc){
                    inverseBlock(rr,cc) = fallbackBlockSolution(
                      static_cast<Eigen::Index>(rowStart + rr - rowOffset),
                      static_cast<Eigen::Index>(cc));
                  }
                }
                // trace(Ai*inverseBlock) via a full dense product is O(q^3)
                // (q=blockHeight), dominant when Ai is a dense relationship
                // matrix. Ai is symmetric, so trace(Ai*B) == accu(Ai % B);
                // accumulate directly over Ai's stored entries instead,
                // which is O(nnz(Ai)) (== O(q^2) even when Ai is dense, but
                // avoids the extra BLAS-3 pass entirely).
                double tsum = 0.0;
                for(arma::sp_mat::const_iterator ait2 = Ai(iR).begin();
                    ait2 != Ai(iR).end(); ++ait2){
                  tsum += (*ait2) * inverseBlock(ait2.row(), ait2.col());
                }
                trAiCuu = tsum;
              }
            }

            traces(iRow, iCol) = trAiCuu;
          }
        }
        }
        if(std::find(traceBatchFailed.begin(), traceBatchFailed.end(), 1) != traceBatchFailed.end()){
          Rcpp::stop("Non-finite batched exact trace solution.");
        }
        }

        traces =
          arma::symmatu(traces);
        if(saddleActive){
          saddleInversePriorTrace += arma::accu(arma::mat(lambda(iR)) % traces);
        }

        {

          arma::vec structuredScore(
            covPar(iR).n_elem,
            arma::fill::zeros
          );

          const DescriptorPrecisionState & precisionState =
            randomPrecisionState[static_cast<std::size_t>(iR)];

          const arma::mat scoreQuadratic =
            randomQuadraticBase[static_cast<std::size_t>(iR)]
            +
            traces;

          for(arma::uword k = 0; k < covPar(iR).n_elem; ++k){
            double precisionContraction = 0.0;
            if(precisionState.factorWise){
              std::vector<arma::mat> derivativeFactors;
              double derivativeScale = 1.0;
              factorEngine.descriptorPrecisionD1Factors(
                precisionState,
                k,
                derivativeFactors,
                derivativeScale
              );
              precisionContraction =
                factorEngine.contractKronecker(
                  scoreQuadratic,
                  derivativeFactors,
                  derivativeScale
                );
              if(profileCommonScale){
                commonScaleQuadratic(nVcStart(iR) - 1 + k) =
                  -factorEngine.contractKronecker(randomQuadraticBase[iR],
                    derivativeFactors, derivativeScale);
              }
              if(solverName == "pcg"){
                for(int probe = 0; probe < pcgTraceProbes; ++probe){
                  pcgScoreProbeValues(static_cast<Eigen::Index>(nVcStart(iR) - 1 + k), probe) =
                    factorEngine.contractKronecker(pcgRandomProbeTraces[probe],
                      derivativeFactors, derivativeScale);
                }
              }
            }else{
              precisionContraction =
                arma::accu(
                  scoreQuadratic
                  %
                  precisionState.precisionD1[static_cast<std::size_t>(k)]
                );
              if(profileCommonScale){
                commonScaleQuadratic(nVcStart(iR) - 1 + k) =
                  -arma::accu(randomQuadraticBase[iR] %
                    precisionState.precisionD1[k]);
              }
              if(solverName == "pcg"){
                for(int probe = 0; probe < pcgTraceProbes; ++probe){
                  pcgScoreProbeValues(static_cast<Eigen::Index>(nVcStart(iR) - 1 + k), probe) =
                    arma::accu(pcgRandomProbeTraces[probe] %
                      precisionState.precisionD1[static_cast<std::size_t>(k)]);
                }
              }
            }

            structuredScore(k) =
              arma::as_scalar(nUsTotal(iR))
              *
              precisionState.logDetD1[static_cast<std::size_t>(k)]
              +
              precisionContraction;
          }

          dLu =
            arma::join_cols(
              dLu,
              structuredScore
            );

          // A true EM update in nonlinear AR1 coordinates is not the same
          // diagonal closed form used for variance-cell parameters.  Use the
          // positive diagonal of the exact AI block as the stabilizing
          // information contribution.
          if(!emInformationOnly){
          const arma::uword g0 =
            static_cast<arma::uword>(
              nVcStart(iR) - 1
            );

          const arma::uword g1 =
            static_cast<arma::uword>(
              nVcEnd(iR) - 1
            );

          arma::vec aiDiag =
            arma::diagvec(
              avInf.submat(
                g0,
                g0,
                g1,
                g1
              )
            );

          for(arma::uword k = 0; k < aiDiag.n_elem; ++k){
            if(!std::isfinite(aiDiag(k)) || aiDiag(k) <= tolParInv){
              aiDiag(k) = 1.0;
            }
          }

          emInfList(iR) =
            arma::diagmat(
              aiDiag
            );
          }
        }
      }
    }

    // ============================================================
    // Residual score traces in generic working-parameter coordinates
    // ============================================================
    arma::vec dLe(
      nResidualPar,
      arma::fill::zeros
    );

    // RiWsp is constant across the iP loop below; converting it to
    // Armadillo once here avoids repeating a full triplet-based
    // Eigen->Armadillo conversion for every residual parameter.
    const arma::sp_mat RiWspArmaShared =
      RiWisSparse && !compressedData ? eigenSparseToArmaGlobal(RiWsp) : arma::sp_mat();

    // With a diagonal effective R^{-1} and diagonal dR/dphi, the C-trace
    // operand is W' diag(s .* r^2) W; a row-major W touches only rows with s != 0.
    const bool diagonalRiW =
      effectiveRiDiag.n_elem == static_cast<arma::uword>(nR);
    Eigen::SparseMatrix<double, Eigen::RowMajor> WRowMajor;
    if(diagonalRiW && !compressedData){
      WRowMajor = W;
      WRowMajor.makeCompressed();
    }
    auto isDiagonalSparse = [](const arma::sp_mat & S) -> bool {
      for(arma::sp_mat::const_iterator it = S.begin(); it != S.end(); ++it){
        if(it.row() != it.col()){ return false; }
      }
      return true;
    };

    Eigen::MatrixXd pcgRiWZ;
    Eigen::MatrixXd pcgRiWX;
    if(solverName == "pcg" && RiWisSparse){
      if(compressedData){
        pcgRiWZ.noalias() = sufficient.design * pcgTraceZ;
        pcgRiWX.noalias() = sufficient.design * pcgTraceX;
      }else{
        pcgRiWZ.noalias() = RiWsp * pcgTraceZ;
        pcgRiWX.noalias() = RiWsp * pcgTraceX;
      }
    }

    Eigen::VectorXd saddleResidualTraces;
    if(saddleActive){
      const Eigen::VectorXd shapeTraces = sommer::exactInverseTraces(
        nR, nResidualPar - 1, sommer::WorkspacePolicy::panelColumns(),
        [&](const Eigen::MatrixXd & rhs) -> Eigen::MatrixXd {
          Rcpp::checkUserInterrupt();
          if(reml){ return saddleJoint.solveResidual(rhs, "REML residual score traces"); }
          if(Nu > 0){ return saddleRandom.solveResidual(rhs, "ML residual score traces"); }
          const arma::mat values(rhs.data(), rhs.rows(), rhs.cols());
          const arma::mat transformed = useH ? arma::mat(HsInvT * values) : values;
          arma::mat applied = applyRiDense(transformed);
          if(useH){ applied = HsInvT.t() * applied; }
          return Eigen::Map<const Eigen::MatrixXd>(applied.memptr(), nR, rhs.cols());
        },
        [&](Eigen::Index shapeParameter, const Eigen::MatrixXd & values) -> Eigen::MatrixXd {
          const Eigen::Index parameter = shapeParameter + 1;
          const arma::mat input(values.data(), values.rows(), values.cols());
          arma::mat applied(nR, values.cols(), arma::fill::zeros);
          if(residualDerivativesBlockLocal){
            for(int block = 0; block < residualKronNBlocks; ++block){
              const arma::uvec & rows = residualKronRows[static_cast<std::size_t>(block)];
              if(residualPrecisionState.factorWise){
                applied.rows(rows) = factorEngine.applyDescriptorCovarianceD1Right(
                  input.rows(rows).t(), residualPrecisionState, parameter).t();
              }else{
                applied.rows(rows) = residualLocalD1(parameter) * input.rows(rows);
              }
            }
          }else if(useH){
            applied = Hs.t() * (residualDerivativeBasis(parameter) * (Hs * input));
          }else{
            applied = residualDerivativeBasis(parameter) * input;
          }
          return Eigen::Map<const Eigen::MatrixXd>(applied.memptr(), nR, values.cols());
        });
      saddleResidualTraces.resize(nResidualPar);
      // tr(P V) = n-p (REML), or tr(V^-1 V) = n (ML).
      saddleResidualTraces(0) = nR - (reml ? nX : 0) - Nu + saddleInversePriorTrace;
      saddleResidualTraces.tail(nResidualPar - 1) = shapeTraces;
    }

    for(arma::uword iP = 0; iP < nResidualPar; ++iP){

      const arma::sp_mat & Sprov =
        residualDerivativeBasis(iP);

      arma::sp_mat traceBasisDeriv;

      if(Rgrid || compressedData){
        // traceBasisDeriv is unused on the grid path.
      }else if(useH){
        traceBasisDeriv =
          Hs.t()
          *
          Sprov
          *
          Hs;
      }else{
        traceBasisDeriv =
          Sprov;
      }

      // First trace: tr((dR/dphi_i) * effective R^{-1}).
      double traceSRi = 0.0;
      bool usedCTraceFallback = false;
      double traceCorrection = 0.0;
      Eigen::VectorXd pcgResidualSamples;
      if(solverName == "pcg"){ pcgResidualSamples = Eigen::VectorXd::Zero(pcgTraceProbes); }

      if(compressedData){
        for(std::size_t group = 0; group < sufficient.groups.size(); ++group){
          const auto & statistics = sufficient.groups[group];
          const int category = statistics.residualClass;
          const double variance = compressedVariance(category);
          const double derivative = compressedDerivative(category, iP);
          traceSRi += statistics.count * derivative / variance;
          const double weight = statistics.weight * derivative / (variance * variance);
          if(solverName == "pcg"){
            const Eigen::ArrayXd samples = weight *
              pcgRiWZ.row(group).transpose().array() * pcgRiWX.row(group).transpose().array();
            traceCorrection += samples.sum();
            pcgResidualSamples.array() += samples.head(pcgTraceProbes);
          }
        }
        if(solverName == "pcg"){
          traceCorrection /= static_cast<double>(pcgTraceProbes);
          pcgScoreProbeValues.row(nVcStart(residualStruct) - 1 + iP) = -pcgResidualSamples.transpose();
        }else{
          const Eigen::VectorXd scales = compressedDerivative.col(iP).array() /
            compressedVariance.array().square();
          const EigenSpMat derivativeMatrix = sufficient.derivativeCross(scales);
          const arma::sp_mat Btrace = eigenSparseToArmaGlobal(derivativeMatrix);
          if(reml){
            traceCorrection = sparseTraceInverseTimes(
              Btrace, Cfactor, CselectedTopology,
              "compressed residual trace tr(C^{-1} W'Ri(dR/dphi)RiW)",
              usedCTraceFallback,
              solverName == "cholmod"
                ? std::function<Eigen::MatrixXd(const Eigen::Ref<const Eigen::MatrixXd> &, const std::string &)>(solveCMatrix)
                : nullptr,
              solverName == "cholmod");
          }else if(Nu > 0){
            traceCorrection = sparseTraceInverseTimes(
              arma::sp_mat(Btrace.submat(nX, nX, nEffects - 1, nEffects - 1)),
              Dfactor, DselectedTopology,
              "compressed residual trace tr(D^{-1} Z'Ri(dR/dphi)RiZ)",
              usedCTraceFallback,
              solverName == "cholmod"
                ? std::function<Eigen::MatrixXd(const Eigen::Ref<const Eigen::MatrixXd> &, const std::string &)>(solveDMatrixCholmod)
                : nullptr);
          }
        }
      }else if(saddleActive){
        traceSRi = saddleResidualTraces(iP);
      }else{
      if(Rdiag){

        for(arma::sp_mat::const_iterator it = traceBasisDeriv.begin();
            it != traceBasisDeriv.end();
            ++it){

          if(it.row() == it.col()){
            traceSRi +=
              (*it)
              *
              RdiagInv(
                static_cast<arma::uword>(
                  it.row()
                )
              );
          }
        }

      }else if(Rgrid){

        traceSRi = gridTraceDR(iP);

      }else if(Rkron){

        if(!useH && residualPrecisionState.factorWise){
          traceSRi =
            static_cast<double>(residualKronNBlocks)
            *
            residualPrecisionState.logDetD1[
              static_cast<std::size_t>(iP)
            ];
        }else if(residualDerivativesBlockLocal){
          // The record-space basis is not built in block-local mode; all
          // complete blocks share the same q x q derivative and precision.
          traceSRi =
            static_cast<double>(residualKronNBlocks)
            *
            arma::accu(
              residualLocalD1(iP)
              %
              RkronInv
            );
        }else{

        // Independent per-block extraction/accumulation: parallelize over
        // blocks, and use trace(A*B) = accu(A % B) (valid since RkronInv is
        // symmetric) instead of forming the full q x q product.
        #ifdef _OPENMP
          #pragma omp parallel for if(residualKronNBlocks > 1) schedule(static) reduction(+:traceSRi)
        #endif
        for(int b = 0; b < residualKronNBlocks; ++b){

          const arma::uvec & rows =
            residualKronRows[
              static_cast<std::size_t>(b)
            ];

          arma::mat local(
            residualKronBlockSize,
            residualKronBlockSize,
            arma::fill::zeros
          );

          for(arma::uword aa = 0; aa < rows.n_elem; ++aa){
            for(arma::uword bb = 0; bb < rows.n_elem; ++bb){
              local(aa,bb) =
                traceBasisDeriv(
                  rows(aa),
                  rows(bb)
                );
            }
          }

          traceSRi +=
            arma::accu(
              local
              %
              RkronInv
            );
        }

        }

      }else{

        bool usedRTraceFallback = false;

        traceSRi =
          sparseTraceInverseTimes(
            traceBasisDeriv,
            Rfactor,
            Rselected,
            "residual trace tr((dR/dphi) R^{-1})",
            usedRTraceFallback,
            solverName == "cholmod" ? std::function<Eigen::MatrixXd(const Eigen::Ref<const Eigen::MatrixXd> &, const std::string &)>(solveRMatrixCholmod) : nullptr
          );

        if(verbose && usedRTraceFallback && iIter == 0){
          Rcpp::Rcout
            << "Residual R trace required exact sparse-solve fallback for parameter "
            << iP + 1
            << arma::endl;
        }
      }

      if(solverName == "pcg" && Rdiag && RiWisSparse && isDiagonalSparse(Sprov)){
        for(arma::sp_mat::const_iterator it = Sprov.begin();
            it != Sprov.end(); ++it){
          if(it.row() == it.col()){
            traceCorrection +=
              (*it)
              *
              (
                pcgRiWZ.row(static_cast<Eigen::Index>(it.row())).array()
                *
                pcgRiWX.row(static_cast<Eigen::Index>(it.row())).array()
              ).sum();
            for(int probe = 0; probe < pcgTraceProbes; ++probe){
              pcgResidualSamples(probe) += (*it) *
                pcgRiWZ(it.row(), probe) * pcgRiWX(it.row(), probe);
            }
          }
        }
        traceCorrection /= static_cast<double>(pcgTraceProbes);
      }else if(solverName == "pcg" && RiWisSparse){
        const arma::mat probeZ(pcgRiWZ.data(), nR, pcgTraceZ.cols(), false, true);
        const arma::mat probeX(pcgRiWX.data(), nR, pcgTraceX.cols(), false, true);
        if(residualDerivativesBlockLocal){
          for(int block = 0; block < residualKronNBlocks; ++block){
            const arma::uvec & rows = residualKronRows[static_cast<std::size_t>(block)];
            const arma::mat localX = probeX.rows(rows);
            arma::mat derivativeApplied;
            if(Rgrid){
              derivativeApplied = gridApplyDR(block, iP, localX);
            }else if(residualPrecisionState.factorWise){
              derivativeApplied = factorEngine.applyDescriptorCovarianceD1Right(
                localX.t(), residualPrecisionState, iP).t();
            }else{
              derivativeApplied = residualLocalD1(iP) * localX;
            }
            const arma::mat localZ = probeZ.rows(rows);
            traceCorrection += arma::accu(localZ % derivativeApplied);
            for(int probe = 0; probe < pcgTraceProbes; ++probe){
              pcgResidualSamples(probe) += arma::dot(localZ.col(probe),
                derivativeApplied.col(probe));
            }
          }
        }else{
          const arma::mat derivativeApplied = arma::mat(Sprov * probeX);
          traceCorrection = arma::accu(probeZ % derivativeApplied);
          for(int probe = 0; probe < pcgTraceProbes; ++probe){
            pcgResidualSamples(probe) = arma::dot(probeZ.col(probe), derivativeApplied.col(probe));
          }
        }
        traceCorrection /= static_cast<double>(pcgTraceProbes);
        ++pcgResidualTraceContractions;
      }else{
        arma::sp_mat Btrace;

        if(residualDerivativesBlockLocal){
          std::vector<EigenTriplet> derivativeTrips;
          for(int block = 0; block < residualKronNBlocks; ++block){
            const arma::uvec & columns =
              residualActiveColumns[static_cast<std::size_t>(block)];
            if(columns.n_elem == 0){
              continue;
            }
            const arma::mat & blockRiW =
              residualKronBlockRiW[static_cast<std::size_t>(block)];
            arma::mat localBasis;
            if(Rgrid){
              localBasis = blockRiW.t() * gridApplyDR(block, iP, blockRiW);
            }else if(residualPrecisionState.factorWise){
              localBasis = factorEngine.applyDescriptorCovarianceD1Right(
                blockRiW.t(), residualPrecisionState, iP
              ) * blockRiW;
            }else{
              localBasis = blockRiW.t() * residualLocalD1(iP) * blockRiW;
            }

            for(arma::uword localCol = 0; localCol < columns.n_elem; ++localCol){
              for(arma::uword localRow = 0; localRow < columns.n_elem; ++localRow){
                const double value = localBasis(localRow, localCol);
                if(value != 0.0){
                  derivativeTrips.emplace_back(
                    static_cast<int>(columns(localRow)),
                    static_cast<int>(columns(localCol)),
                    value
                  );
                }
              }
            }
          }

          EigenSpMat localBtrace(nEffects, nEffects);
          localBtrace.setFromTriplets(
            derivativeTrips.begin(),
            derivativeTrips.end()
          );
          localBtrace.makeCompressed();
          Btrace = eigenSparseToArmaGlobal(localBtrace);

        }else if(diagonalRiW && isDiagonalSparse(Sprov)){
          std::vector<arma::uword> locRows;
          std::vector<arma::uword> locCols;
          std::vector<double> locValues;
          for(arma::sp_mat::const_iterator it = Sprov.begin(); it != Sprov.end(); ++it){
            const arma::uword obs = it.row();
            const double ri = effectiveRiDiag(obs);
            const double weight = (*it) * ri * ri;
            if(weight == 0.0){ continue; }
            const int row = static_cast<int>(obs);
            for(Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator a(WRowMajor, row); a; ++a){
              for(Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator c(WRowMajor, row); c; ++c){
                locRows.push_back(static_cast<arma::uword>(a.col()));
                locCols.push_back(static_cast<arma::uword>(c.col()));
                locValues.push_back(weight * a.value() * c.value());
              }
            }
          }
          arma::umat locations(2, locValues.size());
          for(std::size_t q = 0; q < locValues.size(); ++q){
            locations(0, q) = locRows[q];
            locations(1, q) = locCols[q];
          }
          Btrace = arma::sp_mat(true, locations, arma::vec(locValues),
                                static_cast<arma::uword>(nEffects),
                                static_cast<arma::uword>(nEffects));
        }else{
          Btrace =
            RiWspArmaShared.t()
            *
            Sprov
            *
            RiWspArmaShared;
        }

        traceCorrection =
          solverName == "pcg"
          ? pcgTraceCInverseTimesSparse(Btrace)
          : (reml
              ? sparseTraceInverseTimes(
                  Btrace,
                  Cfactor,
                  CselectedTopology,
                  "residual trace tr(C^{-1} W'Ri(dR/dphi)RiW)",
                  usedCTraceFallback,
                  solverName == "cholmod" ? std::function<Eigen::MatrixXd(const Eigen::Ref<const Eigen::MatrixXd> &, const std::string &)>(solveCMatrix) : nullptr,
                  solverName == "cholmod"
                )
              : (Nu > 0
                  ? sparseTraceInverseTimes(
                      arma::sp_mat(Btrace.submat(
                        static_cast<arma::uword>(nX),
                        static_cast<arma::uword>(nX),
                        static_cast<arma::uword>(nEffects - 1),
                        static_cast<arma::uword>(nEffects - 1)
                      )),
                      Dfactor,
                      DselectedTopology,
                      "residual trace tr(D^{-1} Z'Ri(dR/dphi)RiZ)",
                      usedCTraceFallback,
                      solverName == "cholmod" ? std::function<Eigen::MatrixXd(const Eigen::Ref<const Eigen::MatrixXd> &, const std::string &)>(solveDMatrixCholmod) : nullptr
                    )
                  : 0.0)
            );
      }
      if(solverName == "pcg"){
        if(!RiWisSparse){ pcgResidualSamples = pcgLastSparseTraceSamples; }
        pcgScoreProbeValues.row(nVcStart(residualStruct) - 1 + iP) = -pcgResidualSamples.transpose();
      }

      if(verbose && solverName == "ldlt" && usedCTraceFallback && iIter == 0){
        Rcpp::Rcout
          << "Residual C trace required exact sparse-solve fallback for parameter "
          << iP + 1
          << arma::endl;
      }
      }

      double residualQuadratic = 0.0;
      if(compressedData){
        for(std::size_t group = 0; group < sufficient.groups.size(); ++group){
          const int category = sufficient.groups[group].residualClass;
          const double variance = compressedVariance(category);
          residualQuadratic += compressedSquares(group) * compressedDerivative(category, iP) /
            (variance * variance);
        }
      }else if(residualDerivativesBlockLocal){
        for(int block = 0; block < residualKronNBlocks; ++block){
          const arma::uvec & rows =
            residualKronRows[static_cast<std::size_t>(block)];
          const arma::vec localRie = Rie.elem(rows);
          arma::vec derivativeApplied;
          if(Rgrid){
            derivativeApplied = gridApplyDR(block, iP, localRie);
          }else if(residualPrecisionState.factorWise){
            derivativeApplied = factorEngine.applyDescriptorCovarianceD1Right(
              localRie.t(), residualPrecisionState, iP
            ).t();
          }else{
            derivativeApplied = residualLocalD1(iP) * localRie;
          }
          residualQuadratic +=
            arma::dot(
              localRie,
              derivativeApplied
            );
        }
      }else{
        residualQuadratic =
          arma::dot(
            Rie,
            arma::vec(
              Sprov
              *
              Rie
            )
          );
      }

      dLe(iP) =
        traceSRi
        -
        traceCorrection
        -
        residualQuadratic;
      if(profileCommonScale){
        commonScaleQuadratic(residualVcStart + iP) = residualQuadratic;
      }
    }

    if(!emInformationOnly){
      const arma::uword g0 =
        static_cast<arma::uword>(
          nVcStart(residualStruct) - 1
        );

      const arma::uword g1 =
        static_cast<arma::uword>(
          nVcEnd(residualStruct) - 1
        );

      arma::vec aiDiag =
        arma::diagvec(
          avInf.submat(
            g0,
            g0,
            g1,
            g1
          )
        );

      for(arma::uword k = 0; k < aiDiag.n_elem; ++k){
        if(!std::isfinite(aiDiag(k)) || aiDiag(k) <= tolParInv){
          aiDiag(k) = 1.0;
        }
      }

      emInfList(residualStruct) =
        arma::diagmat(
          aiDiag
        );
    }

    phaseProfile.enter(sommer::FitPhase::Update);
    dLu = join_cols( dLu, dLe );// join the random and residual first derivatives in a single vector
    if(profileCommonScale){
      const arma::vec shapeQuadratic = aiTangent.t() * commonScaleQuadratic;
      avInf = (avInf - shapeQuadratic * shapeQuadratic.t() /
        commonScaleQuadraticTotal) / commonScale;
      avInf = 0.5 * (avInf + avInf.t());
      dLu += commonScaleQuadratic * (1.0 - 1.0 / commonScale);
      if(!avInf.is_finite() || !dLu.is_finite()){
        Rcpp::stop("Non-finite profiled common-scale score or information.");
      }
    }
    if(solverName == "pcg"){
      pcgScoreMaxSE = 0.0;
      for(int parameter = 0; parameter < nVcTotal; ++parameter){
        const double mean = pcgScoreProbeValues.row(parameter).mean();
        pcgScoreProbeValues.row(parameter).array() += dLu(parameter) - mean;
        pcgScoreProbeValues.row(parameter) *= -0.5;
      }
      if(modelFreeCount > 0 && pcgTraceProbes < 2){
        pcgScoreMaxSE = std::numeric_limits<double>::quiet_NaN();
      }else{
        Eigen::MatrixXd freeScores(modelFreeCount, pcgTraceProbes);
        for(arma::uword row = 0; row < modelFreeCount; ++row){
          freeScores.row(row) = pcgScoreProbeValues.row(hybridFree(row));
        }
        if(vcMap.active){
          const arma::mat map = vcMap.T.rows(hybridFree);
          const Eigen::Map<const Eigen::MatrixXd> mapping(map.memptr(), map.n_rows, map.n_cols);
          freeScores = (mapping.transpose() * freeScores).eval();
        }
        for(Eigen::Index parameter = 0; parameter < freeScores.rows(); ++parameter){
          const Eigen::ArrayXd centered = freeScores.row(parameter).transpose().array() -
            freeScores.row(parameter).mean();
          pcgScoreMaxSE = std::max(pcgScoreMaxSE,
            std::sqrt(centered.square().sum() / (static_cast<double>(pcgTraceProbes - 1) * pcgTraceProbes)));
        }
      }
      const bool logDetReady = std::isfinite(pcgLogDetProbeStandardError) &&
        pcgLogDetProbeStandardError <= pcgLogDetTol &&
        std::isfinite(pcgLogDetDepthDifference) && pcgLogDetDepthDifference <= pcgLogDetTol;
      const bool scoreReady = std::isfinite(pcgScoreMaxSE) && pcgScoreMaxSE <= pcgScoreTol;
      if(pcgAccuracy == "adaptive" && iIter == 0 &&
         ((!logDetReady && (pcgTraceProbes < pcgMaxTraceProbes || pcgLanczosSteps < pcgMaxLanczosSteps)) ||
          (!scoreReady && pcgTraceProbes < pcgMaxTraceProbes))){
        pcgTraceProbes = std::min<long long>(pcgMaxTraceProbes, 2LL * pcgTraceProbes);
        pcgLanczosSteps = std::min<long long>(pcgMaxLanczosSteps, 2LL * pcgLanczosSteps);
        ++pcgAccuracyRefinements;
        // Calibrate once, then freeze common probes/depth before accepting any step.
        haveAcceptedLikelihood = false;
        pcgEvaluationAccepted = false;
        pcgCandidateSolutions.clear();
        --iIter;
        continue;
      }
    }
    
    if(krylovPcg && pcgKrylov.candidateCount() > 0){
      const auto selectionStarted = std::chrono::steady_clock::now();
      pcgPendingKrylov = pcgKrylov.selectCandidates(
        [&](const Eigen::VectorXd & value) -> Eigen::VectorXd {
          Rcpp::checkUserInterrupt();
          return applyCurrentPCG(Eigen::MatrixXd(value)).col(0);
        }, pcgPreconditionerDiagonal, krylovRank);
      pcgPendingKrylovReady = true;
      pcgSetupSeconds += std::chrono::duration<double>(
        std::chrono::steady_clock::now() - selectionStarted).count();
    }
    for(int i = 0; !emInformationOnly && i < nRRe; ++i){
      emInf.submat(nVcStart(i)-1, nVcStart(i)-1, nVcEnd(i)-1, nVcEnd(i)-1 ) = emInfList(i);
    }
    
    // ###########################
    // # 6) update the variance paramters using the Newton method
    // # PAPER FORMULA (Lee and Van der Werf, 2006)
    // # theta.n+1 = theta.n + (AInfi * dL/ds2)
    // ###########################
    
    arma::vec thetaUnlisted, thetaCUnlisted;
    arma::vec parLowerUnlisted, parUpperUnlisted, parScaleUnlisted;

    for(int i = 0; i < nRRe; ++i){

      thetaUnlisted =
        arma::join_cols(
          thetaUnlisted,
          covPar(i)
        );

      thetaCUnlisted =
        arma::join_cols(
          thetaCUnlisted,
          covConstraint(i)
        );

      parLowerUnlisted =
        arma::join_cols(
          parLowerUnlisted,
          covLower(i)
        );

      parUpperUnlisted =
        arma::join_cols(
          parUpperUnlisted,
          covUpper(i)
        );

      parScaleUnlisted =
        arma::join_cols(
          parScaleUnlisted,
          covScale(i)
        );
    }
    // Rcpp::Rcout << "thetaCUnlisted" << thetaCUnlisted << arma::endl;
    if(useHybridAI){
      if(!dLu.is_finite()){
        Rcpp::stop("Non-finite variance-component score in hybrid AI.");
      }
      if(!exactAIThisIteration && !hybridState.update(thetaUnlisted, dLu, hybridFree)){
        ++hybridAIResets;
        sensitivityRHSColumnsSkipped -= static_cast<std::size_t>(nVcTotal);
        computeAverageInformation(true, false);
        exactAIThisIteration = true;
        completeAIThisIteration = true;
      }
      if(exactAIThisIteration){
        hybridState.initialize(avInf, thetaUnlisted, dLu, hybridFree);
        lastExactAIIteration = iIter;
      }else{
        avInf = hybridState.metric;
        ++hybridAIUpdates;
      }
      // Between refreshes, stabilization follows the current BFGS metric,
      // not stale exact-AI diagonals. Final uncertainty never uses this metric.
      arma::vec stabilizer = avInf.diag();
      for(arma::uword k = 0; k < stabilizer.n_elem; ++k){
        if(!std::isfinite(stabilizer(k)) || stabilizer(k) <= tolParInv){
          stabilizer(k) = 1.0;
        }
      }
      emInf = arma::diagmat(stabilizer);
      phaseProfile.enter(sommer::FitPhase::Update);
    }
    finalAIExact = exactAIThisIteration && completeAIThisIteration;
    if(emInformationOnly){
      InfMat = avInf;
      monitor.col(iIter) = thetaUnlisted;
      normMonitor.zeros();
      normMonitor(1, iIter) = arma::norm(dLu, 1);
      dLuOut = dLu;
      invokeConditionalCallback();
      break;
    }
    // create the 'weight' EM information matrix (TO BE USED LATER WITHIN THE OPTIMIZATION)
    arma::vec v2(nVcTotal, arma::fill::ones) ;
    arma::mat weightEmInfMat = diagmat(v2) * arma::as_scalar(weightEmInf(iIter));
    arma::mat weightAiInfMat = diagmat(v2) * (1 - arma::as_scalar(weightEmInf(iIter)));
    // Joint information matrix and update
    //                  AVERAGE INFORMATION                         +       EXPECTATION MAXIMIZATION
    InfMat = (weightAiInfMat * avInf) + (weightEmInfMat * emInf);

    arma::vec weightedScore =
      arma::as_scalar(weightInf(iIter))
      *
      dLu;

    bool informationSolveOK =
      solveInformationSystem(
        delta,
        InfMat,
        weightedScore,
        "main variance-component update"
      );

    if(!informationSolveOK){
      Rcpp::stop(
        "Unable to solve the variance-component information system."
      );
    }

    if(useAitken && iIter % 3 == 2 && lineSearchHalvings == 0 &&
       previousAcceptedParameters.n_elem == thetaUnlisted.n_elem){
      const arma::vec previousStep = thetaUnlisted - previousAcceptedParameters;
      const arma::vec proposedStep = -delta;
      const arma::vec difference = proposedStep - previousStep;
      const double denominator = arma::dot(difference, difference);
      if(denominator > std::numeric_limits<double>::epsilon() &&
         arma::dot(previousStep, proposedStep) > 0.0){
        const double multiplier = -arma::dot(previousStep, difference) / denominator;
        if(std::isfinite(multiplier) && multiplier > 1.0){
          delta *= std::min(2.0, multiplier);
          ++acceleratedProposals;
        }
      }
    }
    previousAcceptedParameters = thetaUnlisted;

    arma::vec expectedNewTheta =
      thetaUnlisted
      -
      delta;
    
    phaseProfile.enter(sommer::FitPhase::Constraints);
    // #######################
    // # 7) APPLY CONSTRAINTS to VC
    // # suggestions from Madsen and Jensen (1997) and Gilmour (2019)
    // #######################
    //
    // Robust boundary / positive-definiteness handling:
    //
    //  * Use one numerical floor consistently for positive variance
    //    components and for covariance-matrix PD checks.
    //  * Touching a boundary does NOT immediately remove a parameter
    //    from the information-system update.
    //  * A positive parameter is frozen at its current boundary value
    //    only after 3 boundary hits.
    //  * PD repairs are local to the failing covariance structure.
    //  * If every proposed repair fails, reject only that structure's
    //    update and retain its previous covariance matrix rather than
    //    aborting the whole REML fit.
    // #######################

    const double vcFloor =
      std::max(
        1.0e-8,
        tolParInv
      );

    // A covariance matrix only needs to be numerically inside the open
    // positive-definite cone.  This threshold is deliberately below the
    // enforced scalar variance floor.
    const double pdCheckFloor =
      std::max(
        1.0e-12,
        0.1 * vcFloor
      );

    // Helper to impose descriptor-defined scalar search-space bounds.
    // Legacy structures reproduce the previous behaviour; structured
    // models such as AR1 can supply their own bounds (e.g. |rho| < 1).
    auto applyVarianceBounds =
      [&](arma::vec & candidate,
          const arma::uvec & indices,
          const bool markBoundary){

        for(arma::uword k = 0; k < indices.n_elem; ++k){

          const arma::uword g =
            indices(k);

          if(thetaCUnlisted(g) == 3){
            continue;
          }

          double lower =
            parLowerUnlisted(g);

          double upper =
            parUpperUnlisted(g);

          if(thetaCUnlisted(g) == 1){
            lower =
              std::max(
                lower,
                vcFloor
              );
          }

          if(
              std::isfinite(lower)
              &&
              candidate(g) < lower
          ){
            candidate(g) =
              lower;

            if(markBoundary){
              toBoundary(iIter,g) =
                1;
            }
          }

          if(
              std::isfinite(upper)
              &&
              candidate(g) > upper
          ){
            candidate(g) =
              upper;

            if(markBoundary){
              toBoundary(iIter,g) =
                1;
            }
          }
        }
      };

    // A) Apply ordinary scalar bounds to the main proposal.
    arma::uvec allParameterIndices =
      arma::regspace<arma::uvec>(
        static_cast<arma::uword>(0),
        static_cast<arma::uword>(1),
        static_cast<arma::uword>(nVcTotal - 1)
      );

    applyVarianceBounds(
      expectedNewTheta,
      allParameterIndices,
      true
    );

    // Recompute cumulative boundary-hit counts once per iteration.
    // The previous implementation accidentally treated a single boundary
    // hit as constrained even though the comments intended a 3-hit rule.
    for(int j = 0; j < toBoundary.n_cols; ++j){
      sumToBoundary(j) =
        arma::accu(
          toBoundary.col(j)
        );
    }

    // Parameters explicitly fixed by the model specification.
    arma::uvec modelFixed =
      arma::find(
        thetaCUnlisted == 3
      );

    // Positive variance parameters are frozen only after 3 boundary hits.
    // Keep this state separate from thetaC so dynamically boundary-frozen
    // parameters are not accidentally interpreted through thetaF.
    arma::uvec boundaryForced =
      arma::find(
        (thetaCUnlisted == 1)
        %
        (sumToBoundary >= 3)
      );

    arma::uvec constrained;

    if(modelFixed.n_elem > 0 && boundaryForced.n_elem > 0){
      constrained =
        arma::unique(
          arma::join_cols(
            modelFixed,
            boundaryForced
          )
        );
    }else if(modelFixed.n_elem > 0){
      constrained =
        modelFixed;
    }else{
      constrained =
        boundaryForced;
    }

    // Complement of constrained parameters.
    arma::uvec constrainedMask(
      nVcTotal,
      arma::fill::zeros
    );

    if(constrained.n_elem > 0){
      constrainedMask(constrained).ones();
    }

    arma::uvec unconstrained =
      arma::find(
        constrainedMask == 0
      );

    // Model-defined fixed working parameters simply retain their current
    // accepted value.  thetaF/addScaleParam are no longer part of the
    // Henderson covariance interface.
    if(modelFixed.n_elem > 0){
      expectedNewTheta(modelFixed) =
        thetaUnlisted(modelFixed);
    }

    // Boundary-forced parameters stay exactly at their current accepted
    // value.  They are no longer allowed to move after the third hit.
    if(boundaryForced.n_elem > 0){
      expectedNewTheta(boundaryForced) =
        thetaUnlisted(boundaryForced);
    }

    // B) If there are constrained parameters, re-solve only the genuinely
    // free variance-component block.  Merely touching a boundary once or
    // twice no longer removes a parameter from this solve.
    arma::mat InfMat_uu;
    arma::vec dLu_uu, delta_uu;

    if(constrained.n_elem > 0){

      delta(constrained).zeros();

      if(unconstrained.n_elem > 0){

        InfMat_uu =
          InfMat(
            unconstrained,
            unconstrained
          );

        dLu_uu =
          dLu(unconstrained);

        arma::vec weightedScore_uu =
          arma::as_scalar(
            weightInf(iIter)
          )
          *
          dLu_uu;

        bool constrainedSolveOK =
          solveInformationSystem(
            delta_uu,
            InfMat_uu,
            weightedScore_uu,
            "free variance-component block"
          );

        if(!constrainedSolveOK){
          Rcpp::stop(
            "Unable to solve the free variance-component information block."
          );
        }

        delta(unconstrained) =
          delta_uu;

        expectedNewTheta(unconstrained) =
          thetaUnlisted(unconstrained)
          -
          delta_uu;

        applyVarianceBounds(
          expectedNewTheta,
          unconstrained,
          true
        );
      }
    }

    if(vcMap.active){
      expectedNewTheta =
        projectVarianceParameters(vcMap, expectedNewTheta, unconstrained, InfMat);
      applyVarianceBounds(expectedNewTheta, unconstrained, true);
      expectedNewTheta =
        projectVarianceParameters(vcMap, expectedNewTheta, unconstrained,
                                  arma::eye<arma::mat>(nVcTotal, nVcTotal));
    }

    // C) Delta-change bookkeeping is intentionally delayed until AFTER
    //    local PD repair and global working-coordinate trust scaling, so the
    //    recorded delta is the actual proposal sent to the likelihood line
    //    search rather than the raw unconstrained AI proposal.

    // D) Positive-definiteness fallback PER covariance structure.
    //    A failing structure is repaired locally; valid covariance
    //    structures keep their original AI/EM update unchanged.
    std::string notPD = "-";

    for(int iStruct = 0; iStruct < nRRe; ++iStruct){

      arma::uvec structIdx =
        arma::regspace<arma::uvec>(
          static_cast<arma::uword>(
            nVcStart(iStruct) - 1
          ),
          static_cast<arma::uword>(1),
          static_cast<arma::uword>(
            nVcEnd(iStruct) - 1
          )
        );

      auto structureIsPD =
        [&](const arma::vec & candidate) -> bool {

          const arma::vec localPar =
            candidate(structIdx);

          // Respect descriptor bounds before even constructing Sigma.
          for(arma::uword k = 0; k < localPar.n_elem; ++k){

            const double lo =
              covLower(iStruct)(k);

            const double hi =
              covUpper(iStruct)(k);

            if(
                std::isfinite(lo)
                &&
                localPar(k) < lo
            ){
              return false;
            }

            if(
                std::isfinite(hi)
                &&
                localPar(k) > hi
            ){
              return false;
            }
          }

          return
            descriptorIsPositiveDefinite(
              covDescriptor[static_cast<std::size_t>(iStruct)],
              localPar,
              pdCheckFloor
            );
        };

      if(structureIsPD(expectedNewTheta)){
        continue;
      }

      if(notPD == "-"){
        notPD = std::to_string(iStruct + 1);
      }else{
        notPD += "," + std::to_string(iStruct + 1);
      }

      // Identify genuinely free parameters belonging only to this
      // covariance structure.
      std::vector<arma::uword> localFreeStd;

      for(arma::uword a = 0; a < structIdx.n_elem; ++a){

        const arma::uword g =
          structIdx(a);

        if(constrainedMask(g) == 0){
          localFreeStd.push_back(g);
        }
      }

      arma::uvec localFree(
        localFreeStd.size()
      );

      for(std::size_t k = 0; k < localFreeStd.size(); ++k){
        localFree(k) =
          localFreeStd[k];
      }

      bool repaired =
        false;

      arma::vec bestCandidate =
        expectedNewTheta;

      // ------------------------------------------------------------
      // One-dimensional descriptor structures still use log(sigma2), so the
      // parameter itself is not a variance directly. Reverting unconditionally
      // to the previous accepted value produces zero forward progress every
      // iteration for a boundary-approaching variance component (true value
      // at/near zero): the identical oversized step gets proposed and
      // rejected again next iteration, freezing the parameter permanently.
      // Instead, backtrack along the proposed step, trying progressively
      // smaller fractions of it until the structure is PD again - the same
      // step-halving idea used below for multi-parameter structures. Falls
      // through to the shared final fallback (full revert) further below
      // only if no fraction of the step is repairable (e.g. a genuine
      // non-finite proposal).
      // ------------------------------------------------------------
      if(covarianceDim[static_cast<std::size_t>(iStruct)] == 1 &&
        covPar(iStruct).n_elem == 1){

        const arma::uword g0 = structIdx(0);
        const double previousValue = thetaUnlisted(g0);
        const double proposedValue = expectedNewTheta(g0);

        if(std::isfinite(proposedValue)){

          arma::vec candidate = expectedNewTheta;
          double stepFraction = 1.0;

          for(int halvingAttempt = 0; halvingAttempt < 10 && !repaired; ++halvingAttempt){
            stepFraction *= 0.5;
            candidate(g0) = previousValue + stepFraction * (proposedValue - previousValue);
            repaired = structureIsPD(candidate);
          }

          if(repaired){
            bestCandidate = candidate;
          }
        }
      }

      if(
          !repaired
          &&
          localFree.n_elem > 0
          &&
          covarianceDim[static_cast<std::size_t>(iStruct)] > 1
      ){

        // ----------------------------------------------------------
        // First local fallback: 50/50 AI-EM information.
        // ----------------------------------------------------------
        arma::mat localInfo =
          0.5
          *
          avInf(
            localFree,
            localFree
          )
          +
          0.5
          *
          emInf(
            localFree,
            localFree
          );

        arma::vec localScore =
          arma::as_scalar(
            weightInf(iIter)
          )
          *
          dLu(localFree);

        arma::vec localDelta;

        bool localOK =
          solveInformationSystem(
            localDelta,
            localInfo,
            localScore,
            "per-structure 50/50 AI-EM PD fallback"
          );

        if(localOK){

          bestCandidate =
            expectedNewTheta;

          bestCandidate(localFree) =
            thetaUnlisted(localFree)
            -
            localDelta;

          applyVarianceBounds(
            bestCandidate,
            localFree,
            true
          );

          repaired =
            structureIsPD(
              bestCandidate
            );
        }

        // ----------------------------------------------------------
        // Second local fallback: pure EM block.
        // ----------------------------------------------------------
        if(!repaired){

          localInfo =
            emInf(
              localFree,
              localFree
            );

          localOK =
            solveInformationSystem(
              localDelta,
              localInfo,
              localScore,
              "per-structure pure-EM PD fallback"
            );

          if(localOK){

            bestCandidate =
              expectedNewTheta;

            bestCandidate(localFree) =
              thetaUnlisted(localFree)
              -
              localDelta;

            applyVarianceBounds(
              bestCandidate,
              localFree,
              true
            );

            repaired =
              structureIsPD(
                bestCandidate
              );
          }
        }

        // ----------------------------------------------------------
        // Third local fallback: step-halving from the current accepted
        // covariance structure toward the best local proposal.
        //
        // The current point is the safest anchor.  Thirty halvings are
        // cheap because covariance structures are tiny and provide a
        // much more robust approach to the open PD boundary.
        // ----------------------------------------------------------
        if(!repaired){

          const arma::vec target =
            bestCandidate;

          arma::vec base =
            expectedNewTheta;

          base(structIdx) =
            thetaUnlisted(structIdx);

          // Preserve any model-defined fixed values in the candidate.
          if(modelFixed.n_elem > 0){
            for(arma::uword kk = 0; kk < modelFixed.n_elem; ++kk){
              const arma::uword g = modelFixed(kk);

              if(
                  g >= structIdx.min()
                  &&
                  g <= structIdx.max()
              ){
                base(g) =
                  expectedNewTheta(g);
              }
            }
          }

          for(
              int half = 0;
              half < 30 && !repaired;
              ++half
          ){

            const double alpha =
              std::pow(
                0.5,
                static_cast<double>(half)
              );

            arma::vec trial =
              base;

            trial(localFree) =
              base(localFree)
              +
              alpha
              *
              (
                target(localFree)
                -
                base(localFree)
              );

            applyVarianceBounds(
              trial,
              localFree,
              false
            );

            if(structureIsPD(trial)){
              bestCandidate =
                trial;

              repaired =
                true;
            }
          }
        }

      }

      // ------------------------------------------------------------
      // Final safe fallback:
      //
      // Reject ONLY this covariance structure's proposed step and
      // retain its previous accepted value.  This allows the other
      // covariance structures to continue updating instead of
      // terminating an otherwise well-behaved REML fit.
      // ------------------------------------------------------------
      if(!repaired){

        arma::vec rejectedCandidate =
          expectedNewTheta;

        rejectedCandidate(structIdx) =
          thetaUnlisted(structIdx);

        if(structureIsPD(rejectedCandidate)){

          bestCandidate =
            rejectedCandidate;

          repaired =
            true;
        }
      }

      // If the previous accepted covariance structure itself is not PD,
      // we have an internal/model-specification inconsistency rather than
      // merely a bad proposed update.
      if(!repaired){

        Rcpp::stop(
          "Covariance structure "
          +
          std::to_string(iStruct + 1)
          +
          " is not positive definite even at its previous accepted value."
        );
      }

      expectedNewTheta(structIdx) =
        bestCandidate(structIdx);

      delta(structIdx) =
        thetaUnlisted(structIdx)
        -
        expectedNewTheta(structIdx);
    }

    // ==========================================================
    // GLOBAL WORKING-COORDINATE TRUST SCALING
    // ==========================================================
    //
    // expectedNewTheta is already model-fixed / boundary / PD repaired.
    // Limit the largest transformed-coordinate move while preserving the
    // complete joint AI direction.  The exact likelihood line search above
    // remains the final acceptance criterion.
    // ==========================================================

    arma::vec proposalStep =
      expectedNewTheta
      -
      thetaUnlisted;

    double trustScale = 1.0;

    for(arma::uword kk = 0; kk < unconstrained.n_elem; ++kk){

      const arma::uword g =
        unconstrained(kk);

      const double a =
        std::abs(
          proposalStep(g)
        );

      const double cap =
        workingStepCap(g);

      if(
          std::isfinite(a)
          &&
          std::isfinite(cap)
          &&
          cap > 0.0
          &&
          a > cap
      ){
        trustScale =
          std::min(
            trustScale,
            cap / a
          );
      }
    }

    if(trustScale < 1.0){

      expectedNewTheta =
        thetaUnlisted
        +
        trustScale
        *
        (
          expectedNewTheta
          -
          thetaUnlisted
        );

      // Fixed / dynamically constrained parameters must remain exact.
      if(constrained.n_elem > 0){
        expectedNewTheta(constrained) =
          thetaUnlisted(constrained);
      }

      delta =
        thetaUnlisted
        -
        expectedNewTheta;
    }

    // C) Quantify changes in the ACTUAL safeguarded proposal.
    const arma::vec previousDelta = useHybridAI ? delta_minus1 : arma::vec();
    if(iIter == 0){

      delta_minus1 =
        delta;

    }else{

      percDelta.col(iIter) =
        delta
        /
        delta_minus1;

      delta_minus1 =
        delta;
    }

    // #######################
    // # 8) Save CURRENT ACCEPTED parameter vector in the monitor
    // #######################
    //
    // llik(iIter), C, b, u, thetaUnlisted and this monitor column now all
    // describe the same accepted parameter point.  The new proposal is not
    // committed until after the stopping decision below.
    monitor.col(iIter) =
      thetaUnlisted;

    // #######################
    // # 9) Stopping criteria
    // #######################


    normMonitor(0,iIter) =
      arma::norm(
        delta(unconstrained),
        1
      );

    normMonitor(1,iIter) =
      arma::norm(
        dLu(unconstrained),
        1
      );

    arma::vec stopDelta =
      delta(unconstrained);

    arma::vec stopTheta =
      thetaUnlisted(unconstrained);

    normMonitor(2,iIter) =
      arma::norm(
        stopDelta,
        2
      )
      /
      (
        arma::norm(
          stopTheta,
          2
        )
        +
        tolParInv
      );

    if(useHybridAI){
      const bool activeSetChanged = hybridFree.n_elem != unconstrained.n_elem ||
        (hybridFree.n_elem == unconstrained.n_elem && arma::any(hybridFree != unconstrained));
      hybridFree = unconstrained;
      const bool stopHint = normMonitor(2,iIter) < tolParConvNorm ||
        (iIter > 0 && std::abs(llik(iIter) - llik(iIter-1)) < tolParConvLL);
      if(!finalAIExact && (stopHint || lineSearchStalled || notPD != "-" || activeSetChanged)){
        hybridForceExact = true;
        ++hybridAICertifications;
        toBoundary.row(iIter).zeros();
        delta_minus1 = previousDelta;
        --iIter;
        continue;
      }
      if(activeSetChanged){ hybridState.ready = false; }
      finalProjectedScoreNorm = 0.0;
      for(arma::uword k = 0; k < unconstrained.n_elem; ++k){
        const arma::uword g = unconstrained(k);
        const double point = thetaUnlisted(g);
        const double lower = parLowerUnlisted(g);
        const double upper = parUpperUnlisted(g);
        if((std::isfinite(lower) && point <= lower + 1.0e-8 * (1.0 + std::abs(lower)) && dLu(g) >= 0.0) ||
           (std::isfinite(upper) && point >= upper - 1.0e-8 * (1.0 + std::abs(upper)) && dLu(g) <= 0.0)){
          continue;
        }
        finalProjectedScoreNorm = std::max(finalProjectedScoreNorm,
          std::abs(dLu(g)) / std::sqrt(std::max(tolParInv, std::abs(avInf(g,g)))));
      }
    }

    if(verbose == true){

      if(iIter == 0){
        std::string solverDesc;
        if(saddleActive){
          solverDesc = "saddle-point SparseLU";
        }else if(blockEngine.active){
          if(blockEngine.latent){
            solverDesc = "latent Woodbury Schur";
          }else if(blockEngine.chain){
            solverDesc = "block-chain Schur";
          }else if(blockEngine.banded){
            solverDesc = "block-banded Schur";
          }else{
            solverDesc = "dense block-Schur";
          }
        }else if(dFirstMLActive){
          solverDesc = "D-first ML Schur";
        }else if(matrixFreePCGReady){
          solverDesc = "matrix-free PCG";
        }else if(solverName == "cholmod"){
          solverDesc = "sparse CHOLMOD";
        }else if(solverName == "ldlt"){
          solverDesc = "sparse LDLT";
        }else if(solverName == "pcg"){
          solverDesc = "sparse PCG";
        }else{
          solverDesc = solverName;
        }

        std::string residDesc = residualPlan.diagonal() ? "diagonal"
          : residualPlan.grid() ? "grid"
          : residualPlan.kronecker() ? "Kronecker"
          : "sparse";

        std::string traceDesc = saddleActive ? "exact augmented panels"
          : (solverName == "ldlt") ? "Takahashi"
          : (solverName == "pcg") ? "Hutchinson"
          : reml ? (blockEngine.active ? "Schur-blocks" : "CHOLMOD-blocks")
          : "D-trace";

        Rcpp::Rcout << "MME path: " << solverDesc << " + " << residDesc
                    << " residual (" << traceDesc << " traces) ["
                    << sommer::mmePathStatus(solverName, rrConditional,
                      blockEngine.active, blockEngine.active && blockEngine.chain,
                      blockEngine.active && blockEngine.banded,
                      blockEngine.active && blockEngine.latent, blockSchurMode,
                      Rgrid ? "grid" : Rdiag ? "diagonal" : Rkron ? "kronecker" : "sparse")
                    << "]" << arma::endl;

        progress.header();
      }

      progress.row(iIter + 1, llik(iIter), iIter > 0 ? llik(iIter - 1) : NA_REAL,
        useHybridAI ? (finalAIExact ? "AI" : "hybrid-AI")
          : sommer::MMEProgress::update(arma::as_scalar(weightEmInf(iIter))),
        "off", trustScale, lineSearchAlpha, notPD);
    }

    dLuOut =
      dLu;

    // A failed global line search means the optimizer could not identify an
    // improving step from this accepted point.  Return the accepted solution
    // rather than cycling indefinitely.  This is a stall, not a likelihood
    // decrease falsely labelled as convergence.
    if(lineSearchStalled){

      convergence = false;

      monitor = monitor.cols(0,iIter);
      normMonitor = normMonitor.cols(0,iIter);
      percDelta = percDelta.cols(0,iIter);
      llik = llik.cols(0,iIter);

      if(verbose){
        Rcpp::Rcout
          << "Optimization stopped at the last accepted REML point because the global line search stalled."
          << arma::endl;
      }

      break;
    }

    bool shouldStop = false;

    if(iIter > 0){

      const double delta_llik =
        llik(iIter)
        -
        llik(iIter-1);

      const bool likelihoodConverged =
        std::abs(
          delta_llik
        )
        <
        tolParConvLL;

      const bool parameterConverged =
        normMonitor(2,iIter)
        <
        tolParConvNorm;

      if(
          likelihoodConverged
          ||
          parameterConverged
      ){
        const double scoreTolerance =
          std::sqrt(std::max(tolParConvLL, tolParConvNorm));
        const bool certified = !useHybridAI ||
          (finalAIExact && parameterConverged && finalProjectedScoreNorm <= scoreTolerance);
        convergence = certified;
        shouldStop = certified;
      }
    }

    // Never return an unevaluated proposal simply because nIters was reached.
    if(iIter == nIters - 1){
      shouldStop = true;
    }

    if(shouldStop){

      monitor = monitor.cols(0,iIter);
      normMonitor = normMonitor.cols(0,iIter);
      percDelta = percDelta.cols(0,iIter);
      llik = llik.cols(0,iIter);

      break;
    }

    // ==========================================================
    // Queue the new joint proposal for exact-likelihood checking.
    // ==========================================================

    if(vcMap.active){
      expectedNewTheta =
        projectVarianceParameters(vcMap, expectedNewTheta, unconstrained,
                                  arma::eye<arma::mat>(nVcTotal, nVcTotal));
    }

    lineSearchBase =
      thetaUnlisted;

    lineSearchTarget =
      expectedNewTheta;

    lineSearchAlpha =
      1.0;

    lineSearchHalvings =
      0;

    pendingLineSearch =
      true;
    pendingHybridProposal = useHybridAI && !exactAIThisIteration;

    bool anyStructureReverted = false;

    phaseProfile.enter(sommer::FitPhase::Proposal);
    for(int i = 0; i < nRRe; ++i){

      arma::uvec toFill =
        arma::regspace<arma::uvec>(
          nVcStart(i)-1,
          1,
          nVcEnd(i)-1
        );

      const arma::vec proposedParameters =
        expectedNewTheta(toFill);

      arma::mat proposedTheta;

      if(tryEvaluateStructure(i, proposedParameters, proposedTheta)){
        covPar(i) =
          proposedParameters;

        theta(i) =
          proposedTheta;
      }else{
        anyStructureReverted = true;
        covPar(i) =
          lineSearchBase(toFill);

        lineSearchTarget(toFill) =
          lineSearchBase(toFill);

        if(verbose){
          Rcpp::Rcout
            << "Covariance structure "
            << i + 1
            << " became invalid after joint trust scaling; retaining its previous accepted value."
            << arma::endl;
        }
      }

    }

    // A partial revert would break parameter constraints; revert every structure.
    if(vcMap.active && anyStructureReverted){
      lineSearchTarget = lineSearchBase;
      for(int i = 0; i < nRRe; ++i){
        const arma::vec baseParameters =
          lineSearchBase.subvec(nVcStart(i)-1, nVcEnd(i)-1);
        arma::mat baseTheta;
        covPar(i) = baseParameters;
        if(tryEvaluateStructure(i, baseParameters, baseTheta)){
          theta(i) = baseTheta;
        }
      }
    }
  }// end of iterative optimization
  phaseProfile.enter(sommer::FitPhase::Final);

  // Matrix-free PCG omits random-prior entries from C during REML
  // iterations. Materialize them once at convergence because C is part of
  // the fitted-object contract used by prediction and PEV methods.
  if(matrixFreePCGReady && nZs > 0){
    std::vector<EigenTriplet> finalPriorTriplets;

    for(int iR = 0; iR < nRe; ++iR){
      const arma::mat & finalPrecision =
        lastEvaluatedRandomPrecision[static_cast<std::size_t>(iR)];
      if(finalPrecision.n_rows == 0 || !finalPrecision.is_finite()){
        Rcpp::stop("Unable to materialize the final random-effect precision after matrix-free PCG.");
      }

      const std::size_t iCache = static_cast<std::size_t>(iR);
      for(arma::uword iRow = 0; iRow < finalPrecision.n_rows; ++iRow){
        const arma::uword rowStart =
          partitionStartCache[iCache][static_cast<std::size_t>(iRow)];
        for(arma::uword iCol = 0; iCol < finalPrecision.n_cols; ++iCol){
          const double coefficient = finalPrecision(iRow,iCol);
          if(coefficient == 0.0){ continue; }
          const arma::uword colStart =
            partitionStartCache[iCache][static_cast<std::size_t>(iCol)];
          for(arma::sp_mat::const_iterator relEntry = Ai(iR).begin();
              relEntry != Ai(iR).end(); ++relEntry){
            finalPriorTriplets.emplace_back(
              static_cast<int>(rowStart + relEntry.row()),
              static_cast<int>(colStart + relEntry.col()),
              coefficient * (*relEntry)
            );
          }
        }
      }
    }

    EigenSpMat finalPrior(nEffects, nEffects);
    finalPrior.setFromTriplets(
      finalPriorTriplets.begin(),
      finalPriorTriplets.end()
    );
    finalPrior.makeCompressed();
    C = C + finalPrior;
    C.makeCompressed();
  }
  
  
  double AIC = (-2 * llik((llik.n_cols-1))) + (2 * nX);
  double BIC = (-2 * llik((llik.n_cols-1))) + (log(nR) * nX);
  // move constraints to vector form binding the columns
  arma::vec thetaCUnlistedFinal;
  for (int i = 0; i < nRRe; ++i) {
    thetaCUnlistedFinal =
      arma::join_cols(
        thetaCUnlistedFinal,
        covConstraint(i)
      );
  }
  // ============================================================
  // FINAL-OUTPUT COMPUTATIONS ONLY
  // ============================================================

  // Final inverse information matrix for theta_se.
  // It is formed once AFTER optimisation and is never used for updates.
  // avInf is the AI matrix of -2 logL, so Var(theta) = 2 * inverse(avInf);
  // the EM-blended update matrix InfMat must not be used here.
  if(retainCovarianceInformation){
    arma::uvec freeIdx =
      arma::find(thetaCUnlistedFinal != 3);

    InfMatInv.zeros(
      nVcTotal,
      nVcTotal
    );

    if(freeIdx.n_elem > 0){
      arma::mat aiFree =
        avInf.submat(
          freeIdx,
          freeIdx
        );
      aiFree = 0.5 * (aiFree + aiFree.t());

      arma::mat aiFreeInv;
      bool finalInfoInverseOK = false;
      // Components at a boundary make AI near-singular; use a pseudo-inverse then.
      if(arma::rcond(aiFree) > 1.0e-12){
        finalInfoInverseOK =
          arma::inv_sympd(
            aiFreeInv,
            aiFree
          );
      }

      if(!finalInfoInverseOK){
        aiFreeInv = arma::pinv(aiFree);
      }

      if(vcMap.active){
        arma::mat Ts = vcMap.T.rows(freeIdx);
        arma::uvec used = arma::find(arma::sum(arma::abs(Ts), 0).t() > 0);
        Ts = Ts.cols(used);
        aiFreeInv = Ts * arma::pinv(Ts.t() * aiFree * Ts) * Ts.t();
      }

      InfMatInv.submat(
        freeIdx,
        freeIdx
      ) = 2.0 * aiFreeInv;
    }
  }else{
    InfMatInv.reset();
  }

  // ------------------------------------------------------------
  // Optional final C-inverse work controlled by computeCi:
  //
  //   0 = no extra C-related work after REML convergence
  //   1 = final Takahashi sparse inverse subset
  //   2 = final full C inverse
  //
  // In mode 0, Ci and uPevList are returned empty and C is not
  // factorised again after the iterative optimisation.
  // ------------------------------------------------------------

  if(saddleActive && retainCoefficientMatrix){
    SommerSpatialLDLT finalResidualFactor;
    finalResidualFactor.compute(saddleResidual);
    if(finalResidualFactor.info() != Eigen::Success ||
       !finalResidualFactor.vectorD().allFinite() ||
       finalResidualFactor.vectorD().minCoeff() <= 0.0){
      Rcpp::stop("Residual covariance factorization failed for final compatibility C.");
    }
    C = sommer::assembleNormalEquations(saddleDesign, saddlePrior,
      [&](const Eigen::MatrixXd & rhs) -> Eigen::MatrixXd {
        Rcpp::checkUserInterrupt();
        const Eigen::MatrixXd solved = finalResidualFactor.solve(rhs);
        if(finalResidualFactor.info() != Eigen::Success || !solved.allFinite()){
          Rcpp::stop("Residual precision solve failed for final compatibility C.");
        }
        return solved;
      }, sommer::WorkspacePolicy::panelColumns(), maxFallbackEntries);
    ++saddleCompatibilityAssemblies;
  }

  // Pre-existing bug fix: arma::sp_mat::reset() clears BOTH content and
  // dimensions (to 0x0), so it must be followed by set_size() before the
  // computeCi>0 fill loops below index into Ci.
  Ci.reset();
  Ci.set_size(nEffects, nEffects);
  SelectedInverseSubset & finalCselected = CselectedTopology;

  if(computeCi > 0){

    if(solverName == "ldlt" && !CnumericReady){
      Rcpp::stop("No final sparse LDLT factorisation of C is available.");
    }
    if(solverName == "pcg" && !CpcgReady){
      Rcpp::stop("No final PCG coefficient operator is available.");
    }
    if(solverName == "cholmod" && !CnumericReady){
      Rcpp::stop("No final CHOLMOD factorisation of C is available.");
    }

    // Reuse the numerical factorisation from the final REML iteration.
    // C has not changed since that factorisation.

    if(computeCi == 1){

      bool usedCholmodSelectedInverse = false;
      if(!saddleActive && solverName == "cholmod" && cholmodFinalSelectedMode != "off" && cholmodCFactorCurrent &&
         cholmodState.factor != nullptr &&
         cholmodState.factor->is_super && cholmodState.factor->is_ll &&
         cholmodState.factor->xtype == CHOLMOD_REAL &&
         cholmodState.factor->dtype == CHOLMOD_DOUBLE &&
         cholmodState.factor->itype == CHOLMOD_INT &&
         static_cast<int>(cholmodState.factor->n) == nEffects &&
         static_cast<double>(cholmodState.factor->xsize) <= cholmodInvCacheMaxDoubles){
        cholmodInvCache.ready = false;
        cholmodInvCache.Cxx.reset();
        cholmodInvCache.Cuu.clear();
        cholmodInvCache.Cxu.clear();
        cholmodInvCache.FS.clear();
        cholmodInvCache.crossBlocks.clear();
        const cholmod_factor * factor = cholmodState.factor;
        const int n = static_cast<int>(factor->n);
        const int * permutation = static_cast<const int *>(factor->Perm);
        std::vector<int> originalToPermuted(static_cast<std::size_t>(n));
        for(int permuted = 0; permuted < n; ++permuted){
          originalToPermuted[static_cast<std::size_t>(permutation[permuted])] = permuted;
        }
        if(buildCholmodSelectedInverse(factor, originalToPermuted, cholmodInvCache)){
          const int * super = static_cast<const int *>(factor->super);
          const int * pi = static_cast<const int *>(factor->pi);
          const int * px = static_cast<const int *>(factor->px);
          const int * rows = static_cast<const int *>(factor->s);
          for(std::size_t node = 0; node < factor->nsuper; ++node){
            const int first = super[node];
            const int columns = super[node + 1] - first;
            const int rowOffset = pi[node];
            const int rowCount = pi[node + 1] - rowOffset;
            for(int localColumn = 0; localColumn < columns; ++localColumn){
              const int permutedColumn = first + localColumn;
              const int originalColumn =
                permutation[permutedColumn];
              for(int localRow = 0; localRow < rowCount; ++localRow){
                const int permutedRow = rows[rowOffset + localRow];
                if(permutedRow < permutedColumn){ continue; }
                const int originalRow = permutation[permutedRow];
                const double value = cholmodInvCache.selZ[
                  static_cast<std::size_t>(px[node]) +
                  static_cast<std::size_t>(localColumn) * rowCount + localRow];
                Ci(static_cast<arma::uword>(originalRow),
                   static_cast<arma::uword>(originalColumn)) = value;
                if(originalRow != originalColumn){
                  Ci(static_cast<arma::uword>(originalColumn),
                     static_cast<arma::uword>(originalRow)) = value;
                }
              }
            }
          }
          usedCholmodSelectedInverse = true;
          cholmodFinalSelectedInverse = true;
        }
      }

      // Schur engines have no marginal CHOLMOD factor, and unsupported or
      // over-budget selected workspaces retain the exact Eigen fallback.
      if((saddleActive || solverName == "cholmod") && !usedCholmodSelectedInverse){
        Cfactor.analyzePattern(C);
        if(Cfactor.info() != Eigen::Success){
          Rcpp::stop("Sparse LDLT symbolic analysis failed while preparing the final Takahashi extraction.");
        }
        Cfactor.factorize(C);
        if(Cfactor.info() != Eigen::Success){
          Rcpp::stop("Sparse LDLT factorisation failed while preparing the final Takahashi extraction.");
        }
      }

      if(!usedCholmodSelectedInverse){
        // Takahashi sparse inverse subset in the filled LDLT pattern.
        if(solverName != "ldlt" || !CselectedTopologyReady){
        buildSelectedInverseSubset(
          Cfactor,
          "final C",
          finalCselected,
          false
        );
        CselectedTopologyReady = true;
        }

        // Reconstruct the selected inverse subset in the ORIGINAL ordering.
        std::vector<int> permutedToOriginal(nEffects, -1);

        for(int original = 0; original < nEffects; ++original){
        const int permuted =
          finalCselected.originalToPermuted[original];

        if(permuted < 0 || permuted >= nEffects){
          Rcpp::stop(
            "Invalid permutation while reconstructing the final Takahashi inverse subset."
          );
        }

        permutedToOriginal[permuted] =
          original;
        }

        for(int pcol = 0; pcol < nEffects; ++pcol){

        const int originalCol =
          permutedToOriginal[pcol];

        if(originalCol < 0){
          Rcpp::stop(
            "Unable to invert the final LDLT permutation while reconstructing Ci."
          );
        }

        const std::vector<int> & rr =
          finalCselected.rows[pcol];

        const std::vector<double> & vv =
          finalCselected.values[pcol];

        for(std::size_t k = 0; k < rr.size(); ++k){

          const int originalRow =
            permutedToOriginal[rr[k]];

          if(originalRow < 0){
            Rcpp::stop(
              "Unable to map a Takahashi inverse entry back to original MME ordering."
            );
          }

          const double value =
            vv[k];

          Ci(
            static_cast<arma::uword>(originalRow),
            static_cast<arma::uword>(originalCol)
          ) = value;

          if(originalRow != originalCol){
            Ci(
              static_cast<arma::uword>(originalCol),
              static_cast<arma::uword>(originalRow)
            ) = value;
          }
        }
        }
      }

    }else if(computeCi == 2){

      arma::mat finalCiDense(
        static_cast<arma::uword>(nEffects),
        static_cast<arma::uword>(nEffects)
      );
      const Eigen::Index inverseSize = static_cast<Eigen::Index>(nEffects);
      const Eigen::Index inverseBatch = static_cast<Eigen::Index>(panelCi2Batch);
      for(Eigen::Index first = 0; first < inverseSize; first += inverseBatch){
        const Eigen::Index count = std::min(inverseBatch, inverseSize - first);
        Eigen::MatrixXd identityPanel = Eigen::MatrixXd::Zero(inverseSize, count);
        for(Eigen::Index column = 0; column < count; ++column){
          identityPanel(first + column, column) = 1.0;
        }
        const Eigen::MatrixXd inversePanel =
          solveCMatrix(identityPanel, "batched final full inverse of C");
        for(Eigen::Index column = 0; column < count; ++column){
          std::copy(inversePanel.col(column).data(),
            inversePanel.col(column).data() + inverseSize,
            finalCiDense.colptr(static_cast<arma::uword>(first + column)));
        }
      }

      finalCiDense =
        0.5
        *
        (
          finalCiDense
          +
          finalCiDense.t()
        );

      Ci =
        arma::sp_mat(
          finalCiDense
        );
    }
  }

  // bring back to original scale the variance components
  for (int i = 0; i < nRRe; ++i) {
    theta(i) = theta(i)*vary;
  }
  bu = bu*stdy;
  b = b*stdy;
  if(nZs > 0){ u = bu(uInd); }
  if(hasIntercept && muy != 0.0){
    for(arma::uword j = 0; j < b.n_elem; ++j){
      b(j) += muy * cShift(j);
      bu(j) += muy * cShift(j);
    }
  }
  if(computeCi > 0 && Ci.n_elem > 0){ Ci = Ci*vary; }
  // Transform covariance-parameter uncertainty and monitor output from
  // unconstrained working coordinates to reported natural coordinates.
  auto reportParameters =
    [&](const int iStruct,
        const arma::vec & work,
        arma::vec * jacobianOut = nullptr) -> arma::vec {

      arma::vec out = work;
      arma::vec jac(work.n_elem, arma::fill::ones);

      // Product-level scale is the only covariance convention intentionally
      // owned by the solver: internal tau = log(sigma2 / var(y)).
      out(0) = std::exp(work(0)) * vary;
      jac(0) = out(0);

      Rcpp::List cs =
        covDescriptor[static_cast<std::size_t>(iStruct)];
      Rcpp::List factors = cs["factors"];

      for(int fidx = 0; fidx < factors.size(); ++fidx){
        Rcpp::List f = Rcpp::as<Rcpp::List>(factors[fidx]);
        const int start1 = Rcpp::as<int>(f["par_start"]);
        const int end1 = Rcpp::as<int>(f["par_end"]);
        if(end1 < start1){
          continue;
        }

        if(!f.containsElementNamed("report")){
          Rcpp::stop("CovarianceFactor is missing report specification.");
        }
        Rcpp::List report = Rcpp::as<Rcpp::List>(f["report"]);
        const std::string backend =
          Rcpp::as<std::string>(report["backend"]);

        const int nFactorPar = end1 - start1 + 1;

        if(backend == "builtin"){
          Rcpp::CharacterVector tr = report["transform"];
          Rcpp::NumericVector lo = report["lower"];
          Rcpp::NumericVector hi = report["upper"];
          if(tr.size() != nFactorPar || lo.size() != nFactorPar ||
             hi.size() != nFactorPar){
            Rcpp::stop("CovarianceFactor report specification has incompatible length.");
          }

          for(int local = 0; local < nFactorPar; ++local){
            const arma::uword k =
              static_cast<arma::uword>(start1 - 1 + local);
            const std::string code = Rcpp::as<std::string>(tr[local]);
            const double eta = work(k);

            if(code == "identity"){
              out(k) = eta;
              jac(k) = 1.0;

            }else if(code == "exp"){
              out(k) = std::exp(eta);
              jac(k) = out(k);

            }else if(code == "tanh"){
              out(k) = std::tanh(eta);
              jac(k) = 1.0 - out(k)*out(k);

            }else if(code == "bounded_logit"){
              const double lower = lo[local];
              const double upper = hi[local];
              if(!std::isfinite(lower) || !std::isfinite(upper) ||
                 !(lower < upper)){
                Rcpp::stop("Invalid bounded-logit reporting interval.");
              }
              const double logistic =
                eta >= 0.0
                ? 1.0 / (1.0 + std::exp(-eta))
                : std::exp(eta) / (1.0 + std::exp(eta));
              out(k) = lower + (upper-lower)*logistic;
              jac(k) = (upper-lower)*logistic*(1.0-logistic);

            }else{
              Rcpp::stop("Unknown CovarianceFactor reporting transform: " + code);
            }
          }

        }else{
          Rcpp::stop("Unknown CovarianceFactor reporting backend: " + backend);
        }
      }

      if(jacobianOut != nullptr){
        (*jacobianOut) = jac;
      }
      return out;
    };

  arma::field<arma::vec> covParOut(nRRe);
  arma::field<arma::vec> covParWorking(nRRe);

  arma::vec finalJacobian(
    nVcTotal,
    arma::fill::ones
  );

  arma::mat monitorOut =
    monitor;

  arma::uword offset =
    0;

  for(int i = 0; i < nRRe; ++i){

    arma::vec localJac;

    covParOut(i) =
      reportParameters(
        i,
        covPar(i),
        &localJac
      );
    covParWorking(i) = covPar(i);
    covParWorking(i)(0) += std::log(vary);

    const arma::uword localN =
      covPar(i).n_elem;

    finalJacobian.subvec(
      offset,
      offset + localN - 1
    ) =
      localJac;

    for(arma::uword iter = 0; iter < monitor.n_cols; ++iter){

      if(
          !std::isfinite(
            monitor(
              offset,
              iter
            )
          )
      ){
        continue;
      }

      arma::vec localWork =
        monitor.submat(
          offset,
          iter,
          offset + localN - 1,
          iter
        );

      monitorOut.submat(
        offset,
        static_cast<arma::uword>(iter),
        offset + localN - 1,
        static_cast<arma::uword>(iter)
      ) =
        reportParameters(
          i,
          localWork
        );
    }

    offset +=
      localN;
  }

  if(
      InfMatInv.n_rows == finalJacobian.n_elem
      &&
      InfMatInv.n_cols == finalJacobian.n_elem
  ){

    InfMatInv %=
      finalJacobian
      *
      finalJacobian.t();
  }

  // dLuOut=dLuOut/vary;
  // move the effects from vector to a field with matrices
  arma::field<arma::mat> uList(nRe), uPevList(nRe); // store indices of the random effects

  for (int i = 0; i < nRe; ++i) {

    arma::mat partitionsP = partitions(i);

    arma::mat uMat(
      arma::as_scalar(partitionsP(0,1) - partitionsP(0,0) + 1),
      partitionsP.n_rows
    );

    arma::mat eMat;

    if(computeCi > 0){
      eMat.set_size(
        arma::as_scalar(partitionsP(0,1) - partitionsP(0,0) + 1),
        partitionsP.n_rows
      );
    }

    for (int j = 0; j < partitionsP.n_rows; ++j) {

      const arma::uword firstIndex =
        static_cast<arma::uword>(partitionsP(j,0)-1);

      const arma::uword lastIndex =
        static_cast<arma::uword>(partitionsP(j,1)-1);

      uMat.col(j) =
        bu.submat(
          firstIndex,
          0,
          lastIndex,
          0
        );

      if(computeCi == 1){

        // Diagonal entries are always represented by either selected backend.
        arma::vec pevDiag(lastIndex - firstIndex + 1);

        for(arma::uword k = firstIndex; k <= lastIndex; ++k){
          const double cii = Ci(k, k);
          if(!std::isfinite(cii) || cii <= 0.0){
            Rcpp::stop(
              "Final selected inverse is missing a valid diagonal entry required for uPevList."
            );
          }

          pevDiag(k - firstIndex) =
            cii;
        }

        eMat.col(j) =
          pevDiag;

      }else if(computeCi == 2){

        // Ci has already been rescaled to the original response scale.
        eMat.col(j) =
          arma::diagvec(
            Ci.submat(
              firstIndex,
              firstIndex,
              lastIndex,
              lastIndex
            )
          );
      }
    }

    uList(i) = uMat;

    if(computeCi > 0){
      uPevList(i) = eMat;
    }else{
      // Mode 0: list shape is retained, but no PEV/inverse work is done.
      uPevList(i) = arma::mat();
    }
  }

  // return results in a list form
  // W/C are native Eigen sparse matrices internally; bridge to arma::sp_mat
  // only here, for the R-facing Matrix output.
  // Intermediate EM-blend evaluations never become the returned fit; skip
  // exporting their (possibly dense) coefficient matrix.
  const bool omitCoefficientMatrix = covStructI.hasAttribute("omitCoefficientMatrix") &&
    Rcpp::as<bool>(covStructI.attr("omitCoefficientMatrix"));
  const arma::sp_mat Wout = eigenSparseToArmaGlobal(W);
  const arma::sp_mat Cout = omitCoefficientMatrix
    ? arma::sp_mat(static_cast<arma::uword>(C.rows()), static_cast<arma::uword>(C.cols()))
    : eigenSparseToArmaGlobal(C);
  phaseProfile.finish();
  Rcpp::RObject profileOut = R_NilValue;
  if(phaseProfile.enabled()){
    const std::vector<sommer::PhaseProfile::Row> & profileRows = phaseProfile.iterations();
    if(selectedInverseDemandTrackingComplete){
      for(auto & demand : selectedInverseDemandKeys){
        const SelectedInverseSubset & subset = *demand.first;
        finalizeSelectedInverseDemand(subset, demand.second);
      }
    }
    auto selectedPatternEntries = [](const SelectedInverseSubset & subset) -> std::size_t {
      std::size_t entries = 0;
      for(const std::vector<int> & columnRows : subset.rows){ entries += columnRows.size(); }
      return entries;
    };
    Rcpp::CharacterVector phaseNames(sommer::kFitPhaseCount);
    Rcpp::NumericVector phaseTotals(sommer::kFitPhaseCount);
    Rcpp::NumericMatrix phaseByIteration(static_cast<int>(profileRows.size()),
      static_cast<int>(sommer::kFitPhaseCount));
    for(std::size_t phase = 0; phase < sommer::kFitPhaseCount; ++phase){
      phaseNames[phase] = sommer::fitPhaseName(phase);
      phaseTotals[phase] = phaseProfile.total()[phase];
      for(std::size_t row = 0; row < profileRows.size(); ++row){
        phaseByIteration(static_cast<int>(row), static_cast<int>(phase)) = profileRows[row][phase];
      }
    }
    phaseTotals.names() = phaseNames;
    Rcpp::colnames(phaseByIteration) = phaseNames;
    profileOut = Rcpp::List::create(
      Rcpp::Named("seconds") = phaseTotals,
      Rcpp::Named("secondsByIteration") = phaseByIteration,
      Rcpp::Named("counters") = Rcpp::List::create(
        Rcpp::Named("iterations") = static_cast<int>(profileRows.size()),
        Rcpp::Named("nEffects") = static_cast<int>(nEffects),
        Rcpp::Named("nnzC") = static_cast<double>(C.nonZeros()),
        Rcpp::Named("freeCovarianceParameters") = static_cast<int>(modelFreeCount),
        Rcpp::Named("exactAIRefreshes") = exactAIRefreshes,
        Rcpp::Named("hybridAIUpdates") = hybridAIUpdates,
        Rcpp::Named("sensitivityRHSColumns") = static_cast<double>(sensitivityRHSColumns),
        Rcpp::Named("sensitivityRHSColumnsSkipped") = static_cast<double>(sensitivityRHSColumnsSkipped),
        Rcpp::Named("sensitivityRHSNonzeros") = sensitivityRHSNonzeros,
        Rcpp::Named("sensitivitySolveSeconds") = sensitivitySolveSeconds,
        Rcpp::Named("sensitivitySolveReachEntries") = sensitivitySolveReachEntries,
        Rcpp::Named("sensitivitySolveReachPossibleEntries") =
          sensitivitySolveReachPossibleEntries,
        Rcpp::Named("traceDerivativePreparationSeconds") = traceDerivativePreparationSeconds,
        Rcpp::Named("traceRandomPrecisionDerivativeSeconds") = traceRandomPrecisionDerivativeSeconds,
        Rcpp::Named("traceFactorCovarianceDerivativeSeconds") = traceFactorCovarianceDerivativeSeconds,
        Rcpp::Named("traceFactorPrecisionDerivativeSeconds") = traceFactorPrecisionDerivativeSeconds,
        Rcpp::Named("traceRandomPrecisionOtherSeconds") = traceRandomPrecisionOtherSeconds,
        Rcpp::Named("traceResidualPrecisionDerivativeSeconds") = traceResidualPrecisionDerivativeSeconds,
        Rcpp::Named("traceGridDerivativeBasisSeconds") = traceGridDerivativeBasisSeconds,
        Rcpp::Named("traceResidualDerivativeBasisSeconds") = traceResidualDerivativeBasisSeconds,
        Rcpp::Named("traceDerivativeOtherSeconds") = traceDerivativeOtherSeconds,
        Rcpp::Named("traceSelectedInverseSeconds") = traceSelectedInverseSeconds,
        Rcpp::Named("selectedInverseLookupRequests") = selectedInverseLookupRequests,
        Rcpp::Named("selectedInverseLookupHits") = selectedInverseLookupHits,
        Rcpp::Named("selectedInverseLookupMisses") = selectedInverseLookupMisses,
        Rcpp::Named("selectedInverseDemandTrackingComplete") =
          selectedInverseDemandTrackingComplete,
        Rcpp::Named("selectedInverseDemandUniqueEntries") =
          static_cast<double>(selectedInverseDemandUniqueEntries),
        Rcpp::Named("selectedInverseDemandClosureEntries") =
          static_cast<double>(selectedInverseDemandClosureEntries),
        Rcpp::Named("selectedInverseDemandClosureOperations") =
          static_cast<double>(selectedInverseDemandClosureOperations),
        Rcpp::Named("cSelectedPatternEntries") =
          static_cast<double>(selectedPatternEntries(CselectedTopology)),
        Rcpp::Named("dSelectedPatternEntries") =
          static_cast<double>(selectedPatternEntries(DselectedTopology)),
        Rcpp::Named("scoreTraceFallbackSeconds") = scoreTraceFallbackSeconds,
        Rcpp::Named("aiAssemblyWallSeconds") = aiAssemblyWallSeconds,
        Rcpp::Named("aiRandomFirstTermSeconds") = aiRandomFirstTermSeconds,
        Rcpp::Named("aiRandomCurvatureSeconds") = aiRandomCurvatureSeconds,
        Rcpp::Named("aiResidualBlockSeconds") = aiResidualBlockSeconds,
        Rcpp::Named("solver") = solverName,
        Rcpp::Named("numericFactorCalls") = cholmodTiming.calls + ldltTiming.calls + blockSchurTiming.calls +
          static_cast<int>(saddleJoint.numericFactorizations() + saddleRandom.numericFactorizations()),
        Rcpp::Named("numericFactorSeconds") = cholmodTiming.elapsedSeconds +
          ldltTiming.elapsedSeconds + blockSchurTiming.elapsedSeconds + saddleFactorSeconds,
        Rcpp::Named("selectedTopologyReuses") = DselectedTopologyReuses,
        Rcpp::Named("selectedTopologyBuilds") = DselectedTopologyBuilds,
        Rcpp::Named("factorSchurFallbacks") = factorSchurFallbacks,
        Rcpp::Named("patternInverseCacheBuilds") = patternInverseCacheBuilds,
        Rcpp::Named("patternInverseCacheFallbacks") = patternInverseCacheFallbacks,
        Rcpp::Named("traceFallbackSolves") = traceFallbackSolves,
        Rcpp::Named("traceBatchSolves") = traceBatchSolves,
        Rcpp::Named("traceBatchWorkspaceBytes") = traceBatchWorkspaceBytes,
        Rcpp::Named("sparseFactorNonzeros") = static_cast<double>(blockEngine.sparseFactorNonzeros),
        Rcpp::Named("sparseFactorFlops") = blockEngine.sparseFactorFlops,
        Rcpp::Named("schurCostFlops") = blockEngine.schurCostFlops,
        Rcpp::Named("pcgBatchSolves") = pcgBatchSolves,
        Rcpp::Named("pcgBatchIterations") = pcgBatchIterations,
        Rcpp::Named("cSelectedDependencyEntries") = static_cast<double>(CselectedTopology.dependencyEntries),
        Rcpp::Named("dSelectedDependencyEntries") = static_cast<double>(DselectedTopology.dependencyEntries),
        Rcpp::Named("contractionDigitEntries") = static_cast<double>(factorEngine.contractionDigitsEntries),
        Rcpp::Named("maxIrregularMMEEntries") = static_cast<double>(maxIrregularMMEEntries)
      )
    );
  }
  Rcpp::List saddleDiagnostics;
  if(saddleActive){
    saddleDiagnostics = Rcpp::List::create(
        Rcpp::Named("active") = true,
        Rcpp::Named("backend") = "SparseLU",
        Rcpp::Named("traceMethod") = "exact-augmented-panels",
        Rcpp::Named("jointDimension") = nR + nEffects,
        Rcpp::Named("randomDimension") = !reml && Nu > 0 ? nR + Nu : 0,
        Rcpp::Named("jointNonzeros") = static_cast<double>(saddleJoint.nonzeros()),
        Rcpp::Named("jointFactorNonzeros") = static_cast<double>(saddleJoint.factorNonzeros()),
        Rcpp::Named("randomFactorNonzeros") = static_cast<double>(saddleRandom.factorNonzeros()),
        Rcpp::Named("symbolicAnalyses") = static_cast<double>(
          saddleJoint.symbolicAnalyses() + saddleRandom.symbolicAnalyses()),
        Rcpp::Named("numericFactorizations") = static_cast<double>(
          saddleJoint.numericFactorizations() + saddleRandom.numericFactorizations()),
        Rcpp::Named("solveColumns") = static_cast<double>(
          saddleJoint.solveColumns() + saddleRandom.solveColumns()),
        Rcpp::Named("maxRelativeResidual") = std::max(
          saddleJoint.maxRelativeResidual(), saddleRandom.maxRelativeResidual()),
        Rcpp::Named("factorSeconds") = saddleFactorSeconds,
        Rcpp::Named("normalEquationAssemblies") = normalEquationAssemblies,
        Rcpp::Named("compatibilityAssemblies") = saddleCompatibilityAssemblies,
        Rcpp::Named("compatibilityCMaterialized") = retainCoefficientMatrix);
  }
  auto timingList = [](const FactorizationTiming & timing, const char * measurement) -> Rcpp::List {
    return Rcpp::List::create(
      Rcpp::Named("measurement") = measurement,
      Rcpp::Named("calls") = timing.calls,
      Rcpp::Named("elapsedSeconds") = timing.elapsedSeconds,
      Rcpp::Named("cpuSeconds") = timing.cpuSeconds,
      Rcpp::Named("averageCPUs") = timing.elapsedSeconds > 0.0
        ? timing.cpuSeconds / timing.elapsedSeconds : NA_REAL);
  };
#ifdef _OPENMP
  const int openmpMaxThreads = omp_get_max_threads();
#else
  const int openmpMaxThreads = NA_INTEGER;
#endif
  const bool pcgAccuracyPassed = solverName == "pcg" &&
    std::isfinite(pcgLogDetProbeStandardError) && pcgLogDetProbeStandardError <= pcgLogDetTol &&
    std::isfinite(pcgLogDetDepthDifference) && pcgLogDetDepthDifference <= pcgLogDetTol &&
    (modelFreeCount == 0 || (std::isfinite(pcgScoreMaxSE) && pcgScoreMaxSE <= pcgScoreTol));
  if(solverName == "pcg" && pcgAccuracy == "adaptive" && !pcgAccuracyPassed){
    convergence = false;
    Rcpp::warning("PCG final empirical accuracy checks failed. Increase maxTraceProbes/maxLanczosSteps or inspect the reported log-determinant and score uncertainty; convergence is FALSE.");
  }
  Rcpp::List pcgDiagnostics = Rcpp::List::create();
  const Rcpp::List compressionDiagnostics = Rcpp::List::create(
    Rcpp::Named("method") = compressionMethod,
    Rcpp::Named("observations") = nR,
    Rcpp::Named("groups") = compressedData ? static_cast<double>(sufficient.groups.size()) : nR,
    Rcpp::Named("residualClasses") = compressedData ? residualDim : 0,
    Rcpp::Named("retainedBytes") = static_cast<double>(sufficient.retainedBytes),
    Rcpp::Named("setupSeconds") = sufficient.setupSeconds);
  if(solverName == "pcg"){
    pcgDiagnostics = Rcpp::List::create(
      Rcpp::Named("preconditioner") = pcgPreconditioner,
      Rcpp::Named("rhsMethod") = pcgRhsMethod,
      Rcpp::Named("operatorCalls") = static_cast<double>(pcgOperatorCalls.load()),
      Rcpp::Named("operatorVectorEquivalent") = static_cast<double>(pcgVectorEquivalent.load()),
      Rcpp::Named("matrixPanelApplications") = static_cast<double>(pcgMetrics.matrixPanelApplications),
      Rcpp::Named("preconditionerApplications") = static_cast<double>(pcgMetrics.preconditionerApplications),
      Rcpp::Named("trueResidualChecks") = static_cast<double>(pcgMetrics.trueResidualChecks),
      Rcpp::Named("residualReplacements") = static_cast<double>(pcgMetrics.residualReplacements),
      Rcpp::Named("deflatedColumns") = static_cast<double>(pcgMetrics.deflatedColumns),
      Rcpp::Named("lockedColumns") = static_cast<double>(pcgMetrics.lockedColumns),
      Rcpp::Named("blockFallbacks") = static_cast<double>(pcgMetrics.blockFallbacks),
      Rcpp::Named("blockFallbackReason") = pcgMetrics.fallbackReason,
      Rcpp::Named("warmStarts") = pcgWarmStarts,
      Rcpp::Named("coarseCorrections") = pcgCoarseCorrections,
      Rcpp::Named("compression") = compressionDiagnostics,
      Rcpp::Named("recycleMethod") = pcgRecycleMethod,
      Rcpp::Named("importedRecycleRank") = pcgImportedRank,
      Rcpp::Named("recycleRank") = static_cast<int>(krylovPcg ? pcgKrylov.basis().cols() : pcgRecycleBasis.cols()),
      Rcpp::Named("balancedPreconditionerApplications") = static_cast<double>(pcgBalancedApplications),
      Rcpp::Named("recycleSelectionOperatorApplications") = static_cast<double>(pcgKrylov.selectionOperatorApplications()),
      Rcpp::Named("recyclePreparationOperatorApplications") = static_cast<double>(pcgKrylov.preparationOperatorApplications()),
      Rcpp::Named("traceProbes") = pcgTraceProbes,
      Rcpp::Named("lanczosSteps") = pcgLanczosSteps,
      Rcpp::Named("logDetScaling") = pcgLogDetScaling,
      Rcpp::Named("logDetProbeEstimates") = pcgLogDetProbeValues,
      Rcpp::Named("logDetProbeSE") = pcgLogDetProbeStandardError,
      Rcpp::Named("logDetDepthDifference") = pcgLogDetDepthDifference,
      Rcpp::Named("pairedLogLikProbeSE") = pcgPairedLogLikSE,
      Rcpp::Named("scoreProbeEstimates") = arma::mat(pcgScoreProbeValues.data(),
        pcgScoreProbeValues.rows(), pcgScoreProbeValues.cols()),
      Rcpp::Named("scoreMaxSE") = pcgScoreMaxSE,
      Rcpp::Named("traceEstimator") = pcgTraceEstimator,
      Rcpp::Named("traceRank") = pcgTraceEstimator == "control-variate" ? std::min(pcgTraceRank, nEffects) : 0,
      Rcpp::Named("seed") = pcgSeed,
      Rcpp::Named("accuracy") = pcgAccuracy,
      Rcpp::Named("accuracyRefinements") = pcgAccuracyRefinements,
      Rcpp::Named("empiricalAccuracyPassed") = pcgAccuracy == "adaptive" && pcgAccuracyPassed,
      Rcpp::Named("seconds") = Rcpp::List::create(
        Rcpp::Named("preconditionerAndCoarseSetup") = pcgSetupSeconds,
        Rcpp::Named("solves") = pcgSolveSeconds,
        Rcpp::Named("logDeterminant") = pcgLogDetSeconds,
        Rcpp::Named("tracePreparationIncludingSolves") = pcgTraceSetupSeconds));
  }
  Rcpp::List result = Rcpp::List::create(
    Rcpp::Named("llik") = llik,
    // Rcpp::Named("M") = M,
    Rcpp::Named("W") = Wout,
    Rcpp::Named("C") = Cout,
    Rcpp::Named("Cscale") = vary,
    Rcpp::Named("yWorking") = yDense,
    Rcpp::Named("b") = b,
    Rcpp::Named("u") = u,
    Rcpp::Named("bu") = bu,
    Rcpp::Named("Ci") = Ci,
    // Only mode 2 (full explicit inverse) supports arbitrary D %*% Ci %*% t(D)
    // queries; mode 1 fills only the Takahashi selected-inverse pattern, same
    // as post_mme_Cinverse_cpp()'s CiComputed = (mode == 2) for consistency.
    Rcpp::Named("CiComputed") = (computeCi == 2),
    Rcpp::Named("CiMode") = computeCi,
    Rcpp::Named("solver") = solverName,
    Rcpp::Named("gpuDiagnostics") = Rcpp::List::create(
      Rcpp::Named("enabled") = gpuBackend.enabled(),
      Rcpp::Named("minimumSize") = gpuMinSize,
      Rcpp::Named("choleskyCalls") = static_cast<double>(gpuBackend.choleskyCalls),
      Rcpp::Named("triangularSolveCalls") = static_cast<double>(gpuBackend.triangularSolveCalls),
      Rcpp::Named("productCalls") = static_cast<double>(gpuBackend.productCalls)),
    Rcpp::Named("pcgMatrixFree") = matrixFreePCGReady,
    Rcpp::Named("pcgLogDetProbeStandardError") = pcgLogDetProbeStandardError,
    Rcpp::Named("pcgSolveMaxRelativeResidual") = pcgSolveMaxRelativeResidual,
    Rcpp::Named(".ldltCache") = retainLDLTCache ? static_cast<SEXP>(ldltCache) : R_NilValue,
    Rcpp::Named(".cholmodCache") = retainCholmodCache ? static_cast<SEXP>(cholmodCache) : R_NilValue,
    Rcpp::Named("factorizationDiagnostics") = Rcpp::List::create(
      Rcpp::Named("openmpMaxThreads") = openmpMaxThreads,
      Rcpp::Named("timings") = Rcpp::List::create(
        Rcpp::Named("cholmod") = timingList(cholmodTiming, "numeric factorization"),
        Rcpp::Named("ldlt") = timingList(ldltTiming, "numeric factorization"),
        Rcpp::Named("blockSchur") = timingList(blockSchurTiming,
          "block extraction, factorization, Schur assembly and log determinant"))),
    Rcpp::Named("engineDiagnostics") = Rcpp::List::create(
      Rcpp::Named("compression") = compressionDiagnostics,
      Rcpp::Named("coefficientEvaluations") = coefficientEvaluations,
      Rcpp::Named("emInverseTraceCalls") = emInverseTraceCalls,
      Rcpp::Named("diagonalResidualTraceCalls") = diagonalResidualTraceCalls,
      Rcpp::Named("diagonalResidualTraceSolveColumns") = diagonalResidualTraceSolveColumns,
      Rcpp::Named("emInverseTraceFallbackPanels") = emInverseTraceFallbackPanels,
      Rcpp::Named("emDesignTraceCalls") = emDesignTraceCalls,
      Rcpp::Named("nativeAISensitivityColumns") = static_cast<double>(sensitivityRHSColumns),
      Rcpp::Named("nativeAIAssemblyDimension") = static_cast<double>(aiDimension),
      Rcpp::Named("nativeAISensitivityEntries") =
        static_cast<double>(nEffects) * aiDimension,
      Rcpp::Named("nativeAIAssemblyEntries") = static_cast<double>(aiDimension) * aiDimension,
      Rcpp::Named("dFirstMLActive") = dFirstMLActive,
      Rcpp::Named("dFirstMLEvaluations") = dFirstMLEvaluations,
      Rcpp::Named("cholmodDSelectedBuilds") = cholmodDSelectedBuilds,
      Rcpp::Named("cholmodRSelectedBuilds") = cholmodRSelectedBuilds,
      Rcpp::Named("cholmodFinalSelectedInverse") = cholmodFinalSelectedInverse,
      Rcpp::Named("lineSearchHalvings") = totalLineSearchHalvings,
      Rcpp::Named("acceleratedProposals") = acceleratedProposals,
      Rcpp::Named("hybridAI") = Rcpp::List::create(
        Rcpp::Named("requested") = hybridRequested,
        Rcpp::Named("active") = useHybridAI,
        Rcpp::Named("fallbackReason") = hybridRequested && vcMap.active ? "constraint-map" : "",
        Rcpp::Named("refreshInterval") = hybridRefreshInterval,
        Rcpp::Named("exactRefreshes") = exactAIRefreshes,
        Rcpp::Named("bfgsUpdates") = hybridAIUpdates,
        Rcpp::Named("resets") = hybridAIResets,
        Rcpp::Named("certifications") = hybridAICertifications,
        Rcpp::Named("sensitivityRHSColumns") = static_cast<double>(sensitivityRHSColumns),
        Rcpp::Named("sensitivityRHSColumnsSkipped") = static_cast<double>(sensitivityRHSColumnsSkipped),
        Rcpp::Named("finalInformationExact") = finalAIExact && retainCovarianceInformation,
        Rcpp::Named("finalProjectedScoreNorm") = finalProjectedScoreNorm),
      Rcpp::Named("CsymbolicAnalyses") = CsymbolicAnalyses,
      Rcpp::Named("DsymbolicAnalyses") = DsymbolicAnalyses,
      Rcpp::Named("DselectedTopologyBuilds") = DselectedTopologyBuilds,
      Rcpp::Named("DselectedTopologyReuses") = DselectedTopologyReuses,
      Rcpp::Named("CselectedDependencyEntries") = static_cast<double>(CselectedTopology.dependencyEntries),
      Rcpp::Named("DselectedDependencyEntries") = static_cast<double>(DselectedTopology.dependencyEntries),
      Rcpp::Named("contractionDigitEntries") = static_cast<double>(factorEngine.contractionDigitsEntries),
      Rcpp::Named("factorEvaluationReuses") = factorEngine.factorEvaluationReuses,
      Rcpp::Named("sparseSpatialFactorizationReuses") = factorEngine.sparseSpatialFactorizationReuses,
      Rcpp::Named("pdCertificates") = pdCertificates,
      Rcpp::Named("blockSchurActive") = blockEngine.active,
      Rcpp::Named("blockChainActive") = blockEngine.active && blockEngine.chain,
      Rcpp::Named("blockBandedActive") = blockEngine.active && blockEngine.banded,
      Rcpp::Named("blockBandedWidth") = blockEngine.active && blockEngine.banded
        ? static_cast<int>(blockEngine.bandWidth) : 0,
      Rcpp::Named("factorSchurActive") = blockEngine.active && blockEngine.latent,
      Rcpp::Named("factorSchurFallbacks") = factorSchurFallbacks,
      Rcpp::Named("patternInverseCacheBuilds") = patternInverseCacheBuilds,
      Rcpp::Named("patternInverseCacheFallbacks") = patternInverseCacheFallbacks,
      Rcpp::Named("factorSchurPlanning") = Rcpp::List::create(
        Rcpp::Named("admission") = blockEngine.latentAdmission,
        Rcpp::Named("numericalFallbacks") = factorSchurNumericalFallbacks,
        Rcpp::Named("latentWorkFlops") = blockEngine.latentWorkFlops,
        Rcpp::Named("marginalWorkFlops") = blockEngine.marginalWorkFlops,
        Rcpp::Named("requiredWorkRatio") = 0.8),
      Rcpp::Named("residualEngine") = residualPlan.grid() ? "grid"
        : residualPlan.diagonal() ? "diagonal" : residualPlan.kronecker() ? "kronecker" : "sparse",
      Rcpp::Named("gridFactorizations") = gridFactorizations,
      Rcpp::Named("gridFactorReuses") = gridFactorReuses,
      Rcpp::Named("pcgResidualTraceContractions") = pcgResidualTraceContractions,
      Rcpp::Named("blockSchurGroups") = blockEngine.active ? static_cast<int>(blockEngine.groups.size()) : 0,
      Rcpp::Named("blockSchurBorder") = blockEngine.active ? static_cast<int>(blockEngine.borderCols.size()) : 0,
      Rcpp::Named("denseMemoryMB") = blockSchurPolicy.maxDoubles * sizeof(double) / 1.0e6,
      Rcpp::Named("blockSchurEstimatedDenseMB") = blockEngine.active
        ? (blockEngine.chain || blockEngine.banded || blockEngine.latent ? blockEngine.chainStorageDoubles
            : blockSchurPolicy.storageDoubles(blockEngine.groups, blockEngine.borderCols.size())) * sizeof(double) / 1.0e6
        : 0.0,
      Rcpp::Named("control") = Rcpp::List::create(
        Rcpp::Named("solver") = solverName,
        Rcpp::Named("denseMemoryMB") = blockSchurPolicy.maxDoubles * sizeof(double) / 1.0e6,
        Rcpp::Named("blockBanded") = blockBandedMode,
        Rcpp::Named("factorSchurMultiFactor") = factorSchurMultiFactorMode,
        Rcpp::Named("dFirstML") = dFirstMLMode,
        Rcpp::Named("matrixFreePcg") = matrixFreePcgMode,
        Rcpp::Named("cholmodFinalSelected") = cholmodFinalSelectedMode,
        Rcpp::Named("pdCertificates") = pdCertificatesMode,
        Rcpp::Named("panelCi2") = panelCi2Batch,
        Rcpp::Named("profile") = profileMode,
        Rcpp::Named("blockSchur") = blockSchurMode,
        Rcpp::Named("traceBatch") = traceBatchMode,
        Rcpp::Named("sparseFactors") = sommer::control::modeName(
          parseExpMode("sparseFactors", "sommer.mme.sparseFactors", sommer::control::Mode::Off)),
        Rcpp::Named("memoryPolicy") = memoryPolicyMode
      ),
      Rcpp::Named("profile") = profileOut,
      Rcpp::Named("weightBlockCount") = weightBlockCount,
      Rcpp::Named("pcgBatchSolves") = pcgBatchSolves,
      Rcpp::Named("pcgBatchIterations") = pcgBatchIterations,
      Rcpp::Named("pcgPreconditioner") = pcgPreconditioner,
      Rcpp::Named("pcgNystromRank") = pcgNystromRankUsed,
      Rcpp::Named("maxIrregularDesignEntries") = static_cast<double>(maxIrregularDesignEntries),
      Rcpp::Named("maxIrregularMMEEntries") = static_cast<double>(maxIrregularMMEEntries)
    ),
    Rcpp::Named("residualCovarianceCompact") = compactResidualCovariance,
    Rcpp::Named("pcgTol") = pcgTol,
    Rcpp::Named("pcgMaxIters") = pcgMaxIters,
    Rcpp::Named("theta") = theta,
    Rcpp::Named("covPar") = covParOut,
    Rcpp::Named("covParWorking") = covParWorking,
    Rcpp::Named("covType") = Rcpp::wrap(covType),
    Rcpp::Named("theta_se") = InfMatInv,
    Rcpp::Named("InfMat") = InfMat, //InfMat,
    Rcpp::Named("monitor") = monitorOut,
    // Rcpp::Named("constraints") = thetaCUnlistedFinal,
    Rcpp::Named("uList") = uList,
    Rcpp::Named("uPevList") = uPevList,
    Rcpp::Named("AIC") = AIC,
    Rcpp::Named("BIC") = BIC,
    Rcpp::Named("convergence") = convergence,
    Rcpp::Named("partitions") = partitions,
    Rcpp::Named("percDelta") = percDelta,
    Rcpp::Named("normMonitor") = normMonitor,
    Rcpp::Named("toBoundary") = toBoundary,
    Rcpp::Named("dLu") = dLuOut,
    Rcpp::Named("Cchol") = arma::mat()
    // Rcpp::Named("PMWu_mat") = PMWu_mat, 
    // Rcpp::Named("MWuchol") = MWuchol
  );
  if(solverName == "pcg"){
    Rcpp::List diagnostics = result["engineDiagnostics"];
    diagnostics["pcg"] = pcgDiagnostics;
    result["engineDiagnostics"] = diagnostics;
    if(pcgExportRecycle){
      const Eigen::MatrixXd & exportedBasis =
        pcgEvaluationAccepted && pcgPendingKrylovReady
          ? pcgPendingKrylov.basis : pcgKrylov.basis();
      Rcpp::NumericMatrix exported(exportedBasis.rows(), exportedBasis.cols());
      if(exportedBasis.size() > 0){
        std::copy(exportedBasis.data(), exportedBasis.data() + exportedBasis.size(),
          exported.begin());
      }
      result["pcgRecycleBasis"] = exported;
    }
  }
  if(saddleActive){
    Rcpp::List diagnostics = result["engineDiagnostics"];
    diagnostics["saddlePoint"] = saddleDiagnostics;
    result["engineDiagnostics"] = diagnostics;
  }
  if(!retainCovarianceInformation){
    if(!emInformationOnly){ result["InfMat"] = R_NilValue; }
    result["theta_se"] = R_NilValue;
    Rcpp::List diagnostics = result["engineDiagnostics"];
    diagnostics["covarianceInformation"] = covarianceInformation;
    result["engineDiagnostics"] = diagnostics;
  }
  if(reducedAIAssembly){ result[".aiTangent"] = aiTangent; }
  if(profileCommonScale){
    result[".commonScale"] = commonScale;
    result[".commonScaleDF"] = commonScaleDF;
    result[".commonScaleQuadratic"] = commonScaleQuadraticTotal;
    result[".commonScaleQuadraticVector"] = commonScaleQuadratic;
  }
  if(!retainCoefficientMatrix){
    result["C"] = R_NilValue;
    result["retainCoefficientMatrix"] = false;
  }
  return result;
  
}

// [[Rcpp::export]]
Rcpp::List ai_mme_sp2(const arma::sp_mat & X, const Rcpp::List & ZI,
                     const arma::vec & Zind, const Rcpp::List & AiI,
                     const arma::sp_mat & y0, const arma::sp_mat & H,
                     const bool & useH, const arma::uvec & residualBlockI,
                     const arma::uvec & residualIndexI, int nIters,
                     double tolParConvLL, double tolParConvNorm, double tolParInv,
                     const Rcpp::List & covStructI, const arma::vec & weightEmInf,
                     const arma::vec & weightInf, const bool & verbose,
                     const int & computeCi = 0, const std::string & solver = "ldlt",
                     const double & pcgTol = 1.0e-8, const int & pcgMaxIters = 0,
                     const int & pcgTraceProbes = 8, const int & pcgLanczosSteps = 20,
                     const bool & reml = true, const bool & responsePrepared = false,
                     const double & preparedMean = 0.0, const double & preparedSd = 1.0,
                     const bool & preparedIntercept = false, const bool gpu = false,
                     const int gpuMinSize = 256) {
  return ai_mme_sp2_impl(X, ZI, Zind, AiI, y0, H, useH, residualBlockI,
    residualIndexI, nIters, tolParConvLL, tolParConvNorm, tolParInv, covStructI,
    weightEmInf, weightInf, verbose, computeCi, solver, pcgTol, pcgMaxIters,
    pcgTraceProbes, pcgLanczosSteps, reml, responsePrepared, preparedMean,
    preparedSd, preparedIntercept, gpu, gpuMinSize);
}
