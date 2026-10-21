#ifndef SOMMER_MME_COVARIANCE_H
#define SOMMER_MME_COVARIANCE_H

#include "mme_factors.h"
#include <chrono>
#include <mutex>

namespace sommer {

inline arma::vec covarianceVarianceLowerBounds(
    const Rcpp::List & descriptors, double inverseTolerance, double logResponseVariance = 0.0) {
  std::vector<double> bounds;
  const double floor = std::log(std::max(1.0e-8, inverseTolerance));
  for(int term = 0; term < descriptors.size(); ++term) {
    const Rcpp::List descriptor = descriptors[term];
    const arma::vec parameters = Rcpp::as<arma::vec>(descriptor["par"]);
    const Rcpp::LogicalVector free = descriptor["free"];
    const arma::uword offset = bounds.size();
    bounds.resize(offset + parameters.n_elem, -std::numeric_limits<double>::infinity());
    if(parameters.n_elem == 0) { continue; }
    if(free[0] == TRUE) { bounds[offset] = floor - logResponseVariance; }
    const Rcpp::List factors = descriptor["factors"];
    for(int axis = 0; axis < factors.size(); ++axis) {
      const Rcpp::List factor = factors[axis];
      const int first = Rcpp::as<int>(factor["par_start"]) - 1;
      const Rcpp::CharacterVector names = factor["par_names"];
      const Rcpp::List report = factor["report"];
      if(Rcpp::as<std::string>(report["backend"]) != "builtin") { continue; }
      const Rcpp::CharacterVector transforms = report["transform"];
      for(int j = 0; j < names.size(); ++j) {
        const std::string name = Rcpp::as<std::string>(names[j]);
        const bool variance = name.compare(0, 15, "variance_ratio[") == 0 ||
          name.compare(0, 15, "specific_ratio[") == 0 ||
          name.compare(0, 17, "innovation_ratio[") == 0;
        if(variance && Rcpp::as<std::string>(transforms[j]) == "exp" && free[first + j] == TRUE) {
          bounds[offset + first + j] = floor;
        }
      }
    }
  }
  return arma::vec(bounds);
}

inline arma::vec covarianceWorkingStepCaps(const Rcpp::List & descriptors) {
  std::vector<double> values;
  for(int term = 0; term < descriptors.size(); ++term) {
    const Rcpp::List descriptor = descriptors[term];
    const arma::uword offset = values.size();
    const arma::uword count = Rcpp::as<arma::vec>(descriptor["par"]).n_elem;
    if(count == 0) { continue; }
    values.resize(offset + count, 2.0);
    values[offset] = 1.0;
    const Rcpp::List factors = descriptor["factors"];
    for(int axis = 0; axis < factors.size(); ++axis) {
      const Rcpp::List factor = factors[axis];
      const int first = Rcpp::as<int>(factor["par_start"]);
      const int last = Rcpp::as<int>(factor["par_end"]);
      if(last < first) { continue; }
      const arma::vec caps = Rcpp::as<arma::vec>(factor["trust_cap"]);
      if(caps.n_elem != static_cast<arma::uword>(last - first + 1)) {
        Rcpp::stop("CovarianceFactor trust_cap has incompatible length.");
      }
      for(arma::uword j = 0; j < caps.n_elem; ++j) {
        if(!std::isfinite(caps(j)) || caps(j) <= 0.0) {
          Rcpp::stop("CovarianceFactor trust caps must be positive and finite.");
        }
        values[offset + first - 1 + j] = caps(j);
      }
    }
  }
  return arma::vec(values);
}

class CovarianceFactorEngine {
public:
  typedef FactorPrecisionState DescriptorPrecisionState;

  struct FactorEvaluationCache {
    Rcpp::List owner;
    arma::vec parameters;
    arma::mat covariance;
    std::shared_ptr<SparseSpatialFactor> sparseSpatial;
  };
  std::map<SEXP, FactorEvaluationCache, std::less<SEXP>> factorEvaluationCache;
  std::size_t factorEvaluationEntries = 0;
  std::size_t factorEvaluationMaxEntries = WorkspacePolicy::cacheEntries();
  int factorEvaluationReuses = 0;
  int sparseSpatialFactorizationReuses = 0;

  std::map<std::vector<arma::uword>, arma::umat> contractionDigitsCache;
  std::mutex contractionDigitsMutex;
  std::size_t contractionDigitsEntries = 0;
  std::size_t contractionDigitsMaxEntries = WorkspacePolicy::cacheEntries();

  arma::mat evalNativeFactor(const Rcpp::List & f,
        const arma::vec & localPar,
        const std::string & op) {

      const arma::uword q =
        static_cast<arma::uword>(
          Rcpp::as<int>(f["dim"])
        );

      if(op == "identity"){
        return arma::eye<arma::mat>(q,q);
      }

      if(op == "diag"){
        if(localPar.n_elem + 1 != q){
          Rcpp::stop("Malformed diagonal covariance factor.");
        }
        arma::vec d(q, arma::fill::ones);
        for(arma::uword k = 0; k < localPar.n_elem; ++k){
          d(k+1) = std::exp(localPar(k));
        }
        return arma::diagmat(d);
      }

      if(op == "ar1"){
        if(localPar.n_elem != 1){
          Rcpp::stop("Malformed AR1 covariance factor.");
        }
        const double rho = std::tanh(localPar(0));
        arma::mat K(q,q,arma::fill::zeros);
        for(arma::uword i = 0; i < q; ++i){
          for(arma::uword j = 0; j < q; ++j){
            const arma::uword d = (i > j ? i-j : j-i);
            K(i,j) = std::pow(rho, static_cast<double>(d));
          }
        }
        return K;
      }

      if(op == "us"){
        Rcpp::IntegerVector rr = f["us_row"];
        Rcpp::IntegerVector cc = f["us_col"];
        Rcpp::LogicalVector dd = f["us_diag"];

        if(localPar.n_elem != static_cast<arma::uword>(rr.size()) ||
           rr.size() != cc.size() || rr.size() != dd.size()){
          Rcpp::stop("Malformed unstructured covariance factor.");
        }

        arma::mat L(q,q,arma::fill::zeros);
        L(0,0) = 1.0;

        for(arma::uword k = 0; k < localPar.n_elem; ++k){
          const arma::uword i =
            static_cast<arma::uword>(rr[static_cast<int>(k)] - 1);
          const arma::uword j =
            static_cast<arma::uword>(cc[static_cast<int>(k)] - 1);
          if(dd[static_cast<int>(k)]){
            L(i,j) = std::exp(localPar(k));
          }else{
            L(i,j) = localPar(k);
          }
        }

        return L * L.t();
      }


      if(op == "cor_uniform"){
        if(q < 2 || localPar.n_elem != 1){
          Rcpp::stop("Malformed compound-symmetry/uniform-correlation factor.");
        }

        const double lo =
          -1.0
          /
          static_cast<double>(q - 1);

        const double eta = localPar(0);

        const double s =
          eta >= 0.0
          ?
          1.0 / (1.0 + std::exp(-eta))
          :
          std::exp(eta) / (1.0 + std::exp(eta));

        const double rho =
          lo
          +
          (1.0 - lo) * s;

        arma::mat K(
          q,
          q,
          arma::fill::value(rho)
        );

        K.diag().ones();

        return K;
      }

      if(op == "corh"){
        if(q < 2 || localPar.n_elem != q){
          Rcpp::stop("Malformed heterogeneous uniform-correlation factor.");
        }

        const double lo =
          -1.0
          /
          static_cast<double>(q - 1);

        const double eta = localPar(0);

        const double s =
          eta >= 0.0
          ?
          1.0 / (1.0 + std::exp(-eta))
          :
          std::exp(eta) / (1.0 + std::exp(eta));

        const double rho =
          lo
          +
          (1.0 - lo) * s;

        arma::mat C(
          q,
          q,
          arma::fill::value(rho)
        );

        C.diag().ones();

        arma::vec variances(
          q,
          arma::fill::ones
        );

        for(arma::uword k = 1; k < q; ++k){
          variances(k) =
            std::exp(
              localPar(k)
            );
        }

        arma::vec sd =
          arma::sqrt(
            variances
          );

        return
          arma::diagmat(sd)
          *
          C
          *
          arma::diagmat(sd);
      }

      if(op == "arp"){
        const int order =
          Rcpp::as<int>(
            f["order"]
          );

        if(order < 1 ||
           localPar.n_elem != static_cast<arma::uword>(order) ||
           q <= static_cast<arma::uword>(order)){
          Rcpp::stop("Malformed AR(p) covariance factor.");
        }

        // Reflection/PACF coefficients -> stable AR coefficients.
        arma::vec pacf(
          static_cast<arma::uword>(order),
          arma::fill::zeros
        );

        for(int j = 0; j < order; ++j){
          pacf(static_cast<arma::uword>(j)) =
            std::tanh(
              localPar(static_cast<arma::uword>(j))
            );
        }

        arma::vec phi;

        for(int m = 1; m <= order; ++m){

          arma::vec next(
            static_cast<arma::uword>(m),
            arma::fill::zeros
          );

          next(static_cast<arma::uword>(m-1)) =
            pacf(static_cast<arma::uword>(m-1));

          if(m > 1){

            for(int j = 0; j < m-1; ++j){

              next(static_cast<arma::uword>(j)) =
                phi(static_cast<arma::uword>(j))
                -
                pacf(static_cast<arma::uword>(m-1))
                *
                phi(static_cast<arma::uword>(m-2-j));
            }
          }

          phi = next;
        }

        // Solve Yule-Walker equations for rho_1,...,rho_p with rho_0=1.
        arma::mat A(
          static_cast<arma::uword>(order),
          static_cast<arma::uword>(order),
          arma::fill::zeros
        );

        arma::vec b(
          static_cast<arma::uword>(order),
          arma::fill::zeros
        );

        for(int kk = 1; kk <= order; ++kk){

          A(
            static_cast<arma::uword>(kk-1),
            static_cast<arma::uword>(kk-1)
          ) += 1.0;

          for(int jj = 1; jj <= order; ++jj){

            const int d =
              std::abs(
                kk - jj
              );

            const double pj =
              phi(
                static_cast<arma::uword>(jj-1)
              );

            if(d == 0){

              b(
                static_cast<arma::uword>(kk-1)
              ) += pj;

            }else{

              A(
                static_cast<arma::uword>(kk-1),
                static_cast<arma::uword>(d-1)
              ) -= pj;
            }
          }
        }

        arma::vec rInitial;

        bool ok =
          arma::solve(
            rInitial,
            A,
            b
          );

        if(!ok || !rInitial.is_finite()){
          Rcpp::stop("Unable to solve Yule-Walker equations for AR(p) factor.");
        }

        arma::vec rho(
          q,
          arma::fill::zeros
        );

        rho(0) = 1.0;

        for(int h = 1; h <= order; ++h){
          rho(static_cast<arma::uword>(h)) =
            rInitial(static_cast<arma::uword>(h-1));
        }

        for(arma::uword h = static_cast<arma::uword>(order+1);
            h < q;
            ++h){

          double value = 0.0;

          for(int jj = 1; jj <= order; ++jj){

            value +=
              phi(static_cast<arma::uword>(jj-1))
              *
              rho(
                h
                -
                static_cast<arma::uword>(jj)
              );
          }

          rho(h) = value;
        }

        arma::mat K(
          q,
          q,
          arma::fill::zeros
        );

        for(arma::uword i = 0; i < q; ++i){
          for(arma::uword j = 0; j < q; ++j){
            const arma::uword d =
              i > j
              ?
              i-j
              :
              j-i;

            K(i,j) =
              rho(d);
          }
        }

        return K;
      }

      if(op == "ma"){
        const int order =
          Rcpp::as<int>(
            f["order"]
          );

        if(order < 1 ||
           localPar.n_elem != static_cast<arma::uword>(order) ||
           q <= static_cast<arma::uword>(order)){
          Rcpp::stop("Malformed MA(q) covariance factor.");
        }

        arma::vec coef(
          static_cast<arma::uword>(order + 1),
          arma::fill::zeros
        );

        coef(0) = 1.0;

        for(int j = 1; j <= order; ++j){
          coef(static_cast<arma::uword>(j)) =
            localPar(static_cast<arma::uword>(j-1));
        }

        arma::vec rho(
          q,
          arma::fill::zeros
        );

        double gamma0 = 0.0;

        for(int j = 0; j <= order; ++j){
          const double c =
            coef(static_cast<arma::uword>(j));
          gamma0 += c*c;
        }

        if(!std::isfinite(gamma0) || gamma0 <= 0.0){
          Rcpp::stop("Invalid MA covariance normalization.");
        }

        rho(0) = 1.0;

        for(int h = 1; h <= order; ++h){

          double gamma = 0.0;

          for(int j = 0; j <= order-h; ++j){

            gamma +=
              coef(static_cast<arma::uword>(j))
              *
              coef(static_cast<arma::uword>(j+h));
          }

          rho(static_cast<arma::uword>(h)) =
            gamma / gamma0;
        }

        arma::mat K(
          q,
          q,
          arma::fill::zeros
        );

        for(arma::uword i = 0; i < q; ++i){
          for(arma::uword j = 0; j < q; ++j){

            const arma::uword d =
              i > j
              ?
              i-j
              :
              j-i;

            K(i,j) =
              d <= static_cast<arma::uword>(order)
              ?
              rho(d)
              :
              0.0;
          }
        }

        return K;
      }

      if(op == "corg"){
        Rcpp::IntegerVector rr =
          f["corg_row"];

        Rcpp::IntegerVector cc =
          f["corg_col"];

        if(localPar.n_elem != static_cast<arma::uword>(rr.size()) ||
           rr.size() != cc.size()){
          Rcpp::stop("Malformed general-correlation factor.");
        }

        arma::mat A(
          q,
          q,
          arma::fill::eye
        );

        for(arma::uword k = 0; k < localPar.n_elem; ++k){

          const arma::uword i =
            static_cast<arma::uword>(
              rr[static_cast<int>(k)] - 1
            );

          const arma::uword j =
            static_cast<arma::uword>(
              cc[static_cast<int>(k)] - 1
            );

          A(i,j) =
            localPar(k);
        }

        arma::mat S =
          A * A.t();

        arma::vec sd =
          arma::sqrt(
            S.diag()
          );

        arma::mat denom = sd * sd.t();
        arma::mat K = S % arma::pow(denom, -1.0);

        K.diag().ones();

        return
          0.5
          *
          (
            K
            +
            K.t()
          );
      }


      if(op == "fa" || op == "rr"){
        return sommer::loadingFactorCovariance(f, localPar);
      }

      if(op == "ante"){
        const int ncoef =
          Rcpp::as<int>(
            f["ante_ncoef"]
          );

        Rcpp::IntegerVector rr =
          f["ante_row"];

        Rcpp::IntegerVector cc =
          f["ante_col"];

        if(ncoef < 1 ||
           rr.size() != ncoef ||
           cc.size() != ncoef ||
           localPar.n_elem !=
             static_cast<arma::uword>(
               ncoef
               +
               static_cast<int>(q)
               -
               1
             )){
          Rcpp::stop("Malformed antedependence covariance factor.");
        }

        return anteInnovationCovariance(q, rr, cc, localPar, sparseFactorKernelsEnabled(f));
      }


      Rcpp::stop("Unsupported native covariance evaluator opcode: " + op);
      return arma::mat();
    }

  std::size_t factorCacheEntries(const FactorEvaluationCache & entry) {
    if(entry.sparseSpatial){ return entry.sparseSpatial->covariance.n_elem; }
    return entry.covariance.n_elem;
  }

  SparseSpatialFactor evaluateSparseSpatialCached(const Rcpp::List & f,
                                         const arma::vec & localPar,
                                         const bool derivatives) {
    const SEXP cacheKey = static_cast<SEXP>(f);
    auto cached = factorEvaluationCache.find(cacheKey);
    const bool sameParameters = cached != factorEvaluationCache.end() &&
      cached->second.parameters.n_elem == localPar.n_elem &&
      arma::all(cached->second.parameters == localPar);
    if(sameParameters && cached->second.sparseSpatial){
      ++factorEvaluationReuses;
      ++sparseSpatialFactorizationReuses;
      return evaluateSparseSpatialFactor(f, localPar, derivatives, cached->second.sparseSpatial.get());
    }
    SparseSpatialFactor spatial = evaluateSparseSpatialFactor(f, localPar, false);
    if(cached != factorEvaluationCache.end()){
      factorEvaluationEntries -= factorCacheEntries(cached->second);
      factorEvaluationCache.erase(cached);
    }
    if(WorkspacePolicy::cacheFits(factorEvaluationMaxEntries, factorEvaluationEntries,
                                spatial.covariance.n_elem)){
      FactorEvaluationCache entry;
      entry.owner = f;
      entry.parameters = localPar;
      entry.sparseSpatial = std::make_shared<SparseSpatialFactor>(spatial);
      factorEvaluationEntries += entry.sparseSpatial->covariance.n_elem;
      factorEvaluationCache.emplace(cacheKey, std::move(entry));
    }
    if(derivatives){ return evaluateSparseSpatialFactor(f, localPar, true, &spatial); }
    return spatial;
  }

  // Universal factor evaluator. Backends are deliberately small and stable:
  //   native : optimized built-in primitive selected by evaluator$op
  //   fixed  : a supplied fixed covariance matrix
  //   R      : arbitrary callback fun(par) -> covariance matrix
  // Adding a future covariance model therefore does not require changing
  // ai_mme_sp2(): it can be expressed through the R backend immediately and
  // promoted to a native opcode later only if performance justifies it.
  arma::mat evalFactor(const Rcpp::List & f,
        const arma::vec & localPar) {

      if(!f.containsElementNamed("evaluator")){
        Rcpp::stop("CovarianceFactor is missing evaluator specification.");
      }

      Rcpp::List spec =
        Rcpp::as<Rcpp::List>(f["evaluator"]);

      const std::string backend =
        Rcpp::as<std::string>(spec["backend"]);

      const arma::uword q =
        static_cast<arma::uword>(Rcpp::as<int>(f["dim"]));

      const bool sparseSpatial = sparseFactorKernelsEnabled(f) &&
        sommer::FactorPolicy::spatial(sommer::precisionFamily(f));
      const bool cacheable = backend == "native" || backend == "fixed" ||
        sparseSpatial || (spec.containsElementNamed("cacheable") && Rcpp::as<bool>(spec["cacheable"]));
      const SEXP cacheKey = static_cast<SEXP>(f);
      auto cached = factorEvaluationCache.find(cacheKey);
      const bool cachedCovarianceAvailable = cached != factorEvaluationCache.end() &&
        (cached->second.sparseSpatial
          ? cached->second.sparseSpatial->covariance.n_rows == q
          : cached->second.covariance.n_rows == q);
      if(cacheable && cachedCovarianceAvailable && cached->second.parameters.n_elem == localPar.n_elem &&
         arma::all(cached->second.parameters == localPar)){
        ++factorEvaluationReuses;
        return cached->second.sparseSpatial
          ? cached->second.sparseSpatial->covariance : cached->second.covariance;
      }
      arma::mat K;

      if(sparseSpatial){
        K = evaluateSparseSpatialCached(f, localPar, false).covariance;
      }else if(backend == "native"){
        if(!spec.containsElementNamed("op")){
          Rcpp::stop("Native CovarianceFactor evaluator is missing op.");
        }
        const std::string op =
          Rcpp::as<std::string>(spec["op"]);
        K = evalNativeFactor(f, localPar, op);

      }else if(backend == "fixed"){
        if(!spec.containsElementNamed("matrix")){
          Rcpp::stop("Fixed CovarianceFactor evaluator is missing matrix.");
        }
        K = Rcpp::as<arma::mat>(spec["matrix"]);

      }else if(backend == "R"){
        if(!spec.containsElementNamed("fun") || Rf_isNull(spec["fun"])){
          Rcpp::stop("R CovarianceFactor evaluator is missing fun(par).");
        }
        Rcpp::Function fun = spec["fun"];
        SEXP ans = fun(Rcpp::wrap(localPar));
        K = Rcpp::as<arma::mat>(ans);

      }else{
        Rcpp::stop("Unknown CovarianceFactor evaluator backend: " + backend);
      }

      if(K.n_rows != q || K.n_cols != q || !K.is_finite()){
        Rcpp::stop("CovarianceFactor evaluator returned an invalid matrix.");
      }

      K = 0.5 * (K + K.t());
      if(cacheable && !sparseSpatial){
        if(cached != factorEvaluationCache.end()){
          factorEvaluationEntries -= factorCacheEntries(cached->second);
          factorEvaluationCache.erase(cached);
        }
        if(WorkspacePolicy::cacheFits(factorEvaluationMaxEntries, factorEvaluationEntries, K.n_elem)){
          FactorEvaluationCache entry;
          entry.owner = f;
          entry.parameters = localPar;
          entry.covariance = K;
          factorEvaluationEntries += K.n_elem;
          factorEvaluationCache.emplace(cacheKey, std::move(entry));
        }
      }
      return K;
    }

  // Universal first-derivative evaluator.  The AI algorithm only needs first
  // covariance derivatives.  Native analytic derivatives are used for cheap
  // high-frequency primitives; arbitrary callbacks may provide dfun(par,k);
  // otherwise the covariance engine differentiates the small factor matrix by
  // a central finite difference.  The likelihood itself is never numerically
  // differentiated.
  arma::mat factorD1(const Rcpp::List & f,
        const arma::vec & localPar,
        const arma::uword k) {

      if(k >= localPar.n_elem){
        Rcpp::stop("Invalid CovarianceFactor derivative parameter index.");
      }
      if(!f.containsElementNamed("derivative")){
        Rcpp::stop("CovarianceFactor is missing derivative specification.");
      }
      if(sparseFactorKernelsEnabled(f) && sommer::FactorPolicy::spatial(sommer::precisionFamily(f))){
        return evaluateSparseSpatialCached(f, localPar, true).covarianceD1;
      }

      Rcpp::List spec =
        Rcpp::as<Rcpp::List>(f["derivative"]);
      const std::string backend =
        Rcpp::as<std::string>(spec["backend"]);
      const arma::uword q =
        static_cast<arma::uword>(Rcpp::as<int>(f["dim"]));

      if(backend == "native"){
        if(!spec.containsElementNamed("op")){
          Rcpp::stop("Native derivative specification is missing op.");
        }
        const std::string op =
          Rcpp::as<std::string>(spec["op"]);

        if(op == "diag"){
          arma::mat D(q,q,arma::fill::zeros);
          D(k+1,k+1) = std::exp(localPar(k));
          return D;
        }

        if(op == "ar1"){
          if(k != 0 || localPar.n_elem != 1){
            Rcpp::stop("Malformed native AR1 derivative request.");
          }
          const double rho = std::tanh(localPar(0));
          const double drho = 1.0 - rho*rho;
          arma::mat D(q,q,arma::fill::zeros);
          for(arma::uword i = 0; i < q; ++i){
            for(arma::uword j = 0; j < q; ++j){
              const arma::uword d = (i > j ? i-j : j-i);
              if(d > 0){
                D(i,j) =
                  static_cast<double>(d)
                  * std::pow(rho, static_cast<double>(d-1))
                  * drho;
              }
            }
          }
          return D;
        }

        if(op == "us"){
          Rcpp::IntegerVector rr = f["us_row"];
          Rcpp::IntegerVector cc = f["us_col"];
          Rcpp::LogicalVector dd = f["us_diag"];
          if(localPar.n_elem != static_cast<arma::uword>(rr.size()) ||
             rr.size() != cc.size() || rr.size() != dd.size()){
            Rcpp::stop("Malformed native US derivative metadata.");
          }
          arma::mat L(q,q,arma::fill::zeros);
          L(0,0) = 1.0;
          for(arma::uword a = 0; a < localPar.n_elem; ++a){
            const arma::uword i =
              static_cast<arma::uword>(rr[static_cast<int>(a)] - 1);
            const arma::uword j =
              static_cast<arma::uword>(cc[static_cast<int>(a)] - 1);
            L(i,j) = dd[static_cast<int>(a)]
              ? std::exp(localPar(a)) : localPar(a);
          }
          arma::mat dL(q,q,arma::fill::zeros);
          const arma::uword i =
            static_cast<arma::uword>(rr[static_cast<int>(k)] - 1);
          const arma::uword j =
            static_cast<arma::uword>(cc[static_cast<int>(k)] - 1);
          dL(i,j) = dd[static_cast<int>(k)]
            ? std::exp(localPar(k)) : 1.0;
          return dL * L.t() + L * dL.t();
        }

        if(op == "fa" || op == "rr"){
          const bool isFA = !sommer::zeroSpecificFactor(f);
          const bool heterogeneous = isFA && sommer::heterogeneousSpecificFactor(f);
          const bool rrMetadata = f.containsElementNamed("rr_nload");
          const int order = Rcpp::as<int>(f["order"]);
          const int nload = Rcpp::as<int>(f[rrMetadata ? "rr_nload" : "fa_nload"]);
          const Rcpp::IntegerVector rows = f[rrMetadata ? "rr_row" : "fa_row"];
          const Rcpp::IntegerVector cols = f[rrMetadata ? "rr_col" : "fa_col"];
          const Rcpp::LogicalVector diagonal = f[rrMetadata ? "rr_diag" : "fa_diag"];
          const arma::uword expectedParameters = static_cast<arma::uword>(
            nload + (heterogeneous ? static_cast<int>(q) - 1 : 0));
          if(order < 1 || nload < 1 || rows.size() != nload ||
             cols.size() != nload || diagonal.size() != nload ||
             localPar.n_elem != expectedParameters){
            Rcpp::stop("Malformed native FA/RR derivative metadata.");
          }

          arma::mat loadings(q, static_cast<arma::uword>(order), arma::fill::zeros);
          arma::mat loadingDerivative(q, static_cast<arma::uword>(order), arma::fill::zeros);
          arma::vec specificRaw(q, arma::fill::ones);
          if(!isFA){ specificRaw.zeros(); }
          double normalization = isFA ? 1.0 : 0.0;
          double normalizationDerivative = 0.0;
          for(int parameter = 0; parameter < nload; ++parameter){
            const arma::uword row = static_cast<arma::uword>(rows[parameter] - 1);
            const arma::uword col = static_cast<arma::uword>(cols[parameter] - 1);
            const bool isDiagonal = diagonal[parameter];
            const double value = isDiagonal
              ? std::exp(localPar(static_cast<arma::uword>(parameter)))
              : localPar(static_cast<arma::uword>(parameter));
            loadings(row, col) = value;
            if(k == static_cast<arma::uword>(parameter)){
              loadingDerivative(row, col) = isDiagonal ? value : 1.0;
              if(row == 0 && col == 0 && isDiagonal){
                normalizationDerivative = 2.0 * value * value;
              }
            }
            if(row == 0 && col == 0){ normalization += value * value; }
          }
          if(heterogeneous){
            for(arma::uword level = 1; level < q; ++level){
              specificRaw(level) = std::exp(localPar(
                static_cast<arma::uword>(nload) + level - 1));
            }
          }
          arma::mat covarianceNumerator = loadings * loadings.t() + arma::diagmat(specificRaw);
          arma::mat numeratorDerivative =
            loadingDerivative * loadings.t() + loadings * loadingDerivative.t();
          if(heterogeneous && k >= static_cast<arma::uword>(nload)){
            const arma::uword level = k - static_cast<arma::uword>(nload) + 1;
            numeratorDerivative(level, level) = specificRaw(level);
          }
          arma::mat derivative =
            numeratorDerivative / normalization -
            (normalizationDerivative / (normalization * normalization)) * covarianceNumerator;
          if(!derivative.is_finite()){
            Rcpp::stop("Native FA/RR derivative produced non-finite values.");
          }
          return 0.5 * (derivative + derivative.t());
        }

        Rcpp::stop("Unsupported native derivative opcode: " + op);
      }

      if(backend == "R"){
        if(!spec.containsElementNamed("fun") || Rf_isNull(spec["fun"])){
          Rcpp::stop("R derivative specification is missing fun(par,k).");
        }
        Rcpp::Function dfun = spec["fun"];
        SEXP ans = dfun(Rcpp::wrap(localPar), static_cast<int>(k + 1));
        arma::mat D = Rcpp::as<arma::mat>(ans);
        if(D.n_rows != q || D.n_cols != q || !D.is_finite()){
          Rcpp::stop("CovarianceFactor derivative callback returned an invalid matrix.");
        }
        return 0.5 * (D + D.t());
      }

      if(backend == "numeric"){
        double relStep = 1.0e-6;
        if(spec.containsElementNamed("rel_step")){
          relStep = Rcpp::as<double>(spec["rel_step"]);
        }
        if(!std::isfinite(relStep) || relStep <= 0.0){
          Rcpp::stop("Numerical derivative rel_step must be positive and finite.");
        }
        const double h = relStep * (1.0 + std::abs(localPar(k)));
        arma::vec plus = localPar;
        arma::vec minus = localPar;
        plus(k) += h;
        minus(k) -= h;
        return (evalFactor(f, plus) - evalFactor(f, minus)) / (2.0*h);
      }

      if(backend == "none"){
        Rcpp::stop("Derivative requested for a CovarianceFactor with no parameters.");
      }

      Rcpp::stop("Unknown CovarianceFactor derivative backend: " + backend);
      return arma::mat();
    }

  int factorSections(const Rcpp::List & f) {
      return f.containsElementNamed("section_levels")
        ? Rf_length(f["section_levels"]) : 0;
    }

  // Parameter slice owned by section t, or -1 for an iid (identity) section.
  int sectionListed(const Rcpp::List & f, const int t) {
      if(f.containsElementNamed("section_index")){
        Rcpp::IntegerVector index = f["section_index"];
        return index[t];
      }
      return t;
    }

  arma::vec sectionSlice(const Rcpp::List & f, const arma::vec & localPar, const int t) {
      const arma::uword npar = static_cast<arma::uword>(Rcpp::as<int>(f["section_npar"]));
      const int listed = sectionListed(f, t);
      if(npar == 0 || listed < 0){
        return arma::vec();
      }
      const arma::uword first = static_cast<arma::uword>(listed) * npar;
      return localPar.subvec(first, first + npar - 1);
    }

  arma::mat evalSectionFactor(const Rcpp::List & f, const arma::vec & localPar, const int t) {
      if(sectionListed(f, t) < 0){
        const arma::uword q = static_cast<arma::uword>(Rcpp::as<int>(f["dim"]));
        return arma::eye<arma::mat>(q, q);
      }
      return evalFactor(f, sectionSlice(f, localPar, t));
    }

  arma::mat evaluateObservedDescriptor(const Rcpp::List & descriptor,
        const arma::vec & par, const arma::uvec & index,
        const arma::uword parameter) {
      const Rcpp::List factors = descriptor["factors"];
      const int dimension = Rcpp::as<int>(descriptor["dim"]);
      if(par.is_empty() || dimension < 1 || index.is_empty() ||
         index.max() >= static_cast<arma::uword>(dimension) || parameter > par.n_elem){
        Rcpp::stop("Observed covariance indices or parameter are out of range.");
      }
      // par.n_elem requests covariance; zero differentiates its log scale.
      arma::umat digits(index.n_elem, factors.size());
      arma::uword trailing = dimension;
      int derivativeAxis = -1;
      for(int axis = 0; axis < factors.size(); ++axis){
        const Rcpp::List factor = factors[axis];
        const int size = Rcpp::as<int>(factor["dim"]);
        if(size < 1 || trailing % size != 0){
          Rcpp::stop("Observed covariance factor dimensions are inconsistent.");
        }
        trailing /= size;
        for(arma::uword row = 0; row < index.n_elem; ++row){
          digits(row, axis) = (index(row) / trailing) % size;
        }
        const int first = Rcpp::as<int>(factor["par_start"]) - 1;
        const int last = Rcpp::as<int>(factor["par_end"]) - 1;
        if(parameter > 0 && parameter < par.n_elem &&
           static_cast<int>(parameter) >= first && static_cast<int>(parameter) <= last){
          derivativeAxis = axis;
        }
      }
      if(trailing != 1 || (parameter > 0 && parameter < par.n_elem && derivativeAxis < 0)){
        Rcpp::stop("Observed covariance parameter does not belong to a valid factor.");
      }
      arma::mat result(index.n_elem, index.n_elem, arma::fill::ones);
      for(int axis = 0; axis < factors.size(); ++axis){
        const Rcpp::List factor = factors[axis];
        const int first = Rcpp::as<int>(factor["par_start"]) - 1;
        const int last = Rcpp::as<int>(factor["par_end"]) - 1;
        arma::vec local;
        if(last >= first){ local = par.subvec(first, last); }
        const int sections = factorSections(factor);
        const int owner = sections ? Rcpp::as<int>(factor["section_factor"]) - 1 : -1;
        if(sections && (owner < 0 || owner >= factors.size())){
          Rcpp::stop("A sectioned covariance factor has no valid section owner.");
        }
        std::vector<arma::mat> pieces(sections ? sections : 1);
        for(std::size_t section = 0; section < pieces.size(); ++section){
          if(axis != derivativeAxis){
            pieces[section] = sections ? evalSectionFactor(factor, local, section)
                                      : evalFactor(factor, local);
          }else if(!sections){
            pieces[section] = factorD1(factor, local, parameter - first);
          }else{
            const int listed = sectionListed(factor, section);
            const int count = Rcpp::as<int>(factor["section_npar"]);
            const int offset = static_cast<int>(parameter) - first;
            if(listed >= 0 && count > 0 && offset / count == listed){
              pieces[section] = factorD1(factor, sectionSlice(factor, local, section),
                                         offset % count);
            }
          }
          const arma::mat & piece = pieces[section];
          const arma::uword size = Rcpp::as<int>(factor["dim"]);
          if(!piece.is_empty() && (piece.n_rows != size || piece.n_cols != size)){
            Rcpp::stop("Observed covariance factor returned an invalid matrix.");
          }
        }
        for(arma::uword column = 0; column < index.n_elem; ++column){
          for(arma::uword row = 0; row < index.n_elem; ++row){
            const arma::uword section = sections ? digits(row, owner) : 0;
            if(section >= pieces.size()){
              Rcpp::stop("Observed covariance section index is out of range.");
            }
            const arma::mat & piece = pieces[section];
            result(row, column) *= piece.is_empty() ? 0.0
              : piece(digits(row, axis), digits(column, axis));
          }
        }
      }
      return std::exp(par(0)) * result;
    }

  arma::mat evaluateDescriptor(const Rcpp::List & cs,
        const arma::vec & par) {

      if(par.n_elem < 1){
        Rcpp::stop("Covariance descriptor must contain log_sigma2.");
      }

      Rcpp::List factors = cs["factors"];
      arma::mat K(1,1,arma::fill::ones);

      for(int fidx = 0; fidx < factors.size(); ++fidx){
        Rcpp::List f = Rcpp::as<Rcpp::List>(factors[fidx]);

        const int start1 =
          f.containsElementNamed("par_start")
          ?
          Rcpp::as<int>(f["par_start"])
          :
          1;

        const int end1 =
          f.containsElementNamed("par_end")
          ?
          Rcpp::as<int>(f["par_end"])
          :
          0;

        arma::vec localPar;
        if(end1 >= start1){
          localPar =
            par.subvec(
              static_cast<arma::uword>(start1 - 1),
              static_cast<arma::uword>(end1 - 1)
            );
        }

        K = arma::kron(K, evalFactor(f, localPar));
      }

      return std::exp(par(0)) * K;
    }

  arma::mat descriptorD1(const Rcpp::List & cs,
        const arma::vec & par,
        const arma::uword k) {

      if(k == 0){
        return evaluateDescriptor(cs, par);
      }

      Rcpp::List factors = cs["factors"];
      arma::mat K(1,1,arma::fill::ones);
      bool found = false;

      for(int fidx = 0; fidx < factors.size(); ++fidx){
        Rcpp::List f = Rcpp::as<Rcpp::List>(factors[fidx]);

        const int start1 =
          f.containsElementNamed("par_start")
          ?
          Rcpp::as<int>(f["par_start"])
          :
          1;

        const int end1 =
          f.containsElementNamed("par_end")
          ?
          Rcpp::as<int>(f["par_end"])
          :
          0;

        arma::vec localPar;
        if(end1 >= start1){
          localPar =
            par.subvec(
              static_cast<arma::uword>(start1 - 1),
              static_cast<arma::uword>(end1 - 1)
            );
        }

        arma::mat piece;

        const int k1 = static_cast<int>(k) + 1;
        if(end1 >= start1 && k1 >= start1 && k1 <= end1){
          piece =
            factorD1(
              f,
              localPar,
              static_cast<arma::uword>(k1 - start1)
            );
          found = true;
        }else{
          piece = evalFactor(f, localPar);
        }

        K = arma::kron(K, piece);
      }

      if(!found){
        Rcpp::stop("Covariance derivative parameter does not belong to any factor.");
      }

      return std::exp(par(0)) * K;
    }

  arma::mat factorPrecisionDense(const DescriptorPrecisionState & state,
                                 const std::size_t factor) {
    if(factor >= state.factorPrecision.size()){
      Rcpp::stop("Covariance precision factor index is out of range.");
    }
    if(!state.factorPrecision[factor].is_empty()){
      return state.factorPrecision[factor];
    }
    if(factor < state.sparseFactorPrecision.size() &&
       !state.sparseFactorPrecision[factor].is_empty()){
      return arma::mat(state.sparseFactorPrecision[factor]);
    }
    Rcpp::stop("Covariance precision factor has no stored representation.");
    return arma::mat();
  }

  arma::uword factorPrecisionDimension(const DescriptorPrecisionState & state,
                                     const std::size_t factor) {
    if(factor >= state.factorCovariance.size()){
      Rcpp::stop("Covariance precision factor index is out of range.");
    }
    return state.factorCovariance[factor].n_rows;
  }

  double factorPrecisionValue(const DescriptorPrecisionState & state,
                                 const std::size_t factor,
                                 const arma::uword row,
                                 const arma::uword column) {
    if(factor < state.factorPrecision.size() && !state.factorPrecision[factor].is_empty()){
      return state.factorPrecision[factor](row, column);
    }
    if(factor < state.sparseFactorPrecision.size() &&
       !state.sparseFactorPrecision[factor].is_empty()){
      return state.sparseFactorPrecision[factor](row, column);
    }
    Rcpp::stop("Covariance precision factor has no stored representation.");
    return 0.0;
  }

  double factorPrecisionCovarianceTrace(const DescriptorPrecisionState & state,
                                            const std::size_t factor,
                                            const arma::mat & covarianceDerivative) {
    if(factor < state.sparseFactorPrecision.size() &&
       !state.sparseFactorPrecision[factor].is_empty()){
      double trace = 0.0;
      for(arma::sp_mat::const_iterator entry = state.sparseFactorPrecision[factor].begin();
          entry != state.sparseFactorPrecision[factor].end(); ++entry){
        trace += (*entry) * covarianceDerivative(entry.col(), entry.row());
      }
      return trace;
    }
    return arma::accu(state.factorPrecision[factor] % covarianceDerivative.t());
  }

  std::vector<arma::mat> factorPrecisionFactorsDense(const DescriptorPrecisionState & state) {
    std::vector<arma::mat> factors;
    factors.reserve(state.factorPrecision.size());
    for(std::size_t factor = 0; factor < state.factorPrecision.size(); ++factor){
      factors.push_back(factorPrecisionDense(state, factor));
    }
    return factors;
  }

  arma::mat localPrecisionD1Dense(const DescriptorPrecisionState & state,
                                  const std::size_t parameter) {
    if(parameter < state.localPrecisionD1.size() &&
       !state.localPrecisionD1[parameter].is_empty()){
      return state.localPrecisionD1[parameter];
    }
    if(parameter < state.sparseLocalPrecisionD1.size() &&
       !state.sparseLocalPrecisionD1[parameter].is_empty()){
      return arma::mat(state.sparseLocalPrecisionD1[parameter]);
    }
    Rcpp::stop("Local covariance precision derivative is unavailable.");
    return arma::mat();
  }

  arma::mat localCovarianceD1Dense(const DescriptorPrecisionState & state,
                                    const std::size_t parameter) {
    if(parameter < state.localCovarianceD1.size() &&
       !state.localCovarianceD1[parameter].is_empty()){
      return state.localCovarianceD1[parameter];
    }
    if(parameter >= state.derivativeFactor.size()){
      Rcpp::stop("Local covariance derivative is unavailable.");
    }
    const int factor = state.derivativeFactor[parameter];
    if(factor < 0){ return arma::mat(); }
    const std::size_t axis = static_cast<std::size_t>(factor);
    if(axis >= state.factorDescriptor.size() || axis >= state.factorStart1.size() ||
       axis >= state.factorPar.size() || parameter + 1 < static_cast<std::size_t>(state.factorStart1[axis])){
      Rcpp::stop("Local covariance derivative factor is unavailable.");
    }
    const arma::uword localParameter = static_cast<arma::uword>(parameter + 1) -
      static_cast<arma::uword>(state.factorStart1[axis]);
    return factorD1(state.factorDescriptor[axis], state.factorPar[axis], localParameter);
  }

  bool buildStationaryARPrecision(const Rcpp::List & factor,
                                  const arma::vec & shape,
                                  const arma::mat & covariance,
                                  arma::mat & precision,
                                  double & logDet) {
      const int order = Rcpp::as<int>(Rcpp::as<Rcpp::List>(factor["representation"])["order"]);
      const arma::uword q = covariance.n_rows;
      if(order < 1 || q <= static_cast<arma::uword>(order) ||
         shape.n_elem < static_cast<arma::uword>(order)){
        return false;
      }

      arma::vec phi;
      for(int currentOrder = 1; currentOrder <= order; ++currentOrder){
        const double reflection = std::tanh(shape(static_cast<arma::uword>(currentOrder - 1)));
        arma::vec updated(static_cast<arma::uword>(currentOrder), arma::fill::zeros);
        updated(static_cast<arma::uword>(currentOrder - 1)) = reflection;
        if(currentOrder > 1){
          for(int lag = 0; lag < currentOrder - 1; ++lag){
            updated(static_cast<arma::uword>(lag)) =
              phi(static_cast<arma::uword>(lag)) - reflection *
                phi(static_cast<arma::uword>(currentOrder - 2 - lag));
          }
        }
        phi = std::move(updated);
      }

      const arma::vec standardDeviation = arma::sqrt(covariance.diag());
      if(!standardDeviation.is_finite() || standardDeviation.min() <= 0.0){
        return false;
      }
      const arma::vec inverseSd = 1.0 / standardDeviation;
      const arma::uword p = static_cast<arma::uword>(order);
      arma::mat initial(p, p);
      for(arma::uword column = 0; column < p; ++column){
        for(arma::uword row = 0; row < p; ++row){
          initial(row, column) = inverseSd(row) * covariance(row, column) *
            inverseSd(column);
        }
      }
      arma::mat initialPrecision;
      double initialLogDet = 0.0;
      if(!eigenSpdInverse(initial, initialPrecision, &initialLogDet)){
        return false;
      }

      double innovationVariance = 1.0;
      for(arma::uword lag = 0; lag < p; ++lag){
        innovationVariance -= phi(lag) *
          (inverseSd(0) * covariance(0, lag + 1) * inverseSd(lag + 1));
      }
      if(!std::isfinite(innovationVariance) || innovationVariance <= 0.0){
        return false;
      }

      arma::mat direct(q, q, arma::fill::zeros);
      for(arma::uword column = 0; column < p; ++column){
        for(arma::uword row = 0; row < p; ++row){
          direct(row, column) = initialPrecision(row, column) *
            inverseSd(row) * inverseSd(column);
        }
      }
      arma::vec innovation(p + 1, arma::fill::ones);
      innovation.subvec(1, p) = -phi;
      for(arma::uword record = p; record < q; ++record){
        for(arma::uword column = 0; column <= p; ++column){
          for(arma::uword row = 0; row <= p; ++row){
            direct(record - row, record - column) +=
              innovation(row) * innovation(column) / innovationVariance *
              inverseSd(record - row) * inverseSd(record - column);
          }
        }
      }
      precision = std::move(direct);
      logDet = initialLogDet + static_cast<double>(q - p) *
        std::log(innovationVariance) + 2.0 * arma::accu(arma::log(standardDeviation));
      return precision.is_finite() && precision.diag().min() > 0.0 &&
        std::isfinite(logDet);
    }

  bool buildDescriptorPrecision(const Rcpp::List & cs,
        const arma::vec & par,
        DescriptorPrecisionState & out) {

      out = DescriptorPrecisionState();

      if(par.n_elem < 1){
        Rcpp::stop("Covariance descriptor must contain log_sigma2.");
      }

      Rcpp::List factors = cs["factors"];
      const std::size_t nFactors =
        static_cast<std::size_t>(factors.size());

      std::vector<arma::mat> factorPrecision(nFactors);
      std::vector<arma::sp_mat> sparseFactorPrecision(nFactors);
      std::vector<arma::mat> factorCovariance(nFactors);
      std::vector<arma::vec> factorPar(nFactors);
      std::vector<int> factorStart1(nFactors, 1);
      std::vector<int> factorEnd1(nFactors, 0);
      std::vector<double> factorLogDet(nFactors, 0.0);
      arma::mat lowRankLoadings;
      arma::vec lowRankSpecific;
      arma::vec lowRankPrecisionDiagonal;
      int lowRankFactorIndex = -1;
      double lowRankCoreLogDet = 0.0;

      arma::uword totalDim = 1;

      for(int fidx = 0; fidx < factors.size(); ++fidx){
        const std::size_t factorOffset =
          static_cast<std::size_t>(fidx);
        Rcpp::List f = Rcpp::as<Rcpp::List>(factors[fidx]);

        const int start1 =
          f.containsElementNamed("par_start")
          ? Rcpp::as<int>(f["par_start"])
          : 1;
        const int end1 =
          f.containsElementNamed("par_end")
          ? Rcpp::as<int>(f["par_end"])
          : 0;

        factorStart1[factorOffset] = start1;
        factorEnd1[factorOffset] = end1;

        if(end1 >= start1){
          factorPar[factorOffset] =
            par.subvec(
              static_cast<arma::uword>(start1 - 1),
              static_cast<arma::uword>(end1 - 1)
            );
        }

        const arma::mat covariance =
          evalFactor(f, factorPar[factorOffset]);
        factorCovariance[factorOffset] = covariance;

        bool inverseOK = false;
        bool bandedPrecision = false;
        const sommer::PrecisionFamily family = sommer::precisionFamily(f);
        if(sparseFactorKernelsEnabled(f) && sommer::FactorPolicy::spatial(family)){
          const SparseSpatialFactor spatial = evaluateSparseSpatialCached(f, factorPar[factorOffset], false);
          factorPrecision[factorOffset] = arma::mat(spatial.precision);
          factorLogDet[factorOffset] = spatial.logDet;
          inverseOK = true;
          bandedPrecision = true;
        }
        if(family == sommer::PrecisionFamily::LowRankDiagonal){
          if(sommer::zeroSpecificFactor(f)){
            Rcpp::stop("Zero-specific reduced-rank covariance has no full-space precision; use the latent Henderson path.");
          }
          const Rcpp::List capability = sommer::factorCapability(f, "loadings");
          const bool heterogeneousSpecific = Rcpp::as<bool>(capability["heterogeneous_specific"]);
          const int order = Rcpp::as<int>(capability["rank"]);
          const int nload = Rcpp::as<int>(capability["count"]);
          const Rcpp::IntegerVector rows = capability["rows"];
          const Rcpp::IntegerVector columns = capability["columns"];
          const Rcpp::LogicalVector diagonal = capability["diagonal"];
          if(order < 1 || nload < 1 || rows.size() != nload || columns.size() != nload ||
             diagonal.size() != nload || factorPar[factorOffset].n_elem !=
               static_cast<arma::uword>(nload) + (heterogeneousSpecific ? covariance.n_rows - 1 : 0)){
            Rcpp::stop("Malformed low-rank precision capability.");
          }
          double reference = 0.0;
          for(int loading = 0; loading < nload; ++loading){
            if(rows[loading] == 1 && columns[loading] == 1){
              reference = factorPar[factorOffset](loading);
              break;
            }
          }
          const double twiceReference = 2.0 * reference;
          const double logScale = std::max(0.0, twiceReference) + std::log1p(std::exp(-std::abs(twiceReference)));
          const double inverseSd = std::exp(-0.5 * logScale);
          arma::mat loadings(covariance.n_rows, order, arma::fill::zeros);
          for(int loading = 0; loading < nload; ++loading){
            loadings(rows[loading] - 1, columns[loading] - 1) = diagonal[loading]
              ? std::exp(factorPar[factorOffset](loading) - 0.5 * logScale)
              : factorPar[factorOffset](loading) * inverseSd;
          }
          arma::vec specific(covariance.n_rows, arma::fill::ones);
          specific *= std::exp(-logScale);
          if(heterogeneousSpecific){
            for(arma::uword row = 1; row < specific.n_elem; ++row){
              specific(row) = std::exp(factorPar[factorOffset](nload + row - 1) - logScale);
            }
          }
          if(loadings.is_finite() && specific.is_finite() &&
             specific.min() > 1.0e-8 * covariance.diag().max() &&
             arma::norm(covariance - loadings * loadings.t() - arma::diagmat(specific), "fro") <=
               1.0e-10 * std::max(1.0, arma::norm(covariance, "fro"))){
            arma::mat weightedLoadings = loadings;
            weightedLoadings.each_col() /= specific;
            const arma::mat core = arma::eye<arma::mat>(order, order) + loadings.t() * weightedLoadings;
            arma::mat coreInverse;
            double coreLogDet = 0.0;
            if(eigenSpdInverse(core, coreInverse, &coreLogDet)){
              factorPrecision[factorOffset] = arma::diagmat(1.0 / specific) - weightedLoadings * coreInverse * weightedLoadings.t();
              factorPrecision[factorOffset] = 0.5 * (factorPrecision[factorOffset] + factorPrecision[factorOffset].t());
              factorLogDet[factorOffset] = arma::accu(arma::log(specific)) + coreLogDet;
              inverseOK = factorPrecision[factorOffset].is_finite() &&
                factorPrecision[factorOffset].diag().min() > 0.0 &&
                std::isfinite(factorLogDet[factorOffset]);
              if(inverseOK && lowRankFactorIndex < 0){
                lowRankLoadings = loadings;
                lowRankSpecific = specific;
                lowRankPrecisionDiagonal = factorPrecision[factorOffset].diag();
                lowRankFactorIndex = fidx;
                lowRankCoreLogDet = coreLogDet;
              }
            }
          }
        }else if(family == sommer::PrecisionFamily::CompoundSymmetry &&
                 covariance.n_rows >= 2){
          const arma::vec standardDeviation = arma::sqrt(covariance.diag());
          const double rho = covariance(0,1) /
            (standardDeviation(0) * standardDeviation(1));
          if(standardDeviation.is_finite() &&
             standardDeviation.min() > 0.0 &&
             std::isfinite(rho) && rho >= 0.0 && rho < 1.0){
            arma::mat loadings(covariance.n_rows, 1, arma::fill::zeros);
            loadings.col(0) = standardDeviation * std::sqrt(rho);
            const arma::vec specific =
              covariance.diag() * (1.0 - rho);
            const arma::mat reconstructed =
              loadings * loadings.t() + arma::diagmat(specific);
            if(specific.is_finite() &&
               specific.min() > 1.0e-8 * covariance.diag().max() &&
               arma::norm(covariance - reconstructed, "fro") <=
                 1.0e-10 * std::max(1.0, arma::norm(covariance, "fro"))){
              arma::mat weightedLoadings = loadings;
              weightedLoadings.each_col() /= specific;
              const arma::mat core =
                arma::eye<arma::mat>(1,1) + loadings.t() * weightedLoadings;
              arma::mat coreInverse;
              double coreLogDet = 0.0;
              if(eigenSpdInverse(core, coreInverse, &coreLogDet)){
                factorPrecision[factorOffset] =
                  arma::diagmat(1.0 / specific) -
                  weightedLoadings * coreInverse * weightedLoadings.t();
                factorPrecision[factorOffset] = 0.5 *
                  (factorPrecision[factorOffset] + factorPrecision[factorOffset].t());
                factorLogDet[factorOffset] =
                  arma::accu(arma::log(specific)) + coreLogDet;
                inverseOK = factorPrecision[factorOffset].is_finite() &&
                  factorPrecision[factorOffset].diag().min() > 0.0 &&
                  std::isfinite(factorLogDet[factorOffset]);
                if(inverseOK && lowRankFactorIndex < 0){
                  lowRankLoadings = std::move(loadings);
                  lowRankSpecific = specific;
                  lowRankPrecisionDiagonal = factorPrecision[factorOffset].diag();
                  lowRankFactorIndex = fidx;
                  lowRankCoreLogDet = coreLogDet;
                }
              }
            }
          }
        }
        if(!inverseOK && sommer::FactorPolicy::spatial(family)){
          const bool normalEquations = family == sommer::PrecisionFamily::NormalEquations;
          const Rcpp::List capability = sommer::factorCapability(f, "spatial");
          {
            const arma::uword q = covariance.n_rows;
            const arma::mat W = Rcpp::as<arma::mat>(capability["weights"]);
            const double eta = factorPar[factorOffset](0);
            const double logistic = eta >= 0.0
              ? 1.0 / (1.0 + std::exp(-eta))
              : std::exp(eta) / (1.0 + std::exp(eta));
            double rho = 0.0;
            arma::mat system(q, q, arma::fill::eye);
            const double lower = Rcpp::as<double>(capability["lower"]);
            const double upper = Rcpp::as<double>(capability["upper"]);
            rho = lower + (upper - lower) * logistic;
            if(normalEquations){
              system -= rho * W;
            }else{
              const arma::vec degree = Rcpp::as<arma::vec>(capability["degree"]);
              if(degree.n_elem != q || degree.min() <= 0.0){
                return false;
              }
              system = arma::diagmat(degree) - rho * W;
            }
            double systemLogDet = 0.0;
            double systemSign = 0.0;
            const bool haveLogDet = arma::log_det(systemLogDet, systemSign, system);
            arma::vec firstColumn = arma::zeros<arma::vec>(q);
            firstColumn(0) = 1.0;
            arma::vec normalizationSolve;
            const bool haveSolve = normalEquations
              ? arma::solve(normalizationSolve, system.t(), firstColumn)
              : arma::solve(normalizationSolve, system, firstColumn);
            if(haveLogDet && haveSolve && std::isfinite(systemLogDet) &&
               systemSign != 0.0 && normalizationSolve.is_finite()){
              const double normalization = normalEquations
                ? arma::dot(normalizationSolve, normalizationSolve)
                : normalizationSolve(0);
              if(std::isfinite(normalization) && normalization > 0.0){
                if(normalEquations){
                  factorPrecision[factorOffset] = normalization * system.t() * system;
                  factorLogDet[factorOffset] = -2.0 * systemLogDet -
                    static_cast<double>(q) * std::log(normalization);
                }else{
                  factorPrecision[factorOffset] = normalization * system;
                  factorLogDet[factorOffset] = -systemLogDet -
                    static_cast<double>(q) * std::log(normalization);
                }
                factorPrecision[factorOffset] = 0.5 *
                  (factorPrecision[factorOffset] + factorPrecision[factorOffset].t());
                inverseOK = factorPrecision[factorOffset].is_finite() &&
                  factorPrecision[factorOffset].diag().min() > 0.0 &&
                  std::isfinite(factorLogDet[factorOffset]);
              }
            }
          }
        }
        if(!inverseOK && family == PrecisionFamily::StationaryAR){
          inverseOK = buildStationaryARPrecision(
            f, factorPar[factorOffset], covariance,
            factorPrecision[factorOffset], factorLogDet[factorOffset]);
          bandedPrecision = inverseOK;
        }

        if(!inverseOK){
          {
            const arma::uword q = covariance.n_rows;

            if(family == PrecisionFamily::AR1 && q >= 2){
              const double rho = covariance(0,1);
              const double denominator = 1.0 - rho*rho;
              if(std::isfinite(denominator) && denominator > 0.0){
                arma::mat direct(q,q,arma::fill::zeros);
                direct(0,0) = 1.0 / denominator;
                direct(q-1,q-1) = 1.0 / denominator;
                for(arma::uword row = 1; row + 1 < q; ++row){
                  direct(row,row) = (1.0 + rho*rho) / denominator;
                }
                for(arma::uword row = 0; row + 1 < q; ++row){
                  direct(row,row+1) = -rho / denominator;
                  direct(row+1,row) = direct(row,row+1);
                }
                factorPrecision[factorOffset] = std::move(direct);
                factorLogDet[factorOffset] =
                  static_cast<double>(q - 1) * std::log(denominator);
                inverseOK = true;
                bandedPrecision = true;
              }

            }else if(family == PrecisionFamily::CompoundSymmetry && q >= 2){
              const arma::vec standardDeviation =
                arma::sqrt(covariance.diag());
              const double rho =
                covariance(0,1)
                /
                (standardDeviation(0) * standardDeviation(1));
              const double independent = 1.0 - rho;
              const double common = 1.0 + static_cast<double>(q - 1) * rho;
              if(standardDeviation.is_finite() &&
                 standardDeviation.min() > 0.0 &&
                 independent > 0.0 && common > 0.0){
                arma::mat correlationPrecision =
                  arma::eye<arma::mat>(q,q) / independent;
                correlationPrecision -=
                  arma::ones<arma::mat>(q,q)
                  *
                  (rho / (independent * common));
                const arma::vec inverseSd = 1.0 / standardDeviation;
                factorPrecision[factorOffset] =
                  arma::diagmat(inverseSd)
                  *
                  correlationPrecision
                  *
                  arma::diagmat(inverseSd);
                factorLogDet[factorOffset] =
                  2.0 * arma::accu(arma::log(standardDeviation))
                  +
                  static_cast<double>(q - 1) * std::log(independent)
                  +
                  std::log(common);
                inverseOK = true;
              }

            }else if(family == sommer::PrecisionFamily::ModifiedCholesky){
              const Rcpp::List capability = sommer::factorCapability(f, "innovation");
              const int ncoef = Rcpp::as<int>(capability["coefficient_count"]);
              Rcpp::IntegerVector rr = capability["rows"];
              Rcpp::IntegerVector cc = capability["columns"];
              if(rr.size() != ncoef || cc.size() != ncoef){
                Rcpp::stop("Malformed innovation precision capability.");
              }
              for(int coefficient = 0; coefficient < ncoef; ++coefficient){
                if(cc[coefficient] < 1 || rr[coefficient] > static_cast<int>(q) ||
                   rr[coefficient] <= cc[coefficient]){
                  Rcpp::stop("Invalid innovation precision coordinates.");
                }
              }
              if(ncoef >= 1 &&
                 factorPar[factorOffset].n_elem ==
                   static_cast<arma::uword>(ncoef) + q - 1){
                std::vector<std::vector<std::pair<arma::uword, double> > > innovationRows(q);
                for(arma::uword row = 0; row < q; ++row){
                  innovationRows[row].emplace_back(row, 1.0);
                }
                for(int coefficient = 0; coefficient < ncoef; ++coefficient){
                  const double value = -factorPar[factorOffset](static_cast<arma::uword>(coefficient));
                  if(value != 0.0){
                    innovationRows[static_cast<arma::uword>(rr[coefficient] - 1)].emplace_back(
                      static_cast<arma::uword>(cc[coefficient] - 1), value);
                  }
                }
                arma::vec inverseInnovation(q,arma::fill::ones);
                double innovationLogDet = 0.0;
                for(arma::uword row = 1; row < q; ++row){
                  const double logInnovation =
                    factorPar[factorOffset](
                      static_cast<arma::uword>(ncoef) + row - 1
                    );
                  inverseInnovation(row) = std::exp(-logInnovation);
                  innovationLogDet += logInnovation;
                }
                arma::mat direct(q, q, arma::fill::zeros);
                for(arma::uword record = 0; record < q; ++record){
                  auto & active = innovationRows[record];
                  std::sort(active.begin(), active.end());
                  for(const auto & column : active){
                    for(const auto & row : active){
                      direct(row.first, column.first) +=
                        inverseInnovation(record) * row.second * column.second;
                    }
                  }
                }
                factorPrecision[factorOffset] = std::move(direct);
                factorLogDet[factorOffset] = innovationLogDet;
                inverseOK =
                  factorPrecision[factorOffset].is_finite()
                  &&
                  std::isfinite(innovationLogDet);
                bandedPrecision = inverseOK;
              }
            }
          }
        }

        if(f.containsElementNamed("structurally_diagonal") &&
           !inverseOK &&
           Rcpp::as<bool>(f["structurally_diagonal"])){
          arma::mat offDiagonal = covariance;
          offDiagonal.diag().zeros();
          const arma::vec diagonal = covariance.diag();
          if(arma::accu(arma::abs(offDiagonal)) == 0.0 &&
             diagonal.is_finite() && diagonal.min() > 0.0){
            factorPrecision[factorOffset] =
              arma::diagmat(1.0 / diagonal);
            factorLogDet[factorOffset] =
              arma::accu(arma::log(diagonal));
            inverseOK = true;
          }
        }

        if(!inverseOK){
          inverseOK =
            eigenSpdInverse(
              covariance,
              factorPrecision[factorOffset],
              &factorLogDet[factorOffset]
            );
        }

        if(!inverseOK){
          return false;
        }
        if(bandedPrecision && sparseFactorKernelsEnabled(f)){
          sparseFactorPrecision[factorOffset] = arma::sp_mat(factorPrecision[factorOffset]);
        }

        totalDim *= factorPrecision[factorOffset].n_rows;
      }

      const double inverseScale = std::exp(-par(0));
      if(!std::isfinite(inverseScale)){
        return false;
      }

      double logDet = static_cast<double>(totalDim) * par(0);
      for(std::size_t fidx = 0; fidx < nFactors; ++fidx){
        const arma::uword factorDim = factorPrecision[fidx].n_rows;
        logDet +=
          (
            static_cast<double>(totalDim)
            /
            static_cast<double>(factorDim)
          )
          *
          factorLogDet[fidx];
      }
      if(!std::isfinite(logDet)){
        return false;
      }

      for(std::size_t factor = 0; factor < nFactors; ++factor){
        if(!sparseFactorPrecision[factor].is_empty()){
          factorPrecision[factor].reset();
        }
      }

      out.factorCovariance = std::move(factorCovariance);
      out.factorPrecision = std::move(factorPrecision);
      out.sparseFactorPrecision = std::move(sparseFactorPrecision);
      out.factorDescriptor.clear();
      out.factorDescriptor.reserve(nFactors);
      for(int factor = 0; factor < factors.size(); ++factor){
        out.factorDescriptor.push_back(Rcpp::as<Rcpp::List>(factors[factor]));
      }
      out.factorLogDet = std::move(factorLogDet);
      out.factorPar = std::move(factorPar);
      out.factorStart1 = std::move(factorStart1);
      out.factorEnd1 = std::move(factorEnd1);
      out.lowRankLoadings = std::move(lowRankLoadings);
      out.lowRankSpecific = std::move(lowRankSpecific);
      out.lowRankPrecisionDiagonal = std::move(lowRankPrecisionDiagonal);
      out.lowRankFactorIndex = lowRankFactorIndex;
      out.lowRankCoreLogDet = lowRankCoreLogDet;
      out.inverseScale = inverseScale;
      out.logDet = logDet;
      out.totalDim = totalDim;
      out.factorWise = true;
      return true;
    }

  void buildDescriptorPrecisionDerivatives(const Rcpp::List & cs,
        const arma::vec & par,
        DescriptorPrecisionState & state,
        double * covarianceDerivativeSeconds = nullptr,
        double * precisionDerivativeSeconds = nullptr) {

      if(!state.factorWise){
        return;
      }

      Rcpp::List factors = cs["factors"];
      // Compact residual states never materialize the dense precision.
      const arma::uword totalDim = state.totalDim;

      state.precisionD1.clear();
      state.logDetD1.assign(static_cast<std::size_t>(par.n_elem), 0.0);
      state.localCovarianceD1.resize(static_cast<std::size_t>(par.n_elem));
      state.localPrecisionD1.resize(static_cast<std::size_t>(par.n_elem));
      state.sparseLocalPrecisionD1.resize(static_cast<std::size_t>(par.n_elem));
      state.derivativeFactor.assign(static_cast<std::size_t>(par.n_elem), -2);
      std::vector<bool> derivativeAssigned(
        static_cast<std::size_t>(par.n_elem),
        false
      );

      state.logDetD1[0] = static_cast<double>(totalDim);
      state.derivativeFactor[0] = -1;
      derivativeAssigned[0] = true;

      for(int fidx = 0; fidx < factors.size(); ++fidx){
        const std::size_t factorOffset = static_cast<std::size_t>(fidx);
        const int start1 = state.factorStart1[factorOffset];
        const int end1 = state.factorEnd1[factorOffset];
        if(end1 < start1){
          continue;
        }

        Rcpp::List f = Rcpp::as<Rcpp::List>(factors[fidx]);
        const sommer::PrecisionFamily family = sommer::precisionFamily(f);
        const arma::uword arBandwidth = family == sommer::PrecisionFamily::StationaryAR
          ? static_cast<arma::uword>(Rcpp::as<int>(f["order"]))
          : 0;
        arma::sp_mat stationaryARPrecision;
        if(family == sommer::PrecisionFamily::StationaryAR &&
           factorOffset < state.factorPrecision.size() &&
           !state.factorPrecision[factorOffset].is_empty()){
          const arma::mat & densePrecision = state.factorPrecision[factorOffset];
          const arma::uword capacity = densePrecision.n_rows *
            std::min(densePrecision.n_rows, 2 * arBandwidth + 1);
          arma::umat locations(2, capacity);
          arma::vec values(capacity);
          arma::uword count = 0;
          for(arma::uword column = 0; column < densePrecision.n_cols; ++column){
            const arma::uword first = column > arBandwidth ? column - arBandwidth : 0;
            const arma::uword last = std::min(densePrecision.n_rows - 1, column + arBandwidth);
            for(arma::uword row = first; row <= last; ++row){
              locations(0, count) = row;
              locations(1, count) = column;
              values(count) = densePrecision(row, column);
              ++count;
            }
          }
          locations.resize(2, count);
          values.resize(count);
          stationaryARPrecision = arma::sp_mat(
            locations, values, densePrecision.n_rows, densePrecision.n_cols, true, true);
        }
        for(int k1 = start1; k1 <= end1; ++k1){
          const arma::uword localK =
            static_cast<arma::uword>(k1 - start1);
          arma::mat factorCovarianceD1;
          arma::mat factorPrecisionD1;
          arma::sp_mat sparseFactorPrecisionD1;
          const bool exactInnovationDerivative = sommer::FactorPolicy::innovationDerivative(
            family, sparseFactorKernelsEnabled(f));
          double analyticLogDetD1 = std::numeric_limits<double>::quiet_NaN();
          bool analyticSpatialDerivative = false;
           if(sparseFactorKernelsEnabled(f) && sommer::FactorPolicy::spatial(family)){
            const SparseSpatialFactor spatial = evaluateSparseSpatialCached(f, state.factorPar[factorOffset], true);
            factorCovarianceD1 = spatial.covarianceD1;
            factorPrecisionD1 = arma::mat(spatial.precisionD1);
            analyticLogDetD1 = spatial.logDetD1;
            analyticSpatialDerivative = true;
          }

          if(!analyticSpatialDerivative && sommer::FactorPolicy::spatial(family)){
            const bool normalEquations = family == sommer::PrecisionFamily::NormalEquations;
            const Rcpp::List capability = sommer::factorCapability(f, "spatial");
            if(localK == 0 &&
               state.factorPar[factorOffset].n_elem == 1){
              const arma::uword dimension = state.factorCovariance[factorOffset].n_rows;
              const arma::mat W = Rcpp::as<arma::mat>(capability["weights"]);
              const double eta = state.factorPar[factorOffset](0);
              const double logistic = eta >= 0.0
                ? 1.0 / (1.0 + std::exp(-eta))
                : std::exp(eta) / (1.0 + std::exp(eta));
              double rho = 0.0;
              double drho = 0.0;
              arma::mat system;
              const double lower = Rcpp::as<double>(capability["lower"]);
              const double upper = Rcpp::as<double>(capability["upper"]);
              rho = lower + (upper - lower) * logistic;
              drho = (upper - lower) * logistic * (1.0 - logistic);
              if(normalEquations){
                system.eye(dimension, dimension);
                system -= rho * W;
              }else{
                const arma::vec degree = Rcpp::as<arma::vec>(capability["degree"]);
                system = arma::diagmat(degree) - rho * W;
              }

              arma::vec firstColumn(dimension, arma::fill::zeros);
              firstColumn(0) = 1.0;
              arma::vec normalizationSolve;
              arma::mat inverseSystemTimesW;
              bool solvesOK = false;
              double normalization = 0.0;
              double normalizationD1 = 0.0;
              if(normalEquations){
                solvesOK = arma::solve(normalizationSolve, system.t(), firstColumn) &&
                  arma::solve(inverseSystemTimesW, system, W);
                if(solvesOK){
                  normalization = arma::dot(normalizationSolve, normalizationSolve);
                  arma::vec normalizationSolveD1;
                  solvesOK = arma::solve(normalizationSolveD1, system.t(),
                    drho * W.t() * normalizationSolve);
                  if(solvesOK){
                    normalizationD1 = 2.0 * arma::dot(normalizationSolve, normalizationSolveD1);
                  }
                }
              }else{
                solvesOK = arma::solve(normalizationSolve, system, firstColumn) &&
                  arma::solve(inverseSystemTimesW, system, W);
                if(solvesOK){
                  normalization = normalizationSolve(0);
                  arma::vec normalizationSolveD1;
                  solvesOK = arma::solve(normalizationSolveD1, system,
                    drho * W * normalizationSolve);
                  if(solvesOK){ normalizationD1 = normalizationSolveD1(0); }
                }
              }

              if(solvesOK && std::isfinite(normalization) && normalization > 0.0 &&
                 std::isfinite(normalizationD1) && inverseSystemTimesW.is_finite()){
                const arma::mat systemD1 = -drho * W;
                if(normalEquations){
                  factorPrecisionD1 = normalizationD1 * system.t() * system +
                    normalization * (systemD1.t() * system + system.t() * systemD1);
                  analyticLogDetD1 = -2.0 * drho * arma::trace(inverseSystemTimesW) +
                    static_cast<double>(dimension) * normalizationD1 / normalization;
                }else{
                  factorPrecisionD1 = normalizationD1 * system + normalization * systemD1;
                  analyticLogDetD1 = -drho * arma::trace(inverseSystemTimesW) +
                    static_cast<double>(dimension) * normalizationD1 / normalization;
                }
                factorPrecisionD1 = 0.5 * (factorPrecisionD1 + factorPrecisionD1.t());
                const arma::mat & covariance = state.factorCovariance[factorOffset];
                factorCovarianceD1 = -covariance * factorPrecisionD1 * covariance;
                factorCovarianceD1 = 0.5 * (factorCovarianceD1 + factorCovarianceD1.t());
                analyticLogDetD1 = -analyticLogDetD1;
                analyticSpatialDerivative = factorPrecisionD1.is_finite() &&
                  factorCovarianceD1.is_finite() && std::isfinite(analyticLogDetD1);
              }
            }
          }

          if(exactInnovationDerivative){
            const Rcpp::List capability = sommer::factorCapability(f, "innovation");
            const Rcpp::IntegerVector innovationRows = capability["rows"];
            const Rcpp::IntegerVector innovationColumns = capability["columns"];
            sparseFactorPrecisionD1 = antePrecisionDerivative(
              state.factorCovariance[factorOffset].n_rows,
              innovationRows,
              innovationColumns,
              state.factorPar[factorOffset],
              localK);
          }else if(!analyticSpatialDerivative){
            const auto derivativeStart = covarianceDerivativeSeconds != nullptr
              ? std::chrono::steady_clock::now()
              : std::chrono::steady_clock::time_point();
            factorCovarianceD1 = factorD1(f, state.factorPar[factorOffset], localK);
            if(covarianceDerivativeSeconds != nullptr){
              *covarianceDerivativeSeconds += std::chrono::duration<double>(
                std::chrono::steady_clock::now() - derivativeStart).count();
            }
          }

          if(!exactInnovationDerivative && !analyticSpatialDerivative){
            if(family == PrecisionFamily::AR1 &&
               localK == 0 && state.factorCovariance[factorOffset].n_rows >= 2){
              const arma::mat & covariance = state.factorCovariance[factorOffset];
              const arma::uword dimension = covariance.n_rows;
              const double rho = covariance(0, 1);
              const double denominator = 1.0 - rho * rho;
              factorPrecisionD1.zeros(dimension, dimension);
              factorPrecisionD1(0, 0) = 2.0 * rho / denominator;
              factorPrecisionD1(dimension - 1, dimension - 1) = 2.0 * rho / denominator;
              for(arma::uword row = 1; row + 1 < dimension; ++row){
                factorPrecisionD1(row, row) = 4.0 * rho / denominator;
              }
              for(arma::uword row = 0; row + 1 < dimension; ++row){
                factorPrecisionD1(row, row + 1) = -(1.0 + rho * rho) / denominator;
                factorPrecisionD1(row + 1, row) = factorPrecisionD1(row, row + 1);
              }
            }
          }

          if(!exactInnovationDerivative && !analyticSpatialDerivative && factorPrecisionD1.is_empty()){
            const auto derivativeStart = precisionDerivativeSeconds != nullptr
              ? std::chrono::steady_clock::now()
              : std::chrono::steady_clock::time_point();
            bool stationaryARDerivative = false;
            if(!stationaryARPrecision.is_empty()){
              const double step = 1.0e-6 *
                (1.0 + std::abs(state.factorPar[factorOffset](localK)));
              arma::vec plus = state.factorPar[factorOffset];
              arma::vec minus = state.factorPar[factorOffset];
              plus(localK) += step;
              minus(localK) -= step;
              arma::mat plusPrecision;
              arma::mat minusPrecision;
              double plusLogDet = 0.0;
              double minusLogDet = 0.0;
              if(buildStationaryARPrecision(f, plus, evalFactor(f, plus),
                    plusPrecision, plusLogDet) &&
                 buildStationaryARPrecision(f, minus, evalFactor(f, minus),
                    minusPrecision, minusLogDet)){
                const arma::uword capacity = plusPrecision.n_rows *
                  std::min(plusPrecision.n_rows, 2 * arBandwidth + 1);
                arma::umat locations(2, capacity);
                arma::vec values(capacity);
                arma::uword count = 0;
                for(arma::uword column = 0; column < plusPrecision.n_cols; ++column){
                  const arma::uword first = column > arBandwidth ? column - arBandwidth : 0;
                  const arma::uword last = std::min(plusPrecision.n_rows - 1,
                                                     column + arBandwidth);
                  for(arma::uword row = first; row <= last; ++row){
                    const double value = (plusPrecision(row, column) -
                      minusPrecision(row, column)) / (2.0 * step);
                    if(value == 0.0){ continue; }
                    locations(0, count) = row;
                    locations(1, count) = column;
                    values(count) = value;
                    ++count;
                  }
                }
                locations.resize(2, count);
                values.resize(count);
                sparseFactorPrecisionD1 = arma::sp_mat(
                  locations, values, plusPrecision.n_rows, plusPrecision.n_cols, true, true);
                stationaryARDerivative = true;
              }
            }
            if(!stationaryARDerivative &&
               !state.sparseFactorPrecision[factorOffset].is_empty()){
              factorPrecisionD1 = -state.sparseFactorPrecision[factorOffset] *
                factorCovarianceD1 * state.sparseFactorPrecision[factorOffset];
            }else if(!stationaryARDerivative){
              factorPrecisionD1 = -state.factorPrecision[factorOffset] *
                factorCovarianceD1 * state.factorPrecision[factorOffset];
            }
            if(precisionDerivativeSeconds != nullptr){
              *precisionDerivativeSeconds += std::chrono::duration<double>(
                std::chrono::steady_clock::now() - derivativeStart).count();
            }
          }

          const std::size_t globalK = static_cast<std::size_t>(k1 - 1);
          if(globalK >= state.localPrecisionD1.size() || derivativeAssigned[globalK]){
            Rcpp::stop("CovarianceFactor parameter ranges overlap or are invalid.");
          }
          if(!exactInnovationDerivative){
            state.localCovarianceD1[globalK] = factorCovarianceD1;
            if(!sparseFactorPrecisionD1.is_empty()){
              state.sparseLocalPrecisionD1[globalK] = std::move(sparseFactorPrecisionD1);
            }else{
              state.localPrecisionD1[globalK] = std::move(factorPrecisionD1);
            }
          }else{
            state.sparseLocalPrecisionD1[globalK] = std::move(sparseFactorPrecisionD1);
          }
          state.derivativeFactor[globalK] = fidx;
          if(exactInnovationDerivative){
            double trace = 0.0;
            const arma::mat & covariance = state.factorCovariance[factorOffset];
            for(arma::sp_mat::const_iterator entry = state.sparseLocalPrecisionD1[globalK].begin();
                entry != state.sparseLocalPrecisionD1[globalK].end(); ++entry){
              trace += (*entry) * covariance(entry.col(), entry.row());
            }
            state.logDetD1[globalK] = -
              (static_cast<double>(totalDim) / factorPrecisionDimension(state, factorOffset)) * trace;
          }else{
            state.logDetD1[globalK] = sommer::FactorPolicy::logDetMultiplicity(
              totalDim, factorPrecisionDimension(state, factorOffset)) *
              (analyticSpatialDerivative ? analyticLogDetD1 :
                factorPrecisionCovarianceTrace(state, factorOffset, factorCovarianceD1));
          }
          derivativeAssigned[globalK] = true;
        }
      }

      for(std::size_t k = 0; k < derivativeAssigned.size(); ++k){
        if(!derivativeAssigned[k] ||
           !std::isfinite(state.logDetD1[k])){
          Rcpp::stop("CovarianceFactor parameter ranges do not cover the descriptor parameters.");
        }
      }
    }

  arma::mat applyKroneckerRight(const arma::mat & input,
        const std::vector<arma::mat> & factors,
        const double scale,
        const std::vector<arma::sp_mat> * sparseFactors = nullptr) {

      arma::uword totalDim = 1;
      for(std::size_t factorIndex = 0; factorIndex < factors.size(); ++factorIndex){
        const arma::mat & factor = factors[factorIndex];
        const bool sparse = sparseFactors && !(*sparseFactors)[factorIndex].is_empty();
        const arma::uword rows = factor.is_empty()
          ? (sparse ? (*sparseFactors)[factorIndex].n_rows : 0)
          : factor.n_rows;
        const arma::uword columns = factor.is_empty()
          ? (sparse ? (*sparseFactors)[factorIndex].n_cols : 0)
          : factor.n_cols;
        if(rows != columns || rows == 0){
          Rcpp::stop("Kronecker application requires square factors.");
        }
        totalDim *= rows;
      }
      if(input.n_cols != totalDim){
        Rcpp::stop("Kronecker application received an incompatible input width.");
      }

      arma::mat current = input;
      arma::mat next(input.n_rows, input.n_cols);

      // Column-major: for a fixed leading index the (row, tail) x in slab is
      // contiguous, so each mode product is one GEMM per leading index.
      arma::uword trailing = totalDim;
      std::size_t axis = 0;
      for(std::size_t factorIndex = 0; factorIndex < factors.size(); ++factorIndex){
        const arma::mat & factor = factors[factorIndex];
        const bool sparse = sparseFactors && !(*sparseFactors)[factorIndex].is_empty();
        const arma::uword dim = factor.is_empty()
          ? (*sparseFactors)[factorIndex].n_rows : factor.n_rows;
        trailing /= dim;
        if(dim == 1){
          current *= factor.is_empty()
            ? (*sparseFactors)[factorIndex](0, 0) : factor(0, 0);
          ++axis;
          continue;
        }
        const arma::uword leading = totalDim / (dim * trailing);
        const arma::uword slabRows = current.n_rows * trailing;
        arma::mat factorT;
        arma::sp_mat sparseFactorT;
        if(sparse){ sparseFactorT = (*sparseFactors)[axis].t(); }
        else{ factorT = factor.t(); }
        for(arma::uword lead = 0; lead < leading; ++lead){
          const arma::uword offset = lead * slabRows * dim;
          const arma::mat inSlab(current.memptr() + offset, slabRows, dim, false, true);
          arma::mat outSlab(next.memptr() + offset, slabRows, dim, false, true);
          if(sparse){ outSlab = inSlab * sparseFactorT; }
          else{ outSlab = inSlab * factorT; }
        }
        current.swap(next);
        ++axis;
      }

      current *= scale;
      return current;
    }

  void descriptorPrecisionD1Factors(const DescriptorPrecisionState & state,
        const arma::uword k,
        std::vector<arma::mat> & factors,
        double & scale) {

      if(!state.factorWise || k >= state.derivativeFactor.size()){
        Rcpp::stop("Factor-wise precision derivative is unavailable.");
      }

      factors = state.factorPrecision;
      for(std::size_t factor = 0; factor < factors.size(); ++factor){
        if(factors[factor].is_empty()){
          factors[factor] = factorPrecisionDense(state, factor);
        }
      }
      scale = state.inverseScale;
      const int derivativeFactor = state.derivativeFactor[k];
      if(derivativeFactor < 0){
        scale = -scale;
      }else{
        factors[static_cast<std::size_t>(derivativeFactor)] =
          localPrecisionD1Dense(state, k);
      }
    }

  void descriptorCovarianceD1Factors(const DescriptorPrecisionState & state,
        const arma::uword k,
        std::vector<arma::mat> & factors,
        double & scale) {

      if(!state.factorWise || k >= state.derivativeFactor.size()){
        Rcpp::stop("Factor-wise covariance derivative is unavailable.");
      }

      factors = state.factorCovariance;
      scale = 1.0 / state.inverseScale;
      const int derivativeFactor = state.derivativeFactor[k];
      if(derivativeFactor >= 0){
        factors[static_cast<std::size_t>(derivativeFactor)] =
          localCovarianceD1Dense(state, k);
      }
    }

  double contractKronecker(const arma::mat & matrix,
        const std::vector<arma::mat> & factors,
        const double scale) {

      arma::uword totalDim = 1;
      std::vector<arma::uword> dimensions;
      double supportEntries = 1.0;
      for(const arma::mat & factor : factors){
        if(factor.n_rows == 0 || factor.n_cols != factor.n_rows ||
           totalDim > matrix.n_rows / factor.n_rows){
          Rcpp::stop("Kronecker contraction factor dimensions are incompatible.");
        }
        totalDim *= factor.n_rows;
        dimensions.push_back(factor.n_rows);
        const arma::uword nonzeros = arma::accu(factor != 0.0);
        if(nonzeros == 0){ return 0.0; }
        supportEntries *= static_cast<double>(nonzeros);
      }
      if(matrix.n_rows != totalDim || matrix.n_cols != totalDim){
        Rcpp::stop("Kronecker contraction received incompatible dimensions.");
      }

      if(FactorPolicy::supportContraction(supportEntries, static_cast<double>(totalDim),
                                         contractionDigitsMaxEntries)){
        struct FactorEntry {
          arma::uword row;
          arma::uword column;
          double value;
        };
        std::vector<std::vector<FactorEntry>> support(factors.size());
        for(std::size_t axis = 0; axis < factors.size(); ++axis){
          const arma::mat & factor = factors[axis];
          for(arma::uword row = 0; row < factor.n_rows; ++row){
            for(arma::uword column = 0; column < factor.n_cols; ++column){
              if(factor(row, column) != 0.0){
                support[axis].push_back({row, column, factor(row, column)});
              }
            }
          }
        }
        std::function<double(std::size_t, arma::uword, arma::uword, double)> accumulate;
        accumulate = [&](std::size_t axis, arma::uword row, arma::uword column, double value) -> double {
          if(axis == support.size()){ return matrix(row, column) * value; }
          double sum = 0.0;
          for(const FactorEntry & entry : support[axis]){
            sum += accumulate(axis + 1, row * dimensions[axis] + entry.row,
              column * dimensions[axis] + entry.column, value * entry.value);
          }
          return sum;
        };
        return accumulate(0, 0, 0, scale);
      }

      const arma::umat * cachedDigits = nullptr;
      {
        std::lock_guard<std::mutex> lock(contractionDigitsMutex);
        auto cached = contractionDigitsCache.find(dimensions);
        if(cached == contractionDigitsCache.end() && !factors.empty() &&
           totalDim <= (contractionDigitsMaxEntries - contractionDigitsEntries) / factors.size()){
          arma::umat digits(totalDim, factors.size());
          arma::uword trailing = totalDim;
          for(std::size_t axis = 0; axis < factors.size(); ++axis){
            trailing /= dimensions[axis];
            for(arma::uword coordinate = 0; coordinate < totalDim; ++coordinate){
              digits(coordinate, axis) = (coordinate / trailing) % dimensions[axis];
            }
          }
          cached = contractionDigitsCache.emplace(dimensions, std::move(digits)).first;
          contractionDigitsEntries += cached->second.n_elem;
        }
        if(cached != contractionDigitsCache.end()){
          cachedDigits = &cached->second;
        }
      }
      // Entries are immutable and never erased; map insertion preserves their addresses.
      double result = 0.0;
      for(arma::uword row = 0; row < totalDim; ++row){
        for(arma::uword col = 0; col < totalDim; ++col){
          double value = scale;
          arma::uword trailing = totalDim;
          for(std::size_t axis = 0; axis < factors.size(); ++axis){
            const arma::mat & factor = factors[axis];
            const arma::uword dim = factor.n_rows;
            trailing /= dim;
            const arma::uword factorRow = cachedDigits == nullptr
              ? (row / trailing) % dim : (*cachedDigits)(row, axis);
            const arma::uword factorCol = cachedDigits == nullptr
              ? (col / trailing) % dim : (*cachedDigits)(col, axis);
            value *= factor(factorRow, factorCol);
          }
          result += matrix(row,col) * value;
        }
      }
      return result;
    }

  arma::vec kroneckerDiagonal(const std::vector<arma::mat> & factors,
        const double scale) {

      arma::vec diagonal(1, arma::fill::ones);
      for(const arma::mat & factor : factors){
        diagonal = arma::kron(diagonal, factor.diag());
      }
      diagonal *= scale;
      return diagonal;
    }

  arma::mat applyDescriptorPrecisionD1Right(const arma::mat & input,
        const DescriptorPrecisionState & state,
        const arma::uword k) {

      if(!state.factorWise || k >= state.derivativeFactor.size()){
        Rcpp::stop("Factor-wise precision derivative is unavailable.");
      }

      std::vector<arma::mat> factors = state.factorPrecision;
      double scale = state.inverseScale;
      std::vector<arma::sp_mat> sparseFactors = state.sparseFactorPrecision;
      const int derivativeAxis = state.derivativeFactor[k];
      if(derivativeAxis < 0){
        scale = -scale;
      }else{
        const std::size_t axis = static_cast<std::size_t>(derivativeAxis);
        if(k < state.sparseLocalPrecisionD1.size() && !state.sparseLocalPrecisionD1[k].is_empty()){
          factors[axis].reset();
          sparseFactors[axis] = state.sparseLocalPrecisionD1[k];
        }else{
          factors[axis] = localPrecisionD1Dense(state, k);
          sparseFactors[axis].reset();
        }
      }
      return applyKroneckerRight(input, factors, scale, &sparseFactors);
    }

  arma::mat applyDescriptorCovarianceD1Right(const arma::mat & input,
        const DescriptorPrecisionState & state,
        const arma::uword k) {

      std::vector<arma::mat> factors;
      double scale = 1.0;
      descriptorCovarianceD1Factors(state, k, factors, scale);
      return applyKroneckerRight(input, factors, scale);
    }
};

}

#endif
