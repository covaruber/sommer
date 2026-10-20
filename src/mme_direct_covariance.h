#ifndef SOMMER_MME_DIRECT_COVARIANCE_H
#define SOMMER_MME_DIRECT_COVARIANCE_H

#include "mme_sp2_support.h"

inline arma::mat directEvalNativeFactor(const Rcpp::List & f,
                                        const arma::vec & localPar,
                                        const std::string & op){

  const arma::uword q =
    static_cast<arma::uword>(Rcpp::as<int>(f["dim"]));

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
    const double lo = -1.0 / static_cast<double>(q - 1);
    const double eta = localPar(0);
    const double s = eta >= 0.0 ? 1.0/(1.0+std::exp(-eta)) : std::exp(eta)/(1.0+std::exp(eta));
    const double rho = lo + (1.0-lo)*s;
    arma::mat K(q, q, arma::fill::value(rho));
    K.diag().ones();
    return K;
  }

  if(op == "corh"){
    if(q < 2 || localPar.n_elem != q){
      Rcpp::stop("Malformed heterogeneous uniform-correlation factor.");
    }
    const double lo = -1.0 / static_cast<double>(q - 1);
    const double eta = localPar(0);
    const double s = eta >= 0.0 ? 1.0/(1.0+std::exp(-eta)) : std::exp(eta)/(1.0+std::exp(eta));
    const double rho = lo + (1.0-lo)*s;
    arma::mat C(q, q, arma::fill::value(rho));
    C.diag().ones();
    arma::vec variances(q, arma::fill::ones);
    for(arma::uword k = 1; k < q; ++k){
      variances(k) = std::exp(localPar(k));
    }
    arma::vec sd = arma::sqrt(variances);
    return arma::diagmat(sd) * C * arma::diagmat(sd);
  }

  if(op == "arp"){
    const int order = Rcpp::as<int>(f["order"]);
    if(order < 1 || localPar.n_elem != static_cast<arma::uword>(order) ||
       q <= static_cast<arma::uword>(order)){
      Rcpp::stop("Malformed AR(p) covariance factor.");
    }
    arma::vec pacf(static_cast<arma::uword>(order), arma::fill::zeros);
    for(int j = 0; j < order; ++j){
      pacf(static_cast<arma::uword>(j)) = std::tanh(localPar(static_cast<arma::uword>(j)));
    }
    arma::vec phi;
    for(int m = 1; m <= order; ++m){
      arma::vec next(static_cast<arma::uword>(m), arma::fill::zeros);
      next(static_cast<arma::uword>(m-1)) = pacf(static_cast<arma::uword>(m-1));
      if(m > 1){
        for(int j = 0; j < m-1; ++j){
          next(static_cast<arma::uword>(j)) =
            phi(static_cast<arma::uword>(j)) -
            pacf(static_cast<arma::uword>(m-1)) * phi(static_cast<arma::uword>(m-2-j));
        }
      }
      phi = next;
    }
    arma::mat A(static_cast<arma::uword>(order), static_cast<arma::uword>(order), arma::fill::zeros);
    arma::vec b(static_cast<arma::uword>(order), arma::fill::zeros);
    for(int kk = 1; kk <= order; ++kk){
      A(static_cast<arma::uword>(kk-1), static_cast<arma::uword>(kk-1)) += 1.0;
      for(int jj = 1; jj <= order; ++jj){
        const int d = std::abs(kk - jj);
        const double pj = phi(static_cast<arma::uword>(jj-1));
        if(d == 0){
          b(static_cast<arma::uword>(kk-1)) += pj;
        }else{
          A(static_cast<arma::uword>(kk-1), static_cast<arma::uword>(d-1)) -= pj;
        }
      }
    }
    arma::vec rInitial;
    bool ok = arma::solve(rInitial, A, b);
    if(!ok || !rInitial.is_finite()){
      Rcpp::stop("Unable to solve Yule-Walker equations for AR(p) factor.");
    }
    arma::vec rho(q, arma::fill::zeros);
    rho(0) = 1.0;
    for(int h = 1; h <= order; ++h){
      rho(static_cast<arma::uword>(h)) = rInitial(static_cast<arma::uword>(h-1));
    }
    for(arma::uword h = static_cast<arma::uword>(order+1); h < q; ++h){
      double value = 0.0;
      for(int jj = 1; jj <= order; ++jj){
        value += phi(static_cast<arma::uword>(jj-1)) * rho(h - static_cast<arma::uword>(jj));
      }
      rho(h) = value;
    }
    arma::mat K(q, q, arma::fill::zeros);
    for(arma::uword i = 0; i < q; ++i){
      for(arma::uword j = 0; j < q; ++j){
        const arma::uword d = i > j ? i-j : j-i;
        K(i,j) = rho(d);
      }
    }
    return K;
  }

  if(op == "ma"){
    const int order = Rcpp::as<int>(f["order"]);
    if(order < 1 || localPar.n_elem != static_cast<arma::uword>(order) ||
       q <= static_cast<arma::uword>(order)){
      Rcpp::stop("Malformed MA(q) covariance factor.");
    }
    arma::vec coef(static_cast<arma::uword>(order + 1), arma::fill::zeros);
    coef(0) = 1.0;
    for(int j = 1; j <= order; ++j){
      coef(static_cast<arma::uword>(j)) = localPar(static_cast<arma::uword>(j-1));
    }
    arma::vec rho(q, arma::fill::zeros);
    double gamma0 = 0.0;
    for(int j = 0; j <= order; ++j){
      const double c = coef(static_cast<arma::uword>(j));
      gamma0 += c*c;
    }
    if(!std::isfinite(gamma0) || gamma0 <= 0.0){
      Rcpp::stop("Invalid MA covariance normalization.");
    }
    rho(0) = 1.0;
    for(int h = 1; h <= order; ++h){
      double gamma = 0.0;
      for(int j = 0; j <= order-h; ++j){
        gamma += coef(static_cast<arma::uword>(j)) * coef(static_cast<arma::uword>(j+h));
      }
      rho(static_cast<arma::uword>(h)) = gamma / gamma0;
    }
    arma::mat K(q, q, arma::fill::zeros);
    for(arma::uword i = 0; i < q; ++i){
      for(arma::uword j = 0; j < q; ++j){
        const arma::uword d = i > j ? i-j : j-i;
        K(i,j) = d <= static_cast<arma::uword>(order) ? rho(d) : 0.0;
      }
    }
    return K;
  }

  if(op == "corg"){
    Rcpp::IntegerVector rr = f["corg_row"];
    Rcpp::IntegerVector cc = f["corg_col"];
    if(localPar.n_elem != static_cast<arma::uword>(rr.size()) || rr.size() != cc.size()){
      Rcpp::stop("Malformed general-correlation factor.");
    }
    arma::mat A(q, q, arma::fill::eye);
    for(arma::uword k = 0; k < localPar.n_elem; ++k){
      const arma::uword i = static_cast<arma::uword>(rr[static_cast<int>(k)] - 1);
      const arma::uword j = static_cast<arma::uword>(cc[static_cast<int>(k)] - 1);
      A(i,j) = localPar(k);
    }
    arma::mat S = A * A.t();
    arma::vec sd = arma::sqrt(S.diag());
    arma::mat denom = sd * sd.t();
    arma::mat K = S % arma::pow(denom, -1.0);
    K.diag().ones();
    return 0.5 * (K + K.t());
  }

  if(op == "fa" || op == "rr"){
    return sommer::loadingFactorCovariance(f, localPar);
  }

  if(op == "ante"){
    const int ncoef = Rcpp::as<int>(f["ante_ncoef"]);
    Rcpp::IntegerVector rr = f["ante_row"];
    Rcpp::IntegerVector cc = f["ante_col"];
    if(ncoef < 1 || rr.size() != ncoef || cc.size() != ncoef ||
       localPar.n_elem != static_cast<arma::uword>(ncoef + static_cast<int>(q) - 1)){
      Rcpp::stop("Malformed antedependence covariance factor.");
    }
    return anteInnovationCovariance(q, rr, cc, localPar, sparseFactorKernelsEnabled(f));
  }

  Rcpp::stop("Unsupported native covariance evaluator opcode: " + op);
  return arma::mat();
}

inline arma::mat directEvalFactor(const Rcpp::List & f, const arma::vec & localPar){

  if(!f.containsElementNamed("evaluator")){
    Rcpp::stop("CovarianceFactor is missing evaluator specification.");
  }
  Rcpp::List spec = Rcpp::as<Rcpp::List>(f["evaluator"]);
  const std::string backend = Rcpp::as<std::string>(spec["backend"]);
  const arma::uword q = static_cast<arma::uword>(Rcpp::as<int>(f["dim"]));

  arma::mat K;

  if(sparseFactorKernelsEnabled(f) && f.containsElementNamed("model") &&
     (Rcpp::as<std::string>(f["model"]) == "sar" || Rcpp::as<std::string>(f["model"]) == "car")){
    K = evaluateSparseSpatialFactor(f, localPar, false).covariance;
  }else if(backend == "native"){
    if(!spec.containsElementNamed("op")){
      Rcpp::stop("Native CovarianceFactor evaluator is missing op.");
    }
    K = directEvalNativeFactor(f, localPar, Rcpp::as<std::string>(spec["op"]));
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
  return 0.5 * (K + K.t());
}

// Central finite-difference derivative. Covariance-shaping factors are
// always small (levels/traits/lags, never records or coefficients), so
// re-evaluating the factor twice per free parameter is negligible.
inline arma::mat directFactorD1(const Rcpp::List & f, const arma::vec & localPar,
                                const arma::uword k, const double relStep = 1.0e-6){

  if(k >= localPar.n_elem){
    Rcpp::stop("Invalid CovarianceFactor derivative parameter index.");
  }
  if(!f.containsElementNamed("derivative")){
    Rcpp::stop("CovarianceFactor is missing derivative specification.");
  }
  if(sparseFactorKernelsEnabled(f) && f.containsElementNamed("model") &&
     (Rcpp::as<std::string>(f["model"]) == "sar" || Rcpp::as<std::string>(f["model"]) == "car")){
    return evaluateSparseSpatialFactor(f, localPar, true).covarianceD1;
  }
  Rcpp::List spec = Rcpp::as<Rcpp::List>(f["derivative"]);
  const std::string backend = Rcpp::as<std::string>(spec["backend"]);

  if(backend == "none"){
    Rcpp::stop("Derivative requested for a CovarianceFactor with no parameters.");
  }
  if(backend == "R"){
    if(!spec.containsElementNamed("fun") || Rf_isNull(spec["fun"])){
      Rcpp::stop("R derivative specification is missing fun(par,k).");
    }
    Rcpp::Function dfun = spec["fun"];
    SEXP ans = dfun(Rcpp::wrap(localPar), static_cast<int>(k + 1));
    arma::mat D = Rcpp::as<arma::mat>(ans);
    return 0.5 * (D + D.t());
  }

  const double h = relStep * (1.0 + std::abs(localPar(k)));
  arma::vec plus = localPar;
  arma::vec minus = localPar;
  plus(k) += h;
  minus(k) -= h;
  return (directEvalFactor(f, plus) - directEvalFactor(f, minus)) / (2.0*h);
}

// Residual structures with correlated or section-specific factors are assembled
// directly in record space instead of over the full covariance-product space.
inline bool directResidualCompact(const Rcpp::List & cs){
  Rcpp::List factors = cs["factors"];
  for(int fidx = 0; fidx < factors.size(); ++fidx){
    Rcpp::List f = Rcpp::as<Rcpp::List>(factors[fidx]);
    if(f.containsElementNamed("section_levels") ||
       !f.containsElementNamed("structurally_diagonal") ||
       !Rcpp::as<bool>(f["structurally_diagonal"])){
      return true;
    }
  }
  return false;
}

inline int directSectionListed(const Rcpp::List & f, const int t){
  if(f.containsElementNamed("section_index")){
    Rcpp::IntegerVector index = f["section_index"];
    return index[t];
  }
  return t;
}

inline arma::vec directSectionSlice(const Rcpp::List & f, const arma::vec & localPar, const int listed){
  const arma::uword npar = static_cast<arma::uword>(Rcpp::as<int>(f["section_npar"]));
  if(npar == 0 || listed < 0){
    return arma::vec();
  }
  const arma::uword first = static_cast<arma::uword>(listed) * npar;
  return localPar.subvec(first, first + npar - 1);
}

// Kronecker PD: sigma2 * prod_f min eig(F_f), taking the worst section per factor.
inline bool directDescriptorPD(const Rcpp::List & cs, const arma::vec & par, const double floor){
  if(par.n_elem < 1 || !par.is_finite()){
    return false;
  }
  Rcpp::List factors = cs["factors"];
  double logMinimum = par(0);
  for(int fidx = 0; fidx < factors.size(); ++fidx){
    Rcpp::List f = Rcpp::as<Rcpp::List>(factors[fidx]);
    const int start1 = f.containsElementNamed("par_start") ? Rcpp::as<int>(f["par_start"]) : 1;
    const int end1 = f.containsElementNamed("par_end") ? Rcpp::as<int>(f["par_end"]) : 0;
    arma::vec localPar;
    if(end1 >= start1){
      localPar = par.subvec(static_cast<arma::uword>(start1-1), static_cast<arma::uword>(end1-1));
    }
    const int nSections = f.containsElementNamed("section_levels") ? Rf_length(f["section_levels"]) : 0;
    double factorMinimum = 1.0;
    bool any = false;
    for(int t = 0; t < std::max(nSections, 1); ++t){
      arma::mat K;
      try{
        if(nSections > 0){
          const int listed = directSectionListed(f, t);
          if(listed < 0){
            continue;
          }
          K = directEvalFactor(f, directSectionSlice(f, localPar, listed));
        }else{
          K = directEvalFactor(f, localPar);
        }
      }catch(...){
        return false;
      }
      arma::vec eigenvalues;
      if(!arma::eig_sym(eigenvalues, K) || !eigenvalues.is_finite() || eigenvalues.min() <= 0.0){
        return false;
      }
      factorMinimum = any ? std::min(factorMinimum, eigenvalues.min()) : eigenvalues.min();
      any = true;
    }
    logMinimum += std::log(factorMinimum);
  }
  return std::isfinite(logMinimum) && logMinimum > std::log(floor);
}

inline arma::mat directEvaluateDescriptor(const Rcpp::List & cs, const arma::vec & par){

  if(par.n_elem < 1){
    Rcpp::stop("Covariance descriptor must contain log_sigma2.");
  }
  Rcpp::List factors = cs["factors"];
  arma::mat K(1,1,arma::fill::ones);

  for(int fidx = 0; fidx < factors.size(); ++fidx){
    Rcpp::List f = Rcpp::as<Rcpp::List>(factors[fidx]);
    if(f.containsElementNamed("section_levels")){
      Rcpp::stop("Internal error: dsumm() structures must be assembled in record space.");
    }
    const int start1 = f.containsElementNamed("par_start") ? Rcpp::as<int>(f["par_start"]) : 1;
    const int end1 = f.containsElementNamed("par_end") ? Rcpp::as<int>(f["par_end"]) : 0;
    arma::vec localPar;
    if(end1 >= start1){
      localPar = par.subvec(static_cast<arma::uword>(start1-1), static_cast<arma::uword>(end1-1));
    }
    K = arma::kron(K, directEvalFactor(f, localPar));
  }
  return std::exp(par(0)) * K;
}


inline arma::mat directDescriptorD1(const Rcpp::List & cs, const arma::vec & par,
                                    const arma::uword k){

  if(k == 0){
    return directEvaluateDescriptor(cs, par);
  }
  Rcpp::List factors = cs["factors"];
  arma::mat K(1,1,arma::fill::ones);
  bool found = false;

  for(int fidx = 0; fidx < factors.size(); ++fidx){
    Rcpp::List f = Rcpp::as<Rcpp::List>(factors[fidx]);
    const int start1 = f.containsElementNamed("par_start") ? Rcpp::as<int>(f["par_start"]) : 1;
    const int end1 = f.containsElementNamed("par_end") ? Rcpp::as<int>(f["par_end"]) : 0;
    arma::vec localPar;
    if(end1 >= start1){
      localPar = par.subvec(static_cast<arma::uword>(start1-1), static_cast<arma::uword>(end1-1));
    }
    arma::mat piece;
    const int k1 = static_cast<int>(k) + 1;
    if(end1 >= start1 && k1 >= start1 && k1 <= end1){
      piece = directFactorD1(f, localPar, static_cast<arma::uword>(k1-start1));
      found = true;
    }else{
      piece = directEvalFactor(f, localPar);
    }
    K = arma::kron(K, piece);
  }
  if(!found){
    Rcpp::stop("Covariance derivative parameter does not belong to any factor.");
  }
  return std::exp(par(0)) * K;
}

// Transforms unconstrained working parameters (log_sigma2, atanh(rho), ...)
// to reported natural-scale values, with the elementwise Jacobian used for
// the delta-method transform of the parameter covariance matrix.
inline arma::vec directReportParameters(const Rcpp::List & cs, const arma::vec & work,
                                        const double vary, arma::vec & jacobianOut){

  arma::vec out = work;
  arma::vec jac(work.n_elem, arma::fill::ones);

  out(0) = std::exp(work(0)) * vary;
  jac(0) = out(0);

  Rcpp::List factors = cs["factors"];
  for(int fidx = 0; fidx < factors.size(); ++fidx){
    Rcpp::List f = Rcpp::as<Rcpp::List>(factors[fidx]);
    const int start1 = Rcpp::as<int>(f["par_start"]);
    const int end1 = Rcpp::as<int>(f["par_end"]);
    if(end1 < start1){ continue; }
    if(!f.containsElementNamed("report")){
      Rcpp::stop("CovarianceFactor is missing report specification.");
    }
    Rcpp::List report = Rcpp::as<Rcpp::List>(f["report"]);
    const std::string backend = Rcpp::as<std::string>(report["backend"]);
    const int nFactorPar = end1 - start1 + 1;

    if(backend != "builtin"){
      Rcpp::stop("Unknown CovarianceFactor reporting backend: " + backend);
    }
    Rcpp::CharacterVector tr = report["transform"];
    Rcpp::NumericVector lo = report["lower"];
    Rcpp::NumericVector hi = report["upper"];
    if(tr.size() != nFactorPar || lo.size() != nFactorPar || hi.size() != nFactorPar){
      Rcpp::stop("CovarianceFactor report specification has incompatible length.");
    }
    for(int local = 0; local < nFactorPar; ++local){
      const arma::uword k = static_cast<arma::uword>(start1 - 1 + local);
      const std::string code = Rcpp::as<std::string>(tr[local]);
      const double eta = work(k);
      if(code == "identity"){
        out(k) = eta; jac(k) = 1.0;
      }else if(code == "exp"){
        out(k) = std::exp(eta); jac(k) = out(k);
      }else if(code == "tanh"){
        out(k) = std::tanh(eta); jac(k) = 1.0 - out(k)*out(k);
      }else if(code == "bounded_logit"){
        const double lower = lo[local], upper = hi[local];
        if(!std::isfinite(lower) || !std::isfinite(upper) || !(lower < upper)){
          Rcpp::stop("Invalid bounded-logit reporting interval.");
        }
        const double logistic = eta >= 0.0 ? 1.0/(1.0+std::exp(-eta)) : std::exp(eta)/(1.0+std::exp(eta));
        out(k) = lower + (upper-lower)*logistic;
        jac(k) = (upper-lower)*logistic*(1.0-logistic);
      }else{
        Rcpp::stop("Unknown CovarianceFactor reporting transform: " + code);
      }
    }
  }
  jacobianOut = jac;
  return out;
}


#endif
