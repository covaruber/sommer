#include "mme_direct_residual.h"

namespace {
void directInferenceMemory(const double elements, const double budget){
  if(!std::isfinite(elements) || elements * 8e-6 > budget){
    Rcpp::stop("Direct post-fit covariance workspace exceeds denseMemoryMB/maxMemoryMB; request fewer contrasts or increase the memory budget.");
  }
}

arma::mat directObservationSolve(const arma::mat & upper, const arma::mat & rhs){
  arma::mat intermediate, result;
  if(!arma::solve(intermediate, arma::trimatl(upper.t()), rhs, arma::solve_opts::fast) ||
     !arma::solve(result, arma::trimatu(upper), intermediate, arma::solve_opts::fast)){
    Rcpp::stop("Direct post-fit observation covariance solve failed.");
  }
  return result;
}
}

// [[Rcpp::export]]
Rcpp::List direct_mmes_vcov_factor_cpp(const Rcpp::List & inputs,
                                      const double denseMemoryMB){
  const arma::sp_mat X = Rcpp::as<arma::sp_mat>(inputs["X"]);
  const Rcpp::List Z = inputs["Z"], precision = inputs["Ai"],
    descriptors = inputs["covStruct"], parameters = inputs["parameters"],
    partitions = inputs["partitions"];
  const arma::vec owners = inputs["Zind"];
  const arma::uword n = X.n_rows, p = X.n_cols;
  const int terms = partitions.size();
  if(descriptors.size() != terms + 1 || parameters.size() != terms + 1 ||
     precision.size() < terms || owners.n_elem != static_cast<arma::uword>(Z.size())){
    Rcpp::stop("Direct post-fit covariance inputs are inconsistent; refit the model.");
  }
  double relationshipElements = 0, largestLevels = 0, covarianceElements = 0;
  for(int term = 0; term < terms; ++term){
    const arma::sp_mat A = convertSparse(precision[term]);
    if(!isIdentity_spmat(A)){ relationshipElements += static_cast<double>(A.n_rows) * A.n_cols; }
    largestLevels = std::max(largestLevels, static_cast<double>(A.n_rows));
    const Rcpp::List descriptor = descriptors[term];
    const double dimension = Rcpp::as<int>(descriptor["dim"]);
    covarianceElements += dimension * dimension;
  }
  const double workspace = 10.0*n*n + 5.0*p*p + 5.0*n*p +
    4.0*n*largestLevels + 4.0*relationshipElements + 3.0*covarianceElements;
  directInferenceMemory(workspace, denseMemoryMB);
  const DirectResidualCovariance residual(descriptors[terms],
    Rcpp::as<arma::uvec>(inputs["residualBlock"]),
    Rcpp::as<arma::uvec>(inputs["residualIndex"]),
    Rcpp::as<arma::sp_mat>(inputs["H"]), Rcpp::as<bool>(inputs["useH"]));
  arma::mat V = residual.evaluate(Rcpp::as<arma::vec>(parameters[terms]));
  Rcpp::List blocks(terms), relationships(terms), covariances(terms);
  arma::uword effects = p;
  for(int term = 0; term < terms; ++term){
    const arma::uvec indices = arma::find(owners == term + 1);
    const arma::mat Sigma = directEvaluateDescriptor(descriptors[term],
      Rcpp::as<arma::vec>(parameters[term]));
    if(Sigma.n_rows != indices.n_elem){ Rcpp::stop("Direct post-fit random design dimension mismatch."); }
    const arma::sp_mat A = convertSparse(precision[term]);
    const bool identity = isIdentity_spmat(A);
    arma::mat K;
    if(!identity){
      const arma::mat denseA(A);
      if(!eigenSpdInverse(denseA, K) &&
         !eigenSpdInverse(nearPDcpp(denseA, 100, 1e-6, 1e-7), K)){
        Rcpp::stop("Direct post-fit relationship precision inversion failed.");
      }
      relationships[term] = K;
    }
    Rcpp::List termBlocks(indices.n_elem);
    std::vector<arma::sp_mat> sparseDesigns(indices.n_elem);
    for(arma::uword a = 0; a < indices.n_elem; ++a){
      termBlocks[a] = Z[indices(a)];
      sparseDesigns[a] = Rcpp::as<arma::sp_mat>(Z[indices(a)]);
      const arma::sp_mat & Za = sparseDesigns[a];
      if(Za.n_rows != n || Za.n_cols != A.n_rows){
        Rcpp::stop("Direct post-fit relationship and design dimensions differ.");
      }
      effects += Za.n_cols;
    }
    for(arma::uword a = 0; a < indices.n_elem; ++a){
      const arma::sp_mat & Za = sparseDesigns[a];
      const arma::mat ZK = identity ? arma::mat(Za) : arma::mat(Za * K);
      for(arma::uword b = 0; b < indices.n_elem; ++b){
        if(Sigma(a,b) != 0.0){
          const arma::sp_mat & Zb = sparseDesigns[b];
          V += Sigma(a,b) * (ZK * Zb.t());
        }
      }
    }
    blocks[term] = termBlocks;
    covariances[term] = Sigma;
  }
  arma::mat upper;
  const double observationJitter = Rcpp::as<double>(inputs["observationJitter"]);
  const double fixedInformationJitter = Rcpp::as<double>(inputs["fixedInformationJitter"]);
  if(!std::isfinite(observationJitter) || observationJitter < 0 ||
     !std::isfinite(fixedInformationJitter) || fixedInformationJitter < 0){
    Rcpp::stop("Direct post-fit regularization metadata is invalid; refit the model.");
  }
  V.diag() += observationJitter;
  if(!arma::chol(upper, arma::symmatu(V))){
    Rcpp::stop("Direct post-fit observation covariance is not positive definite.");
  }
  const arma::mat VX = directObservationSolve(upper, arma::mat(X));
  arma::mat fixedCovariance;
  arma::mat information(X.t() * VX);
  information.diag() += fixedInformationJitter;
  if(!arma::inv_sympd(fixedCovariance, information)){
    Rcpp::stop("Direct post-fit fixed-effect information is singular.");
  }
  return Rcpp::List::create(Rcpp::Named("upper")=upper,
    Rcpp::Named("VX")=VX, Rcpp::Named("fixedCovariance")=fixedCovariance,
    Rcpp::Named("blocks")=blocks, Rcpp::Named("relationships")=relationships,
    Rcpp::Named("covariances")=covariances, Rcpp::Named("partitions")=partitions,
    Rcpp::Named("fixedColumns")=p, Rcpp::Named("effects")=effects,
    Rcpp::Named("workspace")=workspace, Rcpp::Named("budget")=denseMemoryMB,
    Rcpp::Named("largestLevels")=largestLevels);
}

// [[Rcpp::export]]
arma::mat direct_mmes_vcov_solve_cpp(const Rcpp::List & factor,
                                    const arma::sp_mat & D,
                                    const bool priorOnly = false){
  const arma::uword p = Rcpp::as<arma::uword>(factor["fixedColumns"]);
  const arma::uword effects = Rcpp::as<arma::uword>(factor["effects"]);
  const arma::mat upper = factor["upper"], VX = factor["VX"],
    fixedCovariance = factor["fixedCovariance"];
  const arma::uword k = D.n_rows, n = upper.n_rows;
  if(D.n_cols != effects || !D.is_finite()){
    Rcpp::stop("Prediction contrasts must be finite with one column per fitted coefficient.");
  }
  directInferenceMemory(Rcpp::as<double>(factor["workspace"]) +
    8.0*k*k + 6.0*k*n + 4.0*k*Rcpp::as<double>(factor["largestLevels"]),
    Rcpp::as<double>(factor["budget"]));
  if(!k){ return arma::mat(0,0); }
  const Rcpp::List blocks = factor["blocks"], relationships = factor["relationships"],
    covariances = factor["covariances"], partitions = factor["partitions"];
  arma::mat cross(k, n, arma::fill::zeros), prior(k, k, arma::fill::zeros);
  for(int term = 0; term < blocks.size(); ++term){
    const Rcpp::List designs = blocks[term];
    const arma::mat ranges = partitions[term], Sigma = covariances[term];
    const bool identity = Rf_isNull(relationships[term]);
    arma::mat K;
    if(!identity){ K = Rcpp::as<arma::mat>(relationships[term]); }
    for(int a = 0; a < designs.size(); ++a){
      const arma::mat Da(D.cols(ranges(a,0)-1, ranges(a,1)-1));
      const arma::mat DK = identity ? Da : arma::mat(Da * K);
      for(int b = 0; b < designs.size(); ++b){
        if(Sigma(a,b) == 0.0){ continue; }
        const arma::sp_mat Zb = Rcpp::as<arma::sp_mat>(designs[b]);
        const arma::sp_mat Db = D.cols(ranges(b,0)-1, ranges(b,1)-1);
        cross += Sigma(a,b) * (DK * Zb.t());
        prior += Sigma(a,b) * (DK * Db.t());
      }
    }
  }
  if(priorOnly){ return 0.5 * (prior + prior.t()); }
  // Error covariance of L beta_hat + M u_hat - (L beta + M u).
  const arma::mat adjusted = arma::mat(D.cols(0,p-1)) - cross * VX;
  const arma::mat result = prior - cross * directObservationSolve(upper, cross.t()) +
    adjusted * fixedCovariance * adjusted.t();
  if(!result.is_finite()){ Rcpp::stop("Direct post-fit covariance is non-finite."); }
  return 0.5 * (result + result.t());
}
