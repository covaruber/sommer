#ifndef SOMMER_MME_GLMM_SYSTEM_H
#define SOMMER_MME_GLMM_SYSTEM_H

#include "mme_glmm_precision.h"
#include "mme_glmm_acceleration.h"
#include <memory>

namespace sommer {

struct LaplaceSystemPlan {
  std::string requested;
  std::string selected;
  std::string reason;
  double workspace;
};

inline LaplaceSystemPlan planLaplaceSystem(const std::string & requested,
    const arma::sp_mat & design, const arma::sp_mat & precision, double memoryMB){
  if(requested != "auto" && requested != "random" && requested != "observation")
    Rcpp::stop("Laplace backend must be 'auto', 'random', or 'observation'.");
  if(!std::isfinite(memoryMB) || memoryMB <= 0)
    Rcpp::stop("Laplace system memory budget must be positive and finite.");
  const double n = design.n_rows, q = design.n_cols;
  const double workspace = 8.0 * (20.0*n*n + 16.0*q +
    8.0*design.n_nonzero + 8.0*precision.n_nonzero);
  const bool fits = workspace + 8.0*(6.0*q+8.0*n) <= memoryMB*1048576.0;
  if(requested == "observation" && !fits)
    Rcpp::stop("Observation-space Laplace workspace exceeds execution$denseMemoryMB.");
  if(requested != "auto")
    return {requested, requested, "explicit backend request", workspace};
  // Dense row coupling, rather than q alone, predicts expensive random-space fill.
  arma::vec support(design.n_rows, arma::fill::zeros);
  for(auto cell = design.begin(); cell != design.end(); ++cell) ++support(cell.row());
  const double pairs = precision.n_nonzero + arma::dot(support,support);
  const bool coupled = pairs > 0.25*q*q;
  if(q >= 2.0*n && coupled && fits)
    return {requested, "observation", "smaller observation system with dense random coupling", workspace};
  return {requested, "random", !fits ? "observation workspace exceeds budget" :
    (q < 2.0*n ? "random system not substantially larger" :
      "sparse random coupling favors random-space factorization"), workspace};
}

class LaplaceSystem {
public:
  int symbolicAnalyses = 0;
  int factorizations = 0;
  int priorFactorizations = 0;
  int covarianceBuilds = 0;
  int selectedBuilds = 0;
  int contractionSolves = 0;
  std::string traceChoice = "panels";
  double traceMemory = 1000;
  double logDet = 0.0;
  double logDetRatio = 0.0;
  virtual ~LaplaceSystem() = default;
  virtual void prepare(const arma::sp_mat & precision, const arma::sp_mat & design,
    double logDetPrecision) = 0;
  virtual bool factorize(const arma::sp_mat & curvature) = 0;
  virtual bool determinant() = 0;
  virtual arma::mat solve(const arma::mat & rhs) = 0;
  virtual bool observationSpace() const = 0;
  virtual bool selectedEntry(int,int,double &){ return false; }
  virtual int blockCount() const { return 1; }
  void configureTrace(const std::string & choice,double memory){
    laplaceChoice(choice,{"panels","selected"},"trace");
    if(choice=="selected" && observationSpace())
      Rcpp::stop("Selected-inverse Laplace traces require the random-space backend.");
    traceChoice=choice;
    traceMemory=memory;
  }
};

class RandomLaplaceSystem : public LaplaceSystem {
  const arma::sp_mat * precision = nullptr;
  const arma::sp_mat * design = nullptr;
  CholmodState factor;
  double priorLogDet = 0.0;
  LaplaceSelectedInverse selected;
  bool selectedReady = false;
  bool useSelected = false;
public:
  void prepare(const arma::sp_mat & suppliedPrecision, const arma::sp_mat & suppliedDesign,
      double logDetPrecision) override {
    precision = &suppliedPrecision;
    design = &suppliedDesign;
    priorLogDet = logDetPrecision;
  }
  bool factorize(const arma::sp_mat & curvature) override {
    const arma::sp_mat hessian = *precision + design->t()*curvature*(*design);
    const auto eigen = laplaceEigenSparse(hessian);
    ++factorizations;
    selectedReady=false;
    useSelected=traceChoice=="selected";
    return factor.analyzePattern(eigen,&symbolicAnalyses) && factor.factorize(eigen);
  }
  bool determinant() override {
    if(!laplaceFactorLogDeterminant(factor,logDet)) return false;
    logDetRatio = logDet-priorLogDet;
    return true;
  }
  arma::mat solve(const arma::mat & rhs) override {
    Eigen::Map<const Eigen::MatrixXd> mapped(rhs.memptr(),rhs.n_rows,rhs.n_cols);
    const Eigen::MatrixXd result = factor.solve(mapped,"Native Laplace system solve failed.");
    return arma::mat(result.data(),result.rows(),result.cols());
  }
  bool observationSpace() const override { return false; }
  bool selectedEntry(int a,int b,double & value) override {
    if(!useSelected) return false;
    if(!selectedReady){
      selected.build(factor,traceMemory);
      selectedReady=true;
      ++selectedBuilds;
    }
    return selected.get(a,b,value);
  }
};

class BlockLaplaceSystem : public LaplaceSystem {
  const arma::sp_mat * precision=nullptr;
  const arma::sp_mat * design=nullptr;
  double priorLogDet=0;
  std::vector<arma::uvec> groups;
  std::vector<int> groupIndex,localIndex;
  std::vector<arma::sp_mat> blocks;
  std::vector<std::unique_ptr<CholmodState>> factors;
  std::vector<LaplaceSelectedInverse> selected;
  std::vector<bool> selectedReady;
  double selectedBytes=0;
public:
  void prepare(const arma::sp_mat & p,const arma::sp_mat & z,double determinant) override {
    precision=&p; design=&z; priorLogDet=determinant;
  }
  bool factorize(const arma::sp_mat & curvature) override {
    const arma::sp_mat h=*precision+design->t()*curvature*(*design);
    const arma::uword q=h.n_rows;
    std::vector<int> parent(q);
    for(arma::uword k=0;k<q;++k) parent[k]=k;
    auto root=[&](int a){
      while(parent[a]!=a){ parent[a]=parent[parent[a]]; a=parent[a]; }
      return a;
    };
    for(auto cell=h.begin();cell!=h.end();++cell){
      const int a=root(cell.row()),b=root(cell.col());
      parent[b]=a;
    }
    std::vector<int> ids(q,-1),nextGroup(q),nextLocal(q);
    std::vector<std::vector<arma::uword>> members;
    for(arma::uword k=0;k<q;++k){
      const int r=root(k);
      if(ids[r]<0){ ids[r]=members.size(); members.emplace_back(); }
      nextGroup[k]=ids[r]; nextLocal[k]=members[ids[r]].size();
      members[ids[r]].push_back(k);
    }
    if(64.0*h.n_nonzero+128.0*q+4096.0*members.size()>traceMemory*1048576.0)
      Rcpp::stop("Laplace block workspace exceeds execution$denseMemoryMB.");
    if(nextGroup!=groupIndex){
      groupIndex=nextGroup; localIndex=nextLocal;
      groups.clear(); factors.clear(); selected.clear();
      for(const auto & member:members){
        groups.push_back(arma::uvec(member));
        factors.emplace_back(new CholmodState());
        selected.emplace_back();
      }
    }
    blocks.clear();
    for(const auto & group:groups)
      blocks.emplace_back(group.n_elem,group.n_elem);
    for(auto cell=h.begin();cell!=h.end();++cell)
      blocks[groupIndex[cell.row()]](localIndex[cell.row()],localIndex[cell.col()])=*cell;
    selectedReady.assign(groups.size(),false);
    selectedBytes=0;
    factorizations+=blocks.size();
    for(std::size_t k=0;k<blocks.size();++k){
      const auto eigen=laplaceEigenSparse(blocks[k]);
      if(!factors[k]->analyzePattern(eigen,&symbolicAnalyses) ||
         !factors[k]->factorize(eigen)) return false;
    }
    return true;
  }
  bool determinant() override {
    logDet=0;
    for(auto & factor:factors){
      double value;
      if(!laplaceFactorLogDeterminant(*factor,value)) return false;
      logDet+=value;
    }
    logDetRatio=logDet-priorLogDet;
    return std::isfinite(logDet);
  }
  arma::mat solve(const arma::mat & rhs) override {
    arma::mat result(rhs.n_rows,rhs.n_cols);
    for(std::size_t k=0;k<groups.size();++k){
      const arma::mat panel=rhs.rows(groups[k]);
      Eigen::Map<const Eigen::MatrixXd> mapped(panel.memptr(),panel.n_rows,panel.n_cols);
      const Eigen::MatrixXd value=factors[k]->solve(mapped,"Laplace block solve failed.");
      result.rows(groups[k])=arma::mat(value.data(),value.rows(),value.cols());
    }
    return result;
  }
  bool observationSpace() const override { return false; }
  int blockCount() const override { return groups.size(); }
  bool selectedEntry(int a,int b,double & value) override {
    if(groupIndex[a]!=groupIndex[b]){ value=0; return true; }
    const int k=groupIndex[a];
    if(traceChoice=="panels") return false;
    if(!selectedReady[k]){
      selected[k].build(*factors[k],traceMemory-selectedBytes/1048576.0);
      selectedBytes+=LaplaceSelectedInverse::workspace(*factors[k]);
      selectedReady[k]=true; ++selectedBuilds;
    }
    return selected[k].get(localIndex[a],localIndex[b],value);
  }
};

class ObservationLaplaceSystem : public LaplaceSystem {
  arma::sp_mat precision;
  arma::sp_mat design;
  CholmodState prior;
  arma::mat root;
  arma::mat lower;
  arma::sp_mat curvature;
  double memoryMB;
  double workspace = 0.0;
  double priorLogDet = 0.0;
  bool ready = false;

  arma::mat priorSolve(const arma::mat & rhs){
    Eigen::Map<const Eigen::MatrixXd> mapped(rhs.memptr(),rhs.n_rows,rhs.n_cols);
    const Eigen::MatrixXd result = prior.solve(mapped,"Laplace prior covariance solve failed.");
    return arma::mat(result.data(),result.rows(),result.cols());
  }
  void checkPanel(arma::uword columns){
    const double bytes = workspace+8.0*(6.0*design.n_cols+8.0*design.n_rows)*columns;
    if(!std::isfinite(bytes) || bytes > memoryMB*1048576.0)
      Rcpp::stop("Observation-space Laplace solve panel exceeds execution$denseMemoryMB.");
  }
public:
  explicit ObservationLaplaceSystem(double memory) : memoryMB(memory) {}
  void prepare(const arma::sp_mat & suppliedPrecision, const arma::sp_mat & suppliedDesign,
      double logDetPrecision) override {
    if(ready && suppliedDesign.n_rows == design.n_rows &&
       suppliedDesign.n_cols == design.n_cols &&
       (suppliedPrecision-precision).eval().n_nonzero == 0 &&
       (suppliedDesign-design).eval().n_nonzero == 0) return;
    ready = false;
    const auto plan = planLaplaceSystem("observation",suppliedDesign,suppliedPrecision,memoryMB);
    workspace = plan.workspace;
    design = suppliedDesign;
    precision = suppliedPrecision;
    priorLogDet = logDetPrecision;
    const auto eigen = laplaceEigenSparse(precision);
    if(!prior.analyzePattern(eigen) || !prior.factorize(eigen))
      Rcpp::stop("Observation-space Laplace requires positive-definite random precision.");
    ++priorFactorizations;
    const arma::uword n = design.n_rows;
    const double available = memoryMB*1048576.0-workspace;
    const arma::uword panel = static_cast<arma::uword>(std::min(32.0,
      std::floor(available/(8.0*(6.0*design.n_cols+8.0*n)))));
    arma::mat covariance(n,n,arma::fill::zeros);
    for(arma::uword first=0; first<n; first+=panel){
      Rcpp::checkUserInterrupt();
      const arma::uword last = std::min(n-1,first+panel-1);
      covariance.cols(first,last) = design*priorSolve(arma::mat(design.rows(first,last).t()));
    }
    covariance = arma::symmatu(0.5*(covariance+covariance.t()));
    arma::vec values;
    arma::mat vectors;
    if(!arma::eig_sym(values,vectors,covariance) || !values.is_finite())
      Rcpp::stop("Observation-space Laplace latent covariance eigendecomposition failed.");
    const double roundoff = 128.0*std::numeric_limits<double>::epsilon()*
      std::max(1.0,arma::abs(values).max())*n;
    if(values.min() < -roundoff)
      Rcpp::stop("Observation-space Laplace latent covariance is not positive semidefinite.");
    values.transform([](double value){ return std::sqrt(std::max(0.0,value)); });
    root = vectors*arma::diagmat(values);
    ++covarianceBuilds;
    ready = true;
  }
  bool factorize(const arma::sp_mat & suppliedCurvature) override {
    curvature = suppliedCurvature;
    const arma::mat system = arma::eye<arma::mat>(root.n_cols,root.n_cols)+
      root.t()*curvature*root;
    ++factorizations;
    if(!system.is_finite() || !arma::chol(lower,arma::symmatu(system),"lower")) return false;
    logDetRatio = 2.0*arma::accu(arma::log(lower.diag()));
    logDet = priorLogDet+logDetRatio;
    return std::isfinite(logDet);
  }
  bool determinant() override { return std::isfinite(logDet); }
  arma::mat solve(const arma::mat & rhs) override {
    checkPanel(rhs.n_cols);
    const arma::mat base = priorSolve(rhs);
    const arma::mat weighted = curvature*(design*base);
    arma::mat intermediate, corrected;
    if(!arma::solve(intermediate,arma::trimatl(lower),root.t()*weighted,
        arma::solve_opts::no_approx) ||
       !arma::solve(corrected,arma::trimatu(lower.t()),intermediate,
        arma::solve_opts::no_approx))
      Rcpp::stop("Observation-space Laplace triangular solve failed.");
    const arma::mat projected = weighted-curvature*(root*corrected);
    arma::mat result = base-priorSolve(arma::mat(design.t()*projected));
    if(!result.is_finite()) Rcpp::stop("Observation-space Laplace solve is non-finite.");
    return result;
  }
  bool observationSpace() const override { return true; }
};

inline std::unique_ptr<LaplaceSystem> makeLaplaceSystem(const LaplaceSystemPlan & plan,
    double memoryMB,const std::string & blocks="off"){
  laplaceChoice(blocks,{"off","on"},"blocks");
  if(plan.selected=="observation" && blocks=="on")
    Rcpp::stop("Explicit Laplace block decomposition requires the random-space backend.");
  if(plan.selected == "observation")
    return std::unique_ptr<LaplaceSystem>(new ObservationLaplaceSystem(memoryMB));
  if(blocks=="on")
    return std::unique_ptr<LaplaceSystem>(new BlockLaplaceSystem());
  return std::unique_ptr<LaplaceSystem>(new RandomLaplaceSystem());
}

inline Rcpp::List laplaceSystemDiagnostics(const LaplaceSystemPlan & plan,
    const LaplaceSystem & system){
  return Rcpp::List::create(Rcpp::_["requested"]=plan.requested,
    Rcpp::_["selected"]=plan.selected,Rcpp::_["reason"]=plan.reason,
    Rcpp::_["workspaceBytes"]=system.observationSpace() ? plan.workspace : 0.0,
    Rcpp::_["factorizations"]=system.factorizations,
    Rcpp::_["priorFactorizations"]=system.priorFactorizations,
    Rcpp::_["covarianceBuilds"]=system.covarianceBuilds,
    Rcpp::_["blocks"]=system.blockCount(),
    Rcpp::_["trace"]=system.traceChoice,
    Rcpp::_["selectedBuilds"]=system.selectedBuilds,
    Rcpp::_["contractionSolves"]=system.contractionSolves);
}

template<class System>
inline bool laplacePredictorVariance(System & system,
    const arma::sp_mat & transposed,arma::uword column,double & result){
  result=0;
  for(auto a=transposed.begin_col(column);a!=transposed.end_col(column);++a)
    for(auto b=transposed.begin_col(column);b!=transposed.end_col(column);++b){
      double value;
      if(!system.selectedEntry(a.row(),b.row(),value)) return false;
      result+=(*a)*(*b)*value;
    }
  return true;
}

template<class System>
inline arma::mat laplacePredictorCovariance(System & system,
    const arma::sp_mat & transposed,arma::uword panelSize){
  const arma::uword d=transposed.n_cols;
  arma::mat result(d,d,arma::fill::zeros);
  bool available=system.traceChoice!="panels" && !system.observationSpace();
  for(arma::uword j=0;available && j<d;++j)
    for(arma::uword k=0;available && k<d;++k){
      double sum=0;
      for(auto a=transposed.begin_col(j);available && a!=transposed.end_col(j);++a)
        for(auto b=transposed.begin_col(k);b!=transposed.end_col(k);++b){
          double value;
          if(!system.selectedEntry(a.row(),b.row(),value)){ available=false; break; }
          sum+=(*a)*(*b)*value;
        }
      result(j,k)=sum;
    }
  if(available) return result;
  for(arma::uword first=0;first<d;first+=panelSize){
    const arma::uword last=std::min(d-1,first+panelSize-1);
    ++system.contractionSolves;
    result.cols(first,last)=transposed.t()*system.solve(arma::mat(transposed.cols(first,last)));
  }
  return result;
}

template<class System>
inline bool laplaceSelectedTrace(System & system,
    const LaplacePrecisionDerivative & derivative,double & result){
  if(system.traceChoice=="panels" || system.observationSpace()) return false;
  result=0;
  for(auto cell=derivative.block.begin();cell!=derivative.block.end();++cell){
    double value;
    if(!system.selectedEntry(derivative.start+cell.row(),
        derivative.start+cell.col(),value)) return false;
    result+=(*cell)*value;
  }
  return true;
}

} // namespace sommer
#endif
