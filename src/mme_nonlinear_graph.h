#ifndef SOMMER_MME_NONLINEAR_GRAPH_H
#define SOMMER_MME_NONLINEAR_GRAPH_H

#include <RcppArmadillo.h>
#include <numeric>
#include <vector>
#include <string>

namespace sommer {

struct GraphJet {
  double value=0;
  double first[3]={};
  double second[3][3]={};
  double third[3][3][3]={};
};

struct GraphNode {
  std::string op;
  arma::uvec inputs;
  arma::uvec inputWidths;
  arma::umat weights;
  arma::uvec bias;
  arma::uword width=0, column=0;
  double center=0, scale=1;
};

inline double graphNaturalParameter(double value, bool logScale){
  if(!logScale) return value;
  const double result=std::exp(value);
  if(!std::isfinite(result) || result<=0)
    Rcpp::stop("Log-transformed graph parameter exceeds numeric range.");
  return result;
}

inline double graphTransformDerivative(double natural, bool logScale){
  return logScale ? natural : 1;
}

inline GraphJet graphParameterJet(double value, bool logScale,
    const arma::mat & directions, arma::uword index){
  GraphJet result;
  const double natural=graphNaturalParameter(value,logScale);
  result.value=natural;
  for(arma::uword a=0;a<3;++a){
    const double da=directions(index,a);
    result.first[a]=logScale ? natural*da : da;
    for(arma::uword b=0;b<3;++b){
      const double db=directions(index,b);
      result.second[a][b]=logScale ? natural*da*db : 0;
      for(arma::uword c=0;c<3;++c){
        const double dc=directions(index,c);
        result.third[a][b][c]=logScale ? natural*da*db*dc : 0;
      }
    }
  }
  return result;
}

inline GraphJet graphProduct(const GraphJet & a,const GraphJet & b){
  GraphJet result;
  result.value=a.value*b.value;
  for(int i=0;i<3;++i){
    result.first[i]=a.first[i]*b.value+a.value*b.first[i];
    for(int j=0;j<3;++j){
      result.second[i][j]=a.second[i][j]*b.value+
        a.first[i]*b.first[j]+a.first[j]*b.first[i]+a.value*b.second[i][j];
      for(int k=0;k<3;++k){
        result.third[i][j][k]=a.third[i][j][k]*b.value+
          a.second[i][j]*b.first[k]+a.second[i][k]*b.first[j]+
          a.second[j][k]*b.first[i]+a.first[i]*b.second[j][k]+
          a.first[j]*b.second[i][k]+a.first[k]*b.second[i][j]+
          a.value*b.third[i][j][k];
      }
    }
  }
  return result;
}

inline GraphJet graphSum(const GraphJet & a,const GraphJet & b){
  GraphJet result;
  result.value=a.value+b.value;
  for(int i=0;i<3;++i){
    result.first[i]=a.first[i]+b.first[i];
    for(int j=0;j<3;++j){
      result.second[i][j]=a.second[i][j]+b.second[i][j];
      for(int k=0;k<3;++k)
        result.third[i][j][k]=a.third[i][j][k]+b.third[i][j][k];
    }
  }
  return result;
}

inline GraphJet graphActivation(const GraphJet & x,bool tanhActivation){
  if(!tanhActivation) return x;
  GraphJet result;
  result.value=std::tanh(x.value);
  const double first=1-result.value*result.value;
  const double second=-2*result.value*first;
  const double third=-2*first*(1-3*result.value*result.value);
  for(int i=0;i<3;++i){
    result.first[i]=first*x.first[i];
    for(int j=0;j<3;++j){
      result.second[i][j]=second*x.first[i]*x.first[j]+first*x.second[i][j];
      for(int k=0;k<3;++k)
        result.third[i][j][k]=third*x.first[i]*x.first[j]*x.first[k]+
          second*(x.second[i][j]*x.first[k]+x.second[i][k]*x.first[j]+
            x.second[j][k]*x.first[i])+first*x.third[i][j][k];
    }
  }
  return result;
}

class NonlinearGraph {
  std::vector<GraphNode> nodes;
  arma::uword output, parameterCount;
  std::vector<bool> logScale;

  std::vector<arma::rowvec> forward(const arma::rowvec & data,
      const arma::rowvec & parameters,bool stopOnNonFinite=true) const {
    std::vector<arma::rowvec> values(nodes.size());
    for(arma::uword i=0;i<nodes.size();++i){
      const auto & node=nodes[i];
      if(node.op=="input"){
        values[i]=arma::rowvec(1);
        values[i](0)=(data(node.column)-node.center)/node.scale;
      }else if(node.op=="affine"){
        arma::rowvec input(arma::accu(node.inputWidths));
        arma::uword offset=0;
        for(arma::uword j=0;j<node.inputs.n_elem;++j){
          const arma::uword source=node.inputs(j), width=node.inputWidths(j);
          input.cols(offset,offset+width-1)=values[source];
          offset+=width;
        }
        values[i].zeros(node.width);
        for(arma::uword r=0;r<node.weights.n_rows;++r){
          for(arma::uword c=0;c<node.weights.n_cols;++c){
            const arma::uword parameter=node.weights(r,c)-1;
            values[i](c)+=input(r)*graphNaturalParameter(parameters(parameter),
              logScale[parameter]);
          }
        }
        for(arma::uword c=0;c<node.bias.n_elem;++c) if(node.bias(c)){
          const arma::uword parameter=node.bias(c)-1;
          values[i](c)+=graphNaturalParameter(parameters(parameter),
            logScale[parameter]);
        }
      }else{
        values[i]=values[node.inputs(0)];
        if(node.op=="tanh") values[i]=arma::tanh(values[i]);
      }
      if(!values[i].is_finite()){
        if(stopOnNonFinite)
          Rcpp::stop("Nonlinear graph produced a non-finite intermediate value.");
        return {};
      }
    }
    return values;
  }

  std::vector<std::vector<GraphJet>> forwardJets(const arma::rowvec & data,
      const arma::rowvec & parameters,const arma::mat & directions) const {
    std::vector<std::vector<GraphJet>> values(nodes.size());
    for(arma::uword i=0;i<nodes.size();++i){
      const auto & node=nodes[i];
      values[i].resize(node.width);
      if(node.op=="input"){
        values[i][0].value=(data(node.column)-node.center)/node.scale;
      }else if(node.op=="affine"){
        std::vector<GraphJet> input;
        for(arma::uword j=0;j<node.inputs.n_elem;++j){
          const auto & source=values[node.inputs(j)];
          input.insert(input.end(),source.begin(),source.end());
        }
        for(arma::uword r=0;r<node.weights.n_rows;++r){
          for(arma::uword c=0;c<node.weights.n_cols;++c){
            const arma::uword parameter=node.weights(r,c)-1;
            const GraphJet product=graphProduct(
              graphParameterJet(parameters(parameter),logScale[parameter],
                directions,parameter),input[r]);
            values[i][c]=graphSum(values[i][c],product);
          }
        }
        for(arma::uword c=0;c<node.bias.n_elem;++c) if(node.bias(c)){
          const arma::uword parameter=node.bias(c)-1;
          const GraphJet term=graphParameterJet(parameters(parameter),
            logScale[parameter],directions,parameter);
          values[i][c]=graphSum(values[i][c],term);
        }
      }else{
        for(arma::uword c=0;c<node.width;++c)
          values[i][c]=graphActivation(values[node.inputs(0)][c],node.op=="tanh");
      }
      for(const auto & value:values[i])
        if(!std::isfinite(value.value))
          Rcpp::stop("Nonlinear graph produced a non-finite intermediate value.");
    }
    return values;
  }

  arma::rowvec reverse(const std::vector<arma::rowvec> & values,
      const arma::rowvec & parameters) const {
    std::vector<arma::rowvec> adjoints(nodes.size());
    for(arma::uword i=0;i<nodes.size();++i) adjoints[i].zeros(nodes[i].width);
    adjoints[output](0)=1;
    arma::rowvec gradient(parameterCount,arma::fill::zeros);
    for(arma::sword rev=static_cast<arma::sword>(nodes.size())-1;rev>=0;--rev){
      const arma::uword i=static_cast<arma::uword>(rev);
      const auto & node=nodes[i];
      if(node.op=="identity" || node.op=="tanh"){
        const arma::uword source=node.inputs(0);
        for(arma::uword c=0;c<node.width;++c){
          const double derivative=node.op=="tanh" ?
            1-values[i](c)*values[i](c) : 1;
          adjoints[source](c)+=adjoints[i](c)*derivative;
        }
      }else if(node.op=="affine"){
        std::vector<arma::rowvec> inputs(node.inputs.n_elem);
        for(arma::uword source=0;source<node.inputs.n_elem;++source)
          inputs[source]=values[node.inputs(source)];
        for(arma::uword r=0;r<node.weights.n_rows;++r){
          arma::uword source=0, sourceColumn=r;
          while(sourceColumn>=node.inputWidths(source))
            sourceColumn-=node.inputWidths(source++);
          const arma::uword sourceNode=node.inputs(source);
          for(arma::uword c=0;c<node.weights.n_cols;++c){
            const arma::uword parameter=node.weights(r,c)-1;
            const double weight=graphNaturalParameter(parameters(parameter),
              logScale[parameter]);
            gradient(parameter)+=adjoints[i](c)*inputs[source](sourceColumn)*
              graphTransformDerivative(weight,logScale[parameter]);
            adjoints[sourceNode](sourceColumn)+=adjoints[i](c)*weight;
          }
        }
        for(arma::uword c=0;c<node.bias.n_elem;++c) if(node.bias(c)){
          const arma::uword parameter=node.bias(c)-1;
          const double bias=graphNaturalParameter(parameters(parameter),
            logScale[parameter]);
          gradient(parameter)+=adjoints[i](c)*
            graphTransformDerivative(bias,logScale[parameter]);
        }
      }
    }
    return gradient;
  }

public:
  arma::uword parameterCountValue() const { return parameterCount; }
  double workspaceBytes(arma::uword rows,arma::uword inputColumns) const {
    const double widths=std::accumulate(nodes.begin(),nodes.end(),0.0,
      [](double total,const GraphNode & node){ return total+node.width; });
    return 8.0*rows*(inputColumns+7.0*parameterCount+16.0*widths);
  }
  bool validParameters(const arma::vec & parameters) const {
    if(parameters.n_elem!=parameterCount || !parameters.is_finite()) return false;
    for(arma::uword i=0;i<parameterCount;++i){
      if(logScale[i]){
        const double natural=std::exp(parameters(i));
        if(!std::isfinite(natural) || natural<=0) return false;
      }
    }
    return true;
  }

  NonlinearGraph(const Rcpp::List & specification,arma::uword selectedOutput,
      arma::uword inputCount,arma::uword parameters,
      const Rcpp::IntegerVector & transformed)
    : output(selectedOutput),parameterCount(parameters),logScale(parameters,false){
    if(!specification.size() || output>=static_cast<arma::uword>(specification.size()) ||
       transformed.size()!=static_cast<R_xlen_t>(parameterCount))
      Rcpp::stop("Invalid native nonlinear graph layout.");
    for(arma::uword i=0;i<parameterCount;++i) logScale[i]=transformed(i)!=0;
    nodes.reserve(specification.size());
    for(R_xlen_t i=0;i<specification.size();++i){
      Rcpp::List source=specification[i];
      GraphNode node;
      node.op=Rcpp::as<std::string>(source["op"]);
      node.width=Rcpp::as<arma::uword>(source["width"]);
      if(node.op=="input"){
        node.column=Rcpp::as<arma::uword>(source["columnIndex"]);
        node.center=Rcpp::as<double>(source["center"]);
        node.scale=Rcpp::as<double>(source["scale"]);
        if(node.column>=inputCount || node.width!=1 || !std::isfinite(node.center) ||
           !std::isfinite(node.scale) || node.scale<=0)
          Rcpp::stop("Invalid native graph input transform.");
      }else{
        node.inputs=Rcpp::as<arma::uvec>(source["inputIndices"]);
        node.inputWidths=Rcpp::as<arma::uvec>(source["inputWidths"]);
        if(!node.inputs.n_elem || node.inputs.n_elem!=node.inputWidths.n_elem)
          Rcpp::stop("Invalid native graph node inputs.");
        for(arma::uword j=0;j<node.inputs.n_elem;++j)
          if(node.inputs(j)>=i || !node.inputWidths(j))
            Rcpp::stop("Native graph nodes must be in topological order.");
        if(node.op=="affine"){
          Rcpp::IntegerMatrix weights=source["weights"];
          Rcpp::IntegerVector bias=source["bias"];
          const arma::uword inputWidth=arma::accu(node.inputWidths);
          if(weights.nrow()!=static_cast<int>(inputWidth) ||
             weights.ncol()!=static_cast<int>(node.width) ||
             bias.size()!=static_cast<R_xlen_t>(node.width))
            Rcpp::stop("Invalid native affine dimensions.");
          node.weights.set_size(weights.nrow(),weights.ncol());
          node.bias.set_size(bias.size());
          for(int r=0;r<weights.nrow();++r) for(int c=0;c<weights.ncol();++c){
            const int index=weights(r,c);
            if(index<1 || index>static_cast<int>(parameterCount))
              Rcpp::stop("Invalid native affine parameter index.");
            node.weights(r,c)=static_cast<arma::uword>(index);
          }
          for(R_xlen_t c=0;c<bias.size();++c){
            const int index=bias(c);
            if(index<0 || index>static_cast<int>(parameterCount))
              Rcpp::stop("Invalid native affine bias index.");
            node.bias(c)=static_cast<arma::uword>(index);
          }
        }else if((node.op!="identity" && node.op!="tanh") ||
                 node.inputs.n_elem!=1 || node.width!=node.inputWidths(0)){
          Rcpp::stop("Unsupported or dimensionally invalid native graph operation.");
        }
      }
      nodes.push_back(node);
    }
    if(nodes[output].width!=1)
      Rcpp::stop("Native nonlinear graph output must be scalar.");
  }

  Rcpp::List evaluate(const arma::mat & data,const arma::mat & parameters,
      const arma::mat & directions,double memoryMB) const {
    const arma::uword n=data.n_rows;
    if(!n || parameters.n_rows!=n || parameters.n_cols!=parameterCount ||
       data.n_cols==0 || !data.is_finite() || !parameters.is_finite() ||
       directions.n_rows!=parameterCount || directions.n_cols!=3*n ||
       !directions.is_finite() || !std::isfinite(memoryMB) || memoryMB<=0)
      Rcpp::stop("Native nonlinear graph data or parameter dimensions are invalid.");
    const double bytes=workspaceBytes(n,data.n_cols);
    if(bytes>memoryMB*1048576)
      Rcpp::stop("Nonlinear graph workspace exceeds execution$denseMemoryMB.");
    Rcpp::NumericVector outputValues(n),gradient(static_cast<R_xlen_t>(n)*parameterCount);
    Rcpp::NumericVector hessianVector(static_cast<R_xlen_t>(n)*parameterCount);
    Rcpp::NumericVector thirdGradient(static_cast<R_xlen_t>(n)*parameterCount);
    Rcpp::NumericVector first(static_cast<R_xlen_t>(n)*3);
    Rcpp::NumericVector second(static_cast<R_xlen_t>(n)*9);
    Rcpp::NumericVector third(static_cast<R_xlen_t>(n)*27);
    for(arma::uword row=0;row<n;++row){
      if(row%256==0) Rcpp::checkUserInterrupt();
      const arma::rowvec observation=data.row(row), parameter=parameters.row(row);
      const auto values=forward(observation,parameter);
      const arma::rowvec grad=reverse(values,parameter);
      outputValues(row)=values[output](0);
      for(arma::uword j=0;j<parameterCount;++j)
        gradient[row+n*j]=grad(j);
      const arma::mat rowDirections=directions.cols(3*row,3*row+2);
      const auto jets=forwardJets(observation,parameter,rowDirections);
      const GraphJet & result=jets[output][0];
      for(int a=0;a<3;++a){
        first[row+n*a]=result.first[a];
        for(int b=0;b<3;++b){
          second[row+n*(a+3*b)]=result.second[a][b];
          for(int c=0;c<3;++c)
            third[row+n*(a+3*(b+3*c))]=result.third[a][b][c];
        }
      }
      if(arma::any(arma::vectorise(rowDirections)!=0)){
        for(arma::uword j=0;j<parameterCount;++j){
          arma::mat contraction(parameterCount,3,arma::fill::zeros);
          contraction.col(1)=rowDirections.col(0);
          contraction.col(2)=rowDirections.col(1);
          contraction(j,0)=1;
          const auto contracted=forwardJets(observation,parameter,contraction);
          const GraphJet & secondResult=contracted[output][0];
          hessianVector[row+n*j]=secondResult.second[0][1];
          thirdGradient[row+n*j]=secondResult.third[0][1][2];
        }
      }
    }
    first.attr("dim")=Rcpp::IntegerVector::create(n,3);
    second.attr("dim")=Rcpp::IntegerVector::create(n,3,3);
    third.attr("dim")=Rcpp::IntegerVector::create(n,3,3,3);
    gradient.attr("dim")=Rcpp::IntegerVector::create(n,parameterCount);
    hessianVector.attr("dim")=Rcpp::IntegerVector::create(n,parameterCount);
    thirdGradient.attr("dim")=Rcpp::IntegerVector::create(n,parameterCount);
    return Rcpp::List::create(Rcpp::_["output"]=outputValues,
      Rcpp::_["gradient"]=gradient,Rcpp::_["hessianVector"]=hessianVector,
      Rcpp::_["thirdDirectionalGradient"]=thirdGradient,
      Rcpp::_["directionalFirst"]=first,
      Rcpp::_["directionalSecond"]=second,Rcpp::_["directionalThird"]=third);
  }

  Rcpp::List evaluateFirst(const arma::mat & data,const arma::vec & parameters,
      double memoryMB) const {
    const arma::uword n=data.n_rows;
    if(!n || data.n_cols==0 || parameters.n_elem!=parameterCount ||
       !data.is_finite() || !parameters.is_finite() ||
       !std::isfinite(memoryMB) || memoryMB<=0)
      Rcpp::stop("Native nonlinear graph data or parameter dimensions are invalid.");
    const double bytes=workspaceBytes(n,data.n_cols);
    if(bytes>memoryMB*1048576)
      Rcpp::stop("Nonlinear graph workspace exceeds execution$denseMemoryMB.");
    arma::vec outputValues(n);
    arma::mat gradient(n,parameterCount);
    for(arma::uword row=0;row<n;++row){
      if(row%256==0) Rcpp::checkUserInterrupt();
      const arma::rowvec observation=data.row(row);
      const auto values=forward(observation,parameters.t(),false);
      if(values.empty())
        return Rcpp::List::create(Rcpp::_["valid"]=false,
          Rcpp::_["reason"]="non-finite graph intermediate");
      outputValues(row)=values[output](0);
      gradient.row(row)=reverse(values,parameters.t());
      if(!gradient.row(row).is_finite())
        return Rcpp::List::create(Rcpp::_["valid"]=false,
          Rcpp::_["reason"]="non-finite graph derivative");
    }
    if(!outputValues.is_finite() || !gradient.is_finite())
      Rcpp::stop("Nonlinear graph produced a non-finite value or derivative.");
    return Rcpp::List::create(Rcpp::_["valid"]=true,Rcpp::_["reason"]="evaluated",
      Rcpp::_["output"]=outputValues,Rcpp::_["jacobian"]=gradient);
  }
};

} // namespace sommer
#endif
