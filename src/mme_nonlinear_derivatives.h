#ifndef SOMMER_MME_NONLINEAR_DERIVATIVES_H
#define SOMMER_MME_NONLINEAR_DERIVATIVES_H

#include "mme_nonlinear_predictor.h"

namespace sommer {

struct CurveDerivatives {
  double mean = 0;
  double d[3] = {};
  double d2[3][3] = {};
  double d3[3][3][3] = {};
};

// Closed-form coordinate derivatives of A*g(s), not a layer/AD evaluator.
inline bool nonlinearCurveDerivatives(const std::string & curve, double x,
    const arma::vec & theta, CurveDerivatives & out){
  const arma::uword k = nonlinearParameterCount(curve);
  if(theta.n_elem != k || !theta.is_finite() || !std::isfinite(x)) return false;
  if(curve != "logistic" && x < 0) return false;
  const double a = std::exp(theta(0)), rate = std::exp(theta(1));
  if(!std::isfinite(a) || a <= 0 || !std::isfinite(rate) || rate <= 0) return false;
  double s, si[3] = {}, sij[3][3] = {}, sijk[3][3][3] = {};
  if(curve == "exponential"){
    s = -rate*x;
    si[1] = sij[1][1] = sijk[1][1][1] = s;
  }else if(curve == "logistic"){
    s = rate*(x-theta(2));
    si[1] = sij[1][1] = sijk[1][1][1] = s;
    si[2] = sij[1][2] = sij[2][1] = -rate;
    sijk[1][1][2] = sijk[1][2][1] = sijk[2][1][1] = -rate;
  }else{
    if(x == 0){ out = CurveDerivatives(); return true; }
    s = std::log(x)-theta(1);
    si[1] = -1;
  }
  if(!std::isfinite(s)) return false;
  double g, g1, g2, g3;
  if(curve == "exponential"){
    g = g1 = g2 = g3 = std::exp(s);
  }else{
    const double tail = std::exp(-std::abs(s));
    g = s >= 0 ? 1/(1+tail) : tail/(1+tail);
    g1 = tail/((1+tail)*(1+tail));
    g2 = g1*(1-2*g);
    g3 = g1*(1-6*g+6*g*g);
  }
  double gi[3] = {}, gij[3][3] = {};
  out.mean = a*g;
  for(arma::uword i=0; i<k; ++i){
    gi[i] = g1*si[i];
    out.d[i] = a*(gi[i]+(i==0 ? g : 0));
    for(arma::uword j=0; j<k; ++j)
      gij[i][j] = g2*si[i]*si[j]+g1*sij[i][j];
  }
  for(arma::uword i=0; i<k; ++i){
    for(arma::uword j=0; j<k; ++j){
      out.d2[i][j] = a*(gij[i][j]+(i==0 ? gi[j] : 0)+
        (j==0 ? gi[i] : 0)+(i==0 && j==0 ? g : 0));
      for(arma::uword l=0; l<k; ++l){
        const double third = g3*si[i]*si[j]*si[l] +
          g2*(sij[i][j]*si[l]+sij[i][l]*si[j]+sij[j][l]*si[i]) +
          g1*sijk[i][j][l];
        out.d3[i][j][l] = a*(third+(i==0 ? gij[j][l] : 0)+
          (j==0 ? gij[i][l] : 0)+(l==0 ? gij[i][j] : 0)+
          (i==0 && j==0 ? gi[l] : 0)+(i==0 && l==0 ? gi[j] : 0)+
          (j==0 && l==0 ? gi[i] : 0)+(i==0 && j==0 && l==0 ? g : 0));
        if(!std::isfinite(out.d3[i][j][l])) return false;
      }
      if(!std::isfinite(out.d2[i][j])) return false;
    }
    if(!std::isfinite(out.d[i])) return false;
  }
  return std::isfinite(out.mean);
}

} // namespace sommer
#endif
