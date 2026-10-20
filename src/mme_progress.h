#ifndef SOMMER_MME_PROGRESS_H
#define SOMMER_MME_PROGRESS_H

#include <Rcpp.h>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

namespace sommer {

inline std::string mmePathStatus(const std::string & solver, bool latentRR,
    bool schur, bool chain, bool banded, bool factorSchur,
    const std::string & schurMode, const std::string & residual) {
  const auto state = [](bool active) { return active ? "active" : "inactive"; };
  std::ostringstream text;
  text << "solver=" << solver << "; latent RR=" << state(latentRR)
       << "; Schur=" << state(schur);
  if(!schur) {
    text << (schurMode == "off" ? " (disabled)"
      : solver != "cholmod" ? " (requires CHOLMOD)"
      : schurMode == "force" ? " (forced, no eligible plan)"
      : " (auto: not selected)");
  }
  text << "; chain=" << state(chain) << "; banded=" << state(banded)
       << "; factor-Schur=" << state(factorSchur)
       << "; residual=" << residual
       << "; grid=" << state(residual == "grid")
       << "; repeated-Kronecker=" << state(residual == "kronecker");
  return text.str();
}

inline std::string mmePathStatus(const Rcpp::List & result, bool latentRR) {
  const Rcpp::List diagnostics = result["engineDiagnostics"];
  const Rcpp::List control = diagnostics["control"];
  return mmePathStatus(Rcpp::as<std::string>(result["solver"]), latentRR,
    Rcpp::as<bool>(diagnostics["blockSchurActive"]),
    Rcpp::as<bool>(diagnostics["blockChainActive"]),
    Rcpp::as<bool>(diagnostics["blockBandedActive"]),
    Rcpp::as<bool>(diagnostics["factorSchurActive"]),
    Rcpp::as<std::string>(control["blockSchur"]),
    Rcpp::as<std::string>(diagnostics["residualEngine"]));
}

class ProfileTimer {
  double * elapsedSeconds;
  std::chrono::steady_clock::time_point started;
public:
  ProfileTimer(double & elapsed, bool enabled)
    : elapsedSeconds(enabled ? &elapsed : nullptr),
      started(enabled ? std::chrono::steady_clock::now()
                      : std::chrono::steady_clock::time_point()) {}
  void stop(){
    if(elapsedSeconds != nullptr){
      *elapsedSeconds += std::chrono::duration<double>(
        std::chrono::steady_clock::now() - started).count();
      elapsedSeconds = nullptr;
    }
  }
  ~ProfileTimer(){ stop(); }
};

class MMEProgress {
  std::chrono::steady_clock::time_point started = std::chrono::steady_clock::now();
  std::vector<int> widths;

  void print(const std::vector<std::string> & values) const {
    std::ostringstream row;
    for(std::size_t column = 0; column < values.size(); ++column) {
      if(column > 0) { row << ' '; }
      row << std::setw(widths[column]) << values[column];
    }
    Rcpp::Rcout << row.str() << '\n';
  }

public:
  explicit MMEProgress(int structures)
    : widths{5, 14, 10, 12, 15, 6, 6, 5, 8} {
    int notPDWidth = 0;
    for(int term = 1; term <= structures; ++term) {
      notPDWidth += static_cast<int>(std::to_string(term).size()) + (term > 1);
    }
    widths[7] = std::max(widths[7], notPDWidth);
  }

  static std::string number(double value, int width, int precision) {
    if(!std::isfinite(value)) { return "-"; }
    std::ostringstream text;
    text << std::fixed << std::setprecision(precision) << value;
    if(text.str().size() > static_cast<std::size_t>(width)) {
      text.str("");
      text << std::scientific << std::setprecision(std::max(0, width - 7)) << value;
    }
    return text.str();
  }

  static std::string update(double weight, bool expanded = false) {
    const std::string method = expanded ? "PXEM" : "EM";
    return weight == 0.0 ? "AI" : weight == 1.0 ? method
      : "AI+" + method + ":" + number(weight, 4, 2);
  }

  static std::string policy(bool automatic, bool expanded, double weight,
                            const std::string & reason) {
    if(!automatic) { return weight == 0.0 ? "off" : "scheduled"; }
    if(reason == "repeated-backtracking") { return "backtracking"; }
    if(reason == "negligible-gains") { return "low-gain"; }
    if(reason == "information-solve") { return "info-solve"; }
    if(!reason.empty()) { return reason; }
    return expanded ? "auto:pxem" : "auto:em";
  }

  void header() const {
    print({"iter", "LogLik", "deltaLL", "update", "policy", "trust", "line", "notPD", "time(s)"});
  }

  void row(int iteration, double likelihood, double previous,
           const std::string & update, const std::string & policy,
           double trust = std::numeric_limits<double>::quiet_NaN(),
           double line = std::numeric_limits<double>::quiet_NaN(),
           const std::string & notPD = "-") const {
    const double elapsed = std::chrono::duration<double>(
      std::chrono::steady_clock::now() - started).count();
    print({std::to_string(iteration), number(likelihood, widths[1], 6),
      number(iteration > 1 ? likelihood - previous
        : std::numeric_limits<double>::quiet_NaN(), widths[2], 4),
      update, policy, number(trust, widths[5], 3), number(line, widths[6], 3),
      notPD, number(elapsed, widths[8], 2)});
  }
};

} // namespace sommer

#endif
