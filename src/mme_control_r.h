#ifndef SOMMER_MME_CONTROL_R_H
#define SOMMER_MME_CONTROL_R_H

#include <Rcpp.h>
#include "mme_control.h"

namespace sommer {
namespace control {

// Descriptor attributes remain private transport. Read at the owning engine's
// original stage, and on every recursive fit: wrappers rewrite these attributes.
class Transport {
  const Rcpp::List & descriptors;
public:
  explicit Transport(const Rcpp::List & value) : descriptors(value) {}

  bool logical(const char * name, bool fallback = false) const {
    return descriptors.hasAttribute(name)
      ? Rcpp::as<bool>(descriptors.attr(name)) : fallback;
  }

  std::string string(const char * name, const char * fallback) const {
    return descriptors.hasAttribute(name)
      ? Rcpp::as<std::string>(descriptors.attr(name)) : fallback;
  }

  Method method() const { return control::method(string("acceleration", "none")); }

  Blend blend(const char * error) const {
    Blend result = Blend::EM;
    if(!parseBlend(string("blend", "em"), result)) { Rcpp::stop(error); }
    return result;
  }

  Information information() const {
    Information result = Information::Full;
    if(!parseInformation(string("covarianceInformation", "full"), result)) {
      Rcpp::stop("covarianceInformation must be 'full' or 'none'.");
    }
    return result;
  }

  Memory memory(const Rcpp::Function & getOption) const {
    SEXP value = descriptors.hasAttribute("memoryPolicy")
      ? static_cast<SEXP>(descriptors.attr("memoryPolicy"))
      : static_cast<SEXP>(getOption("sommer.mme.memoryPolicy"));
    Memory result = Memory::Speed;
    if(TYPEOF(value) == STRSXP && Rf_xlength(value) == 1) {
      if(!parseMemory(Rcpp::as<std::string>(value), result)) {
        Rcpp::stop("memoryPolicy must be 'speed' or 'limit'.");
      }
    } else if(value != R_NilValue) {
      Rcpp::stop("memoryPolicy must be 'speed' or 'limit'.");
    }
    return result;
  }

  Mode experimentalMode(const std::string & attribute, const std::string & option,
                        Mode fallback, const Rcpp::Function & getOption) const {
    if(descriptors.hasAttribute(attribute)) {
      SEXP value = static_cast<SEXP>(descriptors.attr(attribute));
      if(TYPEOF(value) == STRSXP && Rf_xlength(value) == 1) {
        Mode result = fallback;
        if(parseMode(Rcpp::as<std::string>(value), result)) { return result; }
        Rcpp::stop(attribute + " must be 'auto', 'off', 'force', TRUE, or FALSE.");
      }
      if(TYPEOF(value) == LGLSXP && Rf_xlength(value) == 1) {
        const int flag = LOGICAL(value)[0];
        if(flag == NA_LOGICAL) { return fallback; }
        return flag ? Mode::Force : Mode::Off;
      }
    }
    const Rcpp::RObject value = getOption(option);
    if(!value.isNULL()) {
      if(TYPEOF(value) != LGLSXP || Rf_xlength(value) != 1) {
        Rcpp::stop(option + " must be TRUE or FALSE.");
      }
      const Rcpp::LogicalVector flags(value);
      if(flags[0] == NA_LOGICAL) { Rcpp::stop(option + " must be TRUE or FALSE."); }
      return flags[0] == TRUE ? Mode::Force : Mode::Off;
    }
    return fallback;
  }
};

} // namespace control
} // namespace sommer

#endif
