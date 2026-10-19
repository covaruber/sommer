#ifndef SOMMER_MME_CONTROL_H
#define SOMMER_MME_CONTROL_H

#include <algorithm>
#include <cctype>
#include <string>

namespace sommer {
namespace control {

enum class Mode { Auto, Off, Force };
enum class Curvature { ExactAI, HybridAI };
enum class Extrapolation { None, Aitken };
enum class Blend { EM, PXEM };
enum class Information { Full, None };
enum class Memory { Speed, Limit };

struct Representation {
  bool commonScale = false;
  bool conditionalReducedRank = false;
};

struct CurvaturePolicy {
  Curvature method = Curvature::ExactAI;
  bool informationOnly = false;
};

struct Proposal {
  Extrapolation extrapolation = Extrapolation::None;
  bool accelerationRequested = false;
  bool directionBlend = false;
  Blend blend = Blend::EM;
  bool adaptiveRecovery = false;
};

struct Outputs {
  Information covarianceInformation = Information::Full;
  bool retainCoefficientMatrix = true;
};

struct Execution {
  Mode profile = Mode::Off;
  Memory memory = Memory::Speed;
};

struct Elimination {
  Mode blockBanded = Mode::Auto;
  Mode factorSchurMultiFactor = Mode::Auto;
  Mode dFirstML = Mode::Auto;
  Mode matrixFreePcg = Mode::Auto;
  Mode cholmodFinalSelected = Mode::Auto;
  Mode pdCertificates = Mode::Auto;
  Mode blockSchur = Mode::Auto;
  Mode saddlePoint = Mode::Off;
};

struct Traces {
  Mode batch = Mode::Off;
};

struct Method {
  CurvaturePolicy curvature;
  Proposal proposal;
};

// Unknown legacy acceleration strings were ignored by the update routines,
// but still excluded from deterministic-only representations.
inline Method method(const std::string & acceleration) {
  Method result;
  result.curvature.method = acceleration == "hybrid-ai"
    ? Curvature::HybridAI : Curvature::ExactAI;
  result.proposal.extrapolation = acceleration == "aitken"
    ? Extrapolation::Aitken : Extrapolation::None;
  result.proposal.accelerationRequested = acceleration != "none";
  return result;
}

inline const char * modeName(Mode mode) {
  return mode == Mode::Force ? "force" : mode == Mode::Off ? "off" : "auto";
}

inline bool parseMode(std::string value, Mode & mode) {
  std::transform(value.begin(), value.end(), value.begin(),
    [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
  if(value == "auto") { mode = Mode::Auto; }
  else if(value == "off" || value == "false") { mode = Mode::Off; }
  else if(value == "force" || value == "true") { mode = Mode::Force; }
  else { return false; }
  return true;
}

inline bool parseBlend(const std::string & value, Blend & blend) {
  if(value == "em") { blend = Blend::EM; }
  else if(value == "pxem") { blend = Blend::PXEM; }
  else { return false; }
  return true;
}

inline const char * blendName(Blend blend) {
  return blend == Blend::PXEM ? "pxem" : "em";
}

inline bool parseInformation(const std::string & value, Information & information) {
  if(value == "full") { information = Information::Full; }
  else if(value == "none") { information = Information::None; }
  else { return false; }
  return true;
}

inline bool parseMemory(const std::string & value, Memory & memory) {
  if(value == "speed") { memory = Memory::Speed; }
  else if(value == "limit") { memory = Memory::Limit; }
  else { return false; }
  return true;
}

inline const char * memoryName(Memory memory) {
  return memory == Memory::Speed ? "speed" : "limit";
}

} // namespace control
} // namespace sommer

#endif
