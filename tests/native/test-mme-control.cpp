#include "../../src/mme_control.h"
#include <cassert>

int main() {
  using namespace sommer::control;
  const Representation representation;
  assert(!representation.commonScale && !representation.conditionalReducedRank);
  const Outputs outputs;
  assert(outputs.covarianceInformation == Information::Full);
  assert(outputs.retainCoefficientMatrix);
  const Proposal proposal;
  assert(proposal.blend == Blend::EM && !proposal.adaptiveRecovery);
  assert(!proposal.directionBlend && !proposal.accelerationRequested);
  assert(proposal.extrapolation == Extrapolation::None);
  assert(CurvaturePolicy().method == Curvature::ExactAI);
  assert(!CurvaturePolicy().informationOnly);
  assert(Execution().profile == Mode::Off);
  assert(Execution().memory == Memory::Speed);
  assert(Traces().batch == Mode::Off);
  const Elimination elimination;
  assert(elimination.blockBanded == Mode::Auto);
  assert(elimination.factorSchurMultiFactor == Mode::Auto);
  assert(elimination.dFirstML == Mode::Auto);
  assert(elimination.matrixFreePcg == Mode::Auto);
  assert(elimination.cholmodFinalSelected == Mode::Auto);
  assert(elimination.pdCertificates == Mode::Auto);
  assert(elimination.blockSchur == Mode::Auto);
  assert(elimination.saddlePoint == Mode::Off);

  const Method exact = method("none");
  assert(exact.curvature.method == Curvature::ExactAI);
  assert(!exact.proposal.accelerationRequested);
  const Method hybrid = method("hybrid-ai");
  assert(hybrid.curvature.method == Curvature::HybridAI);
  assert(hybrid.proposal.extrapolation == Extrapolation::None);
  assert(hybrid.proposal.accelerationRequested);
  const Method aitken = method("aitken");
  assert(aitken.curvature.method == Curvature::ExactAI);
  assert(aitken.proposal.extrapolation == Extrapolation::Aitken);
  assert(aitken.proposal.accelerationRequested);
  for(const std::string & value : {"other", "AITKEN", "", "NONE"}) {
    const Method unknown = method(value);
    assert(unknown.curvature.method == Curvature::ExactAI);
    assert(unknown.proposal.extrapolation == Extrapolation::None);
    assert(unknown.proposal.accelerationRequested);
  }

  for(const Mode expected : {Mode::Auto, Mode::Off, Mode::Force}) {
    Mode parsed = Mode::Auto;
    assert(parseMode(modeName(expected), parsed) && parsed == expected);
  }
  Mode parsed = Mode::Auto;
  assert(parseMode("TRUE", parsed) && parsed == Mode::Force);
  assert(parseMode("False", parsed) && parsed == Mode::Off);
  assert(parseMode("AuTo", parsed) && parsed == Mode::Auto);
  for(const std::string & value : {"on", "1", "", " off", "force "}) {
    assert(!parseMode(value, parsed));
    assert(parsed == Mode::Auto);
  }
  Blend blend = Blend::EM;
  assert(parseBlend("pxem", blend) && blend == Blend::PXEM);
  assert(std::string(blendName(blend)) == "pxem");
  assert(parseBlend("em", blend) && blend == Blend::EM);
  assert(std::string(blendName(blend)) == "em");
  assert(!parseBlend("PXEM", blend) && blend == Blend::EM);
  Information information = Information::Full;
  assert(parseInformation("none", information) && information == Information::None);
  assert(parseInformation("full", information) && information == Information::Full);
  assert(!parseInformation("FULL", information) && information == Information::Full);
  Memory memory = Memory::Speed;
  assert(parseMemory("limit", memory) && memory == Memory::Limit);
  assert(std::string(memoryName(memory)) == "limit");
  assert(parseMemory("speed", memory) && memory == Memory::Speed);
  assert(std::string(memoryName(memory)) == "speed");
  assert(!parseMemory("Speed", memory) && memory == Memory::Speed);
}
