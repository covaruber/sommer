#ifndef SOMMER_MME_RECOVERY_H
#define SOMMER_MME_RECOVERY_H

namespace sommer {

struct AIRecoveryPolicy {
  int difficultSteps = 0;
  int negligibleGains = 0;

  const char * observe(bool hasPrevious, int halvings, double gain,
                       double likelihoodTolerance, bool scoreSmall) {
    if(scoreSmall) {
      reset();
      return nullptr;
    }
    if(hasPrevious) {
      difficultSteps = halvings >= 4 ? difficultSteps + 1 : 0;
      negligibleGains = gain < likelihoodTolerance ? negligibleGains + 1 : 0;
    }
    if(difficultSteps >= 3) { return "repeated-backtracking"; }
    if(negligibleGains >= 3) { return "negligible-gains"; }
    return nullptr;
  }

  void reset() {
    difficultSteps = 0;
    negligibleGains = 0;
  }
};

} // namespace sommer

#endif
