// FastVelocityWeakeningLaw.h

#ifndef SEISSOL_DR_FASTVELOCITYWEAKENINGLAW_H
#define SEISSOL_DR_FASTVELOCITYWEAKENINGLAW_H

#include "FrictionLaw.h"

namespace seissol::dr {

class FastVelocityWeakeningLaw : public FrictionLaw {
public:
  FastVelocityWeakeningLaw(DynamicRupture::DRParameters* drParameters, seissol::initializers::Layer& layer);

  void updateFrictionAndSlipRate(double time, DRGodunovData& godunovData, DRGodunovData& neighborData) override;

private:
  // Parameters for the nucleation region
  double tstart;        // Start time of nucleation
  double tend;        // End time of nucleation
  double R0;        // Initial radius of nucleation region
  double R1;        // Final radius of nucleation region
  double x0, y0;    // Center coordinates of nucleation region

  // Function to compute nucleation friction coefficient
  double computeMuNuc(double x, double y, double time);
};

} // namespace seissol::dr

#endif // SEISSOL_DR_FASTVELOCITYWEAKENINGLAW_H