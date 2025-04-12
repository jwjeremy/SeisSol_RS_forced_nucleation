// FastVelocityForcedNuc.cpp

#include "FastVelocityForcedNuc.h"
#include <cmath>

namespace seissol::dr {

  FastVelocityForcedNuc::FastVelocityForcedNuc(DynamicRupture::DRParameters* drParameters, seissol::initializers::Layer& layer)
    : FrictionLaw(drParameters, layer) {
  // Initialize nucleation parameters
  t0 = 0.0;         // Start time
  t1 = 1.0;         // End time
  R0 = 0.0;         // Initial radius
  R1 = 1000.0;      // Final radius
  x0 = 0.0;         // Center x-coordinate
  y0 = 0.0;         // Center y-coordinate
}

double FastVelocityForcedNuc::computeMuNuc(double x, double y, double time, double rs_f0, double rs_muW) {
  if (time < t0 || time > t1) {
    return rs_f0;
  }

  // Linear interpolation of radius over time
  double R = R0 + (R1 - R0) * ((time - t0) / (t1 - t0));

  // Compute distance from center
  double dx = x - x0;
  double dy = y - y0;
  double distance = std::sqrt(dx * dx + dy * dy);

  if (distance <= R) {
    // Linear reduction of friction coefficient within the nucleation region
    return rs_f0 - (rs_f0 - rs_muW) * ((R - distance) / R);
  } else {
    return rs_f0;
  }
}

void FastVelocityForcedNuc::updateFrictionAndSlipRate(double time, DRGodunovData& godunovData, DRGodunovData& neighborData) {
  // Example loop over fault points (pseudo-code)
  for (size_t i = 0; i < godunovData.numPoints; ++i) {
    double x = godunovData.x[i];
    double y = godunovData.y[i];
    double V = godunovData.slipRate[i];
    double theta = godunovData.stateVariable[i];

    // Rate-and-state friction calculation: rs_muW
    double mu_rs = mu0 + a * std::log(V / V0) + b * std::log(theta * V0 / Dc);

    // Compute target friction coefficient based on nucleation region and rate-and-state friction
    double mu_nuc = computeMuNuc(x, y, time, mu_rs);

    double mu;

    // Only compare mu_nuc and mu_rs if time is within the nucleation window (t0 <= time <= t1)
    if (time >= t0 && time <= t1) {
        mu = std::min(mu_nuc, mu_rs);  // Inside nucleation window, use the minimum
    } else {
        mu = mu_rs;  // Outside the nucleation window, use mu_rs directly
    }

    // Update friction coefficient and other relevant variables
    godunovData.frictionCoefficient[i] = mu;

    // Continue with slip rate and state variable updates...
  }
}

} // namespace seissol::dr
