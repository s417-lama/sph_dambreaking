#pragma once

#include <cmath>

#include "defs.hpp"

void calc_dens(Particle* const ps_i, const int ni,
               const Particle* const ps_j, const int nj);

void calc_hydro(Particle* const ps_i, const int ni,
                const Particle* const ps_j, const int nj);

inline real calc_pressure(const real dens) {
  return std::max(0.0, C_B * (pow(dens / DENS0, 7) - 1));
}
