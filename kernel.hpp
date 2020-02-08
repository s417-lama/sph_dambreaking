#pragma once

#include <cmath>

#include "defs.hpp"

void calc_dens(const Particle* const ps_i, const int ni,
               const Particle* const ps_j, const int nj, Dens* const dens);

void calc_hydro(const Particle* const ps_i, const int ni,
                const Particle* const ps_j, const int nj, Hydro* const hydro);

inline real calc_pressure(const real dens) {
  return std::max(0.0, C_B * (pow(dens / DENS0, 7) - 1));
}
