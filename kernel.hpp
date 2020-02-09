#pragma once

#include <cmath>

#include "config.hpp"
#include "defs.hpp"

SPH_KERNEL
void calc_dens(Particle* const ps_i, const int ni,
               const Particle* const ps_j, const int nj);

SPH_KERNEL
void calc_hydro(Particle* const ps_i, const int ni,
                const Particle* const ps_j, const int nj);

SPH_KERNEL
inline real calc_pressure(const real dens) {
  return max_(0.0, C_B * (pow(dens / DENS0, 7) - 1));
}

#ifdef SPH_CUDA_PARALLEL
extern __device__ calc_kernel_t calc_dens_kernel;
extern __device__ calc_kernel_t calc_hydro_kernel;

inline calc_kernel_t get_calc_kernel(calc_type type) {
  calc_kernel_t h_func;
  switch (type) {
    case CALC_TYPE_DENS:
      cudaMemcpyFromSymbol(&h_func, calc_dens_kernel, sizeof(calc_kernel_t));
      break;
    case CALC_TYPE_HYDRO:
      cudaMemcpyFromSymbol(&h_func, calc_hydro_kernel, sizeof(calc_kernel_t));
      break;
  }
  return h_func;
}
#else
inline calc_kernel_t get_calc_kernel(calc_type type) {
  calc_kernel_t h_func;
  switch (type) {
    case CALC_TYPE_DENS:
      h_func = calc_dens;
      break;
    case CALC_TYPE_HYDRO:
      h_func = calc_hydro;
      break;
  }
  return h_func;
}
#endif
