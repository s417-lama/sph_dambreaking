#include "kernel.hpp"

#ifdef SPH_CUDA_PARALLEL
__device__ calc_kernel_t calc_dens_kernel  = calc_dens;
__device__ calc_kernel_t calc_hydro_kernel = calc_hydro;
#endif

SPH_KERNEL
inline real W(const realvec dr, const real dr2) {
  constexpr real H = SLEN / 2.0;
#if SPH_2D
  /* constexpr real COEF = 10.0 / 7.0 / M_PI / pow(H, 2); */
  constexpr real COEF = 10.0 / 7.0 / M_PI / (H * H);
#else
  /* constexpr real COEF = 1.0 / M_PI / pow(H, 3); */
  constexpr real COEF = 1.0 / M_PI / (H * H * H);
#endif
  const real s = sqrt(dr2) / H;
  real v;
  if (s < 1.0) {
    v = 1.0 - 1.5 * pow(s, 2) + 0.75 * pow(s, 3);
  } else if (s < 2.0) {
    v = 0.25 * pow(2.0 - s, 3);
  }
  return COEF * v;
}

SPH_KERNEL
inline realvec gradW(const realvec dr, const real dr2) {
  constexpr real H = SLEN / 2.0;
#if SPH_2D
  /* constexpr real COEF = 45.0 / 14.0 / M_PI / pow(H, 4); */
  constexpr real COEF = 45.0 / 14.0 / M_PI / (H * H * H * H);
#else
  /* constexpr real COEF = 2.25 / M_PI / pow(H, 5); */
  constexpr real COEF = 2.25 / M_PI / (H * H * H * H * H);
#endif
  const real s = sqrt(dr2) / H;
  realvec v;
  if (s < 1.0) {
    v = (s - 4.0 / 3.0) * dr;
  } else if (s < 2.0) {
    v = - pow(2.0 - s, 2) / 3.0 / s * dr;
  }
  return COEF * v;
}

// calculation of density
SPH_KERNEL
void calc_dens(Particle* const ps_i, const int ni,
               const Particle* const ps_j, const int nj) {
  constexpr real slen2 = SLEN * SLEN;
  for (int i = 0; i < ni; i++) {
    ps_i[i].dens = 0;
    ps_i[i].pres = 0;
    for (int j = 0; j < nj; j++) {
      const realvec dr  = ps_i[i].pos - ps_j[j].pos;
      const real    dr2 = dr * dr;
      if (dr2 >= slen2) continue;
      const real W_ij = W(dr, dr2);
      ps_i[i].dens += ps_j[j].mass * W_ij;
    }
    ps_i[i].pres = calc_pressure(ps_i[i].dens);
  }
}

// calculation of hydro force
SPH_KERNEL
void calc_hydro(Particle* const ps_i, const int ni,
                const Particle* const ps_j, const int nj) {
  constexpr real slen2 = SLEN * SLEN;
#if SPH_2D
  const realvec gravity(0.0, -9.81);
#else
  const realvec gravity(0.0, 0.0, -9.81);
#endif
  for (int i = 0; i < ni; i++) {
    ps_i[i].acc = 0;
#if SPH_CFL_DT
    ps_i[i].f = 0;
#endif
    /* if (ps_i[i].type != FLUID) continue; */
    const real tmp_pd_i = ps_i[i].pres / pow(ps_i[i].dens, 2);
    for (int j = 0; j < nj; j++) {
      const realvec dr  = ps_i[i].pos - ps_j[j].pos;
      const real    dr2 = dr * dr;
      if (dr2 >= slen2) continue;
      const real tmp_pd_j = ps_j[j].pres / pow(ps_j[j].dens, 2);
      const realvec gradW_ij = gradW(dr, dr2);
      const realvec dv = ps_i[i].vel - ps_j[j].vel;
      const real vr = dv * dr;
      const real AV = (vr <= 0) ? 0 : - VISC * vr / (dr2 + 0.01 * slen2);
      ps_i[i].acc -= ps_j[j].mass * (tmp_pd_i + tmp_pd_j + AV) * gradW_ij;
    }
    ps_i[i].acc += gravity;
#if SPH_CFL_DT
    ps_i[i].f = ps_i[i].mass * sqrt(ps_i[i].acc * ps_i[i].acc);
#endif
  }
}
