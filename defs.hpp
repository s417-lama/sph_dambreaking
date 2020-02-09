#pragma once

#include "config.hpp"
#include "vector2.hpp"
#include "bounding_box2.hpp"

#if SPH_DOUBLE
typedef double real;
#else
typedef float  real;
#endif

#if SPH_2D
constexpr int DIM = 2;
typedef Vector2<real>      realvec;
typedef BoundingBox2<real> BoundingBox;
#else
constexpr int DIM = 3;
typedef Vector3<real>      realvec;
typedef BoundingBox3<real> BoundingBox;
#endif

typedef enum {
  FLUID = 1,
  WALL  = 2,
} particle_type;

typedef struct {
  real          mass;
  realvec       pos;
#if SPH_REUSE_TREE
  realvec       prev_pos;
#endif
  realvec       vel;
  realvec       acc;
  real          dens;
  real          pres;
  realvec       vel_half;
#if SPH_CFL_DT
  real          f;
#endif
  particle_type type;
} Particle;

typedef void (* calc_kernel_t)(Particle* const, const int, const Particle* const, const int);

typedef enum {
  CALC_TYPE_DENS,
  CALC_TYPE_HYDRO,
} calc_type;

/* Parameters */
constexpr real END_TIME = 1.5;
constexpr real L0       = 0.55 / 30 / SPH_DATA_SCALE;
constexpr real SLEN     = L0 * 2.1;
#if SPH_REUSE_TREE
constexpr real SKIN     = SLEN * 0.3;
#else
constexpr real SKIN     = 0.0;
#endif
constexpr real DENS0    = 1000.0;
constexpr real SND      = 31.3;
constexpr real C_B      = DENS0 * SND * SND / 7;
constexpr real ALPHA    = 0.1;
constexpr real VISC     = ALPHA * SLEN * SND / DENS0;
constexpr real DT       = 0.4 * SLEN / SND / (1 + 0.6 * ALPHA);
/* constexpr real DT       = 0.0005; */
