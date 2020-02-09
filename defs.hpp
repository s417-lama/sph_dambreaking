#pragma once

#include <cmath>

#include "vector2.hpp"
#include "bounding_box2.hpp"

#ifndef SPH_2D
# define SPH_2D 1
#endif

#ifndef SPH_DOUBLE
# define SPH_DOUBLE 1
#endif

#ifndef SPH_CFL_DT
# define SPH_CFL_DT 0
#endif

#ifndef SPH_REUSE_TREE
# define SPH_REUSE_TREE 0
#endif

#ifndef SPH_DATA_SCALE
# define SPH_DATA_SCALE 1
#endif

#ifndef SPH_MAX_STEP
# define SPH_MAX_STEP 1000
#endif

#ifndef SPH_OUTPUT_INTERVAL
# define SPH_OUTPUT_INTERVAL 0
#endif

#ifndef SPH_PARTICLES_CUTOFF
# define SPH_PARTICLES_CUTOFF 64
#endif

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
constexpr real PI       = atan(1.0) * 4.0;
#if SPH_2D
const realvec gravity(0.0, -9.81);
#else
const realvec gravity(0.0, 0.0, -9.81);
#endif

/* Structs */
typedef struct {
  real dens;
  real pres;
} Dens;

typedef struct {
  realvec acc;
#if SPH_CFL_DT
  real f;
#endif
} Hydro;

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

  inline void apply_from(const Dens& dens) {
    this->dens = dens.dens;
    this->pres = dens.pres;
  }

  inline void apply_from(const Hydro& hydro) {
    this->acc = hydro.acc;
#if SPH_CFL_DT
    this->f = hydro.f;
#endif
  }
} Particle;
