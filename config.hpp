#pragma once

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

#ifndef SPH_RECORD_CPU
# define SPH_RECORD_CPU 0
#endif

#ifdef SPH_CUDA_PARALLEL
# define SPH_KERNEL __host__ __device__
#else
# define SPH_KERNEL
#endif
