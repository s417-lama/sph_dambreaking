#include <cstdio>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>

#include "config.hpp"
#include "util.hpp"
#include "defs.hpp"
#include "particle_tree.hpp"
#include "kernel.hpp"

void setup_particles(std::vector<Particle>& particles, const char* filename) {
  std::ifstream ifs(filename);
  realvec pos;
  int type;

  while (ifs >> pos >> type) {
    Particle p;
    p.pos  = pos;
    p.type = (particle_type)type;
#if SPH_2D
    p.mass = DENS0 * pow(L0, 2);
#else
    p.mass = DENS0 * pow(L0, 3);
#endif
    p.vel  = 0;
    p.acc  = 0;
    p.dens = DENS0;
    p.pres = calc_pressure(DENS0);
    particles.push_back(p);
  }
}

void output_particles(ParticleTree& ptree, const char* filename) {
  std::ofstream ofs(filename);

  ptree.for_particle([&] (Particle& p) {
    ofs << p.pos << " " << p.type << std::endl;
  });
}

inline double get_time_step(ParticleTree& ptree) {
#if SPH_CFL_DT
  real fmax = 0.0;
  ptree.for_particle([&] (Particle& p) {
    fmax = std::max(fmax, p.f);
  });
  if (fmax == 0.0) {
    return DT;
  } else {
    return std::min(0.25 * SLEN / fmax, DT);
  }
#else
  return DT;
#endif
}

void initial_kick(ParticleTree& ptree, const double dt) {
  ptree.pfor_particle([&] (Particle& p) {
    if (p.type == FLUID) {
      p.vel_half = p.vel + 0.5 * dt * p.acc;
    }
  });
}

#if SPH_REUSE_TREE
 bool full_drift(ParticleTree& ptree, const double dt) {
  bool reuse = true;
  // time becomes t + dt;
  ptree.pfor_particle([&] (Particle& p) {
    if (p.type == FLUID) {
      p.pos += dt * p.vel_half;
      // check whether we should reuse the list
      realvec dp = p.pos - p.prev_pos;
      if (sqrt(dp * dp) >= SKIN * 0.5) reuse = false;
    }
  });
  return reuse;
}
#else
void full_drift(ParticleTree& ptree, const double dt) {
  // time becomes t + dt;
  ptree.pfor_particle([&] (Particle& p) {
    if (p.type == FLUID) {
      p.pos += dt * p.vel_half;
    }
  });
}
#endif

void final_kick(ParticleTree& ptree, const double dt) {
  ptree.pfor_particle([&] (Particle& p) {
    if (p.type == FLUID) {
      p.vel = p.vel_half + 0.5 * dt * p.acc;
    }
  });
}

#if SPH_REUSE_TREE
void set_prev_pos(ParticleTree& ptree) {
  ptree.pfor_particle([&] (Particle& p) {
    p.prev_pos = p.pos;
  });
}
#endif

int main(int argc, char* argv[]) {
#if PARTICLE_SIMULATOR_TASK_PARALLEL
#if DISABLE_STEAL
  myth_adws_set_stealable(0);
#endif
#endif

#if SPH_2D
  const char* datafile = "data/data2d.txt";
#else
  const char* datafile = "data/data3d.txt";
#endif

  std::vector<Particle> particles;
  setup_particles(particles, datafile);
  ParticleTree ptree(particles);

#if SPH_REUSE_TREE
  bool reuse = false;
  int reuse_count = 0;
#endif

  calc_kernel_t dens_kernel  = get_calc_kernel(CALC_TYPE_DENS);
  calc_kernel_t hydro_kernel = get_calc_kernel(CALC_TYPE_HYDRO);

  // Main loop for time integration
  double dt = 0;
  int step = 0;
  uint64_t t_all = 0;
  for (double time = 0; time < END_TIME && step < SPH_MAX_STEP; time += dt, step++) {
    uint64_t t1 = gettime_in_nsec();

    if (step > 0) {
      // Leap frog: Initial Kick & Full Drift
      initial_kick(ptree, dt);
#if SPH_REUSE_TREE
      reuse = full_drift(ptree, dt);
#else
      full_drift(ptree, dt);
#endif
    }

#if SPH_REUSE_TREE
    if (!reuse) {
      set_prev_pos(ptree);
      ptree.build();
      reuse_count++;
    }
#else
    ptree.build();
#endif

    // particle interactions
    ptree.calc(dens_kernel);
    ptree.calc(hydro_kernel);

    if (step > 0) {
      // Leap frog: Final Kick
      final_kick(ptree, dt);
    }

    // Get a new timestep
    dt = get_time_step(ptree);

    uint64_t t2 = gettime_in_nsec();
    t_all += t2 - t1;

    // Output result files
#if SPH_OUTPUT_INTERVAL
    if (step % SPH_OUTPUT_INTERVAL == 0) {
      char filename[256];
#if SPH_2D
      sprintf(filename, "result/dambreaking2d.txt.%d", step / SPH_OUTPUT_INTERVAL);
#else
      sprintf(filename, "result/dambreaking3d.txt.%d", step / SPH_OUTPUT_INTERVAL);
#endif
      output_particles(ptree, filename);

      std::cout << "================================" << std::endl;
      std::cout << "output " << filename << "." << std::endl;
      std::cout << "================================" << std::endl;
#ifdef RECORD_SPLIT
      {
        sprintf(filename, "result/split.txt.%d", step / SPH_OUTPUT_INTERVAL);
        std::ofstream fout(filename);
        hydr_tree.dumpSplit(fout);
        fout.close();
      }
#endif
#ifdef RECORD_TIMELINE
      {
        sprintf(filename, "result/timeline.txt.%d", step / SPH_OUTPUT_INTERVAL);
        std::ofstream fout(filename);
        hydr_tree.dumpTimeline(fout);
        fout.close();
      }
#endif
    }
#endif
    // Output information to STDOUT
    std::cout << "time: "    << std::fixed   << std::setprecision(5) << time                           << " [s] "
              << "step: "    << std::setw(5) << std::right           << step                           << " "
              << "elapsed: " << std::setw(7) << std::right           << (double)(t2 - t1) / 1000000000 << " [s] "
#if SPH_REUSE_TREE
              << "reuse: "   << (reuse ? "true" : "false")
#endif
              << std::endl;
  }

#if SPH_REUSE_TREE
  std::cout << "reuse = " << reuse_count << std::endl;
  std::cout << "total time = " << (double)t_all / 1000000000 << " sec" << std::endl;
#endif

  return 0;
}
