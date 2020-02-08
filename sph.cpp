#include <cstdio>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <sys/stat.h>
#include <sys/time.h>

#include "defs.hpp"
#include "kernel.hpp"

static inline uint64_t gettime_in_nsec() {
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  return (uint64_t)ts.tv_sec * 1000000000 + (uint64_t)ts.tv_nsec;
}

template <typename Func>
inline void parallel_for(int begin, int end, const Func body) {
#ifdef SPH_TASK_PARALLEL
  mtbb::parallel_for(begin, end, 1, CUTOFF_PFOR, [&] (int a, int b) {
    for (int i = a; i < b; i++) {
      body(i);
    }
  });
#else
#ifdef SPH_LOOP_PARALLEL
#pragma omp parallel for
#endif
  for (int i = begin; i < end; i++) {
    body(i);
  }
#endif
}

template <typename Func>
inline void pfor_particles(std::vector<Particle>& particles, const Func body) {
  parallel_for(0, particles.size(), [&] (int i) {
    body(particles[i]);
  });
}

template <typename T>
void apply_result(Particle* particles, const int n, const T* results) {
  parallel_for(0, n, [&] (int i) {
    particles[i].apply_from(results[i]);
  });
}

void setup_particles(std::vector<Particle>& particles, const char* filename) {
  std::ifstream ifs(filename);
  realvec pos;
  int type;

  while (ifs >> pos >> type) {
    Particle p;
    p.pos  = pos;
    p.type = (particle_type)type;
#if SPH_2D
    p.mass = DENS0 * pow(l0, 2);
#else
    p.mass = DENS0 * pow(l0, 3);
#endif
    p.vel  = 0;
    p.acc  = 0;
    p.dens = DENS0;
    p.pres = calc_pressure(DENS0);
    particles.push_back(p);
  }
}

void output_particles(std::vector<Particle>& particles, const char* filename) {
  std::ofstream ofs(filename);

  for (int i = 0; i < particles.size(); i++) {
    ofs << particles[i].pos << " " << particles[i].type << std::endl;
  }
}

inline double get_time_step(const std::vector<Particle>& particles) {
#if SPH_CFL_DT
   real fmax = 0.0;
   for (int i = 0; i < particles.size(); i++) {
     fmax = std::max(fmax, particles[i].f);
   }
   if (fmax == 0.0) {
     return DT;
   } else {
     return std::min(0.25 * slen / fmax, DT);
   }
#else
   return DT;
#endif
}

void initial_kick(std::vector<Particle>& particles, const double dt) {
  pfor_particles(particles, [&] (Particle& p) {
    if (p.type == FLUID) {
      p.vel_half = p.vel + 0.5 * dt * p.acc;
    }
  });
}

#if SPH_REUSE_TREE
bool full_drift(std::vector<Particle>& particles, const double dt) {
  bool reuse = true;
  // time becomes t + dt;
  pfor_particles(particles, [&] (Particle& p) {
    if (p.type == FLUID) {
      p.pos += dt * p.vel_half;
      // check whether we should reuse the list
      realvec dp = p.pos - p.prev_pos;
      if (sqrt(dp * dp) >= skin / 2.0) reuse = false;
    }
  });
  return reuse;
}
#else
void full_drift(std::vector<Particle>& particles, const double dt) {
  // time becomes t + dt;
  pfor_particles(particles, [&] (Particle& p) {
    if (p.type == FLUID) {
      p.pos += dt * p.vel_half;
    }
  });
}
#endif

void final_kick(std::vector<Particle>& particles, const double dt) {
  pfor_particles(particles, [&] (Particle& p) {
    if (p.type == FLUID) {
      p.vel = p.vel_half + 0.5 * dt * p.acc;
    }
  });
}

#if SPH_REUSE_TREE
void set_prev_pos(std::vector<Particle>& particles) {
  pfor_particles(particles, [&] (Particle& p) {
    p.prev_pos = p.pos;
  });
}
#endif

int main(int argc, char* argv[]) {
  std::vector<Particle> particles;

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
  setup_particles(particles, datafile);

  Dens*  dens  = new Dens[particles.size()];
  Hydro* hydro = new Hydro[particles.size()];

#if SPH_REUSE_TREE
  bool reuse = false;
  int reuse_count = 0;
#endif

  // Main loop for time integration
  double dt;
  int step = 0;
  uint64_t t_all = 0;
  for (double time = 0; time < END_TIME && step < SPH_MAX_STEP; time += dt, step++) {
    uint64_t t1 = gettime_in_nsec();

    if (step > 0) {
      // Leap frog: Initial Kick & Full Drift
      initial_kick(particles, dt);
#if SPH_REUSE_TREE
      reuse = full_drift(particles, dt);
#else
      full_drift(particles, dt);
#endif
    }

#if SPH_REUSE_TREE
    if (reuse) {
      dens_tree.calcForceAllAndWriteBack(CalcDensity()   , particles, dinfo, false, PS::REUSE_LIST);
      hydr_tree.calcForceAllAndWriteBack(CalcHydroForce(), particles, dinfo, false, PS::REUSE_LIST);
    } else {
      set_prev_pos(particles);
      dens_tree.calcForceAllAndWriteBack(CalcDensity()   , particles, dinfo, false, PS::MAKE_LIST_FOR_REUSE);
      hydr_tree.calcForceAllAndWriteBack(CalcHydroForce(), particles, dinfo, false, PS::MAKE_LIST_FOR_REUSE);
      reuse_count++;
    }
#else
    calc_dens(particles.data(), particles.size(), particles.data(), particles.size(), dens);
    apply_result(particles.data(), particles.size(), dens);
    calc_hydro(particles.data(), particles.size(), particles.data(), particles.size(), hydro);
    apply_result(particles.data(), particles.size(), hydro);
#endif
    if (step > 0) {
      // Leap frog: Final Kick
      final_kick(particles, dt);
    }
    // Get a new timestep
    dt = get_time_step(particles);

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
      output_particles(particles, filename);

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
    std::cout << "time: "    << std::fixed   << std::setprecision(5) << time    << " [s] "
              << "step: "    << std::setw(5) << std::right           << step    << " "
              << "elapsed: " << std::setw(7) << std::right           << t2 - t1 << " [usec] "
#if SPH_REUSE_TREE
              << "reuse: "   << (reuse ? "true" : "false")
#endif
              << std::endl;
  }

#if SPH_REUSE_TREE
  std::cout << "reuse = " << reuse_count << std::endl;
  std::cout << "total time = " << (double)t_all / 1000 / 1000 << " sec" << std::endl;
#endif

  return 0;
}
