#pragma once

#include <ctime>

inline uint64_t gettime_in_nsec() {
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

template <typename T>
inline const T& max_(const T& a, const T& b) {
  return (a > b) ? a : b;
}

template <typename T>
inline const T& min_(const T& a, const T& b) {
  return (b < a) ? b : a;
}
