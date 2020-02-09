#pragma once

#include <vector>

#include "config.hpp"
#include "defs.hpp"

typedef struct ParticleTreeNode {
  Particle*                             particles_i;
  Particle*                             particles_j; // TODO: remove?
  int                                   n_particles;
  bool                                  is_leaf;
  BoundingBox                           bbox;
  BoundingBox                           outer_bbox;
  struct ParticleTreeNode*              children[(1 << DIM)];
  std::vector<struct ParticleTreeNode*> neighbors;
  int                                   n_neighbors;
#if SPH_CUDA_PARALLEL
  int                                   index;
#endif

  ParticleTreeNode(Particle* ps_i, Particle* ps_j, const int n)
    : particles_i(ps_i), particles_j(ps_j), n_particles(n), is_leaf(false) {}
} ParticleTreeNode;

template <typename Func>
inline void pfor_leaf_impl(ParticleTreeNode* node, const Func body) {
  if (node->is_leaf) {
    body(node);
  } else {
    for (int i = 0; i < (1 << DIM); i++) {
      if (ParticleTreeNode* child = node->children[i]) {
        pfor_leaf_impl(child, body);
      }
    }
  }
}

template <typename Func>
inline void for_leaf_impl(ParticleTreeNode* node, const Func body) {
  if (node->is_leaf) {
    body(node);
  } else {
    for (int i = 0; i < (1 << DIM); i++) {
      if (ParticleTreeNode* child = node->children[i]) {
        for_leaf_impl(child, body);
      }
    }
  }
}

#if SPH_CUDA_PARALLEL
#define cudaCheckError(ans) { cudaAssert((ans), __FILE__, __LINE__); }
inline void cudaAssert(cudaError_t code, const char *file, int line) {
   if (code != cudaSuccess) {
      fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
      if (abort) exit(code);
   }
}

template <typename Func>
__global__
static void on_gpu(Particle* ps_i, Particle* ps_j, int np,
                   int* pi_offsets, int* pj_offsets, Func body) {
  int l = blockIdx.x;
  int ibegin = pi_offsets[l];
  int iend   = pi_offsets[l + 1];

  int ni = (iend - ibegin + blockDim.x - 1) / blockDim.x;
  int ibegin_ = ibegin + threadIdx.x * ni;
  if (ibegin_ + ni >= iend) {
    ni = iend - ibegin_;
  }

  int jbegin = pj_offsets[l];
  int jend   = pj_offsets[l + 1];
  int nj = jend - jbegin;

  Particle* ps_i_ = &ps_i[ibegin_];
  Particle* ps_j_ = &ps_j[jbegin];

  body(ps_i_, ni, ps_j_, nj);
}
#endif

class ParticleTree {
  private:
    int               n_particles_;
    Particle*         particles_i_;
    Particle*         particles_j_;
    ParticleTreeNode* root_;
#if SPH_CUDA_PARALLEL
    int               n_leafs_;
    int*              pi_offsets_;
    int*              pj_offsets_;
    int               pj_buf_size_;
    Particle*         pj_buf_;
#endif

  public:
    ParticleTree(const std::vector<Particle>& particles);
    ~ParticleTree();

    void build();

    template <typename Func>
    void calc(const Func body) {
#if SPH_CUDA_PARALLEL
      pfor_leaf_impl(root_, [&] (ParticleTreeNode* leaf) {
        int idx = leaf->index;
        Particle* ps_j = &pj_buf_[pj_offsets_[idx]];
        int c = 0;
        for (const auto& nb : leaf->neighbors) {
          for (int j = 0; j < nb->n_particles; j++) {
            ps_j[c++] = nb->particles_i[j];
          }
        }
      });

      Particle* d_ps_i;
      Particle* d_ps_j;
      int*      d_pi_offsets;
      int*      d_pj_offsets;

      cudaCheckError(cudaMalloc(&d_ps_i, sizeof(Particle) * n_particles_));
      cudaCheckError(cudaMalloc(&d_ps_j, sizeof(Particle) * pj_buf_size_));
      cudaCheckError(cudaMalloc(&d_pi_offsets, sizeof(int) * (n_leafs_ + 1)));
      cudaCheckError(cudaMalloc(&d_pj_offsets, sizeof(int) * (n_leafs_ + 1)));

      cudaCheckError(cudaMemcpy(d_ps_i, particles_i_, sizeof(Particle) * n_particles_, cudaMemcpyHostToDevice));
      cudaCheckError(cudaMemcpy(d_ps_j, pj_buf_     , sizeof(Particle) * pj_buf_size_, cudaMemcpyHostToDevice));
      cudaCheckError(cudaMemcpy(d_pi_offsets, pi_offsets_, sizeof(int) * (n_leafs_ + 1), cudaMemcpyHostToDevice));
      cudaCheckError(cudaMemcpy(d_pj_offsets, pj_offsets_, sizeof(int) * (n_leafs_ + 1), cudaMemcpyHostToDevice));

      on_gpu<<<n_leafs_, 32>>>(d_ps_i, d_ps_j, n_particles_, d_pi_offsets, d_pj_offsets, body);
      cudaCheckError(cudaPeekAtLastError());
      cudaCheckError(cudaDeviceSynchronize());

      cudaCheckError(cudaMemcpy(particles_i_, d_ps_i, sizeof(Particle) * n_particles_, cudaMemcpyDeviceToHost));

      cudaCheckError(cudaFree(d_ps_i));
      cudaCheckError(cudaFree(d_ps_j));
      cudaCheckError(cudaFree(d_pi_offsets));
      cudaCheckError(cudaFree(d_pj_offsets));
#else
      pfor_leaf_impl(root_, [&] (ParticleTreeNode* leaf) {
        int nj = leaf->n_neighbors;
        Particle* ps_j = new Particle[nj];
        int c = 0;
        for (const auto& nb : leaf->neighbors) {
          for (int j = 0; j < nb->n_particles; j++) {
            ps_j[c++] = nb->particles_i[j];
          }
        }
        Particle* ps_i = leaf->particles_i;
        int ni = leaf->n_particles;
        body(ps_i, ni, ps_j, nj);
        delete[] ps_j;
      });
#endif
    }

    template <typename Func>
    inline void for_particle(const Func body) {
      for (int i = 0; i < n_particles_; i++) {
        body(particles_i_[i]);
      }
    }

    template <typename Func>
    inline void pfor_particle(const Func body) {
      parallel_for(0, n_particles_, [&] (int i) {
        body(particles_i_[i]);
      });
    }

    template <typename Func>
    inline void pfor_leaf(const Func body) {
      pfor_leaf_impl(root_, body);
    }

    template <typename Func>
    inline void for_leaf(const Func body) {
      for_leaf_impl(root_, body);
    }

#if SPH_CUDA_PARALLEL
    void setup_global_array();
    void destroy_global_array();
#endif
};
