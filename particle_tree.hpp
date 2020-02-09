#pragma once

#include <vector>

#include "defs.hpp"

typedef struct ParticleTreeNode {
  Particle*                             particles_i;
  Particle*                             particles_j;
  int                                   n_particles;
  bool                                  is_leaf;
  BoundingBox                           bbox;
  BoundingBox                           outer_bbox;
  struct ParticleTreeNode*              children[(1 << DIM)];
  std::vector<struct ParticleTreeNode*> neighbors;

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

class ParticleTree {
  private:
    int               n_particles_;
    Particle*         particles_i_;
    Particle*         particles_j_;
    ParticleTreeNode* root_;
    bool              j_outdated;

    void barrier() {
      parallel_for(0, n_particles_, [&] (int i) {
        particles_j_[i] = particles_i_[i];
      });
    }

  public:
    ParticleTree(const std::vector<Particle>& particles);
    ~ParticleTree();

    void build();

    template <typename Func>
    void calc(const Func body) {
      if (j_outdated) {
        barrier();
      }
      pfor_leaf_impl(root_, [&] (ParticleTreeNode* leaf) {
        int nj = 0;
        for (const auto& nb : leaf->neighbors) {
          nj += nb->n_particles;
        }
        Particle* ps_j = new Particle[nj];
        int c = 0;
        for (const auto& nb : leaf->neighbors) {
          for (int j = 0; j < nb->n_particles; j++) {
            ps_j[c++] = nb->particles_j[j];
          }
        }
        body(leaf->particles_i, leaf->n_particles,
             ps_j, nj);
        delete[] ps_j;
      });
      j_outdated = true;
    }

    template <typename Func>
    inline void for_particle(const Func body) {
      for (int i = 0; i < n_particles_; i++) {
        body(particles_i_[i]);
      }
      j_outdated = true;
    }

    template <typename Func>
    inline void pfor_particle(const Func body) {
      parallel_for(0, n_particles_, [&] (int i) {
        body(particles_i_[i]);
      });
      j_outdated = true;
    }

    template <typename Func>
    inline void pfor_leaf(const Func body) {
      pfor_leaf_impl(root_, body);
    }
};
