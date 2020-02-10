#include "particle_tree.hpp"

inline BoundingBox get_bbox(const Particle *particles, const int n) {
  BoundingBox bbox;
  for (int i = 0; i < n; i++) {
    bbox.merge(particles[i].pos);
  }
  return bbox;
}

typedef struct {
  int count[(1 << DIM)];

  int& operator [] (const int i) {
    return count[i];
  }
  int operator [] (const int i) const {
    return count[i];
  }
} ParticleCounter;

ParticleCounter count_particles(const Particle* particles, const int n,
                                const BoundingBox bbox) {
  ParticleCounter counter;
  for (int i = 0; i < (1 << DIM); i++) {
    counter[i] = 0;
  }
  for (int i = 0; i < n; i++) {
    int orthant = particles[i].pos.orthant(bbox.center());
    counter[orthant]++;
  }
  return counter;
}

inline ParticleCounter prefix_sum(const ParticleCounter& a) {
  ParticleCounter b;
  int acc = 0;
  for (int i = 0; i < (1 << DIM); i++) {
    b[i] = acc;
    acc += a[i];
  }
  return b;
}

void partition(const Particle* ps_from, Particle* ps_to, const int n,
               ParticleCounter offsets, const BoundingBox bbox) {
  for (int i = 0; i < n; i++) {
    int orthant = ps_from[i].pos.orthant(bbox.center());
    ps_to[offsets[orthant]++] = ps_from[i];
  }
}

ParticleTreeNode* build_tree(Particle *particles1, Particle* particles2,
                             const int n, const BoundingBox bbox, bool flip) {
  // create a node
  ParticleTreeNode* node;
  if (flip) {
    node = new ParticleTreeNode(particles2, particles1, n, bbox);
  } else {
    node = new ParticleTreeNode(particles1, particles2, n, bbox);
  }

  if (n <= SPH_PARTICLES_CUTOFF) {
    node->is_leaf = true;
    for (int i = 0; i < n; i++) particles2[i] = particles1[i];
  } else {
    // count particles
    ParticleCounter counter = count_particles(particles1, n, bbox);
    // prefix sum
    ParticleCounter offsets = prefix_sum(counter);
    // partition
    partition(particles1, particles2, n, offsets, bbox);
    // recursive tree build
    for (int i = 0; i < (1 << DIM); i++) {
      if (counter[i] == 0) {
        node->children[i] = NULL;
      } else {
        BoundingBox child_bbox = bbox.orthant(i);
        node->children[i] = build_tree(&particles2[offsets[i]],
                                       &particles1[offsets[i]],
                                       counter[i], child_bbox, !flip);
      }
    }
  }

  return node;
}

void refine_bbox(ParticleTreeNode* node) {
  if (node->is_leaf) {
    BoundingBox bbox = get_bbox(node->particles_i, node->n_particles);
    node->inner_bbox = bbox;
    node->outer_bbox = bbox.expand(SLEN + SKIN);
  } else {
    for (int i = 0; i < (1 << DIM); i++) {
      if (ParticleTreeNode* child = node->children[i]) {
        refine_bbox(child);
        node->inner_bbox.merge(child->inner_bbox);
        node->outer_bbox.merge(child->outer_bbox);
      }
    }
  }
}

void search_neighbors_impl(ParticleTreeNode* node, ParticleTreeNode* target) {
  if (node->is_leaf) {
    target->neighbors.push_back(node);
    target->n_neighbors += node->n_particles;
  } else {
    for (int i = 0; i < (1 << DIM); i++) {
      if (ParticleTreeNode* child = node->children[i]) {
        if (target->outer_bbox.intersect(child->inner_bbox)) {
          search_neighbors_impl(child, target);
        }
      }
    }
  }
}

void search_neighbors(ParticleTreeNode* root) {
  pfor_leaf_impl(root, [&] (ParticleTreeNode* leaf) {
    leaf->neighbors.clear();
    leaf->n_neighbors = 0;
    search_neighbors_impl(root, leaf);
  });
}

void delete_nodes(ParticleTreeNode* node) {
  if (node->is_leaf) {
    delete node;
  } else {
    for (int i = 0; i < (1 << DIM); i++) {
      if (ParticleTreeNode* child = node->children[i]) {
        delete_nodes(child);
      }
    }
    delete node;
  }
}

ParticleTree::ParticleTree(const std::vector<Particle>& particles) {
  n_particles_ = particles.size();
  particles_i_ = new Particle[n_particles_];
  particles_j_ = new Particle[n_particles_];
  for (int i = 0; i < n_particles_; i++) {
    particles_i_[i] = particles[i];
  }
  root_ = NULL;
}

ParticleTree::~ParticleTree() {
  delete_nodes(root_);
#if SPH_LOOP_PARALLEL || SPH_CUDA_PARALLEL
  destroy_global_array();
#endif
  delete[] particles_i_;
  delete[] particles_j_;
}

void ParticleTree::build() {
  if (root_) {
    delete_nodes(root_);
#if SPH_LOOP_PARALLEL || SPH_CUDA_PARALLEL
    destroy_global_array();
#endif
  }
  BoundingBox bbox = get_bbox(particles_i_, n_particles_).square();
  root_ = build_tree(particles_i_, particles_j_, n_particles_, bbox, false);
  refine_bbox(root_);
  search_neighbors(root_);
#if SPH_LOOP_PARALLEL || SPH_CUDA_PARALLEL
  setup_global_array();
#endif
}

#if SPH_LOOP_PARALLEL || SPH_CUDA_PARALLEL
void ParticleTree::setup_global_array() {
  int idx = 0;
  int acc_neighbors = 0;

  for_leaf_impl(root_, [&] (ParticleTreeNode* leaf) {
    leaf->index = idx++;
    acc_neighbors += leaf->n_neighbors;
  });

  n_leafs_ = idx;
  pi_offsets_ = new int[n_leafs_ + 1];
  pj_offsets_ = new int[n_leafs_ + 1];
  leaf_array_ = new ParticleTreeNode*[n_leafs_];

  pj_buf_size_ = acc_neighbors;
  pj_buf_ = new Particle[pj_buf_size_];

  int pi_acc = 0;
  int pj_acc = 0;
  // scan (prefix sum)
  pi_offsets_[0] = 0;
  pj_offsets_[0] = 0;
  for_leaf_impl(root_, [&] (ParticleTreeNode* leaf) {
    int idx = leaf->index;
    pi_acc += leaf->n_particles;
    pi_offsets_[idx + 1] = pi_acc;
    pj_acc += leaf->n_neighbors;
    pj_offsets_[idx + 1] = pj_acc;

    leaf_array_[idx] = leaf;
  });
}

void ParticleTree::destroy_global_array() {
  delete[] pi_offsets_;
  delete[] pj_offsets_;
  delete[] pj_buf_;
  delete[] leaf_array_;
}
#endif
