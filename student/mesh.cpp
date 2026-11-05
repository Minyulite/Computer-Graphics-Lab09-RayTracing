#include "student.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>

#include "core/bvh.h"

// TASK 2a: Moller-Trumbore (DEF-10.01): does r hit the triangle at t > EPS? Write t, u, v.
bool ray_triangle(const Ray& r, const vec3& v0, const vec3& v1, const vec3& v2,
                  float& t, float& u, float& v) {
  vec3 e1 = v1 - v0;
  vec3 e2 = v2 - v0;
  vec3 h = cross(r.d, e2);
  float a = dot(e1, h);
  if(std::fabs(a) < 1e-8f) {
    return false;
  }
  float inv_a = 1.0f / a;
  vec3 s = r.o - v0;
  float hit_u = dot(s, h) * inv_a;
  if(hit_u < 0.0f || hit_u > 1.0f) {
    return false;
  }
  vec3 q = cross(s, e1);
  float hit_v = dot(r.d, q) * inv_a;
  if(hit_v < 0.0f || hit_u + hit_v > 1.0f) {
    return false;
  }
  float hit_t = dot(e2, q) * inv_a;
  if(hit_t <= EPS) {
    return false;
  }
  t = hit_t;
  u = hit_u;
  v = hit_v;
  return true;
}

// Both traversal methods use the same closest-hit and normal calculation.
static bool test_mesh_triangle(const Ray& r, const Mesh& m, int index,
                               float& closest, Hit& hit) {
  count_triangle_test();
  const auto& tri = m.tris[size_t(index)];
  const vec3& v0 = m.positions[size_t(tri[0])];
  const vec3& v1 = m.positions[size_t(tri[1])];
  const vec3& v2 = m.positions[size_t(tri[2])];
  float t, u, v;
  if(!ray_triangle(r, v0, v1, v2, t, u, v) || t >= closest) {
    return false;
  }
  closest = t;
  hit.t = t;
  hit.tri = index;
  if(m.normals.empty()) {
    hit.n = normalize(cross(v1 - v0, v2 - v0));
  } else {
    hit.n = normalize(m.normals[size_t(tri[0])] * (1.0f - u - v) +
                      m.normals[size_t(tri[1])] * u + m.normals[size_t(tri[2])] * v);
  }
  return true;
}

// TASK 2b: the closest triangle of m below t_max, testing every one (count_triangle_test()).
// Fill hit.t, hit.tri and hit.n.
bool intersect_mesh(const Ray& r, const Mesh& m, float t_max, Hit& hit) {
  float closest = t_max;
  bool found = false;
  for(size_t i = 0; i < m.tris.size(); ++i) {
    if(test_mesh_triangle(r, m, int(i), closest, hit)) {
      found = true;
    }
  }
  return found;
}

// TASK 2c (optional): the slab test (EQ-10.03). inv_d is 1 / r.d; write the entry distance.
bool ray_aabb(const Ray& r, const vec3& inv_d, const AABB& box, float t_max, float& t_entry) {
  float entry = -std::numeric_limits<float>::infinity();
  float exit = std::numeric_limits<float>::infinity();
  for(int axis = 0; axis < 3; ++axis) {
    // A parallel ray must lie within the slab; do not multiply zero by infinity.
    if(r.d[axis] == 0.0f) {
      if(r.o[axis] < box.lo[axis] || r.o[axis] > box.hi[axis]) {
        return false;
      }
      continue;
    }
    float t0 = (box.lo[axis] - r.o[axis]) * inv_d[axis];
    float t1 = (box.hi[axis] - r.o[axis]) * inv_d[axis];
    entry = std::max(entry, std::min(t0, t1));
    exit = std::min(exit, std::max(t0, t1));
    if(entry > exit) {
      return false;
    }
  }
  if(exit <= 0.0f || entry >= t_max) {
    return false;
  }
  t_entry = entry;
  return true;
}

// TASK 2d (optional): TASK 2b through m.bvh (count_node_visit() for every node).
bool intersect_bvh(const Ray& r, const Mesh& m, float t_max, Hit& hit) {
  if(m.bvh.nodes.empty()) {
    return false;
  }
  vec3 inv_d{1.0f / r.d.x, 1.0f / r.d.y, 1.0f / r.d.z};
  float entry;
  if(!ray_aabb(r, inv_d, m.bvh.nodes[0].box, t_max, entry)) {
    count_node_visit();
    return false;
  }

  struct Pending {
    int index;
    float entry;
  };
  // The supplied builder limits depth to BVH_MAX_DEPTH.
  std::array<Pending, BVH_MAX_DEPTH + 1> stack;
  int size = 0;
  stack[size++] = Pending{0, entry};
  float closest = t_max;
  bool found = false;
  while(size > 0) {
    Pending next = stack[--size];
    count_node_visit();
    if(next.entry >= closest) {
      continue;
    }
    const BVHNode& node = m.bvh.nodes[size_t(next.index)];
    if(node.count > 0) {
      for(int i = 0; i < node.count; ++i) {
        int tri = m.bvh.order[size_t(node.first + i)];
        if(test_mesh_triangle(r, m, tri, closest, hit)) {
          found = true;
        }
      }
      continue;
    }

    int left = next.index + 1;
    int right = node.right;
    float left_entry, right_entry;
    bool hit_left = ray_aabb(r, inv_d, m.bvh.nodes[size_t(left)].box, closest, left_entry);
    bool hit_right = ray_aabb(r, inv_d, m.bvh.nodes[size_t(right)].box, closest, right_entry);
    if(hit_left && hit_right) {
      // Push the farther child first so the nearer child is processed next.
      if(left_entry <= right_entry) {
        stack[size++] = Pending{right, right_entry};
        stack[size++] = Pending{left, left_entry};
      } else {
        stack[size++] = Pending{left, left_entry};
        stack[size++] = Pending{right, right_entry};
      }
    } else if(hit_left) {
      stack[size++] = Pending{left, left_entry};
    } else if(hit_right) {
      stack[size++] = Pending{right, right_entry};
    }
  }
  return found;
}