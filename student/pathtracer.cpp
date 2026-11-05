#include "student.h"

#include <algorithm>
#include <cmath>
#include <limits>

// TASK 3a: a unit direction around n with density cos(theta) / pi.
vec3 cosine_sample_hemisphere(const vec3& n, float u1, float u2) {
  float radius = std::sqrt(u1);
  float phi = 2.0f * PI * u2;
  float z = std::sqrt(std::max(0.0f, 1.0f - u1));
  vec3 tangent, bitangent;
  make_basis(n, tangent, bitangent);
  return normalize(tangent * (radius * std::cos(phi)) +
                   bitangent * (radius * std::sin(phi)) + n * z);
}

// TASK 3b: the light reaching hit straight from one random point on the lamp.
vec3 sample_light(const Hit& hit, const PtScene& scene, Rng& rng) {
  float u1 = rng.next();
  float u2 = rng.next();
  vec3 point = scene.light.corner + scene.light.e1 * u1 + scene.light.e2 * u2;
  vec3 origin = hit.p + hit.n * EPS;
  vec3 to_light = point - origin;
  float distance_squared = dot(to_light, to_light);
  float area = scene.light.area();
  if(distance_squared <= 0.0f || area <= 0.0f) {
    return vec3{0, 0, 0};
  }
  float distance = std::sqrt(distance_squared);
  vec3 direction = to_light / distance;
  float cos_surface = dot(direction, hit.n);
  float cos_light = -dot(direction, scene.light.normal());
  if(cos_surface <= 0.0f || cos_light <= 0.0f) {
    return vec3{0, 0, 0};
  }
  Hit blocker;
  if(pt_intersect(Ray{origin, direction}, scene, distance - EPS, blocker)) {
    return vec3{0, 0, 0};
  }
  vec3 albedo = pt_material(scene, hit).albedo;
  return albedo * scene.light.emission *
         (cos_surface * cos_light * area / (PI * distance_squared));
}

// TASK 3c: Russian roulette (EQ-10.08): false stops the path.
bool russian_roulette(vec3& throughput, int bounce, Rng& rng) {
  if(bounce < 3) {
    return true;
  }
  float q = std::clamp(1.0f - luminance(throughput), 0.05f, 0.95f);
  if(rng.next() < q) {
    return false;
  }
  throughput /= 1.0f - q;
  return true;
}

// TASK 3d: the light arriving along r, following one path.
vec3 path_trace(const Ray& r, const PtScene& scene, Rng& rng, const PtOptions& opt) {
  vec3 throughput{1, 1, 1};
  vec3 result{0, 0, 0};
  Ray ray = r;
  for(int bounce = 0; bounce < opt.max_depth; ++bounce) {
    Hit hit;
    if(!pt_intersect(ray, scene, std::numeric_limits<float>::infinity(), hit)) {
      break;
    }
    if(dot(hit.n, ray.d) > 0.0f) {
      hit.n = -hit.n;
    }
    const PtMaterial& material = pt_material(scene, hit);
    if(material.emission.x > 0.0f || material.emission.y > 0.0f || material.emission.z > 0.0f) {
      // With NEE, later lamp hits were already accounted for by light sampling.
      if((bounce == 0 || !opt.nee) && dot(ray.d, scene.light.normal()) < 0.0f) {
        result += throughput * material.emission;
      }
      break;
    }
    if(opt.nee) {
      result += throughput * sample_light(hit, scene, rng);
    }
    float u1 = rng.next();
    float u2 = rng.next();
    vec3 direction = cosine_sample_hemisphere(hit.n, u1, u2);
    ray = Ray{hit.p + hit.n * EPS, direction};
    throughput *= material.albedo;
    if(opt.rr && !russian_roulette(throughput, bounce, rng)) {
      break;
    }
  }
  return result;
}