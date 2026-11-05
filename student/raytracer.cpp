#include "student.h"

#include <algorithm>
#include <cmath>

// TASK 1a: the ray from the camera through the centre of pixel (x, y).
Ray camera_ray(int x, int y, int width, int height, float fov) {
    float t = std::tan(fov / 2.0f);
    float aspect = static_cast<float>(width) / height;
    vec3 d = {(2 * (x + 0.5f) / width - 1) * t * aspect,
              (1 - 2 * (y + 0.5f) / height) * t, -1};
    return Ray{vec3{0, 0, 0}, normalize(d)};
}

// TASK 1b: the smallest t > EPS where r meets the sphere (EQ-10.02).
bool ray_sphere(const Ray& r, const Sphere& s, float& t) {
  float a = dot(r.d, r.d);
    float b = 2 * dot(r.d, r.o - s.center);
    float c = dot(r.o - s.center, r.o - s.center) - s.radius * s.radius;
    float discriminant = b * b - 4 * a * c;
    if (discriminant < 0) {
        return false;
    }
    float sqrt_discriminant = std::sqrt(discriminant);
    float t1 = (-b - sqrt_discriminant) / (2 * a);
    float t2 = (-b + sqrt_discriminant) / (2 * a);
    if (t1 > EPS) {
        t = t1;
        return true;
    }
    if (t2 > EPS) {
        t = t2;
        return true;
    }
    return false;
}


// TASK 1c: d reflected about the unit normal n (EQ-10.04).
vec3 reflect(const vec3& d, const vec3& n) {
  return d - 2.0f * dot(d, n) * n;
}

// TASK 1c: d refracted through a surface whose normal n points out of the object.
// false on total internal reflection.
bool refract(const vec3& d, const vec3& n, float ior, vec3& out) {
  float cos_i = std::clamp(-dot(d, n), -1.0f, 1.0f);
  float eta = 1.0f / ior;
  vec3 normal = n;
  if(cos_i < 0.0f) {
    // Leaving the object: face the normal inward and exchange the two media.
    cos_i = -cos_i;
    normal = -n;
    eta = ior;
  }
  float k = 1.0f - eta * eta * (1.0f - cos_i * cos_i);
  if(k < 0.0f) {
    return false;
  }
  out = eta * d + (eta * cos_i - std::sqrt(k)) * normal;
  return true;
}

// TASK 1d: the colour seen along r: Phong, shadows, reflection, refraction (Whitted).
vec3 cast_ray(const Ray& r, const Scene& scene, int depth) {
  if(depth > MAX_DEPTH) {
    return scene.background;
  }
  Hit hit;
  if(!scene_intersect(r, scene, hit)) {
    return scene.background;
  }
  const Material& material = *hit.material;
  if(STAGE < 2) {
    return material.diffuse;
  }

  // Offset every new ray to the side of the surface that its direction enters.
  auto spawn_ray = [&](const vec3& direction) {
    vec3 offset = hit.n * EPS;
    return Ray{dot(direction, hit.n) < 0.0f ? hit.p - offset : hit.p + offset,
               direction};
  };

  float diffuse = 0.0f;
  float specular = 0.0f;
  for(const Light& light : scene.lights) {
    vec3 to_light = light.position - hit.p;
    float light_distance = length(to_light);
    vec3 light_direction = to_light / light_distance;
    if(STAGE >= 4) {
      Hit blocker;
      if(scene_intersect(spawn_ray(light_direction), scene, blocker) &&
         blocker.t < light_distance) {
        continue;
      }
    }
    diffuse += light.intensity * std::max(0.0f, dot(light_direction, hit.n));
    if(STAGE >= 3) {
      float alignment = std::max(0.0f, -dot(reflect(-light_direction, hit.n), r.d));
      specular += light.intensity * std::pow(alignment, material.specular_exponent);
    }
  }

  vec3 color = material.diffuse * diffuse * material.albedo[0] +
               vec3{1, 1, 1} * specular * material.albedo[1];
  if(STAGE >= 5) {
    vec3 direction = normalize(reflect(r.d, hit.n));
    color += cast_ray(spawn_ray(direction), scene, depth + 1) * material.albedo[2];
  }
  if(STAGE >= 6) {
    vec3 direction;
    if(refract(r.d, hit.n, material.ior, direction)) {
      direction = normalize(direction);
      color += cast_ray(spawn_ray(direction), scene, depth + 1) * material.albedo[3];
    }
  }
  return color;
}