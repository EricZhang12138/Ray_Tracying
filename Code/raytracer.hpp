#pragma once

#include "shapes.hpp"
#include "material.hpp"
#include "light.hpp"
#include "acceleration.hpp"
#include <vector>
#include <random>

// Define a maximum recursion depth to prevent stack overflow
const int MAX_RECURSION_DEPTH = 10;

/**
 * @brief Calculates the local Blinn-Phong color at a given hit point.
 * This function also handles shadow ray casting.
 * @param hit The intersection data (point, normal, shape).
 * @param viewRay The ray that caused the hit (to find the view direction).
 * @param lights A list of all lights in the scene.
 * @param bvh The acceleration structure, used for casting shadow rays.
 * @return The calculated Color (ambient + diffuse + specular).
 */
Color shade(
    const Hit& hit, 
    const Ray& viewRay, 
    const std::vector<Light>& lights, 
    BVH& bvh, // Non-const because your BVH::intersect is non-const
    bool use_bvh,
    std::mt19937& gen,                      // <--- ADD THIS
    std::uniform_real_distribution<double>& dist,
    int light_samples
);

/**
 * @brief The main recursive raytracing function.
 * Finds the closest intersection and combines local color (from shade)
 * with reflected and refracted colors.
 * @param ray The ray to trace into the scene.
 * @param bvh The acceleration structure to test against.
 * @param lights A list of all lights (for the shade function).
 * @param depth The current recursion depth.
 * @return The final calculated Color for the ray.
 */
Color Trace(
    const Ray& ray, 
    BVH& bvh, // Non-const to pass to shade()
    const std::vector<Light>& lights, 
    int depth,
    bool use_bvh,
    std::mt19937& gen,                     
    std::uniform_real_distribution<double>& dist,
    int light_samples
);