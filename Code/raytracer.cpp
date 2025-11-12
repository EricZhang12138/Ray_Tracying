#include "raytracer.hpp"
#include "camera.hpp"
#include "image.hpp"
#include "json_loader.hpp"
#include "shapes.hpp"
#include "light.hpp"
#include "material.hpp"
#include "acceleration.hpp"

#include <iostream>
#include <vector>
#include <memory>
#include <cmath>
#include <random>
#include <algorithm>

// === Vector Math Helpers ===
// A small namespace to hold the vector math your shade function will need.
namespace VecMath {
    inline std::array<float, 3> normalize(const std::array<float, 3>& v) {
        float mag = std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
        if (mag == 0.0f) return {0, 0, 0};
        return {v[0] / mag, v[1] / mag, v[2] / mag};
    }

    inline float dot(const std::array<float, 3>& a, const std::array<float, 3>& b) {
        return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    }

    inline std::array<float, 3> sub(const std::array<float, 3>& a, const std::array<float, 3>& b) {
        return {a[0] - b[0], a[1] - b[1], a[2] - b[2]};
    }

    inline std::array<float, 3> add(const std::array<float, 3>& a, const std::array<float, 3>& b) {
        return {a[0] + b[0], a[1] + b[1], a[2] + b[2]};
    }

    /**
     * @brief Multiplies a vector by a scalar value.
     */
    inline std::array<float, 3> mult(const std::array<float, 3>& v, float scalar) {
        return {v[0] * scalar, v[1] * scalar, v[2] * scalar};
    }

    // NEW FUNCTION: Calculates a reflection ray 
    inline Ray createReflectionRay(const Ray& ray, const Hit& hit) {
        // R = I - 2 * (I . N) * N
        // I = incoming ray direction
        // N = surface normal
        std::array<float, 3> I = ray.direction; // Already normalized from camera
        std::array<float, 3> N = hit.normal;
        float I_dot_N = dot(I, N);
        
        std::array<float, 3> R_dir = sub(I, mult(N, 2.0f * I_dot_N));
        
        // Start new ray slightly off the surface to avoid self-intersection
        std::array<float, 3> origin = add(hit.intersection_point, mult(N, 1e-4f));
        
        return {origin, R_dir}; // Direction is already normalized
    }

    // NEW FUNCTION: Calculates a refraction ray (Snell's Law)
    inline Ray createRefractionRay(const Ray& ray, const Hit& hit, float n_out) {
        std::array<float, 3> I = ray.direction; // Already normalized
        std::array<float, 3> N = hit.normal;
        float n_in = 1.0f; // Assuming ray starts in air (n=1.0)
        
        float cos_i = dot(I, N);
        
        // Handle ray hitting from *inside* the object
        if (cos_i > 0) { 
            std::swap(n_in, n_out);
            N = mult(N, -1.0f); // Normal must point against the ray
        }

        float eta = n_in / n_out;
        float cos_i_abs = std::abs(cos_i);
        
        // Check for Total Internal Reflection (TIR)
        float discriminant = 1.0f - eta * eta * (1.0f - cos_i_abs * cos_i_abs);
        if (discriminant < 0) {
            // TIR: No refraction, return an "invalid" ray
            return {{0,0,0}, {0,0,0}};
        }

        float cos_t = std::sqrt(discriminant);
        
        // T = (eta * I) + (eta * |cos_i| - cos_t) * N
        std::array<float, 3> T_dir = add(mult(I, eta), mult(N, (eta * cos_i_abs - cos_t)));
        
        // Start new ray slightly off the surface (on the *other* side)
        std::array<float, 3> origin = add(hit.intersection_point, mult(N, -1e-4f));
        
        return {origin, normalize(T_dir)};
    }

} // namespace VecMath
using namespace VecMath;


/**
 * @brief Calculates the local Blinn-Phong color at a given hit point.
 */
Color shade(
    const Hit& hit, 
    const Ray& viewRay, 
    const std::vector<Light>& lights,
    BVH& bvh // Non-const because your BVH::intersect is non-const
) {
    // 1. Get the material from the hit object
    Material mat = hit.shape->material;
    
    // 2. Get the base diffuse color for this point
    // TODO: For Task 3.3, you will replace (0,0) with hit.u, hit.v
    Color base_diffuse_color = mat.getDiffuseColor(0.0f, 0.0f); 

    // 3. Start with the Ambient term
    Color final_color = base_diffuse_color * mat.k_ambient;

    // 4. Calculate View Vector (V)
    // Vector from the hit point TO the camera (ray origin)
    std::array<float, 3> V = normalize(sub(viewRay.origin, hit.intersection_point));

    // 5. Loop through every light to add its contribution
    for (const Light& light : lights) {
        
        // 6. Calculate Light Vector (L)
        std::array<float, 3> light_vec = sub(light.location, hit.intersection_point);
        float light_distance = std::sqrt(dot(light_vec, light_vec));
        // Avoid division by zero if light is at the hit point
        if (light_distance < 1e-6f) continue; 
        std::array<float, 3> L = normalize(light_vec);

        // 7. --- Check for Shadows ---
        // Create a shadow ray from the hit point towards the light


        /*float is_normalised = hit.normal[0]*hit.normal[0] + hit.normal[1]*hit.normal[1] + hit.normal[2]*hit.normal[2];
        if (is_normalised != 1){
            std::cerr << "Error: normal is not normalised!!!!!"<< "normal is "<< is_normalised << std::endl;
            exit(1);
        }*/

        Ray shadowRay;
        // Start the ray slightly off the surface to avoid "shadow acne"
        shadowRay.origin = add(hit.intersection_point, mult(hit.normal, 1e-4f));
        shadowRay.direction = L;

        // Trace the shadow ray
        // Note: Your BVH::intersect clears its temp_hit_vec, so this is safe
        Hit shadow_hit = bvh.intersect(shadowRay, bvh.root); 
        
        // Check if the ray hit something, and if that hit is *between* us and the light
        if (shadow_hit.shape != nullptr && shadow_hit.t < light_distance && shadow_hit.t > 1e-4f) {
            continue; // This light is blocked, skip to the next light
        }

        // --- If not shadowed, continue with shading ---

        // 8. Calculate Diffuse Term
        float N_dot_L = std::max(0.0f, dot(hit.normal, L));
        Color diffuse = base_diffuse_color * N_dot_L;

        // 9. Calculate Specular Term (Blinn-Phong)
        std::array<float, 3> H = normalize(add(L, V)); // Halfway vector
        float N_dot_H = std::max(0.0f, dot(hit.normal, H));
        float spec_intensity = std::pow(N_dot_H, mat.shininess);
        Color specular = mat.specular_color * spec_intensity;

        // 10. Combine with light color and attenuation (1/d^2)
        //float attenuation = light.intensity / (light_distance * light_distance);
        // More aggressive attenuation
        float attenuation = 0.1f * light.intensity / (1.0f + light_distance + light_distance * light_distance);
        Color light_color = {light.color[0], light.color[1], light.color[2]};

        /*final_color = final_color + (light_color * (diffuse * mat.k_diffuse + specular * mat.k_specular) * attenuation);
        // Clamp each component
        final_color.r = std::min(1.0f, final_color.r);
        final_color.g = std::min(1.0f, final_color.g);
        final_color.b = std::min(1.0f, final_color.b);*/

        Color light_contribution = light_color * (diffuse * mat.k_diffuse + specular * mat.k_specular) * attenuation;

        // Clamp each component
        /*light_contribution.r = std::min(1.0f, light_contribution.r);
        light_contribution.g = std::min(1.0f, light_contribution.g);
        light_contribution.b = std::min(1.0f, light_contribution.b);*/

        final_color = final_color + light_contribution;
    }
    
    return final_color;
}


/**
 * @brief The main recursive raytracing function.
 */
Color Trace(
    const Ray& ray, 
    BVH& bvh, // Non-const to pass to shade()
    const std::vector<Light>& lights, 
    int depth
) {
    // 1. Base Case: Stop recursion
    if (depth > MAX_RECURSION_DEPTH) {
        return {0, 0, 0}; // Return black
    }

    // 2. Find closest intersection
    Hit hit = bvh.intersect(ray, bvh.root);

    // 3. If ray hits nothing, return background color
    if (!hit.shape) {
        return {0.1f, 0.1f, 0.1f}; // Dark grey background
    }

    // 4. If it hits, get the local color first
    Color localColor = shade(hit, ray, lights, bvh);

    // 5. --- Handle Recursion (The Whitted-Style part) ---
    Material mat = hit.shape->material;
    Color reflectedColor = {0, 0, 0};
    Color refractedColor = {0, 0, 0};

    // --- Reflection Calculation ---
    if (mat.reflectivity > 0.0f) {
        // 1. Calculate reflection ray
        Ray reflectionRay = VecMath::createReflectionRay(ray, hit);
        
        // 2. reflectedColor = Trace(...)
        reflectedColor = Trace(reflectionRay, bvh, lights, depth + 1);
    }

    // --- Refraction Calculation ---
    if (mat.transparency > 0.0f) {
        // 1. Calculate refraction ray (using Snell's Law)
        Ray refractionRay = VecMath::createRefractionRay(ray, hit, mat.refractive_index);

        // Check for valid ray (i.e., no Total Internal Reflection)
        // We check direction, as origin might be valid
        if (VecMath::dot(refractionRay.direction, refractionRay.direction) > 1e-6f) {
            // 2. refractedColor = Trace(...)
            refractedColor = Trace(refractionRay, bvh, lights, depth + 1);
        }
    }

        // 6. Return the final combined color
        // Scale localColor by the energy that isn't reflected or refracted
        float local_contribution = 1.0f - mat.reflectivity - mat.transparency;

        // Clamp this to 0 in case reflectivity/transparency add up to > 1
        local_contribution = std::max(0.0f, local_contribution); 

        return (local_contribution * localColor) + 
            (mat.reflectivity * reflectedColor) + 
            (mat.transparency * refractedColor);
}


/**
 * @brief Main function to run the raytracer.
 */
int main() {
    try {
        std::string scene_file = "/Users/ericzhang/Documents/Computer_graphics/Coursework/s2286795/ASCII/scene.json";

        // 1. Load Camera (Module 1)
        Camera camera(scene_file);
        // Note: You must implement getResolution() to return width/height
        // Assuming getResolution() returns std::tuple<int, int>
        auto [width, height] = camera.getResolution(); 
        if (width == 0 || height == 0) {
            std::cerr << "Error: Camera resolution is 0. Check scene.json." << std::endl;
            return 1;
        }
        Image outputImage(width, height);

        // 2. Load Scene Data (Modules 1, 2, 3)
        std::vector<Light> lights = load_lights_from_json(scene_file);
        std::vector<std::unique_ptr<Shapes>> shapes = load_shapes_from_json(scene_file);
        
        // 3. Build Acceleration Structure (Module 2)
        // Create a vector of raw pointers for the BVH
        std::vector<Shapes*> shape_ptrs;
        for (const auto& shape : shapes) {
            shape_ptrs.push_back(shape.get());
        }
        if (shape_ptrs.empty()) {
            std::cerr << "Warning: No shapes loaded to render." << std::endl;
        }
        
        BVH bvh(shape_ptrs);
        // We must manually pass bvh.root to intersect, as per your acceleration.cpp
        // You MUST make 'root' public in BVH or add a getRoot() method.
        // For now, I'll assume it's public.
        std::cout << "BVH built." << std::endl;

        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<double> dist(0.0, 1.0);
        // 4. Main Render Loop
        std::cout << "Rendering " << width << "x" << height << "..." << std::endl;
        for (int y = 0; y < height; ++y) {
            for (int x = 0; x < width; ++x) {
                

            // 1. Define the grid size. (e.g., 3x3 grid = 9 total samples)
            const int n_samples_sqrt = 3; 
            const int total_samples = n_samples_sqrt * n_samples_sqrt;

            // 2. Initialize accumulator color to black
            // (Your struct default {0,0,0} handles this, but being explicit is good)
            Color totalColor = {0.0f, 0.0f, 0.0f};

            // 3. Loop over the grid cells
            for (int j = 0; j < n_samples_sqrt; ++j) { // y-cell
                for (int i = 0; i < n_samples_sqrt; ++i) { // x-cell
                    
                    // Get a random offset *within* this cell
                    double random_offset_x = dist(gen); // random val [0.0, 1.0)
                    double random_offset_y = dist(gen); // random val [0.0, 1.0)

                    // Calculate the sample's position in the pixel
                    // (i + rand_x) / 3.0 gives a jittered value in the i-th column
                    double sample_x = (i + random_offset_x) / n_samples_sqrt;
                    double sample_y = (j + random_offset_y) / n_samples_sqrt;

                    // Generate the ray for this specific sample
                    auto [origin, direction] = camera.pixelToRay({x + sample_x, y + sample_y});
                    Ray ray = {origin, direction};

                    // b. Get Color (Module 3) and accumulate
                    // (Using your Color::operator+)
                    totalColor = totalColor + Trace(ray, bvh, lights, 0);
                }
            }

            // 4. Average the final color
            // (Using your Color::operator/)
            Color Averaged_color = totalColor / (float)total_samples;

            /*// Apply gamma correction (power of 1.0/2.2)
            // This converts the linear color to sRGB space for your monitor.
            Averaged_color.r = std::pow(Averaged_color.r, 1.0f / 2.2f);
            Averaged_color.g = std::pow(Averaged_color.g, 1.0f / 2.2f);
            Averaged_color.b = std::pow(Averaged_color.b, 1.0f / 2.2f);*/

            // 5. Clamp the final averaged color
            Averaged_color.r = std::max(0.0f, std::min(1.0f, Averaged_color.r));
            Averaged_color.g = std::max(0.0f, std::min(1.0f, Averaged_color.g));
            Averaged_color.b = std::max(0.0f, std::min(1.0f, Averaged_color.b));

            // 6. Write to the image
            outputImage.setPixel(x, y, 
                static_cast<int>(Averaged_color.r * 255.999),
                static_cast<int>(Averaged_color.g * 255.999),
                static_cast<int>(Averaged_color.b * 255.999)
            );
            }
            // Print progress
            if (y > 0 && y % 100 == 0) {
                std::cout << "Progress: " << (100 * y / height) << "%\n";
            }
        }
        std::cout << "Rendering complete.\n";

        // 5. Save Final Image (Module 1)
        outputImage.write("output.ppm");

    } catch (const std::exception& e) {
        std::cerr << "An error occurred: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}