#include "raytracer.hpp"
#include "camera.hpp"
#include "image.hpp"
#include "json_loader.hpp"
#include "shapes.hpp"
#include "light.hpp"
#include "material.hpp"
#include "acceleration.hpp"
#include <cstring> // for strcmp
#include <iostream>
#include <vector>
#include <memory>
#include <cmath>
#include <algorithm>



Color compute_pixel_color(
    int x, int y, 
    int samples_sqrt, 
    Camera& camera, 
    BVH& bvh, 
    const std::vector<Light>& lights, 
    bool use_bvh,
    std::mt19937& gen,
    std::uniform_real_distribution<double>& dist,
    int light_samples
) {
    // Optimization: If samples is 1, shoot one ray through the center
    if (samples_sqrt <= 1) {
        // Pixel center is x + 0.5, y + 0.5
        //auto [origin, direction] = camera.pixelToRay({x + 0.5f, y + 0.5f});
        auto [origin, direction] = camera.pixelToRay_thin_lens({x + 0.5f, y + 0.5f}, gen, dist);

        Ray ray = {origin, direction};

        ray.time = (float)dist(gen);

        return Trace(ray, bvh, lights, 0, use_bvh, gen, dist, light_samples);
    }

    Color totalColor = {0.0f, 0.0f, 0.0f};
    int total_samples = samples_sqrt * samples_sqrt;

    // Stratified Sampling Loop (Jittering)
    for (int j = 0; j < samples_sqrt; ++j) { // y-subpixel
        for (int i = 0; i < samples_sqrt; ++i) { // x-subpixel
            
            // Get a random offset [0, 1)
            double random_offset_x = dist(gen); 
            double random_offset_y = dist(gen); 

            // Map to the specific sub-grid
            double sample_x = (i + random_offset_x) / samples_sqrt;
            double sample_y = (j + random_offset_y) / samples_sqrt;

            // Generate ray
            auto [origin, direction] = camera.pixelToRay_thin_lens({x + sample_x, y + sample_y}, gen, dist);
            //auto [origin, direction] = camera.pixelToRay({x + sample_x, y + sample_y});
            Ray ray = {origin, direction};
            ray.time = (float)dist(gen);

            // Accumulate
            totalColor = totalColor + Trace(ray, bvh, lights, 0, use_bvh, gen, dist, light_samples);
        }
    }

    // Average the result
    return totalColor / (float)total_samples;
}

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
    // Used for Light sampling and Glossy Reflection fuzz
    inline std::array<float, 3> random_in_unit_sphere(std::mt19937& gen, std::uniform_real_distribution<double>& dist) {
        while (true) {
            // Generate x, y, z in range [-1, 1]
            float r1 = (float)dist(gen);
            float r2 = (float)dist(gen);
            float r3 = (float)dist(gen);
            
            std::array<float, 3> p = {
                2.0f * r1 - 1.0f, 
                2.0f * r2 - 1.0f, 
                2.0f * r3 - 1.0f
            };

            // Check if point is inside sphere (length squared < 1)
            if (dot(p, p) < 1.0f) {
                return p;
            }
            // If outside, try again (Rejection Sampling)
        }
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
    BVH& bvh, 
    bool use_bvh,
    std::mt19937& gen,                      // Passed by Reference
    std::uniform_real_distribution<double>& dist, // Passed by Reference
    int light_samples
) {
    Material mat = hit.shape->material;
    Color base_diffuse_color = mat.getDiffuseColor(hit.u, hit.v); 

    // Ambient Term
    Color final_color = base_diffuse_color * mat.k_ambient;

    // View Vector (V)
    std::array<float, 3> V = normalize(sub(viewRay.origin, hit.intersection_point));

    for (const Light& light : lights) {
        
        // --- SOFT SHADOW CALCULATION ---
        float visibility = 0.0f;
        
        // Settings: More samples = smoother shadows but slower
        // For 'distributed' raytracing, we usually take multiple samples per light.
        // If radius is 0, we can optimize to 1 sample (hard shadow).
        int shadow_samples = (light.radius > 0.0f) ? light_samples : 1; 

        for (int s = 0; s < shadow_samples; ++s) {
            std::array<float, 3> target_light_pos = light.location;

            // If light has radius, jitter the target position
            if (light.radius > 0.0f) {
                std::array<float, 3> random_offset = VecMath::random_in_unit_sphere(gen, dist);
                // Scale offset by radius
                random_offset = mult(random_offset, light.radius);
                target_light_pos = add(target_light_pos, random_offset);
            }

            // Vector to this specific point on the light
            std::array<float, 3> light_vec = sub(target_light_pos, hit.intersection_point);
            float light_dist = std::sqrt(dot(light_vec, light_vec));
            std::array<float, 3> L_sample = normalize(light_vec);

            // Shadow Ray
            Ray shadowRay;
            shadowRay.origin = add(hit.intersection_point, mult(hit.normal, 1e-4f));
            shadowRay.direction = L_sample;

            Hit shadow_hit = bvh.get_intersection(shadowRay, use_bvh); 
            
            // Check visibility
            if (!shadow_hit.shape || shadow_hit.t > light_dist) {
                visibility += 1.0f;
            }
        }

        // Average the visibility (0.0 to 1.0)
        visibility /= (float)shadow_samples;

        // Optimization: If completely in shadow, skip shading math
        if (visibility <= 0.0f) continue;

        // --- STANDARD BLINN-PHONG (Using center of light) ---
        // We calculate lighting based on the main light direction to keep it stable
        std::array<float, 3> light_vec_center = sub(light.location, hit.intersection_point);
        float dist_sq = dot(light_vec_center, light_vec_center);
        float light_distance = std::sqrt(dist_sq);
        std::array<float, 3> L = normalize(light_vec_center);

        // Diffuse
        float N_dot_L = std::max(0.0f, dot(hit.normal, L));
        Color diffuse = base_diffuse_color * N_dot_L;

        // Specular
        std::array<float, 3> H = normalize(add(L, V)); 
        float N_dot_H = std::max(0.0f, dot(hit.normal, H));
        float spec_intensity = std::pow(N_dot_H, mat.shininess);
        Color specular = mat.specular_color * spec_intensity;

        // Attenuation
        float attenuation = 0.15f * light.intensity / (1.0f + light_distance + dist_sq);
        Color light_color_vec = {light.color[0], light.color[1], light.color[2]};

        // Combine
        Color contribution = light_color_vec * (diffuse * mat.k_diffuse + specular * mat.k_specular) * attenuation;

        // APPLY VISIBILITY
        final_color = final_color + (contribution * visibility);
    }
    
    return final_color;
}


/**
 * @brief The main recursive raytracing function.
 */
Color Trace(
    const Ray& ray, 
    BVH& bvh, 
    const std::vector<Light>& lights, 
    int depth,
    bool use_bvh,
    std::mt19937& gen,                      // Passed by Reference
    std::uniform_real_distribution<double>& dist, // Passed by Reference
    int light_samples
) {
    if (depth > MAX_RECURSION_DEPTH) {
        return {0, 0, 0}; 
    }

    Hit hit = bvh.get_intersection(ray, use_bvh);

    if (!hit.shape) {
        return {0.1f, 0.1f, 0.1f}; // Background
    }

    // 1. Local Shading (includes Soft Shadows now)
    Color localColor = shade(hit, ray, lights, bvh, use_bvh, gen, dist, light_samples);

    Material mat = hit.shape->material;
    Color reflectedColor = {0, 0, 0};
    Color refractedColor = {0, 0, 0};

    // 2. Reflection (Glossy)
    if (mat.reflectivity > 0.0f) {
        Ray reflectionRay = VecMath::createReflectionRay(ray, hit);
        
        // --- GLOSSY PERTURBATION ---
        if (mat.roughness > 0.0f) {
            std::array<float, 3> fuzz = VecMath::random_in_unit_sphere(gen, dist);
            
            // Add fuzz scaled by roughness
            // Result = Normalized(PerfectVector + (RandomSphere * Roughness))
            std::array<float, 3> perturbed = add(reflectionRay.direction, mult(fuzz, mat.roughness));
            reflectionRay.direction = normalize(perturbed);

            // Check: Did we perturb it so much it points INTO the surface?
            // If dot(Dir, Normal) < 0, it's going inside. 
            if (dot(reflectionRay.direction, hit.normal) < 0.0f) {
                // Determine what to do: absorb the ray (black) or re-sample.
                // For simplicity, we just absorb it (return black for this bounce).
                reflectionRay.direction = {0,0,0}; // Invalid ray
            }
        }

        // Only trace if direction is valid
        if (dot(reflectionRay.direction, reflectionRay.direction) > 0.001f) {
            reflectedColor = Trace(reflectionRay, bvh, lights, depth + 1, use_bvh, gen, dist, light_samples);
        }
    }

    // 3. Refraction (Standard - Glossy transmission is harder, keeping it simple)
    if (mat.transparency > 0.0f) {
        Ray refractionRay = VecMath::createRefractionRay(ray, hit, mat.refractive_index);
        
        // Check for TIR (invalid ray)
        if (dot(refractionRay.direction, refractionRay.direction) > 1e-6f) {
             // Pass gen/dist down, even if we don't use them for glossy refraction yet
            refractedColor = Trace(refractionRay, bvh, lights, depth + 1, use_bvh, gen, dist, light_samples);
        }
    }

    float local_contribution = std::max(0.0f, 1.0f - mat.reflectivity - mat.transparency);

    return (local_contribution * localColor) + 
           (mat.reflectivity * reflectedColor) + 
           (mat.transparency * refractedColor);
}




int main(int argc, char* argv[]) {
    try {
        std::string scene_file = "/Users/ericzhang/Documents/Computer_graphics/Coursework/s2286795/ASCII/scene.json";
        
        // --- Command Line Parsing ---
        bool use_bvh = false;    // Default: false (as per your code)
        int n_samples_sqrt = 4;  // Default: 4x4 samples (as per your code)
        int light_samples = 1;

        for (int i = 1; i < argc; ++i) {
            // Check for BVH flag
            if (std::strcmp(argv[i], "-bvh") == 0) {
                use_bvh = true;
            }
            // Check for Samples flag
            else if (std::strcmp(argv[i], "-s") == 0 && i + 1 < argc) {
                n_samples_sqrt = std::atoi(argv[i + 1]);
                i++; // Skip next arg
            }
            else if(std::strcmp(argv[i], "-light_sample") == 0 && i+1 < argc){
                light_samples = std::atoi(argv[i+1]);
                i++;
            }
        }
        
        // 1. Load Camera
        Camera camera(scene_file);
        auto [width, height] = camera.getResolution(); 
        if (width == 0 || height == 0) {
            std::cerr << "Error: Camera resolution is 0. Check scene.json." << std::endl;
            return 1;
        }
        Image outputImage(width, height);

        // 2. Load Scene Data
        std::vector<Light> lights = load_lights_from_json(scene_file);
        std::vector<std::unique_ptr<Shapes>> shapes = load_shapes_from_json(scene_file);
        
        // 3. Build Acceleration Structure
        std::vector<Shapes*> shape_ptrs;
        for (const auto& shape : shapes) {
            shape_ptrs.push_back(shape.get());
        }
        if (shape_ptrs.empty()) {
            std::cerr << "Warning: No shapes loaded to render." << std::endl;
        }
        
        BVH bvh(shape_ptrs);
        std::cout << "BVH built. Mode: " << (use_bvh ? "ON" : "OFF") << std::endl;

        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<double> dist(0.0, 1.0);

        // 4. Main Render Loop
        std::cout << "Rendering " << width << "x" << height << " with " 
                  << n_samples_sqrt << "x" << n_samples_sqrt << " samples..." << std::endl;

        for (int y = 0; y < height; ++y) {
            for (int x = 0; x < width; ++x) {
                
                // CALL NEW FUNCTION HERE
                Color Averaged_color = compute_pixel_color(
                    x, y, 
                    n_samples_sqrt, 
                    camera, bvh, lights, use_bvh, 
                    gen, dist, light_samples
                );

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

        // 5. Save Final Image
        outputImage.write("output.ppm");

    } catch (const std::exception& e) {
        std::cerr << "An error occurred: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}