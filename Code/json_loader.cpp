#include "json_loader.hpp" 

// --- Standard Library Includes ---
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <array>
#include <memory>      // For std::unique_ptr and std::make_shared
#include <stdexcept>   // For std::runtime_error
#include <algorithm>   // For std::max

   
#include "light.hpp"       
#include "shapes.hpp"       
#include "material.hpp"     
#include "image.hpp"        // Image class (for texture loading)

// Use nlohmann::json
using json = nlohmann::json;

/**
 * @brief (HELPER FUNCTION)
 * Parses the JSON material block and returns a C++ Material object.
 * This is called by load_shapes_from_json.
 *
 * @param mat_json The nlohmann::json object for the "material" block.
 * @return A configured Material object.
 */
Material parse_material(const json& mat_json) {
    Material mat; // Starts with default values from material.hpp

    try {
        // --- 1. Blinn-Phong Properties ---
        if (mat_json.contains("diffuse_color")) {
            auto dc = mat_json["diffuse_color"].get<std::array<float, 3>>();
            mat.diffuse_color = {dc[0], dc[1], dc[2]};
        }
        if (mat_json.contains("specular_color")) {
            auto sc = mat_json["specular_color"].get<std::array<float, 3>>();
            mat.specular_color = {sc[0], sc[1], sc[2]};
        }

        // Use .value() to provide a default if the key is missing
        mat.k_ambient = mat_json.value("k_ambient", 0.1f);
        mat.k_diffuse = mat_json.value("k_diffuse", 0.9f);
        mat.k_specular = mat_json.value("k_specular", 0.3f);
        
        // --- Convert Roughness to Shininess ---
        // Your exporter saves 'roughness', but Blinn-Phong needs 'shininess'
        // We use 1/roughness^2 for a simple conversion.
        float roughness = mat_json.value("roughness", 0.5f);
        roughness = std::max(0.01f, roughness); // Avoid divide-by-zero
        mat.shininess = 500000.0f / (roughness * roughness);

        // --- 2. Whitted-Style Properties ---
        mat.reflectivity = mat_json.value("reflectivity", 0.0f);
        mat.transparency = mat_json.value("transparency", 0.0f);
        mat.refractive_index = mat_json.value("refractive_index", 1.0f);

        // --- 3. Texture Loading ---
        if (mat_json.contains("texture_file") && !mat_json["texture_file"].empty()) {
            std::string texture_filename = mat_json["texture_file"];
            if (!texture_filename.empty()) {
                // This is where the Image class is used
                // It creates a new Image by reading the file
                // and stores it in the material's shared_ptr.
                texture_filename = "../../Textures/" + texture_filename;
                mat.texture = std::make_shared<Image>(texture_filename);

                if (mat.texture->getWidth() == 0 || !mat.texture->loaded_successfully) {
                std::cerr << "Warning: Failed to load texture file: " << texture_filename << std::endl;
                mat.texture = nullptr; // Reset to nullptr if loading failed
            }
            }
        }

    } catch (json::exception& e) {
        std::cerr << "Warning: Error parsing material data: " << e.what() << std::endl;
        // If it fails, just return the default material
        return Material();
    }

    return mat;
}


/**
 * @brief Loads all Light objects from the scene JSON file.
 */
std::vector<Light> load_lights_from_json(const std::string& filename) {
    std::vector<Light> loaded_lights;

    std::ifstream json_file(filename);
    if (!json_file.is_open()) {
        throw std::runtime_error("Error: Could not open JSON file: " + filename);
    }

    json scene_data;
    try {
        json_file >> scene_data;
    } catch (json::parse_error& e) {
        throw std::runtime_error(std::string("Error: JSON parsing error: ") + e.what());
    }

    // --- Process Lights ---
    if (scene_data.contains("lights") && scene_data["lights"].is_array()) {
        for (const auto& light_json : scene_data["lights"]) {
             if (!light_json.is_object()) {
                 std::cerr << "Warning: Skipping non-object entry in 'lights' array." << std::endl;
                 continue;
             }
            try {
                if (!light_json.contains("location") || !light_json.contains("color") || !light_json.contains("intensity")) {
                     std::cerr << "Warning: Skipping invalid light definition." << std::endl;
                     continue;
                }
                
                std::array<float, 3> location = light_json["location"].get<std::array<float, 3>>();
                std::array<float, 3> color = light_json["color"].get<std::array<float, 3>>();
                float intensity = light_json["intensity"].get<float>();

                if (intensity <= 0) {
                     std::cerr << "Warning: Skipping light with non-positive intensity." << std::endl;
                     continue;
                }
                
                loaded_lights.emplace_back(location, color, intensity);

            } catch (json::exception& e) {
                std::cerr << "Warning: Error parsing light entry: " << e.what() << std::endl;
            }
        }
    } else if (scene_data.contains("lights")) {
         std::cerr << "Warning: 'lights' key found but is not an array in " << filename << ". No lights loaded." << std::endl;
    }

    if (loaded_lights.empty()) {
        std::cerr << "Warning: No valid lights were loaded from " << filename << "." << std::endl;
    }

    return loaded_lights;
}


/**
 * @brief Loads all Shape objects from the scene JSON file.
 */
std::vector<std::unique_ptr<Shapes>> load_shapes_from_json(const std::string& filename) {
    std::vector<std::unique_ptr<Shapes>> loaded_shapes;

    std::ifstream json_file(filename);
    if (!json_file.is_open()) {
        throw std::runtime_error("Error: Could not open JSON file: " + filename);
    }

    json scene_data;
    try {
        json_file >> scene_data;
    } catch (json::parse_error& e) {
        throw std::runtime_error(std::string("Error: JSON parsing error: ") + e.what());
    }

    // --- Process Spheres ---
    if (scene_data.contains("spheres") && scene_data["spheres"].is_array()) {
        for (const auto& sphere_json : scene_data["spheres"]) {
             if (!sphere_json.is_object()) continue;
            try {
                if (!sphere_json.contains("location") || !sphere_json.contains("radius")) {
                     std::cerr << "Warning: Skipping invalid sphere definition." << std::endl;
                     continue;
                }
                std::array<float, 3> center = sphere_json["location"].get<std::array<float, 3>>();
                float radius = sphere_json["radius"].get<float>();
                if (radius <= 0) {
                     std::cerr << "Warning: Skipping sphere with non-positive radius." << std::endl;
                     continue;
                }

                // --- MODIFIED SECTION ---
                Material mat; // Create a default material
                if (sphere_json.contains("material")) {
                    mat = parse_material(sphere_json["material"]); // Parse the real one
                }
                
                // Pass the material to the constructor
                loaded_shapes.push_back(std::make_unique<Sphere>(center, radius, mat));
                // --- END MODIFIED SECTION ---

            } catch (json::exception& e) {
                std::cerr << "Warning: Error parsing sphere entry: " << e.what() << std::endl;
            }
        }
    }

    // --- Process Cubes ---
    if (scene_data.contains("cubes") && scene_data["cubes"].is_array()) {
        for (const auto& cube_json : scene_data["cubes"]) {
             if (!cube_json.is_object()) continue;
            try {
                 if (!cube_json.contains("translation") || !cube_json.contains("rotation") || !cube_json.contains("scale")) {
                     std::cerr << "Warning: Skipping invalid cube definition." << std::endl;
                     continue;
                 }
                std::array<float, 3> translation = cube_json["translation"].get<std::array<float, 3>>();
                std::array<float, 3> rotation = cube_json["rotation"].get<std::array<float, 3>>();
                float scale = cube_json["scale"].get<float>();
                 if (scale <= 0) {
                     std::cerr << "Warning: Skipping cube with non-positive scale." << std::endl;
                     continue;
                 }

                // --- MODIFIED SECTION ---
                Material mat; // Create a default material
                if (cube_json.contains("material")) {
                    mat = parse_material(cube_json["material"]); // Parse the real one
                }

                // Pass the material to the constructor
                loaded_shapes.push_back(std::make_unique<Cube>(translation, rotation, scale, mat));
                // --- END MODIFIED SECTION ---

            } catch (json::exception& e) {
                 std::cerr << "Warning: Error parsing cube entry: " << e.what() << std::endl;
            }
        }
    }

    // --- Process Planes ---
    if (scene_data.contains("planes") && scene_data["planes"].is_array()) {
        for (const auto& plane_json : scene_data["planes"]) {
             if (!plane_json.is_object()) continue;
            try {
                if (!plane_json.contains("corners") || !plane_json["corners"].is_array() || plane_json["corners"].size() != 4) {
                    std::cerr << "Warning: Skipping invalid plane definition." << std::endl;
                    continue;
                }
                
                std::array<std::array<float, 3>, 4> corners_array;
                for (size_t i = 0; i < 4; ++i) {
                    corners_array[i] = plane_json["corners"][i].get<std::array<float, 3>>();
                }
                 
                // --- MODIFIED SECTION ---
                Material mat; // Create a default material
                if (plane_json.contains("material")) {
                    mat = parse_material(plane_json["material"]); // Parse the real one
                }

                // Pass the material to the constructor
                loaded_shapes.push_back(std::make_unique<Plane>(corners_array.data(), mat));
                // --- END MODIFIED SECTION ---

            } catch (json::exception& e) {
                std::cerr << "Warning: Error parsing plane entry: " << e.what() << std::endl;
            }
        }
    }

    if (loaded_shapes.empty()) {
        std::cerr << "Warning: No valid shapes were loaded from " << filename << "." << std::endl;
    }

    return loaded_shapes;
}