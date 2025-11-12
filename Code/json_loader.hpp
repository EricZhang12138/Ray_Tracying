#ifndef JSON_LOADER_HPP
#define JSON_LOADER_HPP

#include <vector>
#include <string>
#include <memory> // For std::unique_ptr
#include "shapes.hpp" // Needs Shapes definition
#include "light.hpp" 
#include "json.hpp"

using json = nlohmann::json;

// Function declaration
Material parse_material(const json& mat_json);
std::vector<std::unique_ptr<Shapes>> load_shapes_from_json(const std::string& filename);
std::vector<Light> load_lights_from_json(const std::string& filename);
#endif // JSON_LOADER_HPP
