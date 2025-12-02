#pragma once

#include <string>
#include <memory>   // For std::shared_ptr
#include "image.hpp"   // For storing a texture (Task 3.3)
#include <iostream>

/**
 * @brief A simple helper struct for storing and passing colors.
 * We use 0.0 to 1.0 for shading calculations.
 */
struct Color {
    float r = 0.0f, g = 0.0f, b = 0.0f;

    // --- Operator Overloads (makes math *much* easier) ---

    // Add colors (c1 + c2)
    Color operator+(const Color& other) const {
        return {r + other.r, g + other.g, b + other.b};
    }
    
    // Multiply colors (c1 * c2)
    Color operator*(const Color& other) const {
        return {r * other.r, g * other.g, b * other.b};
    }
    
    // Multiply by a scalar (c * s)
    Color operator*(float scalar) const {
        return {r * scalar, g * scalar, b * scalar};
    }
    
    // Divide by a scalar (c / s)
    Color operator/(float scalar) const {
        return {r / scalar, g / scalar, b / scalar};
    }
};

// Also define scalar * color (s * c)
inline Color operator*(float scalar, const Color& c) {
    return c * scalar;
}


/**
 * @brief Holds all shading and interaction properties for a surface.
 */
class Material {
public:
    // === 1. Blinn-Phong Properties (Task 3.1) ===
    // These define the local shading of the object.

    Color diffuse_color = {0.8f, 0.8f, 0.8f};
    Color specular_color = {1.0f, 1.0f, 1.0f};

    float k_ambient = 0.1f;    // Ambient coefficient
    float k_diffuse = 0.9f;    // Diffuse coefficient
    float k_specular = 0.3f;   // Specular coefficient
    float shininess = 20.0f;   // Specular "shininess" exponent


    // --- NEW: Distributed Raytracing Property ---
    // 0.0 = Perfect Mirror (Sharp), 1.0 = Diffuse (Blurry)
    float roughness = 0.0f;

    // === 2. Whitted-Style Properties (Task 3.1) ===
    // These define the recursive ray properties.

    float reflectivity = 0.0f;    // 0 = not reflective, 1 = perfect mirror
    float transparency = 0.0f;    // 0 = opaque, 1 = fully transparent
    float refractive_index = 1.0f;  // 1.0 = air, 1.33 = water, 1.5 = glass

    // === 3. Texture Properties (Task 3.3) ===
    // This will hold a pointer to a loaded image file.
    
    // It's a smart pointer, so it will be 'nullptr' by default.
    // Your json_loader will be responsible for creating and storing the Image here.
    std::shared_ptr<Image> texture = nullptr;

public:
    // --- Constructor ---
    Material() {} // Default constructor is fine.


    /**
     * @brief Checks if this material has a texture loaded.
     */
    bool hasTexture() const {
        return (texture != nullptr && (texture->loaded_successfully));
    }

    /**
     * @brief Gets the diffuse color for shading.
     * If a texture exists, it samples the texture.
     * Otherwise, it returns the base diffuse color.
     * * @param u The U texture coordinate (from 0.0 to 1.0).
     * @param v The V texture coordinate (from 0.0 to 1.0).
     * @return The calculated diffuse color.
     */
    Color getDiffuseColor(float u, float v) const {
        if (!hasTexture()) {
            //No texture, just return the simple diffuse color.
            //std::cout << "No texture file is going to be read" << std::endl;
            if (texture != nullptr && !(texture->loaded_successfully)){
                std::cout << "Loaded texture image is a nullptr" << std::endl;
                exit(1);
            }
            //exit(1);
            return diffuse_color;
        }
        

        // --- Texture Mapping ---
        // We have a texture, so we need to sample it.
        
        // Note: Your Image class must have public 'width' and 'height' members,
        // or getters like getWidth()/getHeight() for this to work.
        // I will assume 'width' and 'height' are public based on image.cpp.
        
        // Convert (u, v) from [0, 1] to pixel (x, y) coordinates
        // We use (1.0 - v) because .ppm files (and many image formats)
        // store the top-left pixel first (v=0 is the top).
        int x = static_cast<int>(u * (texture->getWidth() - 1));
        int y = static_cast<int>((1.0f - v) * (texture->getHeight() - 1));

        int r, g, b;
        texture->getPixel(x, y, r, g, b); // Get the int[0,255] color

        // 1. Normalize texture color to 0.0 - 1.0
        Color tex_color = {r / 255.0f, g / 255.0f, b / 255.0f};

        // 2. Multiply by the material's base diffuse color (The Orange Tint)
        // This acts exactly like the "Multiply" node in Blender.
        return tex_color * diffuse_color;
    }
};