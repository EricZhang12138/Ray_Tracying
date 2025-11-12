#pragma once

#include <array>

struct Light {
    std::array<float, 3> location;
    std::array<float, 3> color;
    float intensity;

    // Constructor for easy creation
    Light(const std::array<float, 3>& loc, const std::array<float, 3>& col, float i)
        : location(loc), color(col), intensity(i) {}
};