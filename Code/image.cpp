// image.cpp
#include "image.hpp"
#include <fstream>
#include <iostream>
#include <sstream>
#include <algorithm>

// Constructor: create blank image
Image::Image(int w, int h) : width(w), height(h) {
    pixels.resize(width * height * 3, 0); // 3 values per pixel (RGB)
}

// Constructor: read from file
Image::Image(const std::string& filename) {
    read(filename);
}

int Image::clamp(int value, int min, int max) {
    return std::max(min, std::min(value, max));
}

// Helper: convert 2D coordinates to 1D array index
int Image::getIndex(int x, int y) const {
    return (y * width + x) * 3;
}

// Set pixel color (with bounds checking)
void Image::setPixel(int x, int y, int r, int g, int b) {
    if (x < 0 || x >= width || y < 0 || y >= height) {
        return; // Out of bounds
    }
    
    int index = getIndex(x, y);
    pixels[index] = clamp(r, 0, 255);
    pixels[index + 1] = clamp(g, 0, 255);
    pixels[index + 2] = clamp(b, 0, 255);
}

// Get pixel color
void Image::getPixel(int x, int y, int& r, int& g, int& b) const {
    if (x < 0 || x >= width || y < 0 || y >= height) {
        r = g = b = 0;
        return;
    }
    
    int index = getIndex(x, y);
    r = pixels[index];
    g = pixels[index + 1];
    b = pixels[index + 2];
}

// Write PPM file (P3 format - ASCII)
void Image::write(const std::string& filename) const {
    std::ofstream file(filename);
    
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << filename << " for writing\n";
        return;
    }
    
    // PPM header
    file << "P3\n";
    file << width << " " << height << "\n";
    file << "255\n";
    
    // Pixel data
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            int index = getIndex(x, y);
            file << (int)pixels[index] << " "
                 << (int)pixels[index + 1] << " "
                 << (int)pixels[index + 2];
            
            if (x < width - 1) {
                file << "  ";
            }
        }
        file << "\n";
    }
    
    file.close();
    std::cout << "Image written to " << filename << "\n";
}

// Read PPM file (P3 format)
void Image::read(const std::string& filename) {
    std::ifstream file(filename);
    
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << filename << " for reading\n";
        return;
    }
    
    std::string line;
    std::string magic;
    
    // Read magic number
    file >> magic;
    if (magic != "P3") {
        std::cerr << "Error: Only P3 PPM format is supported\n";
        return;
    }
    
    // Skip comments
    file >> std::ws;
    while (file.peek() == '#') {
        std::getline(file, line);
        file >> std::ws;
    }
    
    // Read dimensions
    file >> width >> height;
    
    // Read max color value
    int maxColor;
    file >> maxColor;
    
    if (maxColor != 255) {
        std::cerr << "Warning: Max color value is " << maxColor << ", expected 255\n";
    }
    
    // Read pixel data
    pixels.resize(width * height * 3);
    for (int i = 0; i < width * height * 3; i++) {
        int value;
        file >> value;
        pixels[i] = clamp(value, 0, 255);
    }
    
    file.close();
    std::cout << "Image read from " << filename << " (" << width << "x" << height << ")\n";
    loaded_successfully = true; 

}