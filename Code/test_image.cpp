// test_image.cpp
#include "image.hpp"
#include <cassert>
#include <iostream>
#include <fstream>
#include <cstdio>

// Helper function to check if two images are equal
bool imagesEqual(const Image& img1, const Image& img2) {
    if (img1.getWidth() != img2.getWidth() || img1.getHeight() != img2.getHeight()) {
        return false;
    }
    
    for (int y = 0; y < img1.getHeight(); y++) {
        for (int x = 0; x < img1.getWidth(); x++) {
            int r1, g1, b1, r2, g2, b2;
            img1.getPixel(x, y, r1, g1, b1);
            img2.getPixel(x, y, r2, g2, b2);
            if (r1 != r2 || g1 != g2 || b1 != b2) {
                return false;
            }
        }
    }
    return true;
}

// Test 1: Constructor creates blank image
void test_blank_image_creation() {
    std::cout << "Test 1: Blank image creation... ";
    
    Image img(10, 20);
    assert(img.getWidth() == 10);
    assert(img.getHeight() == 20);
    
    // Check that all pixels are initialized to black (0, 0, 0)
    int r, g, b;
    img.getPixel(0, 0, r, g, b);
    assert(r == 0 && g == 0 && b == 0);
    
    img.getPixel(5, 10, r, g, b);
    assert(r == 0 && g == 0 && b == 0);
    
    std::cout << "PASSED\n";
}

// Test 2: Set and get pixel
void test_set_get_pixel() {
    std::cout << "Test 2: Set and get pixel... ";
    
    Image img(10, 10);
    
    img.setPixel(5, 5, 255, 128, 64);
    
    int r, g, b;
    img.getPixel(5, 5, r, g, b);
    assert(r == 255 && g == 128 && b == 64);
    
    std::cout << "PASSED\n";
}

// Test 3: Bounds checking for setPixel
void test_set_pixel_bounds() {
    std::cout << "Test 3: SetPixel bounds checking... ";
    
    Image img(10, 10);
    
    // These should not crash or cause errors
    img.setPixel(-1, 5, 255, 0, 0);    // Negative x
    img.setPixel(5, -1, 255, 0, 0);    // Negative y
    img.setPixel(10, 5, 255, 0, 0);    // x == width
    img.setPixel(5, 10, 255, 0, 0);    // y == height
    img.setPixel(100, 100, 255, 0, 0); // Way out of bounds
    
    std::cout << "PASSED\n";
}

// Test 4: Bounds checking for getPixel
void test_get_pixel_bounds() {
    std::cout << "Test 4: GetPixel bounds checking... ";
    
    Image img(10, 10);
    int r, g, b;
    
    // Out of bounds should return (0, 0, 0)
    img.getPixel(-1, 5, r, g, b);
    assert(r == 0 && g == 0 && b == 0);
    
    img.getPixel(10, 5, r, g, b);
    assert(r == 0 && g == 0 && b == 0);
    
    img.getPixel(5, -1, r, g, b);
    assert(r == 0 && g == 0 && b == 0);
    
    img.getPixel(5, 10, r, g, b);
    assert(r == 0 && g == 0 && b == 0);
    
    std::cout << "PASSED\n";
}

// Test 5: Color clamping
void test_color_clamping() {
    std::cout << "Test 5: Color clamping... ";
    
    Image img(10, 10);
    
    // Values above 255 should be clamped to 255
    img.setPixel(0, 0, 300, 400, 500);
    int r, g, b;
    img.getPixel(0, 0, r, g, b);
    assert(r == 255 && g == 255 && b == 255);
    
    // Values below 0 should be clamped to 0
    img.setPixel(1, 1, -50, -100, -200);
    img.getPixel(1, 1, r, g, b);
    assert(r == 0 && g == 0 && b == 0);
    
    std::cout << "PASSED\n";
}

// Test 6: Write and read PPM file
void test_write_read_ppm() {
    std::cout << "Test 6: Write and read PPM file... ";
    
    const std::string filename = "test_output.ppm";
    
    // Create a test image
    Image img1(3, 2);
    img1.setPixel(0, 0, 255, 0, 0);     // Red
    img1.setPixel(1, 0, 0, 255, 0);     // Green
    img1.setPixel(2, 0, 0, 0, 255);     // Blue
    img1.setPixel(0, 1, 255, 255, 0);   // Yellow
    img1.setPixel(1, 1, 255, 0, 255);   // Magenta
    img1.setPixel(2, 1, 0, 255, 255);   // Cyan
    
    // Write to file
    img1.write(filename);
    
    // Read back from file
    Image img2(filename);
    
    // Compare
    assert(img2.getWidth() == 3);
    assert(img2.getHeight() == 2);
    assert(imagesEqual(img1, img2));
    
    // Clean up
    std::remove(filename.c_str());
    
    std::cout << "PASSED\n";
}

// Test 7: Read from constructor
void test_constructor_with_file() {
    std::cout << "Test 7: Constructor with filename... ";
    
    const std::string filename = "test_constructor.ppm";
    
    // Create and save a test image
    Image img1(5, 5);
    img1.setPixel(2, 2, 100, 150, 200);
    img1.write(filename);
    
    // Load using constructor
    Image img2(filename);
    
    assert(img2.getWidth() == 5);
    assert(img2.getHeight() == 5);
    
    int r, g, b;
    img2.getPixel(2, 2, r, g, b);
    assert(r == 100 && g == 150 && b == 200);
    
    // Clean up
    std::remove(filename.c_str());
    
    std::cout << "PASSED\n";
}

// Test 8: Large image
void test_large_image() {
    std::cout << "Test 8: Large image handling... ";
    
    Image img(1000, 1000);
    
    // Set some pixels
    img.setPixel(0, 0, 255, 0, 0);
    img.setPixel(999, 999, 0, 0, 255);
    img.setPixel(500, 500, 0, 255, 0);
    
    // Verify
    int r, g, b;
    img.getPixel(0, 0, r, g, b);
    assert(r == 255 && g == 0 && b == 0);
    
    img.getPixel(999, 999, r, g, b);
    assert(r == 0 && g == 0 && b == 255);
    
    img.getPixel(500, 500, r, g, b);
    assert(r == 0 && g == 255 && b == 0);
    
    std::cout << "PASSED\n";
}

// Test 9: Edge pixels
void test_edge_pixels() {
    std::cout << "Test 9: Edge pixel handling... ";
    
    Image img(10, 10);
    
    // Set corner pixels
    img.setPixel(0, 0, 255, 0, 0);           // Top-left
    img.setPixel(9, 0, 0, 255, 0);           // Top-right
    img.setPixel(0, 9, 0, 0, 255);           // Bottom-left
    img.setPixel(9, 9, 255, 255, 0);         // Bottom-right
    
    // Verify corners
    int r, g, b;
    img.getPixel(0, 0, r, g, b);
    assert(r == 255 && g == 0 && b == 0);
    
    img.getPixel(9, 0, r, g, b);
    assert(r == 0 && g == 255 && b == 0);
    
    img.getPixel(0, 9, r, g, b);
    assert(r == 0 && g == 0 && b == 255);
    
    img.getPixel(9, 9, r, g, b);
    assert(r == 255 && g == 255 && b == 0);
    
    std::cout << "PASSED\n";
}

// Test 10: PPM file with comments
void test_ppm_with_comments() {
    std::cout << "Test 10: PPM file with comments... ";
    
    const std::string filename = "test_comments.ppm";
    
    // Create a PPM file with comments manually
    std::ofstream file(filename);
    file << "P3\n";
    file << "# This is a comment\n";
    file << "# Another comment\n";
    file << "2 2\n";
    file << "255\n";
    file << "255 0 0  0 255 0\n";
    file << "0 0 255  255 255 255\n";
    file.close();
    
    // Read the file
    Image img(filename);
    
    assert(img.getWidth() == 2);
    assert(img.getHeight() == 2);
    
    int r, g, b;
    img.getPixel(0, 0, r, g, b);
    assert(r == 255 && g == 0 && b == 0);
    
    img.getPixel(1, 1, r, g, b);
    assert(r == 255 && g == 255 && b == 255);
    
    // Clean up
    std::remove(filename.c_str());
    
    std::cout << "PASSED\n";
}

// Test 11: Empty/1x1 image
void test_minimal_image() {
    std::cout << "Test 11: Minimal (1x1) image... ";
    
    Image img(1, 1);
    img.setPixel(0, 0, 123, 45, 67);
    
    int r, g, b;
    img.getPixel(0, 0, r, g, b);
    assert(r == 123 && g == 45 && b == 67);
    
    std::cout << "PASSED\n";
}

// Main test runner
int main() {
    std::cout << "Running Image module tests...\n\n";
    
    try {
        test_blank_image_creation();
        test_set_get_pixel();
        test_set_pixel_bounds();
        test_get_pixel_bounds();
        test_color_clamping();
        test_write_read_ppm();
        test_constructor_with_file();
        test_large_image();
        test_edge_pixels();
        test_ppm_with_comments();
        test_minimal_image();
        
        std::cout << "\n✓ All tests passed!\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "\n✗ Test failed with exception: " << e.what() << "\n";
        return 1;
    }
}