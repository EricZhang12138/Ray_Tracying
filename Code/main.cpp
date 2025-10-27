// test_camera.cpp
#include "camera.hpp"  // Your camera header
#include <iostream>
#include <iomanip>
#include "shapes.hpp"

void printArray(const std::array<float,3>& arr, const std::string& name) {
    std::cout << name << ": [" << arr[0] << ", " << arr[1] << ", " << arr[2] << "]" << std::endl;
}

int main() {
    std::string scene_path = "/Users/ericzhang/Documents/Computer_graphics/Coursework/s2286795/ASCII/scene.json";
    
    // Test 1: Constructor and file reading
    std::cout << "=== Test 1: Loading Camera ===" << std::endl;
    Camera camera(scene_path);
    
    // Test 2: Check if values were read correctly
    std::cout << "\n=== Test 2: Camera Parameters ===" << std::endl;
    auto resolution = camera.getResolution();
    std::cout << "Resolution: " << std::get<0>(resolution) << " x " << std::get<1>(resolution) << std::endl;
    
    auto sensor = camera.getSensorDim();
    std::cout << "Sensor: " << std::get<0>(sensor) << " x " << std::get<1>(sensor) << std::endl;
    
    std::cout << "Focal Length: " << camera.getFocalLength() << std::endl;
    
    printArray(camera.getLocation(), "Location");
    printArray(camera.getGazeVec(), "Gaze Vector");
    printArray(camera.getUpVec(), "Up Vector");
    
    // Test 3: Test pixelToRay
    std::cout << "\n=== Test 3: Pixel to Ray Conversion ===" << std::endl;
    
    // Test center pixel
    int centerX = std::get<0>(resolution) / 2;
    int centerY = std::get<1>(resolution) / 2;
    auto [origin, direction] = camera.pixelToRay(std::make_tuple(centerX, centerY));
    
    std::cout << "Center pixel (" << centerX << ", " << centerY << "):" << std::endl;
    printArray(origin, "  Origin");
    printArray(direction, "  Direction");
    
    // Test corner pixels
    auto [origin_tl, dir_tl] = camera.pixelToRay(std::make_tuple(0, 0));
    std::cout << "\nTop-left corner (0, 0):" << std::endl;
    printArray(dir_tl, "  Direction");
    
    auto [origin_br, dir_br] = camera.pixelToRay(std::make_tuple(std::get<0>(resolution)-1, std::get<1>(resolution)-1));
    std::cout << "\nBottom-right corner:" << std::endl;
    printArray(dir_br, "  Direction");
    
    // Test 4: Verify direction vectors are normalized
    std::cout << "\n=== Test 4: Direction Normalization ===" << std::endl;
    float mag = std::sqrt(direction[0]*direction[0] + direction[1]*direction[1] + direction[2]*direction[2]);
    std::cout << "Center ray magnitude (should be ~1.0): " << mag << std::endl;


    // verifying intersections 

    // test cube
    std::array<float, 3> translation = {
                2.0292587280273438,
                0.37665316462516785,
                -0.17845650017261505
            };
    std::array<float, 3> rotation = {
                0.0,
                0.0,
                0.0
            };
    float scale = 2.8;

    Cube cube(translation, rotation, scale);
    Ray ray = {origin,direction};
    Hit hit;
    hit.shape = &cube;

    bool intersect = cube.intersect(hit, ray);
    if (intersect == false){
        std::cout << "Didn't intersect" << std::endl;
    }

    std::cout << "intersection point:" << " ";
    for (float coor: hit.intersection_point){
        std::cout << coor << ",";
    }

    std::cout << "  intersection normal:" << " ";
    for (float normal: hit.normal){
        std::cout << normal << ",";
    }

    std::cout << "  intersection distance:" << hit.t <<" ";
    
    return 0;
}