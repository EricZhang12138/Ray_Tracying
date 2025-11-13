#ifndef CAMERA_HPP
#define CAMERA_HPP

// These are the only headers needed for the class *declaration*.
// The implementation file (Camera.cpp) will include others like <fstream> and "json.hpp".
#include <string>
#include <tuple>
#include <array>

class Camera {
public:
    /**
     * @brief Constructs a Camera object by loading its configuration from a JSON file.
     * The constructor is marked 'explicit' to prevent unintended implicit conversions.
     * @param filename The path to the scene's JSON file.
     */
    explicit Camera(const std::string& filename);

    /**
     * @brief Generates a ray in world space for a given pixel coordinate.
     * This is the main function for using the camera in a ray tracer.
     * @param pixel A tuple containing the (x, y) integer coordinates of the pixel.
     * @return A tuple containing the ray's origin (std::array<float, 3>) 
     * and normalized direction (std::array<float, 3>).
     */
    std::tuple<std::array<float, 3>, std::array<float, 3>> pixelToRay(std::tuple<float, float> pixel);

    //functions for testing
    std::tuple<int,int> getResolution()const;
    std::tuple<float,float> getSensorDim()const;
    float getFocalLength()const;
    std::array<float,3> getLocation()const;
    std::array<float,3> getGazeVec()const;
    std::array<float,3> getUpVec()const;

private:
    // --- Member Variables ---
    std::tuple<int, int> resolution;
    std::tuple<float, float> sensor_dim;
    float focal_length;
    std::array<float, 3> location;
    std::array<float, 3> gaze_vec;
    std::array<float, 3> up_vec;

    // --- Private Helper Functions (Implementation Details) ---

    // Reads and parses the camera specifications from the JSON file.
    int readCameraSpec(const std::string& filename);

    // Computes the normalized (unit) vector of a given 3D vector.
    std::array<float, 3> normalize(const std::array<float, 3>& v);

    // Computes the cross product of two 3D vectors and stores it in 'result'.
    void cross_product(const std::array<float, 3>& a, const std::array<float, 3>& b, std::array<float, 3>& result);
};

#endif // CAMERA_HPP








