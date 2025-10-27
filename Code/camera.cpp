// read the data from the scene.json
#include <fstream> //file stream operation 
#include <iostream>
#include "json.hpp"
#include <string>
#include <tuple>
#include <map>
#include <array>
#include <cmath>
#include "camera.hpp"

using json = nlohmann::json;

    int Camera::readCameraSpec(const std::string &filename){
        std::ifstream input(filename);
        //check if file is open correctly 
        if (!input.is_open()){
            std::cerr << "Error: Could not open file " << filename << std::endl;
            return 0;
        }

        try {
            json camera_data;
            input >> camera_data;  // json is parsed into my camera data

            if (camera_data.contains("cameras") && camera_data.contains("render")){
                focal_length = camera_data["cameras"][0]["focal_length"].get<float>();
                
                
                location = camera_data["cameras"][0]["location"].get<std::array<float, 3>>();
                gaze_vec = camera_data["cameras"][0]["gaze_vector"].get<std::array<float, 3>>();
                up_vec = camera_data["cameras"][0]["up_vector"].get<std::array<float,3>>();
                
                
                std::get<0>(sensor_dim) = camera_data["cameras"][0]["sensor_width"].get<int>();
                std::get<1>(sensor_dim) = camera_data["cameras"][0]["sensor_height"].get<int>();

                std::get<0>(resolution) = camera_data["render"]["resolution_x"].get<int>();
                std::get<1>(resolution) = camera_data["render"]["resolution_y"].get<int>();
            }
            else{
                std::cerr << "Error: JSON file is missing required keys." << std::endl;
            }
        }   catch (const json::parse_error& e) {
        // Handle parsing errors (e.g., malformed JSON)
                std::cerr << "JSON Parse Error: " << e.what() << std::endl;
                return 0;
        }   catch (const std::exception& e) {
        // Handle other exceptions
                std::cerr << "An unexpected error occurred: " << e.what() << std::endl;
                return 0;
        }
        return 1;
    }

    std::array<float, 3> Camera::normalize(const std::array<float, 3>& v) {
    float magnitude = std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
    
    if (magnitude == 0.0f) {
        return {0.0f, 0.0f, 0.0f};
    }
    
    return {v[0] / magnitude, v[1] / magnitude, v[2] / magnitude};
    }

    void Camera::cross_product(const std::array<float,3>& a,const std::array<float,3>& b, std::array<float,3>& result) {
    // We use temporary variables for the results to ensure that 
    // if 'result' is the same array as 'a' or 'b' (though not recommended), 
    // the calculation uses the original values.

    // x component: a[1] * b[2] - a[2] * b[1]
    float result_x = a[1] * b[2] - a[2] * b[1];

    // y component: a[2] * b[0] - a[0] * b[2]
    float result_y = a[2] * b[0] - a[0] * b[2]; 

    // z component: a[0] * b[1] - a[1] * b[0]
    float result_z = a[0] * b[1] - a[1] * b[0]; 

    result[0] = result_x;
    result[1] = result_y;
    result[2] = result_z;
    }



    std::tuple<std::array<float,3>, std::array<float,3>> Camera::pixelToRay(std::tuple<int,int> pixel){
        // image space ----> Normalised Device Coordinates
        float nx, ny;
        nx = (((float)std::get<0>(pixel))/(float)std::get<0>(resolution)) * 2 - 1 ;
        ny = (((float)std::get<1>(pixel))/(float)std::get<1>(resolution)) * 2 - 1 ;

        // NDC ----> Camera Space
        // actual length on the sensor
        float nx_r = nx * std::get<0>(sensor_dim);
        float ny_r = ny * std::get<1>(sensor_dim);
        
        std::array<float,3> direction = {nx_r, ny_r, focal_length};

        //normalised gaze vector      z direction
        std::array<float,3> z_dir = normalize(gaze_vec);
        //x direction 
        std::array<float,3> x_dir;
        cross_product(up_vec, z_dir, x_dir);
        x_dir = normalize(x_dir);
        //y direction. This is the corrected Vup, ensuring perfect orthogonality
        std::array<float,3> y_dir;
        cross_product(z_dir,x_dir,y_dir);
        y_dir = normalize(y_dir);


        std::array<std::array<float, 4>, 4> M_C2W = {{
            {x_dir[0], y_dir[0], z_dir[0], location[0]},
            {x_dir[1], y_dir[1], z_dir[1], location[1]},
            {x_dir[2], y_dir[2], z_dir[2], location[2]},
            {0.0f,     0.0f,     0.0f,     1.0f}
        }};

        std::array<float, 4> dir_camera = {nx_r, ny_r, focal_length, 0.0f};


        // Transform using matrix multiplication
        std::array<float, 3> direction_world;
        for (int i = 0; i < 3; i++) {
            direction_world[i] = M_C2W[i][0] * dir_camera[0] + 
                                M_C2W[i][1] * dir_camera[1] + 
                                M_C2W[i][2] * dir_camera[2] + 
                                M_C2W[i][3] * dir_camera[3];
        }
           
        direction_world = normalize(direction_world);
        
        return std::make_tuple(location, direction_world);
    }




// need to use a initialiser list so that every member variable has a initial value even when readCameraSpec fails
    Camera::Camera(const std::string &filename) 
    : resolution{0, 0}, 
      sensor_dim{0, 0}, 
      focal_length(0.0f), 
      location{0.0f, 0.0f, 0.0f}, 
      gaze_vec{0.0f, 0.0f, 0.0f}, 
      up_vec{0.0f, 0.0f, 0.0f}
    {
        // Now, readCameraSpec is called AFTER members have safe default values.
        if (readCameraSpec(filename) == 0) {
            std::cerr << "Camera configuration failed to load. Using default values." << std::endl;
        }
    }



    // Getter methods for testing
    std::tuple<int,int> Camera::getResolution() const { return resolution; }
    std::tuple<float,float> Camera::getSensorDim() const { return sensor_dim; }
    float Camera::getFocalLength() const { return focal_length; }
    std::array<float,3> Camera::getLocation() const { return location; }
    std::array<float,3> Camera::getGazeVec() const { return gaze_vec; }
    std::array<float,3> Camera::getUpVec() const { return up_vec; }
    




