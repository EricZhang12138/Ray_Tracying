# Specify the minimum version of CMake required
cmake_minimum_required(VERSION 4.1)

# Define the project name (this will be the name of your executable)
project(Raytracing_Test LANGUAGES CXX)

# Set the C++ standard required (e.g., C++17)
# nlohmann/json requires at least C++11, C++17 is a safe bet for modern features.
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF) # Optional: disable compiler-specific extensions

# List all your C++ source files (.cpp)
set(SOURCE_FILES
    raytracer.cpp
    camera.cpp
    shapes.cpp
    image.cpp
    acceleration.cpp
    json_loader.cpp
)

# Create the executable
# The first argument is the name of the executable file to create
# The following arguments are the source files needed to build it
add_executable(${PROJECT_NAME} ${SOURCE_FILES})

# --- Optional Section: Add include directories if headers are elsewhere ---
# If json.hpp or your own .hpp files were in a subdirectory (e.g., "include"),
# you would uncomment and modify the following line:
# target_include_directories(${PROJECT_NAME} PRIVATE include)

# --- Optional Section: Link libraries if needed ---
# For this project using only standard libraries and a header-only JSON library,
# you likely don't need to link anything explicitly.
# The math library (-lm on Linux/macOS) is usually linked automatically by modern CMake/compilers.
# If you needed another library (e.g., SDL, OpenGL), you'd add lines like:
# find_package(SDL2 REQUIRED)
# target_link_libraries(${PROJECT_NAME} PRIVATE SDL2::SDL2)

# --- Optional: Add compiler flags for optimization or warnings ---
# Example for GCC/Clang:
# target_compile_options(${PROJECT_NAME} PRIVATE -Wall -Wextra -O2)
# Example for MSVC:
# if(MSVC)
#   target_compile_options(${PROJECT_NAME} PRIVATE /W4 /O2)
# endif()

# Print a message when configuration is done
message(STATUS "CMake configuration complete for ${testing_intersection}")
