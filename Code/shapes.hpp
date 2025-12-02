#ifndef SHAPES_HPP
#define SHAPES_HPP

#include <array>
#include <tuple>
#include <string>
#include <cmath>
#include <limits>    // For std::numeric_limits
#include <algorithm> // For std::min, std::max
#include "material.hpp"


class Shapes; // Shapes uses Hit and Hit uses Shape. Forward declaration to avoid compilation error

class Hit{
public:
    std::array<float,3> intersection_point;
    std::array<float,3> normal;
    float t; // distance along the ray
    Shapes *shape;
    float u =0.0f;
    float v =0.0f;
};

struct Ray{
  std::array<float,3> origin;
  std::array<float,3> direction;  
  float time = 0.0f; // NEW: Time stamp for motion blur (0.0 to 1.0)
};

struct AABB {
    std::array<float, 3> min_point;
    std::array<float, 3> max_point;
    // Default constructor for an empty/invalid box
    AABB() {
        float min_val = std::numeric_limits<float>::lowest();
        float max_val = std::numeric_limits<float>::max();
        min_point = {max_val, max_val, max_val};
        max_point = {min_val, min_val, min_val};
    }
    // Constructor with min/max
    AABB(const std::array<float, 3>& min_p, const std::array<float, 3>& max_p)
        : min_point(min_p), max_point(max_p) {}

    // The fast ray-AABB "slab test"
    bool intersect(const Ray& ray) const;

    // merge with another one
    void merge(const AABB& other); 

    // Merges this AABB with a single point
    void merge(const std::array<float, 3>& point);

    // Function to find the longest dimension of the box
    // Returns 0 for X, 1 for Y, 2 for Z
    int get_longest_axis() const;
};

class Shapes{
public:
    virtual ~Shapes() {};
    virtual bool intersect(Hit &hit, const Ray &ray) = 0;
    virtual AABB get_bounding_box() const = 0;
    std::string type = "Shape";
    Material material;
    std::array<float, 3> velocity; 

protected:
    // === Transformation Helpers (Shared by all transformed shapes) ===
    std::array<std::array<float,4>,4> world_to_object;
    std::array<std::array<float,4>,4> object_to_world;

    void buildTransformationMatrices(const std::array<float, 3>& translation,
                                   const std::array<float, 3>& rotation,
                                   const std::array<float, 3>& scale); // Scale is now x,y,z

    std::array<float, 3> transformPoint(const std::array<std::array<float, 4>, 4>& matrix, 
                                      const std::array<float, 3>& point) const;
    std::array<float, 3> transformVector(const std::array<std::array<float, 4>, 4>& matrix, 
                                       const std::array<float, 3>& vector) const;
    std::array<float, 3> transformNormal(const std::array<std::array<float, 4>, 4>& matrix, 
                                       const std::array<float, 3>& normal) const;
    
    // Static helper to multiply 4x4 matrices
    static std::array<std::array<float, 4>, 4> multiplyMatrices(
        const std::array<std::array<float, 4>, 4>& A, 
        const std::array<std::array<float, 4>, 4>& B);
};

class Plane : public Shapes{
public:
    Plane(std::array<float,3> corner[4], const Material& mat);
    bool intersect(Hit &hit, const Ray &ray) override;
    AABB get_bounding_box() const override;
private:
    std::array<float, 3> corners[4]; // four corners, each with 3 coordinates
    bool isPointInQuad(const std::array<float, 3>& point, const std::array<float, 3> corners[4]);
};


// --- Rectangle (Transformed Unit Square) ---
class Rectangle : public Shapes {
public:
    // Constructor takes transform properties
    Rectangle(const std::array<float, 3>& translation, 
              const std::array<float, 3>& rotation, 
              const std::array<float, 3>& scale, 
              const Material& mat);
              
    bool intersect(Hit &hit, const Ray &ray) override;
    AABB get_bounding_box() const override;
};

class Sphere : public Shapes{
// Constructor takes transform properties instead of center/radius
public:
    Sphere(const std::array<float, 3>& translation, 
           const std::array<float, 3>& rotation, 
           const std::array<float, 3>& scale, 
           const Material& mat,
           const std::array<float, 3>& vel);
           
    bool intersect(Hit &hit, const Ray &ray) override;
    AABB get_bounding_box() const override;
};




class Cube : public Shapes{
public:
    Cube(const std::array<float, 3>& translation, 
         const std::array<float, 3>& rotation,  
         const std::array<float, 3>& scale, // Changed float scale to array scale for consistency
         const Material& mat);
         
    bool intersect(Hit &hit, const Ray &ray) override;
    AABB get_bounding_box() const override;
};

#endif