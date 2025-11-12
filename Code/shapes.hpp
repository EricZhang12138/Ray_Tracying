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
};

struct Ray{
  std::array<float,3> origin;
  std::array<float,3> direction;  
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

class Sphere : public Shapes{
private:
    std::array<float, 3> center;  // Sphere center position
    float radius;                  // Sphere radius
    
public:
    // Constructor
    Sphere(const std::array<float, 3>& c, float r, const Material& mat);
    bool intersect(Hit &hit, const Ray &ray) override;
    AABB get_bounding_box() const override;


};

class Cube : public Shapes{
public:
    // Constructor that takes transformation matrices
    Cube(const std::array<float, 3>& translation, const std::array<float, 3>& rotation,  float scale, const Material& mat);
    bool intersect(Hit &hit, const Ray &ray) override;
    AABB get_bounding_box() const override;

private:
    // transformation matrices
    std::array<std::array<float,4>,4> world_to_object;
    std::array<std::array<float,4>,4> object_to_world;
    std::array<float, 3> transformPoint(const std::array<std::array<float, 4>, 4>& matrix, 
                                       const std::array<float, 3>& point) const;
    std::array<float, 3> transformVector(const std::array<std::array<float, 4>, 4>& matrix, 
                                        const std::array<float, 3>& vector);
    std::array<float, 3> transformNormal(const std::array<std::array<float, 4>, 4>& matrix, 
                                        const std::array<float, 3>& normal);
    
    void buildTransformationMatrices(const std::array<float, 3>& translation,
                                       const std::array<float, 3>& rotation,
                                       float scale);
    std::array<std::array<float, 4>, 4> multiplyMatrices(const std::array<std::array<float, 4>, 4>& A, const std::array<std::array<float, 4>, 4>& B);
};

#endif