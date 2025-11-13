#include "shapes.hpp"
#include <cmath>  


namespace {

// --- Vector Math Helpers ---
std::array<float, 3> vecSub(const std::array<float, 3>& a, const std::array<float, 3>& b) {
    return {a[0] - b[0], a[1] - b[1], a[2] - b[2]};
}

std::array<float, 3> vecCross(const std::array<float, 3>& a, const std::array<float, 3>& b) {
    return {
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0]
    };
}

float vecDot(const std::array<float, 3>& a, const std::array<float, 3>& b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

/**
 * @brief Checks if a point P (on the same plane) is inside a triangle (A, B, C).
 * This uses the 3D "same-side" test, which is robust to winding order.
 */
// In shapes.cpp

// In shapes.cpp

// In shapes.cpp

bool isPointInTriangle(const std::array<float, 3>& P,
                       const std::array<float, 3>& A,
                       const std::array<float, 3>& B,
                       const std::array<float, 3>& C,
                       const std::array<float, 3>& quad_normal) // Use the passed-in normal
{
    // Use the quad_normal directly for the side checks.
    // The check will now be: (Cross(Edge, Point-to-Edge) . quad_normal) >= 0

    // Check side of edge AB
    std::array<float, 3> edge1 = vecSub(B, A);
    std::array<float, 3> vp1 = vecSub(P, A);
    std::array<float, 3> C1 = vecCross(edge1, vp1);
    if (vecDot(C1, quad_normal) < -1e-6f) return false; // Use quad_normal

    // Check side of edge BC
    std::array<float, 3> edge2 = vecSub(C, B);
    std::array<float, 3> vp2 = vecSub(P, B);
    std::array<float, 3> C2 = vecCross(edge2, vp2);
    if (vecDot(C2, quad_normal) < -1e-6f) return false; // Use quad_normal

    // Check side of edge CA
    std::array<float, 3> edge3 = vecSub(A, C);
    std::array<float, 3> vp3 = vecSub(P, C);
    std::array<float, 3> C3 = vecCross(edge3, vp3);
    if (vecDot(C3, quad_normal) < -1e-6f) return false; // Use quad_normal

    return true;
}

}


// get_longest_axis() for AABB
int AABB::get_longest_axis() const {
    // Calculate the dimensions of the bounding box
    float x_dim = max_point[0] - min_point[0];
    float y_dim = max_point[1] - min_point[1];
    float z_dim = max_point[2] - min_point[2];

    // Find the largest dimension
    if (x_dim > y_dim && x_dim > z_dim) {
        return 0; // X-axis is longest
    } else if (y_dim > z_dim) {
        return 1; // Y-axis is longest
    } else {
        return 2; // Z-axis is longest
    }
}

bool AABB::intersect(const Ray& ray) const {
    float t_near = std::numeric_limits<float>::lowest();
    float t_far = std::numeric_limits<float>::max();

    for (int i = 0; i < 3; i++) {
        if (fabs(ray.direction[i]) < 1e-6) {
            // Ray is parallel to this pair of planes
            if (ray.origin[i] < min_point[i] || ray.origin[i] > max_point[i]) {
                return false; // Ray origin is outside the slab
            }
        } else {
            // Calculate intersection distances with the two planes
            float t1 = (min_point[i] - ray.origin[i]) / ray.direction[i];
            float t2 = (max_point[i] - ray.origin[i]) / ray.direction[i];

            if (t1 > t2) {
                std::swap(t1, t2); // Ensure t1 is the near intersection
            }

            t_near = std::max(t_near, t1);
            t_far = std::min(t_far, t2);

            if (t_near > t_far || t_far < 0) {
                return false;
            }
        }
    }
    return true; // Ray intersects the box
}

void AABB::merge(const AABB& other) {
    min_point[0] = std::min(min_point[0], other.min_point[0]);
    min_point[1] = std::min(min_point[1], other.min_point[1]);
    min_point[2] = std::min(min_point[2], other.min_point[2]);
    
    max_point[0] = std::max(max_point[0], other.max_point[0]);
    max_point[1] = std::max(max_point[1], other.max_point[1]);
    max_point[2] = std::max(max_point[2], other.max_point[2]);
}

void AABB::merge(const std::array<float, 3>& point) {
    min_point[0] = std::min(min_point[0], point[0]);
    min_point[1] = std::min(min_point[1], point[1]);
    min_point[2] = std::min(min_point[2], point[2]);
    
    max_point[0] = std::max(max_point[0], point[0]);
    max_point[1] = std::max(max_point[1], point[1]);
    max_point[2] = std::max(max_point[2], point[2]);
}



Plane::Plane(std::array<float,3> corner[4], const Material& mat){
        corners[0] = corner[0];
        corners[1] = corner[1];
        corners[2] = corner[2];
        corners[3] = corner[3];
        type = "Plane";
        this -> material = mat;
    };


bool Plane::intersect(Hit &hit, const Ray &ray) {
    // For a plane defined by four corners, we need to:
    // 1. Find the plane equation (point + normal)
    // 2. Check ray-plane intersection
    // 3. Verify the intersection point is within the quad bounds
    
    // Calculate plane normal from three corners (assuming corners[0], corners[1], corners[2])
    std::array<float, 3> edge1 = {
        corners[1][0] - corners[0][0],
        corners[1][1] - corners[0][1],
        corners[1][2] - corners[0][2]
    };
    
    std::array<float, 3> edge2 = {
        corners[2][0] - corners[0][0],
        corners[2][1] - corners[0][1],
        corners[2][2] - corners[0][2]
    };
    
    // Cross product to get normal
    std::array<float, 3> normal = {
        edge1[1] * edge2[2] - edge1[2] * edge2[1],
        edge1[2] * edge2[0] - edge1[0] * edge2[2],
        edge1[0] * edge2[1] - edge1[1] * edge2[0]
    };
    
    // Normalize the normal vector
    float normal_length = sqrt(normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2]);
    normal[0] /= normal_length;
    normal[1] /= normal_length;
    normal[2] /= normal_length;
    
    // Check if ray is parallel to plane
    float denom = normal[0] * ray.direction[0] + 
                  normal[1] * ray.direction[1] + 
                  normal[2] * ray.direction[2];
    
    if (fabs(denom) < 1e-6) {
        return false; // Ray is parallel to plane
    }
    
    // Calculate t value for intersection
    std::array<float, 3> p0_to_origin = {
        corners[0][0] - ray.origin[0],
        corners[0][1] - ray.origin[1],
        corners[0][2] - ray.origin[2]
    };
    
    float t = (p0_to_origin[0] * normal[0] + 
               p0_to_origin[1] * normal[1] + 
               p0_to_origin[2] * normal[2]) / denom;
    
    // Check if intersection is behind the ray origin
    if (t < 0) {
        return false;
    }
    
    // Calculate intersection point
    std::array<float, 3> intersection = {
        ray.origin[0] + t * ray.direction[0],
        ray.origin[1] + t * ray.direction[1],
        ray.origin[2] + t * ray.direction[2]
    };
    
    // Check if intersection point is inside the quad
    // Using the point-in-quad test
    if (!isPointInQuad(intersection, corners)) {
        return false;
    }

    // Calculate u, v coordinates for the plane
    // We project the hit point onto the plane's edge vectors
    std::array<float, 3> vec_u = vecSub(corners[1], corners[0]);
    std::array<float, 3> vec_v = vecSub(corners[3], corners[0]);
    std::array<float, 3> hit_vec = vecSub(intersection, corners[0]);

    // Calculate projection lengths
    float u = vecDot(hit_vec, vec_u) / vecDot(vec_u, vec_u);
    float v = vecDot(hit_vec, vec_v) / vecDot(vec_v, vec_v);

    // Clamp values to [0, 1] just in case of precision issues
    u = std::max(0.0f, std::min(1.0f, u));
    v = std::max(0.0f, std::min(1.0f, v));
    
    // Fill in the hit information
    hit.u = u;
    hit.v = v;
    hit.intersection_point = intersection;
    hit.normal = normal;
    hit.t = t;
    hit.shape = this;
    
    return true;
}

bool Plane::isPointInQuad(const std::array<float, 3>& point, const std::array<float, 3> corners[4]) {
    
    // 1. Recalculate the plane normal (as requested)
    std::array<float, 3> edge1 = vecSub(corners[1], corners[0]);
    std::array<float, 3> edge2 = vecSub(corners[2], corners[0]);
    std::array<float, 3> normal = vecCross(edge1, edge2);
    
    float normal_length = sqrt(vecDot(normal, normal));
    // Avoid division by zero if the plane is degenerate
    if (normal_length < 1e-6) return false; 
    
    normal[0] /= normal_length;
    normal[1] /= normal_length;
    normal[2] /= normal_length;

    // 2. Test if the point is in the first triangle (corners 0, 1, 2)
    if (isPointInTriangle(point, corners[1], corners[3], corners[2], normal)) {
        return true;
    }
    
    // 3. Test if the point is in the second triangle (corners 0, 2, 3)
    if (isPointInTriangle(point, corners[0], corners[1], corners[2], normal)) {
        return true;
    }
    
    // The point is not in either triangle
    return false;
}

AABB Plane::get_bounding_box() const {
    // Find the min and max coordinates of the 4 corners
    AABB box; // Starts as an invalid box
    
    // To avoid issues with infinitely-thin planes, we can add a tiny padding
    const float padding = 1e-4f; 

    for (int i = 0; i < 4; ++i) {
        box.merge(corners[i]);
    }

    // Apply padding
    box.min_point[0] -= padding;
    box.min_point[1] -= padding;
    box.min_point[2] -= padding;
    box.max_point[0] += padding;
    box.max_point[1] += padding;
    box.max_point[2] += padding;

    return box;
}

 Sphere::Sphere(const std::array<float, 3>& c, float r, const Material& mat) : center(c), radius(r) {
    type = "Sphere";
    this -> material = mat;
 }


bool Sphere::intersect(Hit &hit, const Ray &ray) {
    // Sphere equation: |P - C|² = r²
    // Ray equation: P(t) = O + tD
    // Substitute ray into sphere equation and solve for t
    
    // Vector from ray origin to sphere center
    std::array<float, 3> oc = {
        ray.origin[0] - center[0],
        ray.origin[1] - center[1],
        ray.origin[2] - center[2]
    };
    
    // Quadratic equation coefficients: at² + bt + c = 0
    // a = D · D (if direction is normalized, a = 1)
    float a = ray.direction[0] * ray.direction[0] + 
              ray.direction[1] * ray.direction[1] + 
              ray.direction[2] * ray.direction[2];
    
    // b = 2(D · (O - C))
    float b = 2.0f * (oc[0] * ray.direction[0] + 
                      oc[1] * ray.direction[1] + 
                      oc[2] * ray.direction[2]);
    
    // c = |O - C|² - r²
    float c = (oc[0] * oc[0] + oc[1] * oc[1] + oc[2] * oc[2]) - 
              (radius * radius);
    
    // Calculate discriminant
    float discriminant = b * b - 4 * a * c;
    
    // No intersection if discriminant is negative
    if (discriminant < 0) {
        return false;
    }
    
    // Calculate the two possible t values
    float sqrt_discriminant = sqrt(discriminant);
    float t1 = (-b - sqrt_discriminant) / (2.0f * a);
    float t2 = (-b + sqrt_discriminant) / (2.0f * a);
    
    // We want the closest positive t value
    float t;
    if (t1 > 0) {
        t = t1;  // Use the closer intersection
    } else if (t2 > 0) {
        t = t2;  // Ray origin is inside sphere, use far intersection
    } else {
        return false;  // Both intersections are behind the ray
    }
    
    // Calculate intersection point
    std::array<float, 3> intersection = {
        ray.origin[0] + t * ray.direction[0],
        ray.origin[1] + t * ray.direction[1],
        ray.origin[2] + t * ray.direction[2]
    };
    
    // Calculate normal (vector from center to intersection point)
    std::array<float, 3> normal = {
        intersection[0] - center[0],
        intersection[1] - center[1],
        intersection[2] - center[2]
    };
    
    // Normalize the normal vector
    float normal_length = sqrt(normal[0] * normal[0] + 
                              normal[1] * normal[1] + 
                              normal[2] * normal[2]);
    normal[0] /= normal_length;
    normal[1] /= normal_length;
    normal[2] /= normal_length;

    // This maps the normal vector (x,y,z) to (u,v)
    // atan2 handles all quadrants for longitude (u)
    // asin is simple for latitude (v)
    const float PI = 3.1415926535f;
    float phi = atan2(normal[2], normal[0]); // Longitude, from -PI to +PI
    float theta = asin(normal[1]);          // Latitude, from -PI/2 to +PI/2
    
    float u = 1.0f - (phi + PI) / (2.0f * PI); // Map to [0, 1]
    float v = (theta + PI / 2.0f) / PI;      // Map to [0, 1]

    hit.u = u;
    hit.v = v;
    
    // Fill in the hit information
    hit.intersection_point = intersection;
    hit.normal = normal;
    hit.t = t;
    hit.shape = this;
    
    return true;
}


AABB Sphere::get_bounding_box() const {
    std::array<float, 3> min_p = {
        center[0] - radius,
        center[1] - radius,
        center[2] - radius
    };
    std::array<float, 3> max_p = {
        center[0] + radius,
        center[1] + radius,
        center[2] + radius
    };
    return AABB(min_p, max_p);
}


Cube::Cube(const std::array<float, 3>& translation,
           const std::array<float, 3>& rotation,  
           float scale, const Material& mat) {
    buildTransformationMatrices(translation, rotation, scale);
    type = "Cube";
    this -> material = mat;
}


void Cube::buildTransformationMatrices(const std::array<float, 3>& translation,
                                       const std::array<float, 3>& rotation,
                                       float scale) {
    // First build the object-to-world matrix: M = T * R * S
    // Where S is scale, R is rotation, T is translation
    
    // Initialize as identity matrix
    object_to_world = {{
        {1, 0, 0, 0},
        {0, 1, 0, 0},
        {0, 0, 1, 0},
        {0, 0, 0, 1}
    }};
    
    float model_scale = scale * 2.0f; 

    // Step 1: Build scale matrix
    std::array<std::array<float, 4>, 4> scale_matrix = {{
        {model_scale, 0, 0, 0},
        {0, model_scale, 0, 0},
        {0, 0, model_scale, 0},
        {0, 0, 0, 1}
    }};
    
    // Step 2: Build rotation matrix (using Euler angles: X-Y-Z order)
    float cx = cos(rotation[0]), sx = sin(rotation[0]);
    float cy = cos(rotation[1]), sy = sin(rotation[1]);
    float cz = cos(rotation[2]), sz = sin(rotation[2]);
    
    // Combined rotation matrix R = Rz * Ry * Rx
    std::array<std::array<float, 4>, 4> rotation_matrix = {{
        {cy*cz, sx*sy*cz - cx*sz, cx*sy*cz + sx*sz, 0},
        {cy*sz, sx*sy*sz + cx*cz, cx*sy*sz - sx*cz, 0},
        {-sy, sx*cy, cx*cy, 0},
        {0, 0, 0, 1}
    }};
    
    // Step 3: Build translation matrix
    std::array<std::array<float, 4>, 4> translation_matrix = {{
        {1, 0, 0, translation[0]},
        {0, 1, 0, translation[1]},
        {0, 0, 1, translation[2]},
        {0, 0, 0, 1}
    }};
    
    // Step 4: Combine transformations: T * R * S
    // First compute R * S
    std::array<std::array<float, 4>, 4> RS = multiplyMatrices(rotation_matrix, scale_matrix);
    // Then compute T * (R * S)
    object_to_world = multiplyMatrices(translation_matrix, RS);
    
    // Step 5: Compute the inverse transformation (world-to-object)
    // For M = T * R * S, the inverse is M^-1 = S^-1 * R^T * T^-1
    
    // Inverse scale
    float inv_scale = 1.0f / model_scale;
    std::array<std::array<float, 4>, 4> inv_scale_matrix = {{
        {inv_scale, 0, 0, 0},
        {0, inv_scale, 0, 0},
        {0, 0, inv_scale, 0},
        {0, 0, 0, 1}
    }};
    
    // Transpose of rotation matrix (rotation matrices are orthogonal, so transpose = inverse)
    std::array<std::array<float, 4>, 4> inv_rotation_matrix = {{
        {rotation_matrix[0][0], rotation_matrix[1][0], rotation_matrix[2][0], 0},
        {rotation_matrix[0][1], rotation_matrix[1][1], rotation_matrix[2][1], 0},
        {rotation_matrix[0][2], rotation_matrix[1][2], rotation_matrix[2][2], 0},
        {0, 0, 0, 1}
    }};
    
    // Inverse translation
    std::array<std::array<float, 4>, 4> inv_translation_matrix = {{
        {1, 0, 0, -translation[0]},
        {0, 1, 0, -translation[1]},
        {0, 0, 1, -translation[2]},
        {0, 0, 0, 1}
    }};
    
    // Combine: S^-1 * R^T * T^-1
    std::array<std::array<float, 4>, 4> SR = multiplyMatrices(inv_scale_matrix, inv_rotation_matrix);
    world_to_object = multiplyMatrices(SR, inv_translation_matrix);
}

// Helper function to multiply two 4x4 matrices
std::array<std::array<float, 4>, 4> Cube::multiplyMatrices(
    const std::array<std::array<float, 4>, 4>& A,
    const std::array<std::array<float, 4>, 4>& B) {
    
    std::array<std::array<float, 4>, 4> result = {{
        {0, 0, 0, 0},
        {0, 0, 0, 0},
        {0, 0, 0, 0},
        {0, 0, 0, 0}
    }};
    
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            for (int k = 0; k < 4; k++) {
                result[i][j] += A[i][k] * B[k][j];
            }
        }
    }
    
    return result;
}


bool Cube::intersect(Hit &hit, const Ray &ray) {
    // Transform ray from world space to object space
    // In object space, we test against a unit cube centered at origin
    // with bounds from -0.5 to +0.5 in each dimension
    
    Ray object_ray;
    object_ray.origin = transformPoint(world_to_object, ray.origin);
    object_ray.direction = transformVector(world_to_object, ray.direction);
    
    // Define the cube bounds in object space
    float t_min = -0.5f;
    float t_max = 0.5f;
    
    // Variables to track the intersection distances
    float t_near = -INFINITY;
    float t_far = INFINITY;
    
    // Track which axis we hit (for normal calculation)
    int hit_axis = -1;
    int hit_sign = 0;  // -1 for negative face, +1 for positive face
    
    // Test intersection with each pair of parallel planes (x, y, z)
    for (int i = 0; i < 3; i++) {
        if (fabs(object_ray.direction[i]) < 1e-6) {
            // Ray is parallel to this pair of planes
            if (object_ray.origin[i] < t_min || object_ray.origin[i] > t_max) {
                return false;  // Ray origin is outside the slab, no intersection
            }
        } else {
            // Calculate intersection distances with the two planes
            float t1 = (t_min - object_ray.origin[i]) / object_ray.direction[i];
            float t2 = (t_max - object_ray.origin[i]) / object_ray.direction[i];
            
            // Make sure t1 is the near intersection and t2 is the far
            if (t1 > t2) {
                float temp = t1;
                t1 = t2;
                t2 = temp;
            }
            
            // Update the overall near and far intersections
            if (t1 > t_near) {
                t_near = t1;
                hit_axis = i;
                hit_sign = (object_ray.direction[i] > 0) ? -1 : 1;
            }
            if (t2 < t_far) {
                t_far = t2;
            }
            
            // Check if the intersections are invalid
            if (t_near > t_far || t_far < 0) {
                return false;
            }
        }
    }
    
    // Use t_near if it's positive, otherwise use t_far (ray starts inside cube)
    float t_object = (t_near > 0) ? t_near : t_far;
    
    if (t_object < 0) {
        return false;  // Intersection is behind the ray
    }
    
    // Calculate intersection point in object space
    std::array<float, 3> object_intersection = {
        object_ray.origin[0] + t_object * object_ray.direction[0],
        object_ray.origin[1] + t_object * object_ray.direction[1],
        object_ray.origin[2] + t_object * object_ray.direction[2]
    };
    
    // Calculate normal in object space (axis-aligned)
    std::array<float, 3> object_normal = {0.0f, 0.0f, 0.0f};
    if (hit_axis >= 0) {
        object_normal[hit_axis] = hit_sign;
    }

    // Calculate u,v in object space based on the hit axis
    float u, v;
    // Map coords from [-0.5, 0.5] to [0, 1]
    float uc = object_intersection[0] + 0.5f;
    float vc = object_intersection[1] + 0.5f;
    float wc = object_intersection[2] + 0.5f;

    switch (hit_axis) {
        case 0: // Hit X-face (+x or -x)
            u = (hit_sign > 0) ? wc : (1.0f - wc); // Use Z coord
            v = vc;                               // Use Y coord
            break;
        case 1: // Hit Y-face (+y or -y)
            u = uc;                               // Use X coord
            v = (hit_sign > 0) ? wc : (1.0f - wc); // Use Z coord
            break;
        case 2: // Hit Z-face (+z or -z)
            u = (hit_sign > 0) ? uc : (1.0f - uc); // Use X coord
            v = vc;                               // Use Y coord
            break;
        default: // Should not happen
            u = 0; v = 0;
            break;
    }

    hit.u = u;
    hit.v = v;
    
    // Transform intersection point back to world space
    hit.intersection_point = transformPoint(object_to_world, object_intersection);
    
    // Transform normal back to world space (use inverse transpose for normals)
    hit.normal = transformNormal(object_to_world, object_normal);
    
    hit.t = t_object;
    hit.shape = this;
    
    return true;
}

// Transform a point by a 4x4 matrix
std::array<float, 3> Cube::transformPoint(const std::array<std::array<float, 4>, 4>& matrix, 
                                          const std::array<float, 3>& point) const{
    std::array<float, 3> result;
    
    // Homogeneous coordinates: treat point as (x, y, z, 1)
    result[0] = matrix[0][0] * point[0] + matrix[0][1] * point[1] + 
                matrix[0][2] * point[2] + matrix[0][3];
    result[1] = matrix[1][0] * point[0] + matrix[1][1] * point[1] + 
                matrix[1][2] * point[2] + matrix[1][3];
    result[2] = matrix[2][0] * point[0] + matrix[2][1] * point[1] + 
                matrix[2][2] * point[2] + matrix[2][3];
    
    // Handle homogeneous division if needed
    float w = matrix[3][0] * point[0] + matrix[3][1] * point[1] + 
              matrix[3][2] * point[2] + matrix[3][3];
    
    if (fabs(w - 1.0f) > 1e-6) {
        result[0] /= w;
        result[1] /= w;
        result[2] /= w;
    }
    
    return result;
}




// Transform a vector by a 4x4 matrix (no translation)
std::array<float, 3> Cube::transformVector(const std::array<std::array<float, 4>, 4>& matrix, 
                                           const std::array<float, 3>& vector) {
    std::array<float, 3> result;
    
    // Homogeneous coordinates: treat vector as (x, y, z, 0)
    result[0] = matrix[0][0] * vector[0] + matrix[0][1] * vector[1] + 
                matrix[0][2] * vector[2];
    result[1] = matrix[1][0] * vector[0] + matrix[1][1] * vector[1] + 
                matrix[1][2] * vector[2];
    result[2] = matrix[2][0] * vector[0] + matrix[2][1] * vector[1] + 
                matrix[2][2] * vector[2];
    
    return result;
}


AABB Cube::get_bounding_box() const {
    // To find the AABB of a transformed cube, we must transform
    // all 8 of its corners into world space and find the
    // min/max of those new points.
    
    // 1. Define the 8 corners of the local unit cube
    std::array<std::array<float, 3>, 8> local_corners = {{
        {-0.5, -0.5, -0.5},
        { 0.5, -0.5, -0.5},
        { 0.5,  0.5, -0.5},
        {-0.5,  0.5, -0.5},
        {-0.5, -0.5,  0.5},
        { 0.5, -0.5,  0.5},
        { 0.5,  0.5,  0.5},
        {-0.5,  0.5,  0.5}
    }};
    
    // 2. Transform corners and find min/max
    AABB box; // Starts as an invalid box
    for (int i = 0; i < 8; ++i) {
        // We need a non-const version of transformPoint, or make transformPoint const.
        // Let's assume transformPoint can be called from a const method.
        // We need to pass 'this' or make transformPoint static.
        // Easiest is to just call it directly since it's a member.
        std::array<float, 3> world_corner = transformPoint(object_to_world, local_corners[i]);
        box.merge(world_corner);
    }
    
    return box;
}

/// Transform a normal by a 4x4 matrix (use inverse transpose)
std::array<float, 3> Cube::transformNormal(const std::array<std::array<float, 4>, 4>& obj_to_world_matrix, 
                                           const std::array<float, 3>& normal) {
    // For normals, we need to use the inverse transpose
    // Since we're transforming FROM object space TO world space,
    // and we have both object_to_world and world_to_object matrices,
    // we use the TRANSPOSE of world_to_object (which is the inverse of object_to_world)
    
    // Transpose of world_to_object matrix (upper-left 3x3 only)
    // This is equivalent to (M^-1)^T where M is object_to_world
    std::array<float, 3> result;
    
    result[0] = world_to_object[0][0] * normal[0] + 
                world_to_object[1][0] * normal[1] +  // Note: [1][0] not [0][1]
                world_to_object[2][0] * normal[2];
    
    result[1] = world_to_object[0][1] * normal[0] + 
                world_to_object[1][1] * normal[1] + 
                world_to_object[2][1] * normal[2];
    
    result[2] = world_to_object[0][2] * normal[0] + 
                world_to_object[1][2] * normal[1] + 
                world_to_object[2][2] * normal[2];
    
    // Normalize the result
    float length = sqrt(result[0] * result[0] + 
                       result[1] * result[1] + 
                       result[2] * result[2]);
    
    if (length > 1e-6) {
        result[0] /= length;
        result[1] /= length;
        result[2] /= length;
    }
    
    return result;
}




