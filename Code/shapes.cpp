#include "shapes.hpp"
#include <cmath>
#include <iostream>

// ==========================================
// Internal Math Helpers (Anonymous Namespace)
// ==========================================
namespace {
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
    
    // Legacy helper for Plane class
    bool isPointInTriangle(const std::array<float, 3>& P, const std::array<float, 3>& A,
                           const std::array<float, 3>& B, const std::array<float, 3>& C,
                           const std::array<float, 3>& quad_normal) {
        std::array<float, 3> edge1 = vecSub(B, A);
        std::array<float, 3> vp1 = vecSub(P, A);
        if (vecDot(vecCross(edge1, vp1), quad_normal) < -1e-6f) return false;

        std::array<float, 3> edge2 = vecSub(C, B);
        std::array<float, 3> vp2 = vecSub(P, B);
        if (vecDot(vecCross(edge2, vp2), quad_normal) < -1e-6f) return false;

        std::array<float, 3> edge3 = vecSub(A, C);
        std::array<float, 3> vp3 = vecSub(P, C);
        if (vecDot(vecCross(edge3, vp3), quad_normal) < -1e-6f) return false;

        return true;
    }
}

// ==========================================
// AABB Implementation
// ==========================================
int AABB::get_longest_axis() const {
    float x_dim = max_point[0] - min_point[0];
    float y_dim = max_point[1] - min_point[1];
    float z_dim = max_point[2] - min_point[2];
    if (x_dim > y_dim && x_dim > z_dim) return 0;
    else if (y_dim > z_dim) return 1;
    else return 2;
}

bool AABB::intersect(const Ray& ray) const {
    float t_near = std::numeric_limits<float>::lowest();
    float t_far = std::numeric_limits<float>::max();

    for (int i = 0; i < 3; i++) {
        if (fabs(ray.direction[i]) < 1e-6) {
            if (ray.origin[i] < min_point[i] || ray.origin[i] > max_point[i]) return false;
        } else {
            float t1 = (min_point[i] - ray.origin[i]) / ray.direction[i];
            float t2 = (max_point[i] - ray.origin[i]) / ray.direction[i];
            if (t1 > t2) std::swap(t1, t2);
            t_near = std::max(t_near, t1);
            t_far = std::min(t_far, t2);
            if (t_near > t_far || t_far < 0) return false;
        }
    }
    return true;
}

void AABB::merge(const AABB& other) {
    for (int i=0; i<3; ++i) {
        min_point[i] = std::min(min_point[i], other.min_point[i]);
        max_point[i] = std::max(max_point[i], other.max_point[i]);
    }
}

void AABB::merge(const std::array<float, 3>& point) {
    for (int i=0; i<3; ++i) {
        min_point[i] = std::min(min_point[i], point[i]);
        max_point[i] = std::max(max_point[i], point[i]);
    }
}

// ==========================================
// Base Shapes Class (Matrix Logic)
// ==========================================

void Shapes::buildTransformationMatrices(const std::array<float, 3>& t,
                                         const std::array<float, 3>& r,
                                         const std::array<float, 3>& s) {
    // 1. Build Scale Matrix
    std::array<std::array<float, 4>, 4> scale_m = {{
        {s[0], 0, 0, 0}, {0, s[1], 0, 0}, {0, 0, s[2], 0}, {0, 0, 0, 1}
    }};

    // 2. Build Rotation Matrix (Euler X-Y-Z)
    float cx = cos(r[0]), sx = sin(r[0]);
    float cy = cos(r[1]), sy = sin(r[1]);
    float cz = cos(r[2]), sz = sin(r[2]);

    std::array<std::array<float, 4>, 4> rot_m = {{
        {cy*cz, sx*sy*cz - cx*sz, cx*sy*cz + sx*sz, 0},
        {cy*sz, sx*sy*sz + cx*cz, cx*sy*sz - sx*cz, 0},
        {-sy,   sx*cy,            cx*cy,            0},
        {0,     0,                0,                1}
    }};

    // 3. Build Translation Matrix
    std::array<std::array<float, 4>, 4> trans_m = {{
        {1, 0, 0, t[0]}, {0, 1, 0, t[1]}, {0, 0, 1, t[2]}, {0, 0, 0, 1}
    }};

    // 4. Combine: Object -> World = T * R * S
    object_to_world = multiplyMatrices(trans_m, multiplyMatrices(rot_m, scale_m));

    // 5. Inverse (World -> Object) = S^-1 * R^T * T^-1
    std::array<std::array<float, 4>, 4> inv_s = {{
        {1.0f/s[0], 0, 0, 0}, {0, 1.0f/s[1], 0, 0}, {0, 0, 1.0f/s[2], 0}, {0, 0, 0, 1}
    }};
    
    // Transpose of rotation
    std::array<std::array<float, 4>, 4> inv_r = {{
        {rot_m[0][0], rot_m[1][0], rot_m[2][0], 0},
        {rot_m[0][1], rot_m[1][1], rot_m[2][1], 0},
        {rot_m[0][2], rot_m[1][2], rot_m[2][2], 0},
        {0, 0, 0, 1}
    }};
    
    // Inverse translation
    std::array<std::array<float, 4>, 4> inv_t = {{
        {1, 0, 0, -t[0]}, {0, 1, 0, -t[1]}, {0, 0, 1, -t[2]}, {0, 0, 0, 1}
    }};

    world_to_object = multiplyMatrices(multiplyMatrices(inv_s, inv_r), inv_t);
}

std::array<std::array<float, 4>, 4> Shapes::multiplyMatrices(const std::array<std::array<float, 4>, 4>& A, 
                                                           const std::array<std::array<float, 4>, 4>& B) {
    std::array<std::array<float, 4>, 4> R = {{{0}}};
    for(int i=0; i<4; i++)
        for(int j=0; j<4; j++)
            for(int k=0; k<4; k++)
                R[i][j] += A[i][k] * B[k][j];
    return R;
}

std::array<float, 3> Shapes::transformPoint(const std::array<std::array<float, 4>, 4>& m, 
                                          const std::array<float, 3>& p) const {
    std::array<float, 3> res;
    float w = m[3][0]*p[0] + m[3][1]*p[1] + m[3][2]*p[2] + m[3][3];
    for(int i=0; i<3; i++) res[i] = m[i][0]*p[0] + m[i][1]*p[1] + m[i][2]*p[2] + m[i][3];
    if(fabs(w - 1.0f) > 1e-6f && w != 0) { res[0]/=w; res[1]/=w; res[2]/=w; }
    return res;
}

std::array<float, 3> Shapes::transformVector(const std::array<std::array<float, 4>, 4>& m, 
                                           const std::array<float, 3>& v) const {
    std::array<float, 3> res;
    for(int i=0; i<3; i++) res[i] = m[i][0]*v[0] + m[i][1]*v[1] + m[i][2]*v[2];
    return res;
}

std::array<float, 3> Shapes::transformNormal(const std::array<std::array<float, 4>, 4>& m, 
                                           const std::array<float, 3>& n) const {
    // Inverse Transpose for Normals. We use world_to_object transposed if transforming Object->World
    // Here we assume 'm' passed in is object_to_world.
    // However, the cleanest way in this setup is to use the world_to_object matrix directly:
    // Normal_world = Transpose(World_to_Object) * Normal_local
    
    // Note: m passed here is usually OBJECT_TO_WORLD.
    // The implementation in the previous Cube code was slightly mixed. 
    // Correct way: use world_to_object for the math.
    
    std::array<float, 3> res;
    // Using world_to_object TRANSPOSED (columns become rows)
    res[0] = world_to_object[0][0]*n[0] + world_to_object[1][0]*n[1] + world_to_object[2][0]*n[2];
    res[1] = world_to_object[0][1]*n[0] + world_to_object[1][1]*n[1] + world_to_object[2][1]*n[2];
    res[2] = world_to_object[0][2]*n[0] + world_to_object[1][2]*n[1] + world_to_object[2][2]*n[2];
    
    float len = sqrt(res[0]*res[0] + res[1]*res[1] + res[2]*res[2]);
    if(len > 1e-6f) { res[0]/=len; res[1]/=len; res[2]/=len; }
    return res;
}

// ==========================================
// Sphere / Oval Implementation
// ==========================================

Sphere::Sphere(const std::array<float, 3>& t, const std::array<float, 3>& r, const std::array<float, 3>& s, const Material& mat, const std::array<float, 3>& vel) {
    type = "Sphere";
    this->material = mat;
    this->velocity = vel; 
    buildTransformationMatrices(t, r, s);
}

bool Sphere::intersect(Hit &hit, const Ray &ray) {
// --- MOTION BLUR FIX START ---
    // Calculate where the ray assumes the object is at this specific time
    Ray moving_ray = ray;
    
    // Shift the ray origin "backwards" based on velocity and time
    // New Origin = Original Origin - (Velocity * Time)
    moving_ray.origin[0] -= velocity[0] * ray.time;
    moving_ray.origin[1] -= velocity[1] * ray.time;
    moving_ray.origin[2] -= velocity[2] * ray.time;
    // --- MOTION BLUR FIX END ---

    // 1. Ray to Object Space (Use moving_ray instead of ray)
    Ray r_loc;
    r_loc.origin = transformPoint(world_to_object, moving_ray.origin); 
    r_loc.direction = transformVector(world_to_object, moving_ray.direction);
    // Do NOT normalize direction here, or 't' will be wrong for non-uniform scale!

    // 2. Intersect Unit Sphere (Radius 1 at 0,0,0)
    // |O + tD|^2 = 1 => dot(D,D)t^2 + 2dot(O,D)t + dot(O,O) - 1 = 0
    float a = vecDot(r_loc.direction, r_loc.direction);
    float b = 2.0f * vecDot(r_loc.origin, r_loc.direction);
    float c = vecDot(r_loc.origin, r_loc.origin) - 1.0f;

    float disc = b*b - 4*a*c;
    if (disc < 0) return false;

    float sqrt_disc = sqrt(disc);
    float t1 = (-b - sqrt_disc) / (2.0f * a);
    float t2 = (-b + sqrt_disc) / (2.0f * a);

    float t_loc = (t1 > 0.001f) ? t1 : ((t2 > 0.001f) ? t2 : -1.0f);
    if (t_loc < 0) return false;

    // 3. Local Hit Point & Normal
    std::array<float, 3> p_loc = {
        r_loc.origin[0] + t_loc * r_loc.direction[0],
        r_loc.origin[1] + t_loc * r_loc.direction[1],
        r_loc.origin[2] + t_loc * r_loc.direction[2]
    };
    std::array<float, 3> n_loc = p_loc; // For unit sphere, normal is point position

    // 4. Transform to World
    hit.intersection_point = transformPoint(object_to_world, p_loc);
    // The point above is where the sphere was at Time 0.
    // We must move the point to where the sphere IS at Time t.
    hit.intersection_point[0] += velocity[0] * ray.time;
    hit.intersection_point[1] += velocity[1] * ray.time;
    hit.intersection_point[2] += velocity[2] * ray.time;
    hit.normal = transformNormal(object_to_world, n_loc); // Logic handles scale
    
    // Recalculate true World T
    std::array<float, 3> dist_v = vecSub(hit.intersection_point, ray.origin);
    hit.t = sqrt(vecDot(dist_v, dist_v));
    hit.shape = this;

    // UVs
    const float PI = 3.1415926535f;
    hit.u = 0.5f + atan2(n_loc[2], n_loc[0]) / (2.0f * PI);
    hit.v = 0.5f - asin(n_loc[1]) / PI;
    
    return true;
}

AABB Sphere::get_bounding_box() const {
    // 1. Get the standard bounding box at Time = 0
    AABB box;
    std::array<std::array<float, 3>, 8> corners = {{
        {-1,-1,-1}, {1,-1,-1}, {1,1,-1}, {-1,1,-1},
        {-1,-1,1},  {1,-1,1},  {1,1,1},  {-1,1,1}
    }};
    
    for (const auto& c : corners) {
        std::array<float, 3> world_p = transformPoint(object_to_world, c);
        
        // Merge the point at Time = 0
        box.merge(world_p);
        
        // Merge the point at Time = 1 (Time 0 + Velocity)
        std::array<float, 3> world_p_moved = {
            world_p[0] + velocity[0],
            world_p[1] + velocity[1],
            world_p[2] + velocity[2]
        };
        box.merge(world_p_moved);
    }
    return box;
}

// ==========================================
// Rectangle Implementation (Transformed)
// ==========================================

Rectangle::Rectangle(const std::array<float, 3>& t, const std::array<float, 3>& r, const std::array<float, 3>& s, const Material& mat) {
    type = "Rectangle";
    this->material = mat;
    buildTransformationMatrices(t, r, s);
}

bool Rectangle::intersect(Hit &hit, const Ray &ray) {
    // 1. Ray to Object Space
    Ray r_loc;
    r_loc.origin = transformPoint(world_to_object, ray.origin);
    r_loc.direction = transformVector(world_to_object, ray.direction);

    // 2. Intersect Unit Square on Z=0 (from -0.5 to 0.5)
    // Plane: Z=0, Normal=(0,0,1)
    if (fabs(r_loc.direction[2]) < 1e-6f) return false;

    float t_loc = -r_loc.origin[2] / r_loc.direction[2];
    if (t_loc < 0.001f) return false;

    float hit_x = r_loc.origin[0] + t_loc * r_loc.direction[0];
    float hit_y = r_loc.origin[1] + t_loc * r_loc.direction[1];

    if (hit_x < -0.5f || hit_x > 0.5f || hit_y < -0.5f || hit_y > 0.5f) return false;

    // 3. World Props
    std::array<float, 3> p_loc = {hit_x, hit_y, 0.0f};
    std::array<float, 3> n_loc = {0.0f, 0.0f, 1.0f};

    hit.intersection_point = transformPoint(object_to_world, p_loc);
    hit.normal = transformNormal(object_to_world, n_loc);
    
    std::array<float, 3> dist_v = vecSub(hit.intersection_point, ray.origin);
    hit.t = sqrt(vecDot(dist_v, dist_v));
    hit.shape = this;
    
    // UVs: Map -0.5..0.5 to 0..1
    hit.u = hit_x + 0.5f;
    hit.v = hit_y + 0.5f;

    return true;
}

AABB Rectangle::get_bounding_box() const {
    AABB box;
    // Corners of unit square [-0.5, 0.5] on Z=0
    std::array<std::array<float, 3>, 4> corners = {{
        {-0.5, -0.5, 0}, {0.5, -0.5, 0}, {0.5, 0.5, 0}, {-0.5, 0.5, 0}
    }};
    for (const auto& c : corners) box.merge(transformPoint(object_to_world, c));
    return box;
}

// ==========================================
// Cube Implementation
// ==========================================

Cube::Cube(const std::array<float, 3>& t, const std::array<float, 3>& r, const std::array<float, 3>& s, const Material& mat) {
    type = "Cube";
    this->material = mat;
    buildTransformationMatrices(t, r, s);
}

bool Cube::intersect(Hit &hit, const Ray &ray) {
    // 1. Ray to Object Space
    Ray r_loc;
    r_loc.origin = transformPoint(world_to_object, ray.origin);
    r_loc.direction = transformVector(world_to_object, ray.direction);

    float t_min = -0.5f, t_max = 0.5f;
    float t_near = -std::numeric_limits<float>::max();
    float t_far = std::numeric_limits<float>::max();
    int hit_axis = -1;
    int hit_sign = 0;

    // Slab method
    for (int i = 0; i < 3; i++) {
        if (fabs(r_loc.direction[i]) < 1e-6f) {
            if (r_loc.origin[i] < t_min || r_loc.origin[i] > t_max) return false;
        } else {
            float t1 = (t_min - r_loc.origin[i]) / r_loc.direction[i];
            float t2 = (t_max - r_loc.origin[i]) / r_loc.direction[i];
            
            float t_entry = std::min(t1, t2);
            float t_exit = std::max(t1, t2);
            
            if (t_entry > t_near) {
                t_near = t_entry;
                hit_axis = i;
                hit_sign = (r_loc.direction[i] < 0) ? 1 : -1; // Normal points opposite to ray? No, points outward.
                // If ray direction is positive, we hit min plane (normal -1). 
                // If t1 (min) was entry, ray dir must be positive -> normal -1.
                // Wait, simpler: if t1 < t2, we hit min plane first.
                if (t1 < t2) hit_sign = -1; else hit_sign = 1;
            }
            if (t_exit < t_far) t_far = t_exit;
            if (t_near > t_far || t_far < 0) return false;
        }
    }

    float t_loc = (t_near > 0) ? t_near : t_far;
    if (t_loc < 0) return false;

    std::array<float, 3> p_loc = {
        r_loc.origin[0] + t_loc * r_loc.direction[0],
        r_loc.origin[1] + t_loc * r_loc.direction[1],
        r_loc.origin[2] + t_loc * r_loc.direction[2]
    };

    std::array<float, 3> n_loc = {0,0,0};
    if(hit_axis != -1) n_loc[hit_axis] = (float)hit_sign;

    // World Props
    hit.intersection_point = transformPoint(object_to_world, p_loc);
    hit.normal = transformNormal(object_to_world, n_loc);
    
    std::array<float, 3> dist_v = vecSub(hit.intersection_point, ray.origin);
    hit.t = sqrt(vecDot(dist_v, dist_v));
    hit.shape = this;

    // UVs (Cube Mapping)
    // Map -0.5..0.5 to 0..1 based on face
    float uc = p_loc[0] + 0.5f;
    float vc = p_loc[1] + 0.5f;
    float wc = p_loc[2] + 0.5f;
    
    if (hit_axis == 0) { hit.u = (hit_sign>0)?wc:(1.0f-wc); hit.v = vc; }
    else if(hit_axis == 1) { hit.u = uc; hit.v = (hit_sign>0)?wc:(1.0f-wc); }
    else { hit.u = (hit_sign>0)?uc:(1.0f-uc); hit.v = vc; }

    return true;
}

AABB Cube::get_bounding_box() const {
    AABB box;
    std::array<std::array<float, 3>, 8> corners = {{
        {-0.5, -0.5, -0.5}, {0.5, -0.5, -0.5}, {0.5, 0.5, -0.5}, {-0.5, 0.5, -0.5},
        {-0.5, -0.5, 0.5},  {0.5, -0.5, 0.5},  {0.5, 0.5, 0.5},  {-0.5, 0.5, 0.5}
    }};
    for (const auto& c : corners) box.merge(transformPoint(object_to_world, c));
    return box;
}

// ==========================================
// Legacy Plane Implementation
// ==========================================
Plane::Plane(std::array<float,3> corner[4], const Material& mat) {
    for(int i=0; i<4; i++) corners[i] = corner[i];
    type = "Plane";
    this->material = mat;
}

bool Plane::intersect(Hit &hit, const Ray &ray) {
    // 1. Calculate Normal
    std::array<float, 3> edge1 = vecSub(corners[1], corners[0]);
    std::array<float, 3> edge2 = vecSub(corners[2], corners[0]);
    std::array<float, 3> normal = vecCross(edge1, edge2);
    float len = sqrt(vecDot(normal, normal));
    if (len < 1e-6f) return false;
    normal = {normal[0]/len, normal[1]/len, normal[2]/len};

    // 2. Intersect Plane
    float denom = vecDot(normal, ray.direction);
    if (fabs(denom) < 1e-6f) return false;

    std::array<float, 3> p0_to_origin = vecSub(corners[0], ray.origin);
    float t = vecDot(p0_to_origin, normal) / denom;
    if (t < 0) return false;

    std::array<float, 3> intersection = {
        ray.origin[0] + t * ray.direction[0],
        ray.origin[1] + t * ray.direction[1],
        ray.origin[2] + t * ray.direction[2]
    };

    if (!isPointInQuad(intersection, corners)) return false;

    // Simple UV (Projection)
    std::array<float, 3> vec_u = vecSub(corners[1], corners[0]);
    std::array<float, 3> vec_v = vecSub(corners[3], corners[0]);
    std::array<float, 3> hit_vec = vecSub(intersection, corners[0]);
    float u = vecDot(hit_vec, vec_u) / vecDot(vec_u, vec_u);
    float v = vecDot(hit_vec, vec_v) / vecDot(vec_v, vec_v);

    hit.u = std::max(0.0f, std::min(1.0f, u));
    hit.v = std::max(0.0f, std::min(1.0f, v));
    hit.intersection_point = intersection;
    hit.normal = normal;
    hit.t = t;
    hit.shape = this;
    return true;
}

bool Plane::isPointInQuad(const std::array<float, 3>& point, const std::array<float, 3> corners[4]) {
    std::array<float, 3> edge1 = vecSub(corners[1], corners[0]);
    std::array<float, 3> edge2 = vecSub(corners[2], corners[0]);
    std::array<float, 3> normal = vecCross(edge1, edge2); // Use calculated normal
    float l = sqrt(vecDot(normal,normal)); normal[0]/=l; normal[1]/=l; normal[2]/=l;

    if (isPointInTriangle(point, corners[1], corners[3], corners[2], normal)) return true;
    if (isPointInTriangle(point, corners[0], corners[1], corners[2], normal)) return true;
    return false;
}

AABB Plane::get_bounding_box() const {
    AABB box;
    const float padding = 1e-4f;
    for (int i = 0; i < 4; ++i) box.merge(corners[i]);
    box.min_point[0] -= padding; box.min_point[1] -= padding; box.min_point[2] -= padding;
    box.max_point[0] += padding; box.max_point[1] += padding; box.max_point[2] += padding;
    return box;
}