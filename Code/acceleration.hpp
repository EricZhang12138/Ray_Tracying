#ifndef ACCELERATION_HPP
#define ACCELERATION_HPP

#include "shapes.hpp"
#include <vector>

struct node {
    node() = default;
    node(AABB bounding_box);
    AABB bounding_box;
    std::unique_ptr<node> left;
    std::unique_ptr<node> right;
    std::vector<Shapes*> objects;
};

class BVH {
public:
    std::vector<Shapes*> shape_list;
    node root;
    BVH(const std::vector<Shapes*>& shapes_ptrs); 
    Hit intersect(Ray ray, node& root);
private:
    void construct_tree(int start, int end, node& root);
    
    

private:
    std::vector<Hit> temp_hit_vec;    // may not be a very good design here 
    void intersect_helper(Ray ray, node& root);
};

#endif