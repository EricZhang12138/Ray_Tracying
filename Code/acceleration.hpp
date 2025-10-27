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
    BVH (std::vector<Shapes>& shape_list);
    void construct_tree(int start, int end, node& root);
    Hit intersect(Ray ray, node& root);
    

private:
    std::vector<Hit> temp_hit_vec;    // may not be a very good design here 
    void intersect_helper(Ray ray, node& root);
};