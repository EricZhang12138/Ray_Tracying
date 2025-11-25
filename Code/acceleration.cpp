#include "acceleration.hpp"
#include <algorithm>


node::node(AABB bounding_box): bounding_box(bounding_box){};

BVH::BVH(const std::vector<Shapes*>& shapes_ptrs) { // New version takes pointers
    // Directly copy the pointers
    shape_list = shapes_ptrs; // shape_list is std::vector<Shapes*>

    // Start recursive construction
    if (!shape_list.empty()) {
        // Pass a COPY of the root node if construct_tree modifies it directly
        // If construct_tree takes a reference and initializes it, passing root is fine.
        // Assuming your construct_tree works correctly with a reference to root:
        construct_tree(0, shape_list.size(), root);
    }
}

void BVH::construct_tree(int start, int end, node& current_node) {
    AABB node_box;
    for (int i = start; i < end; i++) {
        node_box.merge(shape_list[i]->get_bounding_box());
    }
    current_node.bounding_box = node_box;
    
    int object_count = end - start;
    
    // 4 is the number of shapes a node has 
    if (object_count <= 4) {
        for (int i = start; i < end; i++) {
            current_node.objects.push_back(shape_list[i]);
        }
    
        // A leaf node has NO children.
        current_node.left = nullptr;
        current_node.right = nullptr;
        return; 
    }
    //internal node
    
    // Find the longest axis to split along
    int axis = node_box.get_longest_axis(); 
    
    // Sort the objects in THIS range [start, end) along THAT axis
    auto comparator = [axis](Shapes* a, Shapes* b) {
        AABB box_a = a->get_bounding_box();
        AABB box_b = b->get_bounding_box();
        float center_a = (box_a.min_point[axis] + box_a.max_point[axis]) / 2.0f;
        float center_b = (box_b.min_point[axis] + box_b.max_point[axis]) / 2.0f;
        return center_a < center_b;
    };
    
    // Sort the sub-list we are looking at
    std::sort(shape_list.begin() + start, shape_list.begin() + end, comparator);

    int mid = (start + end) / 2;
    
    current_node.left = std::make_unique<node>();
    current_node.right = std::make_unique<node>();

    construct_tree(start, mid, *current_node.left);  
    construct_tree(mid, end, *current_node.right); 
}


void BVH::intersect_helper(Ray ray, node& root_node){
    AABB area = root_node.bounding_box;
    bool is_intersected = area.intersect(ray);

    if (is_intersected){ // if intersects
 

        // if internal nodes, at least have one root_node.left or root_node.right
        if (root_node.left || root_node.right){
            if (root_node.left){
            node& left_node = *(root_node.left);
            intersect_helper(ray, left_node);
            }
            if (root_node.right){
            node& right_node = *(root_node.right);
            intersect_helper(ray, right_node);
            }

        }else{ // leaf nodes
            for (Shapes* shape: root_node.objects){

                // with polymorphism, I don't have to determine 
                // what object it is before calling intersect
                Hit hit;
                if (shape -> intersect(hit,ray)){
                    temp_hit_vec.push_back(hit); // putting the hit structure in temp_hit_vec
                }
            }
        }

    }else{ // if didn't intersect
        return;  // do nothing
    }
}


Hit BVH::intersect(Ray ray, node& root_node){
    intersect_helper(ray, root_node);
    if (temp_hit_vec.empty()){
        Hit miss;
        miss.t = std::numeric_limits<float>::max();
        miss.shape = nullptr;
        temp_hit_vec.clear();
        return miss;
    }
    auto closest = std::min_element(temp_hit_vec.begin(),temp_hit_vec.end(),[](Hit hit_1, Hit hit_2){
        return hit_1.t<hit_2.t;
    });
    Hit result = *closest;
    temp_hit_vec.clear();
    return result;
}


// ... existing includes ...

// 1. Implement the Brute Force Linear Search
Hit BVH::intersect_linear(const Ray& ray) {
    Hit closest_hit;
    closest_hit.t = std::numeric_limits<float>::max();
    closest_hit.shape = nullptr;

    // Simply loop through the vector of all shapes
    for (Shapes* shape : shape_list) {
        Hit current_hit;
        if (shape->intersect(current_hit, ray)) {
            if (current_hit.t < closest_hit.t) {
                closest_hit = current_hit;
            }
        }
    }
    return closest_hit;
}

// 2. Implement the Master Getter
Hit BVH::get_intersection(const Ray& ray, bool use_acceleration) {
    if (use_acceleration) {
        // Call your existing tree traversal
        // Note: You might need to clear temp_hit_vec here or ensure your existing intersect handles it
        return intersect(ray, root); 
    } else {
        // Call the new linear search
        return intersect_linear(ray);
    }
}