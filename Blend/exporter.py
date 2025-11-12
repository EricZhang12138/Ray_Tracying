import bpy
from collections import defaultdict
from mathutils import Vector
import json 

class Process_objects:
    def __init__(self, objects):
        self.objects = objects
        # Use defaultdict(list) to store a list of items for each object type
        self.dict = defaultdict(list)

    def _get_material_properties(self, blender_object):
        """
        Helper function to extract material properties from an object.
        It looks for Diffuse BSDF and Glossy BSDF nodes.
        It also provides defaults for Blinn-Phong and Whitted-style properties.
        """
        
        # --- Start with default material values ---
        # These will be used if a material or node is not found.
        material_data = {
            # Blinn-Phong (from nodes)
            'diffuse_color': [0.8, 0.8, 0.8], # Default grey
            'specular_color': [1.0, 1.0, 1.0],# Default white
            'roughness': 0.5,                  # C++ side should convert this to shininess
            
            # Blinn-Phong (coefficients)
            # We'll just hard-code these for now.
            # You can add them as Custom Properties in Blender if you want to.
            'k_ambient': 0.1,
            'k_diffuse': 0.9,
            'k_specular': 0.3,
            
            # Whitted-Style (for later)
            'reflectivity': 0.0,
            'transparency': 0.0,
            'refractive_index': 1.0,

            # Texture (for later)
            'texture_file': ""
        }

        # --- Try to get real values from material nodes ---
        try:
            # Check if the object has a material
            if not blender_object.material_slots:
                return material_data

            material = blender_object.material_slots[0].material
            if not material or not material.node_tree:
                return material_data

            nodes = material.node_tree.nodes

            # Find the first Diffuse and Glossy nodes
            # The 'next' function finds the first item or returns None
            diffuse_node = next((n for n in nodes if n.type == 'BSDF_DIFFUSE'), None)
            glossy_node = next((n for n in nodes if n.type == 'BSDF_GLOSSy'), None)

            # Extract data from the Diffuse node
            if diffuse_node:
                color = diffuse_node.inputs['Color'].default_value
                material_data['diffuse_color'] = list(color)[:3] # Get (R, G, B)

            # Extract data from the Glossy node
            if glossy_node:
                color = glossy_node.inputs['Color'].default_value
                material_data['specular_color'] = list(color)[:3] # Get (R, G, B)
                
                roughness = glossy_node.inputs['Roughness'].default_value
                material_data['roughness'] = roughness
                
        except Exception as e:
            # If anything goes wrong (e.g., no material), just print a warning
            # and return the default values.
            print(f"Warning: Could not get material for {blender_object.name}: {e}")
        
        return material_data


    def process(self) -> dict:
        """
        Processes objects and converts their data into JSON-serializable types.
        """
        for obj in self.objects:
            if obj.type == 'MESH':
                
                # Get material for all mesh objects
                material_data = self._get_material_properties(obj)
                
                if 'Sphere' in obj.name:
                    sphere_data = {
                        'location': list(obj.location), 
                        'radius': obj.scale.x,
                        'material': material_data  
                    }
                    self.dict['spheres'].append(sphere_data)
                    
                elif 'Cube' in obj.name:
                    cube_data = {
                        'translation': list(obj.location), 
                        'rotation': list(obj.rotation_euler), 
                        'scale': obj.scale.x,
                        'material': material_data  
                    }
                    self.dict['cubes'].append(cube_data)
                    
                elif 'Plane' in obj.name:
                    plane_data = {
                        'corners': [list(obj.matrix_world @ v.co) for v in obj.data.vertices],
                        'material': material_data  
                    }
                    self.dict['planes'].append(plane_data)
                    
            elif obj.type == 'CAMERA':
                camera_data = {
                    'location': list(obj.location), 
                    'gaze_vector': list(obj.matrix_world.to_quaternion() @ Vector((0.0, 0.0, -1.0))),
                    'focal_length': obj.data.lens,
                    'sensor_width': obj.data.sensor_width,
                    'sensor_height': obj.data.sensor_height,
                    'up_vector': list(obj.matrix_world.to_quaternion() @ Vector((0.0, 1.0, 0.0))),
                }
                self.dict['cameras'].append(camera_data)
                
            elif obj.type == 'LIGHT' and obj.data.type == 'POINT':
                light_data = {
                    'location': list(obj.location), 
                    'intensity': obj.data.energy,
                    'color': list(obj.data.color) 
                }
                self.dict['lights'].append(light_data)

        self.dict['render'] = {
            'resolution_x': bpy.context.scene.render.resolution_x,
            'resolution_y': bpy.context.scene.render.resolution_y
        }
        
        return self.dict
    
    def save_as_json(self, filepath):
        # Using 'with' is safer, automatically handles closing the file
        with open(filepath, 'w') as f:
            json.dump(self.dict, f, indent=4)
        
        print(f"Scene data successfully exported to {filepath}")

        
                    
# only when __name__ == "__main__" will we run the code below
if __name__ == "__main__":
    
    # Set a default output path, preferably relative or in a known location
    # This path is an example; change it to where you want the file saved.
    output_filepath = '/Users/ericzhang/Documents/Computer_graphics/Coursework/s2286795/ASCII/scene.json' # Saves to the same directory as your .blend file
    
    # Try to build a full path from the .blend file path
    if bpy.data.filepath:
        import os
        # Gets the directory where the .blend file is saved
        blend_dir = os.path.dirname(bpy.data.filepath) 
        output_filepath = os.path.join(blend_dir, 'scene.json')
    else:
        print("Warning: .blend file is not saved. Saving scene.json to an unknown location.")

    
    processor = Process_objects(bpy.data.objects)
    processor.process()
    processor.save_as_json(output_filepath)