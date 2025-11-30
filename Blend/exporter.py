import bpy
from collections import defaultdict
from mathutils import Vector
import json 
import os

class Process_objects:
    def __init__(self, objects):
        self.objects = objects
        self.dict = defaultdict(list)

    def _get_material_properties(self, blender_object):
        """
        Extracts material data, handling the modern Principled BSDF node
        to satisfy the Blinn-Phong shading requirement.
        """
        # Defaults
        material_data = {
            'diffuse_color': [0.8, 0.8, 0.8],
            'specular_color': [1.0, 1.0, 1.0],
            'roughness': 0.5, # C++ converts this to shininess
            'k_ambient': 0.1,
            'k_diffuse': 0.9,
            'k_specular': 0.3,
            'reflectivity': 0.0,
            'transparency': 0.0,
            'refractive_index': 1.0,
            'texture_file': ""
        }

        if not blender_object.material_slots:
            return material_data

        material = blender_object.material_slots[0].material
        if not material or not material.node_tree:
            return material_data

        nodes = material.node_tree.nodes
        
        # --- Handle Principled BSDF ---
        principled = next((n for n in nodes if n.type == 'BSDF_PRINCIPLED'), None)
        
        if principled:
            color = principled.inputs['Base Color'].default_value
            material_data['diffuse_color'] = list(color)[:3]
            material_data['roughness'] = principled.inputs['Roughness'].default_value
            material_data['reflectivity'] = principled.inputs['Metallic'].default_value
            material_data['transparency'] = principled.inputs['Transmission Weight'].default_value
            
            if 'IOR' in principled.inputs:
                material_data['refractive_index'] = principled.inputs['IOR'].default_value
            
            # --- Check for Texture Image ---
            if principled.inputs['Base Color'].is_linked:
                link = principled.inputs['Base Color'].links[0]
                node = link.from_node
                if node.type == 'TEX_IMAGE' and node.image:
                    material_data['texture_file'] = os.path.basename(node.image.filepath)

            return material_data

        # --- Fallback: Legacy Nodes ---
        diffuse_node = next((n for n in nodes if n.type == 'BSDF_DIFFUSE'), None)
        glossy_node = next((n for n in nodes if n.type == 'BSDF_GLOSSY'), None)

        if diffuse_node:
            color = diffuse_node.inputs['Color'].default_value
            material_data['diffuse_color'] = list(color)[:3]

        if glossy_node:
            color = glossy_node.inputs['Color'].default_value
            material_data['specular_color'] = list(color)[:3]
            material_data['roughness'] = glossy_node.inputs['Roughness'].default_value
        
        return material_data


    def process(self) -> dict:
        for obj in self.objects:
            if obj.type == 'MESH':
                material_data = self._get_material_properties(obj)
                
                if 'Sphere' in obj.name:
                    sphere_data = {
                        'location': list(obj.location), 
                        'radius': obj.dimensions.x / 2.0,
                        # [UPDATED] Added Velocity for Motion Blur
                        'velocity': list(obj.get("velocity", [0.0, 0.0, 0.0])),
                        'material': material_data  
                    }
                    self.dict['spheres'].append(sphere_data)
                    
                elif 'Cube' in obj.name:
                    cube_data = {
                        'translation': list(obj.location), 
                        'rotation': list(obj.rotation_euler), 
                        'scale': obj.dimensions.x / 2.0,
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
                    # [UPDATED] Added Aperture and Focus Distance for Depth of Field
                    # We check Custom Properties first (obj.get), then fallback to native Blender DoF settings
                    'aperture': obj.get("aperture", 0.0), 
                    'focus_dist': obj.get("focus_dist", obj.data.dof.focus_distance)
                }
                self.dict['cameras'].append(camera_data)
                
            elif obj.type == 'LIGHT' and obj.data.type == 'POINT':
                light_data = {
                    'location': list(obj.location), 
                    'intensity': obj.data.energy,
                    'color': list(obj.data.color),
                    'radius': obj.data.shadow_soft_size 
                }
                self.dict['lights'].append(light_data)

        self.dict['render'] = {
            'resolution_x': bpy.context.scene.render.resolution_x,
            'resolution_y': bpy.context.scene.render.resolution_y
        }
        
        return self.dict
    
    def save_as_json(self, filepath):
        with open(filepath, 'w') as f:
            json.dump(self.dict, f, indent=4)
        print(f"Scene data successfully exported to {filepath}")

if __name__ == "__main__":
    # --- FIX PATH: Use your specific path ---
    output_filepath = '/Users/ericzhang/Documents/Computer_graphics/Coursework/s2286795/ASCII/scene.json'
    
    # Create directory if it doesn't exist
    output_dir = os.path.dirname(output_filepath)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    processor = Process_objects(bpy.data.objects)
    processor.process()
    processor.save_as_json(output_filepath)