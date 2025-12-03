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
        Extracts material data. Handles Principled, Glass, and Mix Shaders.
        CRITICAL UPDATE: Now finds the 'Tint' color from Multiply nodes.
        """
        # --- Defaults ---
        material_data = {
            'diffuse_color': [0.8, 0.8, 0.8], # Default Grey
            'specular_color': [0.0, 0.0, 0.0],
            'roughness': 0.5,
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
        
        # --- Helper 1: Extract color from a node (handles defaults) ---
        def get_input_color(node, input_name):
            if input_name in node.inputs:
                return list(node.inputs[input_name].default_value)[:3]
            return [1.0, 1.0, 1.0] # Return White if not found, to avoid darkening textures

        # --- Helper 2: Recursive Texture Finder ---
        def find_texture_recursive(node_input):
            if not node_input.is_linked:
                return ""
            
            link = node_input.links[0]
            from_node = link.from_node
            
            if from_node.type == 'TEX_IMAGE' and from_node.image:
                return os.path.basename(from_node.image.filepath)
            
            # Also check BUMP nodes
            if from_node.type == 'BUMP':
                # Search the Height input of Bump node
                if from_node.inputs['Height'].is_linked:
                    return find_texture_recursive(from_node.inputs['Height'])
            
            if from_node.type in ['MIX_RGB', 'MATH', 'MIX_SHADER']:
                for input_idx in range(min(2, len(from_node.inputs))):
                    found_tex = find_texture_recursive(from_node.inputs[input_idx])
                    if found_tex:
                        return found_tex
            return ""

        # --- Helper 3: Tint Finder (NEW) ---
        # If we have a Multiply node, find the color that ISN'T the texture
        def find_tint_color(node_input):
            if not node_input.is_linked:
                # If not linked, just return the default color of the socket
                return list(node_input.default_value)[:3]
            
            link = node_input.links[0]
            from_node = link.from_node

            # If it's a Mix RGB (Multiply, Mix, Add, etc.)
            if from_node.type == 'MIX_RGB':
                # Check Input 1 and Input 2
                c1_linked = from_node.inputs[1].is_linked
                c2_linked = from_node.inputs[2].is_linked
                
                # If Input 1 is the texture (linked), return Input 2's color
                if c1_linked and not c2_linked:
                    return list(from_node.inputs[2].default_value)[:3]
                
                # If Input 2 is the texture (linked), return Input 1's color
                if c2_linked and not c1_linked:
                    return list(from_node.inputs[1].default_value)[:3]

            # Default: If we can't find a tint, return white (so we don't change the texture)
            return [1.0, 1.0, 1.0]

        # --- 1. Principled BSDF ---
        principled = next((n for n in nodes if n.type == 'BSDF_PRINCIPLED'), None)
        if principled:
            # Check for tint in Base Color
            material_data['diffuse_color'] = find_tint_color(principled.inputs['Base Color'])
            
            # Fallback: if find_tint returned white, but the socket is NOT linked, use the actual slider color
            if not principled.inputs['Base Color'].is_linked:
                 material_data['diffuse_color'] = list(principled.inputs['Base Color'].default_value)[:3]

            material_data['roughness'] = principled.inputs['Roughness'].default_value
            material_data['reflectivity'] = principled.inputs['Metallic'].default_value
            
            if 'Transmission Weight' in principled.inputs:
                material_data['transparency'] = principled.inputs['Transmission Weight'].default_value
            elif 'Transmission' in principled.inputs:
                material_data['transparency'] = principled.inputs['Transmission'].default_value
            
            if 'IOR' in principled.inputs:
                material_data['refractive_index'] = principled.inputs['IOR'].default_value

            material_data['texture_file'] = find_texture_recursive(principled.inputs['Base Color'])
            return material_data

        # --- 2. Glass BSDF ---
        glass = next((n for n in nodes if n.type == 'BSDF_GLASS'), None)
        if glass:
            material_data['diffuse_color'] = get_input_color(glass, 'Color')
            material_data['specular_color'] = [1.0, 1.0, 1.0]
            material_data['transparency'] = 1.0
            material_data['refractive_index'] = glass.inputs['IOR'].default_value
            material_data['roughness'] = glass.inputs['Roughness'].default_value
            return material_data

        # --- 3. Mix Shader (Screenshot Scenario) ---
        diffuse_node = next((n for n in nodes if n.type == 'BSDF_DIFFUSE'), None)
        glossy_node = next((n for n in nodes if n.type == 'BSDF_GLOSSY'), None)
        mix_node = next((n for n in nodes if n.type == 'MIX_SHADER'), None)

        if diffuse_node:
            # Search Color input first
            material_data['texture_file'] = find_texture_recursive(diffuse_node.inputs['Color'])
            print(f"DEBUG: Texture from Color input: '{material_data['texture_file']}'")
            
            # If not found, also check Normal input (in case texture is in Bump)
            if not material_data['texture_file'] and diffuse_node.inputs['Normal'].is_linked:
                print("DEBUG: Color was empty, checking Normal input...")
                material_data['texture_file'] = find_texture_recursive(diffuse_node.inputs['Normal'])
                print(f"DEBUG: Texture from Normal input: '{material_data['texture_file']}'")
            
            material_data['diffuse_color'] = find_tint_color(diffuse_node.inputs['Color'])
        
        if glossy_node:
            material_data['specular_color'] = get_input_color(glossy_node, 'Color')
            material_data['roughness'] = glossy_node.inputs['Roughness'].default_value

            if mix_node:
                fac = mix_node.inputs['Fac'].default_value
                
                is_glossy_top = False
                if len(mix_node.inputs) > 1:
                     for link in mix_node.inputs[1].links:
                        if link.from_node == glossy_node:
                            is_glossy_top = True
                            break
                
                if is_glossy_top:
                    k_spec = 1.0 - fac
                    k_diff = fac
                else:
                    k_spec = fac
                    k_diff = 1.0 - fac

                material_data['k_specular'] = k_spec
                material_data['k_diffuse'] = k_diff
                material_data['reflectivity'] = k_spec
            
            else:
                material_data['k_specular'] = 1.0
                material_data['k_diffuse'] = 0.0
                material_data['reflectivity'] = 1.0

        return material_data

    def process(self) -> dict:
        for obj in self.objects:
            if obj.type == 'MESH':
                material_data = self._get_material_properties(obj)
                
                # --- 1. SPHERES (and OVALS) ---
                if 'Sphere' in obj.name:
                    # Blender dimensions = Diameter.
                    # C++ Unit Sphere has Radius = 1.
                    # Therefore, Scale = Dimension / 2.0 to match the visual size.
                    scale_vector = [
                        obj.dimensions.x / 2.0,
                        obj.dimensions.y / 2.0,
                        obj.dimensions.z / 2.0
                    ]

                    sphere_data = {
                        'location': list(obj.location), 
                        'rotation': list(obj.rotation_euler), # Now exporting Rotation
                        'scale': scale_vector,                # Now exporting Vector Scale
                        # 'radius' is removed because 'scale' handles it now
                        'velocity': list(obj.get("velocity", [0.0, 0.0, 0.0])),
                        'material': material_data  
                    }
                    self.dict['spheres'].append(sphere_data)
                    
                # --- 2. CUBES ---
                elif 'Cube' in obj.name:
                    # C++ Unit Cube is size 1.0 (-0.5 to 0.5).
                    # Blender dimensions are actual size.
                    # So Scale = Dimensions.
                    scale_vector = [
                        obj.dimensions.x, 
                        obj.dimensions.y, 
                        obj.dimensions.z
                    ]

                    cube_data = {
                        'translation': list(obj.location), 
                        'rotation': list(obj.rotation_euler), 
                        'scale': scale_vector, # Exporting full vector allows non-uniform cubes
                        'material': material_data  
                    }
                    self.dict['cubes'].append(cube_data)
                    
                # --- 3. RECTANGLES (Transformed Planes) ---
                elif 'Plane' in obj.name:
                    # C++ Unit Square is size 1.0 (-0.5 to 0.5).
                    # Blender Plane dimensions are actual size.
                    # So Scale = Dimensions.
                    scale_vector = [
                        obj.dimensions.x,
                        obj.dimensions.y,
                        1.0 # Z scale is irrelevant for a 2D plane
                    ]

                    rect_data = {
                        'translation': list(obj.location),
                        'rotation': list(obj.rotation_euler),
                        'scale': scale_vector,
                        'material': material_data  
                    }
                    
                    # We append to 'rectangles' which matches your new C++ json_loader
                    self.dict['rectangles'].append(rect_data)
                    
            # --- 4. CAMERAS ---
            elif obj.type == 'CAMERA':
                camera_data = {
                    'location': list(obj.location), 
                    'gaze_vector': list(obj.matrix_world.to_quaternion() @ Vector((0.0, 0.0, -1.0))),
                    'focal_length': obj.data.lens,
                    'sensor_width': obj.data.sensor_width,
                    'sensor_height': obj.data.sensor_height,
                    'up_vector': list(obj.matrix_world.to_quaternion() @ Vector((0.0, 1.0, 0.0))),
                    'aperture': obj.get("aperture", 0.0), 
                    'focus_dist': obj.get("focus_dist", obj.data.dof.focus_distance)
                }
                self.dict['cameras'].append(camera_data)
                
            # --- 5. LIGHTS ---
            elif obj.type == 'LIGHT' and obj.data.type == 'POINT':
                light_data = {
                    'location': list(obj.location), 
                    'intensity': obj.data.energy,
                    'color': list(obj.data.color),
                    'radius': obj.data.shadow_soft_size 
                }
                self.dict['lights'].append(light_data)

        # Standard render settings
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
    output_filepath = os.path.join("..","..", "ASCII", "scene.json")
    
    # Create directory if it doesn't exist
    output_dir = os.path.dirname(output_filepath)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    processor = Process_objects(bpy.data.objects)
    processor.process()
    processor.save_as_json(output_filepath)
