import bpy
from collections import defaultdict
from mathutils import Vector
import json 


class Process_objects:
    def __init__(self, objects):
        self.objects = objects
        # Use defaultdict(list) to store a list of items for each object type
        self.dict = defaultdict(list)

    def process(self) -> dict:
        """
        Processes objects and converts their data into JSON-serializable types.
        """
        for obj in self.objects:
            if obj.type == 'MESH':
                if 'Sphere' in obj.name:
                    sphere_data = {
                        'location': list(obj.location), # Convert Vector to list
                        'radius': obj.scale.x
                    }
                    self.dict['spheres'].append(sphere_data)
                elif 'Cube' in obj.name:
                    cube_data = {
                        'translation': list(obj.location), # Convert Vector to list
                        'rotation': list(obj.rotation_euler), # Convert Euler to list
                        'scale': obj.scale.x
                    }
                    self.dict['cubes'].append(cube_data)
                elif 'Plane' in obj.name:
                    plane_data = {
                        # Convert each corner Vector to a list
                        'corners': [list(obj.matrix_world @ v.co) for v in obj.data.vertices]
                    }
                    self.dict['planes'].append(plane_data)
            elif obj.type == 'CAMERA':
                camera_data = {
                    'location': list(obj.location), # Convert Vector to list
                    # Convert the final gaze Vector to a list
                    'gaze_vector': list(obj.matrix_world.to_quaternion() @ Vector((0.0, 0.0, -1.0))),
                    'focal_length': obj.data.lens,
                    'sensor_width': obj.data.sensor_width,
                    'sensor_height': obj.data.sensor_height,
                    'up_vector': list(obj.matrix_world.to_quaternion() @ Vector((0.0, 1.0, 0.0))),
                }
                self.dict['cameras'].append(camera_data)
            elif obj.type == 'LIGHT' and obj.data.type == 'POINT':
                light_data = {
                    'location': list(obj.location), # Convert Vector to list
                    'intensity': obj.data.energy,
                    'color': list(obj.data.color) # Convert Color to list
                }
                self.dict['lights'].append(light_data)

        self.dict['render'] = {
            'resolution_x': bpy.context.scene.render.resolution_x,
            'resolution_y': bpy.context.scene.render.resolution_y
        }
        
        return self.dict
    
    def save_as_json(self, filepath):
        f = open(filepath, 'w')
        try: 
            json.dump(self.dict, f, indent = 4)
        finally:
            f.close()
        
        print(f"Scene data successfully exported to {filepath}")

        
                    
# only when __name__ == "__main__" will we run the code below, which means it is only run when you directly run 
# this program as opposed to importing this file into other modules   
if __name__ == "__main__":
    
    output_filepath = '/Users/ericzhang/Documents/Computer_graphics/Coursework/s2286795/    ASCII/scene.json'
    
    
    processor = Process_objects(bpy.data.objects)
    
    
    processor.process()
    

    processor.save_as_json(output_filepath)
