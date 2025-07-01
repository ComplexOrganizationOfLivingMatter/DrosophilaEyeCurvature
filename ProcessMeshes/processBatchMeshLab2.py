import pymeshlab
import os

def remesh_ply_files(input_folder):
    # Get a list of all .ply files in the folder
    ply_files = [f for f in os.listdir(input_folder) if f.endswith('.stl')]
    
    for ply_file in ply_files:
        # Load the mesh
        ms = pymeshlab.MeshSet()
        input_path = os.path.join(input_folder, ply_file)
        ms.load_new_mesh(input_path)

        # Apply remeshing five times
        for _ in range(25):
            ms.apply_filter('meshing_isotropic_explicit_remeshing')

        # Generate the new filename by adding "_remeshed" before the file extension
        filename, _ = os.path.splitext(ply_file)
        output_filename = f"{filename}_remeshed.stl"
        output_path = os.path.join(input_folder, output_filename)
        
        # Save the remeshed file
        ms.save_current_mesh(output_path)

        print(f"Saved remeshed file: {output_path}")

# Set the folder containing the .ply files
input_folder = '/media/juan/ce04cf3c-0227-4dac-9535-9f014799bf45/juan/Projects/Auxetics/Auxetic_Eyes/Eyes_Expansion_Codes/Analysis/Pressure/D_Mauritiana/mauritaniaTam16 male42h25 retina12 actin/Pupa_Iterations_Inter_DS_LES_stl/'

# Call the function to process all .ply files in the folder
remesh_ply_files(input_folder)

