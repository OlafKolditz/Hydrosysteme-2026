import pyvista as pv
import glob
import os

# ----------------------------------------------------------------------
# 1. Load a single VTK file (static plot)
# ----------------------------------------------------------------------
file_path = "file0.vtk"   # Change to your desired time step
if not os.path.exists(file_path):
    raise FileNotFoundError(f"VTK file not found: {file_path}")

# Read the unstructured grid
mesh = pv.read(file_path)

# Check the available point data arrays
print("Point data arrays:", mesh.point_data.keys())

# Extract the 'HEAD' scalar field (this is u_new)
head = mesh.point_data["HEAD"]   # As written in OutputResultsVTK

# ----------------------------------------------------------------------
# 2. Plot as a 2D colour map (top view)
# ----------------------------------------------------------------------
plotter = pv.Plotter()
plotter.add_mesh(mesh, scalars="HEAD", cmap="viridis", show_edges=False, 
                 point_size=5.0, render_points_as_spheres=False)
plotter.view_xy()   # Top-down view (x-y plane)
plotter.add_axes()
plotter.add_text("Hydraulic Head (m)", font_size=12)
plotter.show()

# ----------------------------------------------------------------------
# 3. Plot as a 3D surface (warped by scalar)
# ----------------------------------------------------------------------
# Warp the mesh by the head value (multiply by a scaling factor for visibility)
warped = mesh.warp_by_scalar(scalars="HEAD", factor=0.01)  # adjust factor
plotter = pv.Plotter()
plotter.add_mesh(warped, scalars="HEAD", cmap="viridis", show_edges=False)
plotter.add_axes()
plotter.add_text("Warped Surface - Hydraulic Head", font_size=12)
plotter.show()

# ----------------------------------------------------------------------
# 4. Animate multiple time steps (if you have many VTK files)
# ----------------------------------------------------------------------
# Find all VTK files (assuming they are named file*.vtk)
vtk_files = sorted(glob.glob("file*.vtk"))
if not vtk_files:
    print("No VTK files found for animation.")
else:
    # Read all meshes
    meshes = [pv.read(f) for f in vtk_files]
    
    # Create plotter and add first mesh
    plotter = pv.Plotter()
    plotter.open_gif("head_animation.gif")   # Save as GIF
    
    # Add initial mesh
    actor = plotter.add_mesh(meshes[0], scalars="HEAD", cmap="viridis", 
                             show_edges=False)
    plotter.view_xy()
    plotter.add_axes()
    
    # Loop over time steps and update scalars
    for i, mesh in enumerate(meshes):
        actor.mapper.dataset.point_data["HEAD"] = mesh.point_data["HEAD"]
        actor.mapper.dataset.Modified()
        plotter.write_frame()   # Write to GIF
    
    plotter.close()
    print("Animation saved as head_animation.gif")