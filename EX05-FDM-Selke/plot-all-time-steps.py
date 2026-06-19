import pyvista as pv
import glob
import os
import numpy as np
from math import ceil

# ----------------------------------------------------------------------
# 1. Find all VTK files and sort them
# ----------------------------------------------------------------------
vtk_files = sorted(glob.glob("file*.vtk"))
if not vtk_files:
    raise FileNotFoundError("No VTK files found (pattern: file*.vtk)")

print(f"Found {len(vtk_files)} time steps.")

# Load all meshes (this can be memory-heavy for many files; use streaming if needed)
meshes = [pv.read(f) for f in vtk_files]

# Check that the 'HEAD' scalar exists in the first mesh
if "HEAD" not in meshes[0].point_data:
    raise KeyError("'HEAD' not found in point data. Check your VTK output.")

# ----------------------------------------------------------------------
# 2. Interactive animation with time slider
# ----------------------------------------------------------------------
def create_slider_animation(meshes, cmap="viridis", title="Hydraulic Head (m)"):
    """Create a PyVista plot with a slider to scroll through time steps."""
    plotter = pv.Plotter()
    
    # Add initial mesh
    actor = plotter.add_mesh(meshes[0], scalars="HEAD", cmap=cmap, show_edges=False)
    plotter.view_xy()
    plotter.add_axes()
    plotter.add_text(f"{title} - Time step 0", font_size=12)
    
    # Define slider callback
    def update_step(value):
        step = int(round(value))
        if step < 0 or step >= len(meshes):
            return
        mesh = meshes[step]
        # Update the scalar data of the existing actor
        actor.mapper.dataset.point_data["HEAD"] = mesh.point_data["HEAD"]
        actor.mapper.dataset.Modified()
        # Update title
        plotter.remove_actor(plotter.renderer._actors.get("text_actor", None))  # simple hack
        plotter.add_text(f"{title} - Time step {step}", font_size=12, name="text_actor")
        plotter.render()
    
    # Add slider
    plotter.add_slider_widget(update_step, rng=[0, len(meshes)-1], value=0,
                              title="Time step", pointa=(0.2, 0.1), pointb=(0.8, 0.1),
                              style="modern", fmt="%3.0f")
    plotter.show()
    
# Uncomment to run the interactive slider:
# create_slider_animation(meshes)

# ----------------------------------------------------------------------
# 3. Static grid of selected time steps (e.g., first 16 steps)
# ----------------------------------------------------------------------
def plot_time_steps_grid(meshes, ncols=4, cmap="viridis", title="Hydraulic Head"):
    """Plot a grid of time steps (first ncols*ncols steps)."""
    nplots = min(len(meshes), ncols * ncols)
    if nplots == 0:
        return
    
    nrows = ceil(nplots / ncols)
    plotter = pv.Plotter(shape=(nrows, ncols))
    
    for idx in range(nplots):
        row = idx // ncols
        col = idx % ncols
        plotter.subplot(row, col)
        plotter.add_mesh(meshes[idx], scalars="HEAD", cmap=cmap, show_edges=False)
        plotter.view_xy()
        plotter.add_text(f"t = {idx}", font_size=8, position='upper_left')
        plotter.add_axes(line_width=1, xlabel='', ylabel='', zlabel='')
    
    plotter.link_views()  # Link camera for consistent zoom/pan
    plotter.show()

# Uncomment to plot the grid:
# plot_time_steps_grid(meshes, ncols=4)

# ----------------------------------------------------------------------
# 4. Save a GIF animation (all time steps)
# ----------------------------------------------------------------------
def save_gif_animation(meshes, filename="head_animation.gif", cmap="viridis", title="Hydraulic Head"):
    """Save all time steps as an animated GIF."""
    if len(meshes) == 0:
        return
    plotter = pv.Plotter()
    plotter.open_gif(filename)
    actor = plotter.add_mesh(meshes[0], scalars="HEAD", cmap=cmap, show_edges=False)
    plotter.view_xy()
    plotter.add_axes()
    plotter.add_text(f"{title} - t=0", font_size=12, name="text_actor")
    plotter.write_frame()  # write initial frame
    
    for i, mesh in enumerate(meshes):
        actor.mapper.dataset.point_data["HEAD"] = mesh.point_data["HEAD"]
        actor.mapper.dataset.Modified()
        # Update text (remove old, add new)
        plotter.remove_actor("text_actor")
        plotter.add_text(f"{title} - t={i}", font_size=12, name="text_actor")
        plotter.write_frame()
    plotter.close()
    print(f"GIF saved as {filename}")

# Uncomment to save the GIF:
# save_gif_animation(meshes)

# ----------------------------------------------------------------------
# 5. (Optional) Use the slider interactively; uncomment the line below
# ----------------------------------------------------------------------
if __name__ == "__main__":
    # By default, run the interactive slider if there are not too many steps
    if len(meshes) <= 100:
        create_slider_animation(meshes)
    else:
        print("Many time steps detected. Consider using save_gif_animation() or plot grid.")
        # Alternatively, plot a grid of first 16 steps
        plot_time_steps_grid(meshes, ncols=4)