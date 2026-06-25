import numpy as np
import pyvista as pv
import warnings

# ------------------------------------------------------------
# Shape function evaluation
# ------------------------------------------------------------
def lagrange_shape_functions(x, nodes):
    x = np.asarray(x)
    nodes = np.asarray(nodes)
    n_nodes = len(nodes)
    N = np.zeros((n_nodes, len(x)))

    for i in range(n_nodes):
        num = np.ones_like(x)
        den = 1.0
        for j in range(n_nodes):
            if j != i:
                num *= (x - nodes[j])
                den *= (nodes[i] - nodes[j])
        N[i, :] = num / den
    return N

# ------------------------------------------------------------
# Helper: create a PolyData line from (x,y) points
# ------------------------------------------------------------
def make_line(x, y):
    points = np.column_stack((x, y, np.zeros_like(x)))
    line = pv.PolyData(points)
    line.lines = np.array([len(x), *range(len(x))])
    return line

# ------------------------------------------------------------
# Helper: add a dashed line (with fallback to solid if not supported)
# ------------------------------------------------------------
def add_dashed_line(plotter, x, y, color='black', line_width=2, label=None):
    line = make_line(x, y)
    actor = plotter.add_mesh(line, color=color, line_width=line_width,
                             label=label, style='surface')
    try:
        prop = actor.GetProperty()
        # Set stipple pattern (dash) – this alone enables stippling in most VTK versions
        prop.SetLineStipplePattern(0xF0F0)
        # Optionally set factor (default 1) – skip if not available
        # prop.SetLineStippleFactor(1)   # may not exist in all VTK builds
    except AttributeError:
        # Fallback: solid line – no stippling
        warnings.warn("Line stippling not supported; using solid line.", UserWarning)
    return actor

# ------------------------------------------------------------
# Plot shape functions
# ------------------------------------------------------------
def plot_shape_functions(p, x_plot=None):
    if x_plot is None:
        x_plot = np.linspace(-1, 1, 100)

    nodes = np.linspace(-1, 1, p + 1)
    N = lagrange_shape_functions(x_plot, nodes)

    plotter = pv.Plotter(window_size=(800, 500))
    plotter.add_title(f"1D Lagrange shape functions (p={p})")

    # Plot each shape function – PyVista picks distinct colours automatically
    for i in range(p + 1):
        line = make_line(x_plot, N[i, :])
        plotter.add_mesh(line, line_width=4, label=f'N_{i}')

    # Partition of unity – dashed black line (falls back to solid if unsupported)
    add_dashed_line(plotter, x_plot, np.sum(N, axis=0),
                    color='black', line_width=2,
                    label='Sum (partition of unity)')

    plotter.add_legend()
    plotter.view_xy()
    plotter.show()

# ------------------------------------------------------------
# Interpolation example
# ------------------------------------------------------------
def interpolate_example(p, f, x_plot=None):
    if x_plot is None:
        x_plot = np.linspace(-1, 1, 200)

    nodes = np.linspace(-1, 1, p + 1)
    f_nodes = f(nodes)
    N = lagrange_shape_functions(x_plot, nodes)
    f_interp = np.dot(f_nodes, N)

    plotter = pv.Plotter(window_size=(800, 500))
    plotter.add_title(f"Interpolation with p={p}")

    # Exact function (solid red)
    exact_line = make_line(x_plot, f(x_plot))
    plotter.add_mesh(exact_line, color='red', line_width=4,
                     label='Exact f(x)')

    # Interpolant (dashed blue, with fallback)
    add_dashed_line(plotter, x_plot, f_interp,
                    color='blue', line_width=3,
                    label='Interpolant')

    # Nodes
    node_points = np.column_stack((nodes, f_nodes, np.zeros_like(nodes)))
    plotter.add_points(node_points, color='black', point_size=12,
                       render_points_as_spheres=True, label='Nodes')

    plotter.add_legend()
    plotter.view_xy()
    plotter.show()

# ------------------------------------------------------------
# Run examples
# ------------------------------------------------------------
if __name__ == "__main__":
    plot_shape_functions(1)
    plot_shape_functions(2)
    plot_shape_functions(3)

    f = lambda x: np.sin(np.pi * x)
    interpolate_example(2, f)