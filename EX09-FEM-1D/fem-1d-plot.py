import math
import sys

try:
    import matplotlib.pyplot as plt
    HAS_PLOT = True
except ImportError:
    HAS_PLOT = False


def gauss(matrix, vecb, vecx, g):
    """
    Gaussian elimination with partial pivoting.
    Modifies matrix and vecb in place, stores solution in vecx.
    matrix is a 2D list (list of lists) of size g x g.
    vecb and vecx are lists of length g.
    """
    s = [0] * (g - 1)

    # LU factorisation
    for k in range(g - 1):
        # Pivot search
        z = 0.0
        sk = k
        for i in range(k, g):
            hilf = abs(matrix[i][k])
            if hilf > z:
                z = hilf
                sk = i
        s[k] = sk

        # Swap rows if needed
        if sk > k:
            matrix[k], matrix[sk] = matrix[sk], matrix[k]

        # Elimination coefficients
        for i in range(k + 1, g):
            matrix[i][k] /= matrix[k][k]

        # Update the remaining submatrix
        for j in range(k + 1, g):
            for i in range(k + 1, g):
                matrix[i][j] -= matrix[i][k] * matrix[k][j]

    # Transform right‑hand side according to pivot swaps
    for k in range(g - 1):
        sk = s[k]
        if sk > k:
            vecb[k], vecb[sk] = vecb[sk], vecb[k]

    # Forward substitution (L y = b)
    for k in range(1, g):
        for j in range(k):
            vecb[k] -= matrix[k][j] * vecb[j]

    # Back substitution (U x = y)
    for k in range(g - 1, -1, -1):
        for j in range(k + 1, g):
            vecb[k] -= matrix[k][j] * vecb[j]
        vecb[k] /= matrix[k][k]

    # Copy result to vecx
    for k in range(g):
        vecx[k] = vecb[k]


class FEM:
    def __init__(self):
        self.node_vector = []          # nodal coordinates (x)
        self.element_vector = []       # list of [node1, node2] (indices)
        self.nn = 0
        self.ne = 0

        self.h_initial = 5.0
        self.h_top = 12.0
        self.h_bottom = 0.0
        self.K = [1e-5, 2e-5, 1e-5, 1e-5]   # element permeabilities

        self.u = []          # current solution
        self.u_new = []      # new solution
        self.bc_nodes = []   # Dirichlet BC node indices

        self.matrix = []     # global stiffness matrix
        self.vecb = []       # RHS
        self.vecx = []       # solution vector

        self.element_matrix_vector = []  # element stiffness matrices
        self.out_file = open("out.txt", "w")

    def read_mesh(self, filename):
        """Reads a mesh file with $NODES and $ELEMENTS sections."""
        with open(filename, 'r') as f:
            lines = f.readlines()

        i = 0
        while i < len(lines):
            line = lines[i].strip()
            if not line:
                i += 1
                continue
            if line.startswith("#STOP"):
                break
            if line.startswith("$NODES"):
                i += 1
                self.nn = int(lines[i].strip())
                i += 1
                self.node_vector = []
                for _ in range(self.nn):
                    x = float(lines[i].strip())
                    self.node_vector.append(x)
                    i += 1
            elif line.startswith("$ELEMENTS"):
                i += 1
                self.ne = int(lines[i].strip())
                i += 1
                self.element_vector = []
                for _ in range(self.ne):
                    parts = lines[i].strip().split()
                    n0 = int(parts[0])
                    n1 = int(parts[1])
                    self.element_vector.append([n0, n1])
                    i += 1
            else:
                i += 1

    def output_mesh(self, filename):
        """Writes the mesh to a file for verification."""
        with open(filename, 'w') as f:
            f.write("#FEM_MSH\n")
            f.write("$NODES\n")
            for x in self.node_vector:
                f.write(f"{x}\n")
            f.write("$ELEMENTS\n")
            for elem in self.element_vector:
                f.write(f"{elem[0]} {elem[1]}\n")
            f.write("#STOP\n")

    def set_initial_conditions(self):
        self.u = [self.h_initial] * self.nn
        self.u_new = [self.h_initial] * self.nn

    def set_boundary_conditions(self):
        self.bc_nodes = [0, self.nn - 1]
        self.u[0] = self.h_top
        self.u_new[0] = self.h_top
        self.u[self.nn - 1] = self.h_bottom
        self.u_new[self.nn - 1] = self.h_bottom

    def calculate_element_matrices(self):
        self.element_matrix_vector = []
        for e in range(self.ne):
            self.calculate_element_matrix(e)

    def calculate_element_matrix(self, e):
        nodes = self.element_vector[e]
        x0 = self.node_vector[nodes[0]]
        x1 = self.node_vector[nodes[1]]
        L = x1 - x0
        k = self.K[e]
        elem_mat = [k / L, -k / L, -k / L, k / L]
        self.element_matrix_vector.append(elem_mat)

    def assemble_equation_system(self):
        self.matrix = [[0.0] * self.nn for _ in range(self.nn)]
        self.vecb = [0.0] * self.nn
        self.vecx = [0.0] * self.nn

        for e in range(self.ne):
            nodes = self.element_vector[e]
            elem_mat = self.element_matrix_vector[e]
            i = nodes[0]
            j = nodes[1]
            self.matrix[i][i] += elem_mat[0]
            self.matrix[i][j] += elem_mat[1]
            self.matrix[j][i] += elem_mat[2]
            self.matrix[j][j] += elem_mat[3]

    def incorporate_boundary_conditions(self):
        for bc_node in self.bc_nodes:
            i_row = bc_node
            for j in range(self.nn):
                if i_row == j:
                    continue
                # Row
                self.matrix[i_row][j] = 0.0
                # Column
                val = self.matrix[j][i_row]
                self.vecb[j] -= val * self.u[i_row]
                self.matrix[j][i_row] = 0.0
            self.vecb[i_row] = self.u[i_row] * self.matrix[i_row][i_row]

    def save_time_step(self):
        for n in range(self.nn):
            self.u_new[n] = self.vecx[n]
            self.u[n] = self.u_new[n]

    def output_results(self, t):
        for val in self.u_new:
            self.out_file.write(f"{val}\n")


def plot_solution(fem):
    """Plot the nodal solution u vs x."""
    if not HAS_PLOT:
        print("matplotlib not installed – skipping plot.")
        return

    x = fem.node_vector
    u = fem.u_new

    plt.figure(figsize=(8, 5))
    plt.plot(x, u, 'bo-', linewidth=2, markersize=8, label='FEM solution')
    plt.xlabel('x')
    plt.ylabel('u')
    plt.title('Steady‑state solution (1D FEM)')
    plt.grid(True)

    # Mark boundary conditions
    for node in fem.bc_nodes:
        plt.axvline(x=x[node], color='gray', linestyle='--', alpha=0.5)

    plt.legend()
    plt.tight_layout()
    plt.show()


def main():
    # Input file
    msh_file = "column.msh"
    try:
        fem = FEM()
        fem.read_mesh(msh_file)
    except FileNotFoundError:
        print(f"Error: file '{msh_file}' not found.")
        print("Create a mesh file with $NODES and $ELEMENTS sections.")
        return

    # Optional output mesh for debugging
    fem.output_mesh("column_test.msh")

    fem.set_initial_conditions()
    fem.set_boundary_conditions()

    # Single time step (steady‑state)
    tn = 1
    for t in range(tn):
        fem.calculate_element_matrices()
        fem.assemble_equation_system()
        fem.incorporate_boundary_conditions()
        gauss(fem.matrix, fem.vecb, fem.vecx, fem.nn)
        fem.save_time_step()
        fem.output_results(t)

    fem.out_file.close()
    print("Solution written to out.txt")

    # Plot results
    plot_solution(fem)


if __name__ == "__main__":
    main()