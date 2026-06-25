import math
import sys

def gauss(matrix, vecb, vecx, g):
    """
    Gaussian elimination with partial pivoting.
    Modifies matrix and vecb in place, stores solution in vecx.
    matrix is a 2D list (list of lists) of size g x g.
    vecb and vecx are lists of length g.
    """
    s = [0] * (g - 1)          # pivot row indices

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
        self.node_vector = []          # list of nodal coordinates (x)
        self.element_vector = []       # list of [node1, node2] (indices into node_vector)
        self.nn = 0
        self.ne = 0

        self.h_initial = 5.0
        self.h_top = 12.0
        self.h_bottom = 0.0
        self.K = [1e-5, 2e-5, 1e-5, 1e-5]   # will be used according to element index

        self.u = []          # current solution (list of length nn)
        self.u_new = []      # new solution
        self.bc_nodes = []   # list of node indices with Dirichlet BCs

        self.matrix = []     # global stiffness matrix (2D list)
        self.vecb = []       # RHS vector
        self.vecx = []       # solution vector

        self.element_matrix_vector = []  # list of element stiffness matrices (each a 4‑element list)
        self.out_file = open("out.txt", "w")

    def read_mesh(self, filename):
        """Reads a mesh file in the format used by the C++ code."""
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
                # next line contains number of nodes
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
        """Sets u and u_new to h_initial for all nodes."""
        self.u = [self.h_initial] * self.nn
        self.u_new = [self.h_initial] * self.nn

    def set_boundary_conditions(self):
        """Sets Dirichlet BCs at nodes 0 and 4 (last node)."""
        self.bc_nodes = [0, self.nn - 1]   # assuming last node index = nn-1
        self.u[0] = self.h_top
        self.u_new[0] = self.h_top
        self.u[self.nn - 1] = self.h_bottom
        self.u_new[self.nn - 1] = self.h_bottom

    def calculate_element_matrices(self):
        """Computes stiffness matrix for each element."""
        self.element_matrix_vector = []
        for e in range(self.ne):
            self.calculate_element_matrix(e)

    def calculate_element_matrix(self, e):
        """Computes element stiffness matrix for element e and stores it."""
        nodes = self.element_vector[e]
        x0 = self.node_vector[nodes[0]]
        x1 = self.node_vector[nodes[1]]
        L = x1 - x0
        k = self.K[e]   # use K[e] (ensure K has enough entries)
        # Element matrix in row‑major order: [k/L, -k/L; -k/L, k/L]
        elem_mat = [k / L, -k / L, -k / L, k / L]
        self.element_matrix_vector.append(elem_mat)

    def assemble_equation_system(self):
        """Assembles the global stiffness matrix and RHS vector."""
        # Initialize matrix and RHS
        self.matrix = [[0.0] * self.nn for _ in range(self.nn)]
        self.vecb = [0.0] * self.nn
        self.vecx = [0.0] * self.nn

        # Loop over elements
        for e in range(self.ne):
            nodes = self.element_vector[e]
            elem_mat = self.element_matrix_vector[e]
            i = nodes[0]
            j = nodes[1]
            # K_ii += elem_mat[0]
            self.matrix[i][i] += elem_mat[0]
            # K_ij += elem_mat[1]
            self.matrix[i][j] += elem_mat[1]
            # K_ji += elem_mat[2]
            self.matrix[j][i] += elem_mat[2]
            # K_jj += elem_mat[3]
            self.matrix[j][j] += elem_mat[3]

    def incorporate_boundary_conditions(self):
        """Modifies matrix and RHS to enforce Dirichlet BCs."""
        for bc_node in self.bc_nodes:
            i_row = bc_node
            # Zero out off‑diagonal entries in row and column, update RHS
            for j in range(self.nn):
                if i_row == j:
                    continue
                # Row
                self.matrix[i_row][j] = 0.0
                # Column
                val = self.matrix[j][i_row]
                # Apply contribution to RHS (vecb[j] -= matrix[j][i_row] * u[i_row])
                self.vecb[j] -= val * self.u[i_row]
                self.matrix[j][i_row] = 0.0
            # Set diagonal and RHS according to BC
            self.vecb[i_row] = self.u[i_row] * self.matrix[i_row][i_row]

    def save_time_step(self):
        """Copies vecx into u_new and u."""
        for n in range(self.nn):
            self.u_new[n] = self.vecx[n]
            self.u[n] = self.u_new[n]

    def output_results(self, t):
        """Writes the current solution to out.txt (appends)."""
        # In the C++ code, they write every time step (here only one step).
        for val in self.u_new:
            self.out_file.write(f"{val}\n")

    def dump_equation_system(self):
        """Writes the matrix and RHS to out.txt for debugging."""
        self.out_file.write("Matrix and RHS:\n")
        for i in range(self.nn):
            row_str = "\t".join(f"{self.matrix[i][j]:.3e}" for j in range(self.nn))
            self.out_file.write(f"{row_str}\tb: {self.vecb[i]:.3e}\n")


def main():
    # File handling
    aux_file = open("cputime.txt", "w")
    msh_file = "column.msh"          # input mesh file
    msh_file_test = "column_test.msh"
    matrix_file_test = "element_matrices.txt"

    # Create FEM object
    fem = FEM()

    # Read mesh
    try:
        fem.read_mesh(msh_file)
    except FileNotFoundError:
        print(f"Error: file '{msh_file}' not found.")
        return

    fem.output_mesh(msh_file_test)

    fem.set_initial_conditions()
    fem.set_boundary_conditions()

    # Time loop (only one step as in original C++ main)
    tn = 1
    for t in range(tn):
        fem.calculate_element_matrices()
        # Dump element matrices to file
        with open(matrix_file_test, 'w') as f:
            for e, em in enumerate(fem.element_matrix_vector):
                f.write("--------------------------\n")
                f.write(f"{em[0]:.3e} {em[1]:.3e}\n")
                f.write(f"{em[2]:.3e} {em[3]:.3e}\n")
        fem.assemble_equation_system()
        fem.incorporate_boundary_conditions()
        # Optionally dump system (commented out to avoid clutter)
        # fem.dump_equation_system()
        # Solve
        gauss(fem.matrix, fem.vecb, fem.vecx, fem.nn)
        fem.save_time_step()
        fem.output_results(t)

    aux_file.write(f"CPU time: {0.0}\n")   # not measuring in Python
    aux_file.close()
    fem.out_file.close()

    print("Solution written to out.txt")
    # Optional: plot results (uncomment if matplotlib is available)
    # import matplotlib.pyplot as plt
    # x = fem.node_vector
    # u = fem.u_new
    # plt.plot(x, u, 'o-')
    # plt.xlabel('x')
    # plt.ylabel('u')
    # plt.title('Steady-state solution')
    # plt.show()


if __name__ == "__main__":
    main()