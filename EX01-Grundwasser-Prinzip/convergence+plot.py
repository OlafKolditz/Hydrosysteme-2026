import math
import matplotlib.pyplot as plt
import numpy as np

def main():
    # 1 Data structures
    n = 6
    solver_iterations = 100
    eps = 1e-3
    x = [20.0] * n
    x0 = [0.0] * n

    # Store history of h1, h2, h3
    history_h1 = []
    history_h2 = []
    history_h3 = []

    # 2 Parameters
    dx1, dx2, dx3, dx4, dx5, dx6 = 100.0, 1000.0, 1000.0, 100.0, 1000.0, 1000.0
    dy1, dy2, dy3, dy4, dy5, dy6 = 500.0, 500.0, 500.0, 500.0, 500.0, 500.0
    T1, T2, T3, T4, T5, T6 = 0.05, 0.05, 0.05, 0.01, 0.01, 0.01

    T12 = (dx1 + dx2) / (dx1 / T1 + dx2 / T2)
    T32 = (dx3 + dx2) / (dx3 / T3 + dx2 / T2)
    T52 = (dy5 + dy2) / (dy5 / T5 + dy2 / T2)

    # Boundary conditions
    h1, h4 = 15.0, 15.0

    QR = 1e-07 * dx2 * dy2
    QP2 = 0.1
    QP3 = 0.0
    QP5 = 0.0
    QP6 = -0.1

    # Node 2
    c121 = dy1 * T12 * 2.0 / (dx1 + dx2)
    c122 = -c121
    c323 = dy1 * T32 * 2.0 / (dx3 + dx2)
    c322 = -c323
    c525 = dx2 * T52 * 2.0 / (dy5 + dy2)
    c522 = -c525

    a21 = c121
    a22 = c122 + c322 + c522
    a23 = c323
    a25 = c525
    a20 = QR + QP2

    b21 = -a21 / a22 * h1
    b23 = -a23 / a22
    b25 = -a25 / a22
    b20 = -a20 / a22

    # Node 3
    dy23 = dy1
    T23 = (dy2 + dy3) / (dy2 / T2 + dy3 / T3)
    c232 = dy23 * T23 / (dx2 / 2.0 + dx3 / 2.0)
    c233 = -c232
    dx63 = dx3
    T63 = (dy2 + dy6) / (dy2 / T6 + dy3 / T3)
    c633 = -dx63 * T63 / (dy6 / 2.0 + dy3 / 2.0)
    c636 = -c633

    a30 = QR + QP3
    a32 = c232
    a33 = c233 + c633
    a36 = c636

    b30 = -a30 / a33
    b32 = -a32 / a33
    b36 = -a36 / a33

    # Node 5
    dx25 = dx2
    T25 = (dy2 + dy5) / (dy2 / T2 + dy5 / T5)
    c252 = dx25 * T25 / (dy2 / 2.0 + dy5 / 2.0)
    c255 = -c252
    dy45 = dy1
    T45 = (dx4 + dx5) / (dx4 / T4 + dx5 / T5)
    c454 = dy45 * T45 / (dx4 / 2.0 + dx5 / 2.0)
    c455 = -c454
    dy65 = dy1
    T65 = (dx6 + dx5) / (dx6 / T6 + dx5 / T5)
    c656 = dy65 * T65 / (dx6 / 2.0 + dx5 / 2.0)
    c655 = -c656

    a50 = QR + QP5
    a52 = c252
    a54 = c454
    a55 = c255 + c455 + c655
    a56 = c656

    b50 = -a50 / a55
    b52 = -a52 / a55
    b54 = -a54 / a55 * h4
    b56 = -a56 / a55

    # Node 6
    dx36 = dx3
    T36 = (dy3 + dy6) / (dy3 / T3 + dy6 / T6)
    c363 = dx36 * T36 / (dy3 / 2.0 + dy6 / 2.0)
    c366 = -c363
    dy56 = dy1
    T56 = (dx5 + dx6) / (dx5 / T5 + dx6 / T6)
    c565 = dy56 * T56 / (dx5 / 2.0 + dx6 / 2.0)
    c566 = -c565

    a60 = QR + QP6
    a63 = c363
    a65 = c565
    a66 = c366 + c566

    b60 = -a60 / a66
    b63 = -a63 / a66
    b65 = -a65 / a66

    # Initial guesses for fixed nodes
    x[0] = 15.0
    x[3] = 15.0

    # Gauss-Seidel iteration
    for k in range(solver_iterations):
        x0[1] = x[1]
        x0[2] = x[2]
        x0[4] = x[4]
        x0[5] = x[5]

        x[1] = b21 + b23 * x[2] + b25 * x[4] + b20
        x[2] = b32 * x[1] + b36 * x[5] + b30
        x[4] = b52 * x[1] + b54 + b56 * x[5] + b50
        x[5] = b63 * x[2] + b65 * x[4] + b60

        # Store history for nodes 1,2,3 (indices 0,1,2)
        history_h1.append(x[0])
        history_h2.append(x[1])
        history_h3.append(x[2])

        # Error calculation for convergence check
        error = (abs(x0[1] - x[1]) + abs(x0[2] - x[2]) +
                 abs(x0[4] - x[4]) + abs(x0[5] - x[5])) / 4.0

        if error < eps:
            break

    # ========== Combine all iterations in ONE graphical plot ==========
    # Node positions (x-axis)
    node_pos = [1, 2, 3]
    node_labels = ['h1', 'h2', 'h3']
    n_iter = len(history_h1)

    # Create a colormap: early = blue, late = red
    cmap = plt.cm.jet  # or 'viridis', 'coolwarm'
    norm = plt.Normalize(vmin=0, vmax=n_iter-1)

    plt.figure(figsize=(8, 6))

    # Plot each iteration as a separate line
    for i in range(n_iter):
        h_vals = [history_h1[i], history_h2[i], history_h3[i]]
        color = cmap(norm(i))
        # Use a thin line, alpha 0.6 for early iterations, higher for later
        alpha = 0.4 + 0.5 * (i / n_iter)  # becomes more opaque towards end
        plt.plot(node_pos, h_vals, color=color, alpha=alpha, linewidth=1.0)

    # Highlight the final iteration (converged solution)
    final_h = [history_h1[-1], history_h2[-1], history_h3[-1]]
    plt.plot(node_pos, final_h, 'k-', linewidth=3, label='Final solution', zorder=10)

    # Add a colorbar to show iteration number
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=plt.gca())
    cbar.set_label('Iteration number')

    plt.xticks(node_pos, node_labels)
    plt.xlabel('Node')
    plt.ylabel('Hydraulic head h (m)')
    plt.title('Evolution of the hydraulic head profile over all iterations\n(blue → red = early → late)')
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()

    plt.savefig("combined_iterations_profile.png", dpi=150)
    plt.show()

if __name__ == "__main__":
    main()