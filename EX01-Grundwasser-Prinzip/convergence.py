import math
import matplotlib.pyplot as plt

def main():
    # 1 Data structures
    n = 6
    solver_iterations = 100
    eps = 1e-3
    x = [10000.0] * n
    x0 = [0.0] * n
    y = 0.0

    # Store history of h1, h2, h3 and the error
    history_h1 = []
    history_h2 = []
    history_h3 = []
    history_error = []

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
    QP2 = -0.1
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

        # Error calculation
        error = (abs(x0[1] - x[1]) + abs(x0[2] - x[2]) +
                 abs(x0[4] - x[4]) + abs(x0[5] - x[5])) / 4.0
        history_error.append(error)

        if error < eps:
            break

    # ========== Graphical plot of time steps ==========
    iterations = range(len(history_h1))

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 10), sharex=True)

    # Upper plot: h1, h2, h3 vs iteration
    ax1.plot(iterations, history_h1, 'o-', label='h1', linewidth=2, markersize=4)
    ax1.plot(iterations, history_h2, 's-', label='h2', linewidth=2, markersize=4)
    ax1.plot(iterations, history_h3, '^-', label='h3', linewidth=2, markersize=4)
    ax1.set_ylabel('Hydraulic head h (m)')
    ax1.set_title('Convergence of hydraulic heads over iterations')
    ax1.grid(True)
    ax1.legend()

    # Lower plot: error vs iteration
    ax2.semilogy(iterations, history_error, 'r-', linewidth=2, marker='.', markersize=4)
    ax2.set_xlabel('Iteration number')
    ax2.set_ylabel('Error (norm of change)')
    ax2.set_title('Error reduction (semi‑log scale)')
    ax2.grid(True)

    plt.tight_layout()
    plt.savefig("convergence_graphical.png", dpi=150)
    plt.show()

if __name__ == "__main__":
    main()