import math

def main():
    # 1 Data structures
    n = 6
    solver_iterations = 100
    eps = 1e-3
    x = [10000.0] * n
    x0 = [0.0] * n
    y = 0.0

    out_file1 = open("out1.txt", "w")
    out_file2 = open("out2.txt", "w")

    # 2 Parameters
    dx1, dx2, dx3, dx4, dx5, dx6 = 100.0, 1000.0, 1000.0, 100.0, 1000.0, 1000.0
    dy1, dy2, dy3, dy4, dy5, dy6 = 500.0, 500.0, 500.0, 500.0, 500.0, 500.0
    T1, T2, T3, T4, T5, T6 = 0.05, 0.05, 0.05, 0.01, 0.01, 0.01

    T12 = (dx1 + dx2) / (dx1 / T1 + dx2 / T2)
    print(f"T12: {T12}")
    T32 = (dx3 + dx2) / (dx3 / T3 + dx2 / T2)
    print(f"T32: {T32}")
    T52 = (dy5 + dy2) / (dy5 / T5 + dy2 / T2)
    print(f"T52: {T52}")

    # Boundary conditions
    h1, h4 = 15.0, 15.0

    QR = 1e-07 * dx2 * dy2
    print(f"QR: {QR}")
    QP2 = -0.1
    QP3 = 0.0
    QP5 = 0.0
    QP6 = -0.1

    # Node 2
    print("Knoten 2:")
    c121 = dy1 * T12 * 2.0 / (dx1 + dx2)
    print(f"c121: {c121}")
    c122 = -c121
    print(f"c122: {c122}")

    c323 = dy1 * T32 * 2.0 / (dx3 + dx2)
    print(f"c323: {c323}")
    c322 = -c323
    print(f"c322: {c322}")

    c525 = dx2 * T52 * 2.0 / (dy5 + dy2)
    print(f"c525: {c525}")
    c522 = -c525
    print(f"c522: {c522}")

    a21 = c121
    a22 = c122 + c322 + c522
    a23 = c323
    a25 = c525
    a20 = QR + QP2

    b21 = -a21 / a22 * h1
    print(f"b21: {b21}")
    b23 = -a23 / a22
    print(f"b23: {b23}")
    b25 = -a25 / a22
    print(f"b25: {b25}")
    b20 = -a20 / a22
    print(f"b20: {b20}")

    # Node 3
    print("Knoten 3:")
    dy23 = dy1
    T23 = (dy2 + dy3) / (dy2 / T2 + dy3 / T3)
    c232 = dy23 * T23 / (dx2 / 2.0 + dx3 / 2.0)
    print(f"c232: {c232}")
    c233 = -c232
    print(f"c233: {c233}")

    dx63 = dx3
    T63 = (dy2 + dy6) / (dy2 / T6 + dy3 / T3)
    c633 = -dx63 * T63 / (dy6 / 2.0 + dy3 / 2.0)
    print(f"c633: {c633}")
    c636 = -c633

    a30 = QR + QP3
    a32 = c232
    a33 = c233 + c633
    a36 = c636

    b30 = -a30 / a33
    print(f"b30: {b30}")
    b32 = -a32 / a33
    print(f"b32: {b32}")
    b36 = -a36 / a33
    print(f"b36: {b36}")

    # Node 5
    print("Knoten 5:")
    dx25 = dx2
    T25 = (dy2 + dy5) / (dy2 / T2 + dy5 / T5)
    c252 = dx25 * T25 / (dy2 / 2.0 + dy5 / 2.0)
    print(f"c252: {c252}")
    c255 = -c252
    print(f"c255: {c255}")

    dy45 = dy1
    T45 = (dx4 + dx5) / (dx4 / T4 + dx5 / T5)
    c454 = dy45 * T45 / (dx4 / 2.0 + dx5 / 2.0)
    print(f"c454: {c454}")
    c455 = -c454
    print(f"c455: {c455}")

    dy65 = dy1
    T65 = (dx6 + dx5) / (dx6 / T6 + dx5 / T5)
    c656 = dy65 * T65 / (dx6 / 2.0 + dx5 / 2.0)
    print(f"c656: {c656}")
    c655 = -c656
    print(f"c655: {c655}")

    a50 = QR + QP5
    a52 = c252
    a54 = c454
    a55 = c255 + c455 + c655
    a56 = c656

    b50 = -a50 / a55
    print(f"b50: {b50}")
    b52 = -a52 / a55
    print(f"b52: {b52}")
    b54 = -a54 / a55 * h4
    print(f"b54: {b54}")
    b56 = -a56 / a55
    print(f"b56: {b56}")

    # Node 6
    print("Knoten 6:")
    dx36 = dx3
    T36 = (dy3 + dy6) / (dy3 / T3 + dy6 / T6)
    c363 = dx36 * T36 / (dy3 / 2.0 + dy6 / 2.0)
    print(f"c363: {c363}")
    c366 = -c363
    print(f"c366: {c366}")

    dy56 = dy1
    T56 = (dx5 + dx6) / (dx5 / T5 + dx6 / T6)
    c565 = dy56 * T56 / (dx5 / 2.0 + dx6 / 2.0)
    print(f"c565: {c565}")
    c566 = -c565

    a60 = QR + QP6
    a63 = c363
    a65 = c565
    a66 = c366 + c566

    b60 = -a60 / a66
    print(f"b60: {b60}")
    b63 = -a63 / a66
    print(f"b63: {b63}")
    b65 = -a65 / a66
    print(f"b65: {b65}")

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

        # Write to files
        for i in range(3):
            out_file1.write(f"h{i+1}, {x[i]}\n")
        for i in range(3, n):
            out_file2.write(f"h{i+1}, {x[i]}\n")

        # Error calculation
        y = (abs(x0[1] - x[1]) + abs(x0[2] - x[2]) +
             abs(x0[4] - x[4]) + abs(x0[5] - x[5])) / 4.0

        if y < eps:
            break

    out_file1.close()
    out_file2.close()

if __name__ == "__main__":
    main()