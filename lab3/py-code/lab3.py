import numpy as np


# =====
# Data
# =====
x_data = np.array([-1.000, -0.010, 0.750, 1.900, 2.815, 3.100, 3.419, 3.875, 4.200])
y_data = np.array([1.3680, 1.9900, 3.1170, 7.6860, 17.6920, 23.1980, 31.5390, 49.1830, 67.6860])
n = len(x_data) - 1  # number of intervals


# =====
# 1. Cubic Spline (Natural / Free-end conditions)
# =====

def build_cubic_spline_natural(x, y):
    """
    Build natural cubic spline (S''(x0) = S''(xn) = 0).
    Returns the slopes m_i at each node.
    """
    n = len(x) - 1
    h = np.diff(x)  # h[i] = x[i+1] - x[i]

    # Build tridiagonal system A * m = d
    # Size: (n+1) x (n+1)
    A = np.zeros((n + 1, n + 1))
    d = np.zeros(n + 1)

    # Interior equations: mu_i * m_{i-1} + 2*m_i + lambda_i * m_{i+1} = d_i
    for i in range(1, n):
        mu_i = h[i - 1] / (h[i - 1] + h[i])
        lambda_i = h[i] / (h[i - 1] + h[i])
        d_i = 3 * (mu_i * (y[i] - y[i - 1]) / h[i - 1] + lambda_i * (y[i + 1] - y[i]) / h[i])

        A[i, i - 1] = mu_i
        A[i, i] = 2.0
        A[i, i + 1] = lambda_i
        d[i] = d_i

    # Natural boundary conditions: S''(x0) = 0, S''(xn) = 0
    # 2*m_0 + m_1 = 3*(y1 - y0)/h0
    # m_{n-1} + 2*m_n = 3*(yn - y_{n-1})/h_{n-1}
    A[0, 0] = 2.0
    A[0, 1] = 1.0
    d[0] = 3.0 * (y[1] - y[0]) / h[0]

    A[n, n - 1] = 1.0
    A[n, n] = 2.0
    d[n] = 3.0 * (y[n] - y[n - 1]) / h[n - 1]

    m = np.linalg.solve(A, d)
    return m


def eval_cubic_spline(x, y, m, x_eval):
    """
    Evaluate cubic spline at given points using Hermite interpolation formula.
    """
    n = len(x) - 1
    result = np.zeros_like(x_eval, dtype=float)

    for j, xv in enumerate(x_eval):
        # Find the interval [x_i, x_{i+1}] containing xv
        if xv <= x[0]:
            i = 0
        elif xv >= x[-1]:
            i = n - 1
        else:
            i = np.searchsorted(x, xv) - 1
            if i < 0:
                i = 0
            if i >= n:
                i = n - 1

        h_i = x[i + 1] - x[i]
        t = (xv - x[i]) / h_i

        # Hermite basis functions
        phi1 = (1 - t) ** 2 * (1 + 2 * t)
        phi2 = t ** 2 * (3 - 2 * t)
        psi1 = t * (1 - t) ** 2
        psi2 = t ** 2 * (t - 1)

        result[j] = y[i] * phi1 + y[i + 1] * phi2 + h_i * (m[i] * psi1 + m[i + 1] * psi2)

    return result


# =====
# 2. Quadratic Spline (with midpoints as join points, natural conditions)
# =====

def build_quadratic_spline(x, y):
    """
    Build quadratic spline with join points at midpoints.
    Each quadratic: S_i(t) = a_i * t^2 + b_i * t + c_i on [xtilde_{i-1}, xtilde_i]
    
    Join points: xtilde_i = (x[i] + x[i+1]) / 2
    Boundary conditions: S'_0(x0) = 0, S'_n(xn) = 0 (natural-like)
    
    Returns list of (a_i, b_i, c_i, left_i, right_i) for each segment.
    """
    n = len(x) - 1
    
    # Join points (midpoints)
    xtilde = np.array([(x[i] + x[i + 1]) / 2 for i in range(n)])
    
    # Extended intervals: [x0, xtilde0], [xtilde0, xtilde1], ..., [xtilde_{n-1}, xn]
    # Total: n+1 quadratic polynomials
    left_bounds = np.concatenate([[x[0]], xtilde])
    right_bounds = np.concatenate([xtilde, [x[-1]]])
    
    # Each quadratic: S_i(t) = a_i*t^2 + b_i*t + c_i
    # Unknowns: 3*(n+1) coefficients
    num_segments = n + 1
    num_unknowns = 3 * num_segments
    
    A = np.zeros((num_unknowns, num_unknowns))
    b_vec = np.zeros(num_unknowns)
    
    row = 0
    
    # 1. Interpolation conditions: S_i(x_i) = y_i
    for i in range(num_segments):
        xi = x[i]
        # S_i(xi) = a_i*xi^2 + b_i*xi + c_i = y_i
        A[row, 3 * i] = xi ** 2
        A[row, 3 * i + 1] = xi
        A[row, 3 * i + 2] = 1.0
        b_vec[row] = y[i]
        row += 1
    
    # 2. Continuity at join points: S_{i-1}(xtilde_i) = S_i(xtilde_i)
    for i in range(1, num_segments):
        xt = xtilde[i - 1]  # join point between segment i-1 and i
        # S_{i-1}(xt) - S_i(xt) = 0
        A[row, 3 * (i - 1)] = xt ** 2
        A[row, 3 * (i - 1) + 1] = xt
        A[row, 3 * (i - 1) + 2] = 1.0
        A[row, 3 * i] = -(xt ** 2)
        A[row, 3 * i + 1] = -xt
        A[row, 3 * i + 2] = -1.0
        b_vec[row] = 0.0
        row += 1
    
    # 3. Smoothness (C1 continuity) at join points: S'_{i-1}(xtilde_i) = S'_i(xtilde_i)
    for i in range(1, num_segments):
        xt = xtilde[i - 1]
        # S'_{i-1}(xt) - S'_i(xt) = 0
        # S'_i(t) = 2*a_i*t + b_i
        A[row, 3 * (i - 1)] = 2 * xt
        A[row, 3 * (i - 1) + 1] = 1.0
        A[row, 3 * i] = -2 * xt
        A[row, 3 * i + 1] = -1.0
        b_vec[row] = 0.0
        row += 1
    
    # 4. Boundary conditions (natural-like): S'_0(x0) = 0, S'_n(xn) = 0
    # S'_0(x0) = 2*a_0*x0 + b_0 = 0
    A[row, 0] = 2 * x[0]
    A[row, 1] = 1.0
    b_vec[row] = 0.0
    row += 1
    
    # S'_n(xn) = 2*a_n*xn + b_n = 0
    A[row, 3 * n] = 2 * x[-1]
    A[row, 3 * n + 1] = 1.0
    b_vec[row] = 0.0
    row += 1
    
    # Solve the system
    coeffs = np.linalg.solve(A, b_vec)
    
    # Extract polynomials
    segments = []
    for i in range(num_segments):
        a_i = coeffs[3 * i]
        b_i = coeffs[3 * i + 1]
        c_i = coeffs[3 * i + 2]
        segments.append((a_i, b_i, c_i, left_bounds[i], right_bounds[i]))
    
    return segments


def eval_quadratic_spline(segments, x_eval):
    """
    Evaluate quadratic spline at given points.
    """
    result = np.zeros_like(x_eval, dtype=float)
    
    for j, xv in enumerate(x_eval):
        for a_i, b_i, c_i, left, right in segments:
            if left <= xv <= right:
                result[j] = a_i * xv ** 2 + b_i * xv + c_i
                break
        else:
            # Handle edge cases with small tolerance
            for a_i, b_i, c_i, left, right in segments:
                if abs(xv - left) < 1e-10 or abs(xv - right) < 1e-10:
                    result[j] = a_i * xv ** 2 + b_i * xv + c_i
                    break
    
    return result


# ====
# Main
# =====

if __name__ == "__main__":
    print("=" * 70)
    print("LAB 3: SPLINE INTERPOLATION")
    print("=" * 70)
    print("\nData points:")
    print(f"  xi = {x_data}")
    print(f"  yi = {y_data}")
    print(f"  Number of intervals: {n}")

    # -----
    # Part 1: Cubic Spline (Natural)
    # -----
    print("\n" + "-" * 70)
    print("PART 1: CUBIC SPLINE (Natural / Free-end)")
    print("-" * 70)

    m = build_cubic_spline_natural(x_data, y_data)
    print("\nComputed slopes (m_i) at nodes:")
    for i, mi in enumerate(m):
        print(f"  m[{i}] = {mi:.6f}")

    # Evaluate at various points
    x_eval = np.linspace(x_data[0], x_data[-1], 20)
    y_eval = eval_cubic_spline(x_data, y_data, m, x_eval)

    print("\nCubic spline evaluation:")
    print(f"  {'x':>10} {'S(x)':>12}")
    for xv, yv in zip(x_eval, y_eval):
        print(f"  {xv:10.4f} {yv:12.6f}")

    # Verification: check interpolation at original nodes
    y_check = eval_cubic_spline(x_data, y_data, m, x_data)
    print("\nVerification (interpolation at original nodes):")
    print(f"  {'x':>10} {'y':>12} {'S(x)':>12} {'Error':>12}")
    for xv, yv, sv in zip(x_data, y_data, y_check):
        print(f"  {xv:10.4f} {yv:12.6f} {sv:12.6f} {abs(yv - sv):12.2e}")

    # -----
    # Part 2: Quadratic Spline
    # -----
    print("\n" + "-" * 70)
    print("PART 2: QUADRATIC SPLINE (midpoint join points)")
    print("-" * 70)

    segments = build_quadratic_spline(x_data, y_data)

    print("\nQuadratic polynomials S_i(x) = a*x^2 + b*x + c:")
    for i, (a_i, b_i, c_i, left, right) in enumerate(segments):
        print(f"  Segment {i}: [{left:.4f}, {right:.4f}]")
        print(f"    a = {a_i:.6f}, b = {b_i:.6f}, c = {c_i:.6f}")

    # Evaluate at various points
    y_quad_eval = eval_quadratic_spline(segments, x_eval)

    print("\nQuadratic spline evaluation:")
    print(f"  {'x':>10} {'S(x)':>12}")
    for xv, yv in zip(x_eval, y_quad_eval):
        print(f"  {xv:10.4f} {yv:12.6f}")

    # Verification: check interpolation at original nodes
    y_quad_check = eval_quadratic_spline(segments, x_data)
    print("\nVerification (interpolation at original nodes):")
    print(f"  {'x':>10} {'y':>12} {'S(x)':>12} {'Error':>12}")
    for xv, yv, sv in zip(x_data, y_data, y_quad_check):
        print(f"  {xv:10.4f} {yv:12.6f} {sv:12.6f} {abs(yv - sv):12.2e}")

    print("\n" + "=" * 70)
    print("DONE")
    print("=" * 70)
