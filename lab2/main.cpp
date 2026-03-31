#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>

// Solve linear system Ax = b using Gaussian elimination with partial pivoting
std::vector<double> solveLinearSystem(std::vector<std::vector<double>> A, std::vector<double> b) {
    int n = b.size();
    
    // Forward elimination with partial pivoting
    for (int k = 0; k < n; k++) {
        // Find pivot
        int maxRow = k;
        for (int i = k + 1; i < n; i++) {
            if (std::abs(A[i][k]) > std::abs(A[maxRow][k])) {
                maxRow = i;
            }
        }
        
        // Swap rows
        std::swap(A[k], A[maxRow]);
        std::swap(b[k], b[maxRow]);
        
        // Eliminate column
        for (int i = k + 1; i < n; i++) {
            double factor = A[i][k] / A[k][k];
            for (int j = k; j < n; j++) {
                A[i][j] -= factor * A[k][j];
            }
            b[i] -= factor * b[k];
        }
    }
    
    // Back substitution
    std::vector<double> x(n);
    for (int i = n - 1; i >= 0; i--) {
        x[i] = b[i];
        for (int j = i + 1; j < n; j++) {
            x[i] -= A[i][j] * x[j];
        }
        x[i] /= A[i][i];
    }
    
    return x;
}

// Evaluate polynomial at point x
double polyval(const std::vector<double>& coeffs, double x) {
    double result = 0.0;
    for (int i = coeffs.size() - 1; i >= 0; i--) {
        result = result * x + coeffs[i];
    }
    return result;
}

// Polynomial approximation using least squares
std::pair<std::vector<double>, double> approximation(
    const std::vector<double>& x, 
    const std::vector<double>& y, 
    int degree) {
    
    int n = x.size();
    int m = degree + 1;
    
    // Build Vandermonde matrix (increasing powers)
    std::vector<std::vector<double>> A(n, std::vector<double>(m));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            A[i][j] = std::pow(x[i], j);
        }
    }
    
    // Solve normal equations: (A^T * A) * coeffs = A^T * y
    std::vector<std::vector<double>> AtA(m, std::vector<double>(m, 0.0));
    std::vector<double> Aty(m, 0.0);
    
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < m; j++) {
            for (int k = 0; k < n; k++) {
                AtA[i][j] += A[k][i] * A[k][j];
            }
        }
        for (int k = 0; k < n; k++) {
            Aty[i] += A[k][i] * y[k];
        }
    }
    
    // Solve for coefficients
    std::vector<double> coeffs = solveLinearSystem(AtA, Aty);
    
    // Calculate fitted values and MSE
    double mse = 0.0;
    for (int i = 0; i < n; i++) {
        double y_fit = polyval(coeffs, x[i]);
        mse += std::pow(y[i] - y_fit, 2);
    }
    mse /= n;
    
    double rmse = std::sqrt(mse);
    
    return {coeffs, rmse};
}

int main() {
    std::vector<double> x = {
        0.1, 0.3, 0.5, 0.7, 0.9, 1.1, 1.3, 1.5, 1.7, 1.9,
        2.1, 2.3, 2.5, 2.7, 2.9, 3.1, 3.3, 3.5, 3.7, 3.9
    };

    std::vector<double> y = {
        -0.86, -0.77, -0.56, -0.46, -0.28, -0.24, -0.36, -0.43, -0.56, -0.59,
        -0.70, -1.01, -1.03, -1.47, -1.68, -1.93, -2.28, -2.53, -2.93, -3.07
    };
    
    // Perform approximations
    auto [coeffs_1, rmse_1] = approximation(x, y, 1);
    auto [coeffs_2, rmse_2] = approximation(x, y, 2);
    
    // Output results
    std::cout << std::fixed << std::setprecision(2);

    std::cout << "\nFirst-degree polynomial (linear approximation):" << std::endl;
    std::cout << "P_1(x) = " << coeffs_1[0] << " + " << coeffs_1[1] << "*x" << std::endl;
    std::cout << "Root mean square error: " << rmse_1 << std::endl;

    std::cout << "\nSecond-degree polynomial (quadratic approximation):" << std::endl;
    std::cout << "P_2(x) = " << coeffs_2[0] << " + " << coeffs_2[1] << "*x + "
              << coeffs_2[2] << "*x^2" << std::endl;
    std::cout << "Root mean square error: " << rmse_2 << std::endl;
    
    return 0;
}
