#include <iostream>
#include <vector>
#include <algorithm>
#include <iomanip>
#include <cmath>

#ifndef M_PI
#define M_PI 3.14
#endif

using Point = std::pair<double, double>;

double lagrange(const std::vector<Point>& points, double x) {
    int n = points.size();
    double res = 0.0;

    for (int i = 0; i < n; ++i) {
        double xi = points[i].first;
        double yi = points[i].second;
        double t = 1.0;
        for (int j = 0; j < n; ++j) {
            if (j != i) {
                double xj = points[j].first;
                t *= (x - xj) / (xi - xj);
            }
        }
        res += yi * t;
    }

    return res;
}

double aitken(const std::vector<Point>& points, double x) {
    int n = points.size();
    std::vector<std::vector<double>> p(n, std::vector<double>(n, 0.0));

    for (int i = 0; i < n; ++i) {
        p[i][i] = points[i].second;
    }

    for (int k = 1; k < n; ++k) {
        for (int i = 0; i < n - k; ++i) {
            int j = i + k;
            double xi = points[i].first;
            double xj = points[j].first;
            p[i][j] = (p[i][j-1] * (xj - x) - p[i+1][j] * (xi - x)) / (xj - xi);
        }
    }

    return p[0][n-1];
}

double newton1(const std::vector<Point>& points, double x) {
    int n = points.size();
    std::vector<Point> sorted_points = points;
    std::sort(sorted_points.begin(), sorted_points.end());

    std::vector<double> work(n);
    for (int i = 0; i < n; ++i) {
        work[i] = sorted_points[i].second;
    }

    double res = work[0];
    double prod = 1.0;

    for (int j = 1; j < n; ++j) {
        for (int i = 0; i < n - j; ++i) {
            work[i] = (work[i+1] - work[i]) / (sorted_points[i+j].first - sorted_points[i].first);
        }
        prod *= (x - sorted_points[j-1].first);
        res += work[0] * prod;
    }

    return res;
}

double newton2(double x, double h, double a, double b) {
    auto f = [](double x) { return std::pow(x, 5.0 / M_PI); };

    int n = static_cast<int>((b - a) / h) + 1;
    std::vector<Point> points(n);
    for (int i = 0; i < n; ++i) {
        points[i] = {a + i * h, f(a + i * h)};
    }

    std::vector<std::vector<double>> diff;
    std::vector<double> y(n);
    for (int i = 0; i < n; ++i) {
        y[i] = points[i].second;
    }
    diff.push_back(y);

    for (int j = 1; j < n; ++j) {
        std::vector<double> row;
        for (int i = 0; i < n - j; ++i) {
            row.push_back(diff[j-1][i+1] - diff[j-1][i]);
        }
        diff.push_back(row);
    }

    double res;
    if (x > (a + b) / 2) {
        double t = (x - points.back().first) / h;
        res = diff[0].back();
        double term = 1;
        double fact = 1;
        for (int j = 1; j < n; ++j) {
            term *= (t + j - 1);
            fact *= j;
            res += diff[j].back() * term / fact;
        }
    } else {
        double t = (x - points[0].first) / h;
        res = diff[0][0];
        double term = 1;
        double fact = 1;
        for (int j = 1; j < n; ++j) {
            term *= (t - j + 1);
            fact *= j;
            res += diff[j][0] * term / fact;
        }
    }

    return res;
}

double bessel(double x, double h, double a, double b) {
    auto f = [](double x) { return std::pow(x, 5.0 / M_PI); };

    int n = static_cast<int>((b - a) / h) + 1;
    std::vector<Point> points(n);
    for (int i = 0; i < n; ++i) {
        points[i] = {a + i * h, f(a + i * h)};
    }

    std::vector<std::vector<double>> diff;
    std::vector<double> y(n);
    for (int i = 0; i < n; ++i) {
        y[i] = points[i].second;
    }
    diff.push_back(y);

    for (int j = 1; j < n; ++j) {
        std::vector<double> row;
        for (int i = 0; i < n - j; ++i) {
            row.push_back(diff[j-1][i+1] - diff[j-1][i]);
        }
        diff.push_back(row);
    }

    int k = 0;
    for (int i = 0; i < n - 1; ++i) {
        if (points[i].first <= x && x < points[i+1].first) {
            k = i;
            break;
        }
    }

    double x0 = points[k].first;
    double q = (x - x0) / h;

    double res = diff[0][k];

    if (n >= 2) {
        res += q * diff[1][k];
    }

    if (n >= 3 && k >= 1) {
        double coeff = q * (q - 1) / 2;
        double avg = (diff[2][k-1] + diff[2][k]) / 2;
        res += coeff * avg;
    }

    if (n >= 4 && k >= 1) {
        double coeff = q * (q - 1) * (q - 0.5) / 6;
        res += coeff * diff[3][k-1];
    }

    if (n >= 5 && k >= 2) {
        double coeff = q * (q * q - 1) * (q - 2) / 24;
        double avg = (diff[4][k-2] + diff[4][k-1]) / 2;
        res += coeff * avg;
    }

    if (n >= 6 && k >= 2) {
        double coeff = q * (q * q - 1) * (q - 0.5) * (q - 2) / 120;
        res += coeff * diff[5][k-2];
    }

    return res;
}

int main() {
    std::vector<Point> points = {
        {-11.02, 21.1230},
        {-10.78, 22.8279},
        {-10.23, 25.0046},
        {-9.89, 27.5928},
        {-9.65, 28.9933},
        {-9.03, 29.1933},
    };
    double h = 0.4;
    double a = 0;
    double b = 2;

    double x = -9.22;
    std::cout << " " << x << ": " << std::fixed << std::setprecision(2) << lagrange(points, x) << std::endl;
    std::cout << " " << x << ": " << std::fixed << std::setprecision(2) << aitken(points, x) << std::endl;
    std::cout << std::endl;

    double x1 = -10.44;
    double x2 = -9.41;
    std::cout << " " << x1 << ": " << std::fixed << std::setprecision(2) << newton1(points, x1) << std::endl;
    std::cout << " " << x2 << ": " << std::fixed << std::setprecision(2) << newton1(points, x2) << std::endl;
    std::cout << std::endl;

    x1 = -0.78;
    x2 = 1.43;
    x = 0.50;
    std::cout << " " << x1 << ": " << std::fixed << std::setprecision(2) << newton2(x1, h, a, b) << std::endl;
    std::cout << " " << x2 << ": " << std::fixed << std::setprecision(2) << newton2(x2, h, a, b) << std::endl;
    std::cout << " " << x << ": " << std::fixed << std::setprecision(2) << bessel(x, h, a, b) << std::endl;
    std::cout << std::endl;

    return 0;
}
