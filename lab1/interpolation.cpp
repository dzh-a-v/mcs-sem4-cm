#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <iomanip>

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
    auto f = [](double x) { return std::exp(1.0) - std::pow(std::log(x), 2); };

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

double stirling(double x, double h, double a, double b) {
    auto f = [](double x) { return std::exp(1.0) - std::pow(std::log(x), 2); };

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

    int center = n / 2;
    double x0 = points[center].first;
    double t = (x - x0) / h;

    double res = diff[0][center];

    if (n > 1) {
        res += t * (diff[1][center-1] + diff[1][center]) / 2;
    }

    if (n > 2) {
        res += t * t / 2 * diff[2][center-1];
    }

    if (n > 3) {
        res += t * (t * t - 1) / 6 * (diff[3][center-2] + diff[3][center-1]) / 2;
    }

    if (n > 4) {
        res += t * t * (t * t - 1) / 24 * diff[4][center-2];
    }

    return res;
}

int main() {
    std::vector<Point> points = {
        {11.153, -3.234},
        {11.454, 5.321},
        {11.673, -1.123},
        {11.879, 0.393},
        {12.009, 8.939},
        {12.231, 141.231},
        {12.549, 15.001}
    };
    double h = 1;
    double a = 1;
    double b = 5;

    double x = 12.776;
    std::cout << " " << x << ": " << std::fixed << std::setprecision(2) << lagrange(points, x) << std::endl;
    std::cout << " " << x << ": " << std::fixed << std::setprecision(2) << aitken(points, x) << std::endl;
    std::cout << std::endl;

    double x1 = 11.515;
    double x2 = 11.995;
    std::cout << " " << x1 << ": " << std::fixed << std::setprecision(2) << newton1(points, x1) << std::endl;
    std::cout << " " << x2 << ": " << std::fixed << std::setprecision(2) << newton1(points, x2) << std::endl;
    std::cout << std::endl;

    x1 = 0.77;
    x2 = 4.82;
    x = 3.2;
    std::cout << " " << x1 << ": " << std::fixed << std::setprecision(2) << newton2(x1, h, a, b) << std::endl;
    std::cout << " " << x2 << ": " << std::fixed << std::setprecision(2) << newton2(x2, h, a, b) << std::endl;
    std::cout << " " << x << ": " << std::fixed << std::setprecision(2) << newton2(x, h, a, b) << std::endl;
    std::cout << std::endl;

    return 0;
}
