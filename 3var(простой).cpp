#include <iostream>
#include <vector>
#include <cmath>
#include <limits>

using namespace std;

struct PortfolioResult {
    vector<double> weights;
    double expected_return;
    double risk; // стандартное отклонение (квадратный корень из дисперсии)
    double sharpe_ratio;
};

// Умножение матрицы и вектора
vector<double> matVecMul(const vector<vector<double>>& mat, const vector<double>& vec) {
    int n = mat.size();
    vector<double> result(n, 0.0);
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            result[i] += mat[i][j] * vec[j];
    return result;
}

// Скалярное произведение
double dot(const vector<double>& a, const vector<double>& b) {
    double res = 0.0;
    for (size_t i = 0; i < a.size(); i++)
        res += a[i] * b[i];
    return res;
}

// Минимизация риска при фиксированной доходности (грубо, перебор с шагом)
PortfolioResult optimizePortfolio(const vector<double>& expected_returns,
                                  const vector<vector<double>>& cov_matrix,
                                  double target_return,
                                  const vector<double>& min_limits,
                                  const vector<double>& max_limits,
                                  double risk_free_rate) {

    int n = expected_returns.size();
    PortfolioResult best;
    best.risk = numeric_limits<double>::max();

    // Для простоты - перебираем веса с шагом 0.05, сумма весов ровна 1
    // и учитываем ограничения
    vector<double> weights(n, 0.0);
    double step = 0.05;

    // Рекурсивный перебор (упрощённый) для 3 активов
    for (int w0 = 0; w0 <= 20; w0++) {
        for (int w1 = 0; w1 <= 20; w1++) {
            double w2 = 20 - w0 - w1;
            if (w2 < 0 || w2 > 20) continue;

            weights[0] = w0 * step;
            weights[1] = w1 * step;
            weights[2] = w2 * step;

            // Проверяем ограничения
            bool valid = true;
            for (int i = 0; i < n; i++) {
                if (weights[i] < min_limits[i] || weights[i] > max_limits[i]) {
                    valid = false;
                    break;
                }
            }
            if (!valid) continue;

            double port_return = dot(weights, expected_returns);
            if (abs(port_return - target_return) > 0.01) continue; // доходность около target_return

            vector<double> temp = matVecMul(cov_matrix, weights);
            double variance = dot(weights, temp);
            double risk = sqrt(variance);

            if (risk < best.risk) {
                best.risk = risk;
                best.expected_return = port_return;
                best.weights = weights;
                best.sharpe_ratio = (port_return - risk_free_rate) / risk;
            }
        }
    }
    return best;
}

int main() {
    // Входные данные (пример 1)
    vector<double> expected_returns = {0.12, 0.08, 0.15};
    vector<vector<double>> cov_matrix = {
        {0.04, 0.01, 0.02},
        {0.01, 0.03, 0.01},
        {0.02, 0.01, 0.06}
    };
    double risk_free_rate = 0.0;
    vector<double> min_limits = {0.0, 0.0, 0.0};
    vector<double> max_limits = {1.0, 1.0, 1.0};

    // Пример: составим эффективную границу, перебрав доходности от 0.08 до 0.15
    cout << "Оптимизация портфеля по Марковицу с построением эффективной границы:\n";
    for (double target_return = 0.08; target_return <= 0.15; target_return += 0.01) {
        PortfolioResult res = optimizePortfolio(expected_returns, cov_matrix, target_return, min_limits, max_limits, risk_free_rate);

        if (res.weights.size() == 0) continue;
        cout << "Целевая доходность: " << target_return * 100 << "%\n";
        cout << "Оптимальные доли: [";
        for (size_t i = 0; i < res.weights.size(); i++) {
            cout << res.weights[i] * 100 << "%";
            if (i != res.weights.size()-1) cout << ", ";
        }
        cout << "]\n";
        cout << "Ожидаемая доходность: " << res.expected_return * 100 << "%\n";
        cout << "Риск портфеля: " << res.risk * 100 << "%\n";
        cout << "Коэффициент Шарпа: " << res.sharpe_ratio << "\n\n";
    }

    return 0;
}
