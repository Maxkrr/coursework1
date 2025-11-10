#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <stdexcept>
#include <limits>
#include <algorithm>

using namespaces std;

// === Функции линейной алгебры ===

// Обратная матрица (метод Гаусса-Жордана)
bool inverseMatrix(const vector<vector<double>>& A, vector<vector<double>>& invA) {
    int n = static_cast<int>(A.size());
    invA.assign(n, vector<double>(n, 0.0));
    vector<vector<double>> temp = A;

    for (int i = 0; i < n; i++) invA[i][i] = 1.0;

    for (int i = 0; i < n; i++) {
        double maxEl = fabs(temp[i][i]);
        int maxRow = i;
        for (int k = i + 1; k < n; k++) {
            if (fabs(temp[k][i]) > maxEl) {
                maxEl = fabs(temp[k][i]);
                maxRow = k;
            }
        }
        if (maxRow != i) {
            std::swap(temp[i], temp[maxRow]);
            std::swap(invA[i], invA[maxRow]);
        }
        if (fabs(temp[i][i]) < 1e-12) return false;

        double divisor = temp[i][i];
        for (int j = 0; j < n; j++) {
            temp[i][j] /= divisor;
            invA[i][j] /= divisor;
        }

        for (int k = 0; k < n; k++) {
            if (k == i) continue;
            double factor = temp[k][i];
            for (int j = 0; j < n; j++) {
                temp[k][j] -= factor * temp[i][j];
                invA[k][j] -= factor * invA[i][j];
            }
        }
    }
    return true;
}

// Умножение матрицы на вектор
vector<double> matVec(const vector<vector<double>>& A, const vector<double>& v) {
    int n = static_cast<int>(A.size());
    vector<double> result(n, 0.0);
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            result[i] += A[i][j] * v[j];
    return result;
}

// Скалярное произведение
double dot(const vector<double>& a, const vector<double>& b) {
    double res = 0.0;
    for (size_t i = 0; i < a.size(); i++) res += a[i] * b[i];
    return res;
}

// Волатильность портфеля
double portfolioRisk(const vector<double>& w, const vector<vector<double>>& cov) {
    double sum = 0.0;
    int n = static_cast<int>(w.size());
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            sum += w[i] * w[j] * cov[i][j];
    return std::sqrt(sum);
}

// Применение ограничений (min/max) и нормировка весов
void applyConstraints(vector<double>& w, const vector<double>& minW, const vector<double>& maxW) {
    int n = w.size();
    for (int i = 0; i < n; i++) {
        if (w[i] < minW[i]) w[i] = minW[i];
        if (w[i] > maxW[i]) w[i] = maxW[i];
    }

    // Нормировка, чтобы сумма = 1
    double sumW = 0.0;
    for (double wi : w) sumW += wi;
    if (fabs(sumW) > 1e-12)
        for (double &wi : w) wi /= sumW;
}

// Оптимальный (касательный) портфель
vector<double> tangencyPortfolio(const vector<double>& mu, const vector<vector<double>>& cov,
                                 double rf, const vector<double>& minW, const vector<double>& maxW) {
    int n = static_cast<int>(mu.size());
    vector<double> excess(n);
    for (int i = 0; i < n; i++) excess[i] = mu[i] - rf;

    vector<vector<double>> invCov;
    if (!inverseMatrix(cov, invCov)) throw std::logic_error("Матрица ковариаций вырождена");

    vector<double> temp = matVec(invCov, excess);
    double sumTemp = 0.0;
    for (double x : temp) sumTemp += x;

    vector<double> w(n);
    for (int i = 0; i < n; i++) w[i] = temp[i] / sumTemp;

    applyConstraints(w, minW, maxW);
    return w;
}

// Эффективная граница (без ограничений)
vector<vector<double>> efficientFrontier(const vector<double>& mu, const vector<vector<double>>& cov, int numP = 20) {
    int n = static_cast<int>(mu.size());
    vector<vector<double>> frontier;

    double minR = *std::min_element(mu.begin(), mu.end());
    double maxR = *std::max_element(mu.begin(), mu.end());

    for (int i = 0; i <= numP; i++) {
        double tr = minR + (maxR - minR) * i / numP;
        vector<double> ones(n, 1.0), invMu, invOnes;
        vector<vector<double>> invCov;
        if (!inverseMatrix(cov, invCov)) continue;

        invMu = matVec(invCov, mu);
        invOnes = matVec(invCov, ones);
        double A = dot(ones, invOnes);
        double B = dot(mu, invOnes);
        double C = dot(mu, invMu);
        double d = A * C - B * B;
        if (fabs(d) < 1e-12) continue;
        double l1 = (C - B * tr) / d;
        double l2 = (A * tr - B) / d;

        vector<double> w(n, 0.0);
        for (int j = 0; j < n; ++j)
            w[j] = l1 * invOnes[j] + l2 * invMu[j];

        double risk = portfolioRisk(w, cov);
        double expected = dot(w, mu);
        frontier.push_back({expected, risk});
    }
    return frontier;
}

// === Главная функция ===
int main() {
    cout << "Введите количество активов: ";
    int n; cin >> n;
    if (n <= 0) {
        cerr << "Ошибка: количество активов должно быть положительным." << endl;
        return 1;
    }

    vector<double> mu(n);
    cout << "Ожидаемые доходности (от -1 до 10):\n";
    for (int i = 0; i < n; i++) {
        cout << "Актив " << (i+1) << ": ";
        cin >> mu[i];
        if (cin.fail() || mu[i] < -1.0 || mu[i] > 10.0) {
            cerr << "Ошибка: доходность актива #" << (i+1)
                 << " должна быть в диапазоне [-1, 10]." << endl;
            return 1;
        }
    }

    vector<vector<double>> cov(n, vector<double>(n));
    cout << "Введите ковариационную матрицу (" << n << "x" << n << "):\n";
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            cout << "Cov[" << (i+1) << "][" << (j+1) << "]: ";
            cin >> cov[i][j];
            if (cin.fail()) {
                cerr << "Ошибка: некорректный ввод элемента матрицы ковариаций." << endl;
                return 1;
            }
        }
    }

    // Проверка симметричности ковариационной матрицы
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (fabs(cov[i][j] - cov[j][i]) > 1e-8) {
                cerr << "Ошибка: матрица ковариаций должна быть симметричной." << endl;
                return 1;
            }
        }
    }

    double rf;
    cout << "Безрисковая ставка (0..1): ";
    cin >> rf;
    if (cin.fail() || rf < 0.0 || rf > 1.0) {
        cerr << "Ошибка: безрисковая ставка должна быть в диапазоне [0, 1]." << endl;
        return 1;
    }

    // Ввод ограничений
    vector<double> minW(n), maxW(n);
    cout << "Введите минимальные доли (в диапазоне -1..1):\n";
    for (int i = 0; i < n; i++) {
        cout << "Мин для актива " << (i+1) << ": ";
        cin >> minW[i];
        if (cin.fail() || minW[i] < -1.0 || minW[i] > 1.0) {
            cerr << "Ошибка: минимальная доля актива #" << (i+1)
                 << " должна быть в диапазоне [-1, 1]." << endl;
            return 1;
        }
    }

    cout << "Введите максимальные доли (в диапазоне -1..1, но больше минимальных):\n";
    for (int i = 0; i < n; i++) {
        cout << "Макс для актива " << (i+1) << ": ";
        cin >> maxW[i];
        if (cin.fail() || maxW[i] < -1.0 || maxW[i] > 1.0 || maxW[i] < minW[i]) {
            cerr << "Ошибка: максимальная доля актива #" << (i+1)
                 << " должна быть ≥ минимальной и в диапазоне [-1, 1]." << endl;
            return 1;
        }
    }

    // === Оптимизация ===
    try {
        vector<double> w = tangencyPortfolio(mu, cov, rf, minW, maxW);
        double portRet = dot(w, mu);
        double portRisk = portfolioRisk(w, cov);
        double sharpe = (portRet - rf) / portRisk;

        auto eff = efficientFrontier(mu, cov);

        cout << fixed << setprecision(2);
        cout << "\n=== РЕЗУЛЬТАТЫ ===\n";
        cout << "Оптимальные доли: [";
        for (size_t i = 0; i < w.size(); i++) {
            cout << w[i]*100 << "%";
            if (i+1 < w.size()) cout << ", ";
        }
        cout << "]\n";
        cout << "Ожидаемая доходность: " << portRet*100 << "%\n";
        cout << "Риск портфеля: " << portRisk*100 << "%\n";
        cout << "Коэффициент Шарпа: " << sharpe << "\n";

        cout << "\nЭффективная граница (первые 5 точек):\nДоходность%\tРиск%\n";
        for (size_t i = 0; i < min<size_t>(5, eff.size()); i++)
            cout << eff[i][0]*100 << "\t\t" << eff[i][1]*100 << "\n";
    }
    catch (const exception& e) {
        cerr << "Ошибка оптимизации: " << e.what() << endl;
        return 1;
    }

    return 0;
}
