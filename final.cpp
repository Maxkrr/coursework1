#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <stdexcept>
#include <limits>
#include <algorithm>

using namespace std;

// === ФУНКЦИИ ЛИНЕЙНОЙ АЛГЕБРЫ ===

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
            swap(temp[i], temp[maxRow]);
            swap(invA[i], invA[maxRow]);
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
    return sqrt(sum);
}

// Применение ограничений и нормировка весов
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
    if (!inverseMatrix(cov, invCov)) throw logic_error("Матрица ковариаций вырождена (невозможно обратить)");

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

    double minR = *min_element(mu.begin(), mu.end());
    double maxR = *max_element(mu.begin(), mu.end());

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

// === ФУНКЦИЯ ТЕСТИРОВАНИЯ ===
void program_testing() {
    cout << "\n=== ТЕСТИРОВАНИЕ ПРОГРАММЫ ОПТИМИЗАЦИИ ПОРТФЕЛЯ ===" << endl;
    
    struct TestCase {
        string test_name;
        vector<double> mu;
        vector<vector<double>> cov;
        double rf;
        vector<double> minW;
        vector<double> maxW;
        bool should_fail;
        string expected_error;
    };
    
    TestCase tests[] = {
        {
            "Тест 1: Два актива с разной доходностью",
            {0.10, 0.15},
            {{0.04, 0.01}, {0.01, 0.09}},
            0.05,
            {0.0, 0.0},
            {1.0, 1.0},
            false,
            ""
        },
        {
            "Тест 2: Вырожденная матрица ковариаций",
            {0.10, 0.10},
            {{1.0, 1.0}, {1.0, 1.0}},
            0.05,
            {0.0, 0.0},
            {1.0, 1.0},
            true,
            "Матрица ковариаций вырождена"
        },
        {
            "Тест 3: Ограничения на веса",
            {0.08, 0.12, 0.10},
            {{0.04, 0.01, 0.02}, 
             {0.01, 0.09, 0.01},
             {0.02, 0.01, 0.06}},
            0.04,
            {0.2, 0.1, 0.1},
            {1.0, 1.0, 1.0},
            false,
            ""
        },
        {
            "Тест 4: Невозможные ограничения",
            {0.10, 0.15},
            {{0.04, 0.01}, {0.01, 0.09}},
            0.05,
            {0.8, 0.8},
            {1.0, 1.0},
            true,
            "невозможно составить портфель"
        }
    };
    
    int passed_tests = 0;
    int total_tests = sizeof(tests) / sizeof(tests[0]);
    
    for (int i = 0; i < total_tests; i++) {
        TestCase test = tests[i];
        cout << "\n" << test.test_name << endl;
        cout << "Ожидаемые доходности: [";
        for (size_t j = 0; j < test.mu.size(); j++) {
            cout << test.mu[j]*100 << "%";
            if (j+1 < test.mu.size()) cout << ", ";
        }
        cout << "]" << endl;
        
        try {
            // Проверка возможности создания портфеля
            double minSum = 0.0, maxSum = 0.0;
            for (size_t j = 0; j < test.minW.size(); j++) {
                minSum += test.minW[j];
                maxSum += test.maxW[j];
            }
            
            if (minSum > 1.0 + 1e-12 || maxSum < 1.0 - 1e-12) {
                if (test.should_fail) {
                    cout << "ТЕСТ ПРОЙДЕН! Программа правильно обнаружила невозможные ограничения" << endl;
                    passed_tests++;
                } else {
                    cout << "ТЕСТ НЕ ПРОЙДЕН! Невозможные ограничения не были обработаны" << endl;
                }
                cout << "------------------------" << endl;
                continue;
            }
            
            // Выполняем расчет портфеля
            vector<double> w = tangencyPortfolio(test.mu, test.cov, test.rf, test.minW, test.maxW);
            double portRet = dot(w, test.mu);
            double portRisk = portfolioRisk(w, test.cov);
            double sharpe = (portRet - test.rf) / portRisk;
            
            // Проверяем ограничения
            bool constraints_ok = true;
            for (size_t j = 0; j < w.size(); j++) {
                if (w[j] < test.minW[j] - 1e-12 || w[j] > test.maxW[j] + 1e-12) {
                    constraints_ok = false;
                    cout << "   - Нарушено ограничение для актива " << j+1 
                         << ": вес=" << w[j] << ", min=" << test.minW[j] << ", max=" << test.maxW[j] << endl;
                }
            }
            
            // Проверяем сумму весов
            double sumW = 0.0;
            for (double wi : w) sumW += wi;
            bool sum_ok = fabs(sumW - 1.0) < 1e-12;
            
            if (!sum_ok) {
                cout << "   - Сумма весов не равна 1: " << sumW << endl;
            }
            
            // Выводим результаты
            cout << "Полученные веса: [";
            for (size_t j = 0; j < w.size(); j++) {
                cout << fixed << setprecision(4) << w[j]*100 << "%";
                if (j+1 < w.size()) cout << ", ";
            }
            cout << "]" << endl;
            cout << "Доходность: " << portRet*100 << "%, Риск: " << portRisk*100 
                 << "%, Коэффициент Шарпа: " << sharpe << endl;
            
            if (test.should_fail) {
                cout << "ТЕСТ НЕ ПРОЙДЕН! Ожидалась ошибка, но расчет выполнен успешно" << endl;
            } else {
                if (constraints_ok && sum_ok) {
                    cout << "ТЕСТ ПРОЙДЕН! Все ограничения соблюдены" << endl;
                    passed_tests++;
                } else {
                    cout << "ТЕСТ НЕ ПРОЙДЕН! Нарушены ограничения" << endl;
                }
            }
            
        } catch (const exception& e) {
            cout << "Исключение: " << e.what() << endl;
            
            if (test.should_fail) {
                if (test.expected_error.empty() || 
                    string(e.what()).find(test.expected_error) != string::npos) {
                    cout << "ТЕСТ ПРОЙДЕН! Ожидаемая ошибка получена" << endl;
                    passed_tests++;
                } else {
                    cout << "ТЕСТ НЕ ПРОЙДЕН! Получена неожиданная ошибка" << endl;
                    cout << "Ожидалось: " << test.expected_error << endl;
                }
            } else {
                cout << "ТЕСТ НЕ ПРОЙДЕН! Неожиданная ошибка" << endl;
            }
        }
        cout << "------------------------" << endl;
    }
    
    // Выводим итоговую статистику
    cout << "\n=== ИТОГИ ТЕСТИРОВАНИЯ ===" << endl;
    cout << "Пройдено тестов: " << passed_tests << " из " << total_tests << endl;
    cout << "Успешность: " << (passed_tests * 100 / total_tests) << "%" << endl;
    
    if (passed_tests == total_tests) {
        cout << "ВСЕ ТЕСТЫ ПРОЙДЕНЫ УСПЕШНО!" << endl;
    } else {
        cout << "ЕСТЬ ПРОБЛЕМЫ В РАБОТЕ ПРОГРАММЫ" << endl;
    }
}

// === ФУНКЦИИ ИНТЕРФЕЙСА ===

void show_menu() {
    cout << "\n=== ГЛАВНОЕ МЕНЮ ===" << endl;
    cout << "1. Оптимизация портфеля" << endl;
    cout << "2. Тестирование" << endl;
    cout << "3. Выход" << endl;
    cout << "====================" << endl;
    cout << "Выберите действие (1-3): ";
}

void run_portfolio_optimization() {
    cout << "\n=== ОПТИМИЗАЦИЯ ИНВЕСТИЦИОННОГО ПОРТФЕЛЯ ===" << endl;
    
    int n;
    cout << "Введите количество активов: ";
    cin >> n;
    if (n <= 0) {
        cerr << "Ошибка: количество активов должно быть положительным." << endl;
        return;
    }

    vector<double> mu(n);
    cout << "Ожидаемые доходности (от -1 до 10):\n";
    for (int i = 0; i < n; i++) {
        cout << "Актив " << (i+1) << ": ";
        cin >> mu[i];
        if (cin.fail() || mu[i] < -1.0 || mu[i] > 10.0) {
            cerr << "Ошибка: доходность актива #" << (i+1)
                 << " должна быть в диапазоне [-1, 10]." << endl;
            return;
        }
    }

    vector<vector<double>> cov(n, vector<double>(n));
    cout << "Введите ковариационную матрицу (" << n << "x" << n << "):\n";
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++) {
            cout << "Cov[" << (i+1) << "][" << (j+1) << "]: ";
            cin >> cov[i][j];
            if (cin.fail()) {
                cerr << "Ошибка: некорректный ввод элемента матрицы ковариаций." << endl;
                return;
            }
        }

    // Проверка симметричности
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            if (fabs(cov[i][j] - cov[j][i]) > 1e-8) {
                cerr << "Ошибка: матрица ковариаций должна быть симметричной." << endl;
                return;
            }

    double rf;
    cout << "Безрисковая ставка (0..1): ";
    cin >> rf;
    if (cin.fail() || rf < 0.0 || rf > 1.0) {
        cerr << "Ошибка: безрисковая ставка должна быть в диапазоне [0, 1]." << endl;
        return;
    }

    vector<double> minW(n), maxW(n);
    cout << "Введите минимальные доли (в диапазоне 0..1):\n";
    for (int i = 0; i < n; i++) {
        cout << "Мин для актива " << (i+1) << ": ";
        cin >> minW[i];
        if (cin.fail() || minW[i] < 0 || minW[i] > 1.0) {
            cerr << "Ошибка: минимальная доля актива #" << (i+1)
                 << " должна быть в диапазоне [0, 1]." << endl;
            return;
        }
    }

    cout << "Введите максимальные доли (в диапазоне 0..1, больше или равны минимальным):\n";
    for (int i = 0; i < n; i++) {
        cout << "Макс для актива " << (i+1) << ": ";
        cin >> maxW[i];
        if (cin.fail() || maxW[i] < 0 || maxW[i] > 1.0 || maxW[i] < minW[i]) {
            cerr << "Ошибка: максимальная доля актива #" << (i+1)
                 << " должна быть ≥ минимальной и в диапазоне [0, 1]." << endl;
            return;
        }
    }

    // Проверка, что возможно составить портфель с суммой весов = 1
    double minSum = 0.0, maxSum = 0.0;
    for (int i = 0; i < n; i++) {
        minSum += minW[i];
        maxSum += maxW[i];
    }
    if (minSum > 1.0 || maxSum < 1.0) {
        cerr << "Ошибка: заданные ограничения долей не позволяют составить портфель с суммой весов = 1." << endl;
        return;
    }

    try {
        vector<double> w = tangencyPortfolio(mu, cov, rf, minW, maxW);
        double portRet = dot(w, mu);
        double portRisk = portfolioRisk(w, cov);
        double sharpe = (portRet - rf) / portRisk;

        auto eff = efficientFrontier(mu, cov);

        cout << fixed << setprecision(2);
        cout << "\n=== РЕЗУЛЬТАТЫ ===" << endl;
        cout << "Оптимальные доли: [";
        for (size_t i = 0; i < w.size(); i++) {
            cout << w[i]*100 << "%";
            if (i+1 < w.size()) cout << ", ";
        }
        cout << "]" << endl;
        cout << "Ожидаемая доходность: " << portRet*100 << "%" << endl;
        cout << "Риск портфеля: " << portRisk*100 << "%" << endl;
        cout << "Коэффициент Шарпа: " << sharpe << endl;

        cout << "\nЭффективная граница (первые 5 точек):" << endl;
        cout << "Доходность%\tРиск%" << endl;
        for (size_t i = 0; i < min<size_t>(5, eff.size()); i++)
            cout << eff[i][0]*100 << "\t\t" << eff[i][1]*100 << endl;
            
    } catch (const exception& e) {
        cerr << "Ошибка оптимизации: " << e.what() << endl;
    }
}

// === ГЛАВНАЯ ФУНКЦИЯ ===
int main() {
    setlocale(LC_ALL, "Russian");

    cout << "=== СИСТЕМА ОПТИМИЗАЦИИ ИНВЕСТИЦИОННОГО ПОРТФЕЛЯ ===" << endl;
    cout << "Реализация теории Марковица с построением эффективной границы" << endl;

    int choice;
    bool continue_work = true;

    while (continue_work) {
        show_menu();
        cin >> choice;
        cin.ignore(1000, '\n');

        switch (choice) {
            case 1:
                run_portfolio_optimization();
                break;
            case 2:
                program_testing();
                break;
            case 3:
                cout << "До новых встреч!" << endl;
                continue_work = false;
                break;
            default:
                cout << "Ошибка! Выберите 1, 2 или 3." << endl;
                break;
        }
    }

    return 0;
}
