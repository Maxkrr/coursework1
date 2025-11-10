#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <iomanip>
#include <limits>

using namespace std;

// Функция для вычисления обратной матрицы методом Гаусса-Жордана
bool inverseMatrix(const vector<vector<double>>& A, vector<vector<double>>& invA) {
    int n = A.size();
    vector<vector<double>> temp = A;
    invA = vector<vector<double>>(n, vector<double>(n, 0.0));
    
    // Инициализация единичной матрицы
    for (int i = 0; i < n; i++) {
        invA[i][i] = 1.0;
    }
    
    // Прямой ход метода Гаусса
    for (int i = 0; i < n; i++) {
        // Поиск главного элемента
        double maxEl = abs(temp[i][i]);
        int maxRow = i;
        for (int k = i + 1; k < n; k++) {
            if (abs(temp[k][i]) > maxEl) {
                maxEl = abs(temp[k][i]);
                maxRow = k;
            }
        }
        
        // Перестановка строк
        if (maxRow != i) {
            swap(temp[i], temp[maxRow]);
            swap(invA[i], invA[maxRow]);
        }
        
        // Проверка на вырожденность
        if (fabs(temp[i][i]) < 1e-10) {
            return false;
        }
        
        // Нормализация текущей строки
        double divisor = temp[i][i];
        for (int j = 0; j < n; j++) {
            temp[i][j] /= divisor;
            invA[i][j] /= divisor;
        }
        
        // Исключение
        for (int k = 0; k < n; k++) {
            if (k != i) {
                double factor = temp[k][i];
                for (int j = 0; j < n; j++) {
                    temp[k][j] -= factor * temp[i][j];
                    invA[k][j] -= factor * invA[i][j];
                }
            }
        }
    }
    
    return true;
}

// Функция для умножения матрицы на вектор
vector<double> matrixVectorMultiply(const vector<vector<double>>& A, const vector<double>& v) {
    int n = A.size();
    vector<double> result(n, 0.0);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            result[i] += A[i][j] * v[j];
        }
    }
    return result;
}

// Функция для вычисления скалярного произведения
double dotProduct(const vector<double>& a, const vector<double>& b) {
    double result = 0.0;
    for (size_t i = 0; i < a.size(); i++) {
        result += a[i] * b[i];
    }
    return result;
}

// Функция для вычисления стандартного отклонения портфеля
double portfolioRisk(const vector<double>& weights, const vector<vector<double>>& covMatrix) {
    double risk = 0.0;
    int n = weights.size();
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            risk += weights[i] * weights[j] * covMatrix[i][j];
        }
    }
    return sqrt(risk);
}

// Функция для вычисления оптимальных весов касательного портфеля
vector<double> calculateTangencyPortfolio(const vector<double>& expectedReturns, 
                                         const vector<vector<double>>& covMatrix,
                                         double riskFreeRate) {
    int n = expectedReturns.size();
    
    // Вычисляем избыточную доходность
    vector<double> excessReturns(n);
    for (int i = 0; i < n; i++) {
        excessReturns[i] = expectedReturns[i] - riskFreeRate;
    }
    
    // Вычисляем обратную ковариационную матрицу
    vector<vector<double>> invCov(n, vector<double>(n));
    if (!inverseMatrix(covMatrix, invCov)) {
        throw runtime_error("Матрица ковариаций вырождена");
    }
    
    // Вычисляем оптимальные веса по формуле касательного портфеля
    vector<double> temp = matrixVectorMultiply(invCov, excessReturns);
    double sumTemp = 0.0;
    for (double val : temp) {
        sumTemp += val;
    }
    
    vector<double> weights(n);
    for (int i = 0; i < n; i++) {
        weights[i] = temp[i] / sumTemp;
    }
    
    return weights;
}

// Функция для построения эффективной границы
vector<vector<double>> buildEfficientFrontier(const vector<double>& expectedReturns,
                                            const vector<vector<double>>& covMatrix,
                                            int numPoints = 20) {
    int n = expectedReturns.size();
    vector<vector<double>> frontier;
    
    // Находим минимальную и максимальную ожидаемую доходность
    double minReturn = *min_element(expectedReturns.begin(), expectedReturns.end());
    double maxReturn = *max_element(expectedReturns.begin(), expectedReturns.end());
    
    // Вычисляем портфели для разных целевых доходностей
    for (int i = 0; i <= numPoints; i++) {
        double targetReturn = minReturn + (maxReturn - minReturn) * i / numPoints;
        
        try {
            // Для каждого уровня доходности находим портфель с минимальным риском
            // Упрощенная версия - используем метод множителей Лагранжа
            vector<double> ones(n, 1.0);
            
            vector<vector<double>> invCov(n, vector<double>(n));
            if (!inverseMatrix(covMatrix, invCov)) {
                continue;
            }
            
            double A = dotProduct(ones, matrixVectorMultiply(invCov, ones));
            double B = dotProduct(expectedReturns, matrixVectorMultiply(invCov, ones));
            double C = dotProduct(expectedReturns, matrixVectorMultiply(invCov, expectedReturns));
            
            double lambda1 = (C - B * targetReturn) / (A * C - B * B);
            double lambda2 = (A * targetReturn - B) / (A * C - B * B);
            
            vector<double> weights(n);
            for (int j = 0; j < n; j++) {
                weights[j] = lambda1;
                for (int k = 0; k < n; k++) {
                    weights[j] += lambda2 * expectedReturns[k] * invCov[j][k];
                }
            }
            
            double risk = portfolioRisk(weights, covMatrix);
            double actualReturn = dotProduct(weights, expectedReturns);
            
            frontier.push_back({actualReturn, risk});
        }
        catch (...) {
            continue;
        }
    }
    
    return frontier;
}

int main() {
    try {
        // Ввод количества активов
        int n;
        cout << "Введите количество активов: ";
        cin >> n;
        
        if (n <= 0) {
            cout << "Ошибка: количество активов должно быть положительным" << endl;
            return 1;
        }
        
        vector<double> expectedReturns(n);
        vector<vector<double>> covMatrix(n, vector<double>(n));
        double riskFreeRate;
        
        // Ввод ожидаемых доходностей
        cout << "Введите ожидаемые доходности для " << n << " активов (от -1 до 10):" << endl;
        for (int i = 0; i < n; i++) {
            cout << "Актив " << i + 1 << ": ";
            cin >> expectedReturns[i];
            while (expectedReturns[i] < -1.0 || expectedReturns[i] > 10.0) {
                cout << "Ошибка: доходность должна быть в диапазоне от -1 до 10. Введите снова: ";
                cin >> expectedReturns[i];
            }
        }
        
        // Ввод ковариационной матрицы
        cout << "Введите ковариационную матрицу " << n << "x" << n << ":" << endl;
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                cout << "Ковариация [" << i + 1 << "][" << j + 1 << "]: ";
                cin >> covMatrix[i][j];
            }
        }
        
        // Ввод безрисковой ставки
        cout << "Введите безрисковую ставку (от 0 до 1): ";
        cin >> riskFreeRate;
        while (riskFreeRate < 0.0 || riskFreeRate > 1.0) {
            cout << "Ошибка: безрисковая ставка должна быть в диапазоне от 0 до 1. Введите снова: ";
            cin >> riskFreeRate;
        }
        
        // Вычисление оптимальных весов касательного портфеля
        vector<double> weights = calculateTangencyPortfolio(expectedReturns, covMatrix, riskFreeRate);
        
        // Расчет характеристик портфеля
        double portfolioReturn = dotProduct(weights, expectedReturns);
        double risk = portfolioRisk(weights, covMatrix);
        double sharpeRatio = (portfolioReturn - riskFreeRate) / risk;
        
        // Построение эффективной границы
        vector<vector<double>> efficientFrontier = buildEfficientFrontier(expectedReturns, covMatrix);
        
        // Вывод результатов
        cout << fixed << setprecision(1);
        cout << "\n=== РЕЗУЛЬТАТЫ ОПТИМИЗАЦИИ ПОРТФЕЛЯ ===" << endl;
        
        cout << "Оптимальные доли: [";
        for (size_t i = 0; i < weights.size(); i++) {
            cout << weights[i] * 100;
            if (i < weights.size() - 1) cout << "%, ";
            else cout << "%]";
        }
        cout << endl;
        
        cout << fixed << setprecision(1);
        cout << "Ожидаемая доходность: " << portfolioReturn * 100 << "%" << endl;
        cout << "Риск портфеля: " << risk * 100 << "%" << endl;
        cout << fixed << setprecision(2);
        cout << "Коэффициент Шарпа: " << sharpeRatio << endl;
        
        // Вывод эффективной границы
        cout << "\nЭффективная граница (первые 5 точек):" << endl;
        cout << "Доходность%\tРиск%" << endl;
        for (size_t i = 0; i < min(efficientFrontier.size(), size_t(5)); i++) {
            cout << efficientFrontier[i][0] * 100 << "\t\t" << efficientFrontier[i][1] * 100 << endl;
        }
        if (efficientFrontier.size() > 5) {
            cout << "... (всего " << efficientFrontier.size() << " точек)" << endl;
        }
        
    } catch (const exception& e) {
        cout << "Ошибка: " << e.what() << endl;
        return 1;
    }
    
    return 0;
}
