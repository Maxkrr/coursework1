#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <algorithm>

using std::cin;
using std::cout;
using std::vector;
using std::string;

// === Улучшенные функции линейной алгебры ===

// Обратная матрица (упрощенный метод)
bool inverseMatrix(const vector<vector<double>>& A, vector<vector<double>>& invA) {
    int n = A.size();
    invA = A;
    vector<vector<double>> temp = A;
    
    // Инициализация единичной матрицы
    for (int i = 0; i < n; i++) {
        invA[i].assign(n, 0.0);
        invA[i][i] = 1.0;
    }
    
    // Прямой ход метода Гаусса
    for (int i = 0; i < n; i++) {
        // Поиск главного элемента
        double maxEl = std::abs(temp[i][i]);
        int maxRow = i;
        for (int k = i + 1; k < n; k++) {
            if (std::abs(temp[k][i]) > maxEl) {
                maxEl = std::abs(temp[k][i]);
                maxRow = k;
            }
        }
        
        // Перестановка строк
        if (maxRow != i) {
            std::swap(temp[i], temp[maxRow]);
            std::swap(invA[i], invA[maxRow]);
        }
        
        // Проверка на вырожденность
        if (std::abs(temp[i][i]) < 1e-10) {
            return false;
        }
        
        // Нормализация текущей строки
        double divisor = temp[i][i];
        for (int j = 0; j < n; j++) {
            temp[i][j] /= divisor;
            invA[i][j] /= divisor;
        }
        
        // Исключение элементов в столбце
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

// Умножение матрицы на вектор
vector<double> multiplyMatrixVector(const vector<vector<double>>& A, const vector<double>& v) {
    int n = A.size();
    vector<double> result(n, 0.0);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            result[i] += A[i][j] * v[j];
        }
    }
    return result;
}

// Скалярное произведение векторов
double dotProduct(const vector<double>& a, const vector<double>& b) {
    double result = 0.0;
    for (size_t i = 0; i < a.size(); i++) {
        result += a[i] * b[i];
    }
    return result;
}

// Расчет риска портфеля
double calculatePortfolioRisk(const vector<double>& weights, const vector<vector<double>>& covariance) {
    double variance = 0.0;
    int n = weights.size();
    
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            variance += weights[i] * weights[j] * covariance[i][j];
        }
    }
    
    return std::sqrt(variance);
}

// Применение ограничений к весам
void applyWeightConstraints(vector<double>& weights, const vector<double>& minWeights, const vector<double>& maxWeights) {
    int n = weights.size();
    
    // Применение минимальных и максимальных ограничений
    for (int i = 0; i < n; i++) {
        if (weights[i] < minWeights[i]) weights[i] = minWeights[i];
        if (weights[i] > maxWeights[i]) weights[i] = maxWeights[i];
    }
    
    // Нормировка суммы весов к 1
    double totalWeight = 0.0;
    for (double w : weights) totalWeight += w;
    
    if (std::abs(totalWeight) > 1e-10) {
        for (int i = 0; i < n; i++) {
            weights[i] /= totalWeight;
        }
    }
}

// Проверка корректности входных данных
bool validateInput(const vector<double>& returns, const vector<vector<double>>& covariance, 
                  double riskFreeRate, const vector<double>& minWeights, const vector<double>& maxWeights) {
    int n = returns.size();
    
    // Проверка размеров
    if (n <= 0) return false;
    if (covariance.size() != n) return false;
    for (const auto& row : covariance) {
        if (row.size() != n) return false;
    }
    if (minWeights.size() != n || maxWeights.size() != n) return false;
    
    // Проверка диапазонов
    if (riskFreeRate < 0 || riskFreeRate > 1) return false;
    
    for (double r : returns) {
        if (r < -1.0 || r > 10.0) return false;
    }
    
    // Проверка ограничений весов
    for (int i = 0; i < n; i++) {
        if (minWeights[i] < 0 || maxWeights[i] > 1 || minWeights[i] > maxWeights[i]) {
            return false;
        }
    }
    
    return true;
}

// Расчет оптимального (касательного) портфеля
vector<double> calculateOptimalPortfolio(const vector<double>& expectedReturns, 
                                        const vector<vector<double>>& covarianceMatrix,
                                        double riskFreeRate,
                                        const vector<double>& minWeights, 
                                        const vector<double>& maxWeights) {
    int n = expectedReturns.size();
    
    // Вектор избыточной доходности
    vector<double> excessReturns(n);
    for (int i = 0; i < n; i++) {
        excessReturns[i] = expectedReturns[i] - riskFreeRate;
    }
    
    // Обратная матрица ковариаций
    vector<vector<double>> inverseCovariance;
    if (!inverseMatrix(covarianceMatrix, inverseCovariance)) {
        throw std::runtime_error("Ошибка: матрица ковариаций не может быть обращена");
    }
    
    // Расчет оптимальных весов по формуле касательного портфеля
    vector<double> tempWeights = multiplyMatrixVector(inverseCovariance, excessReturns);
    
    // Нормировка весов
    double sum = 0.0;
    for (double w : tempWeights) sum += w;
    
    vector<double> optimalWeights(n);
    for (int i = 0; i < n; i++) {
        optimalWeights[i] = tempWeights[i] / sum;
    }
    
    // Применение ограничений
    applyWeightConstraints(optimalWeights, minWeights, maxWeights);
    
    return optimalWeights;
}

// Построение эффективной границы
vector<vector<double>> buildEfficientFrontier(const vector<double>& expectedReturns,
                                             const vector<vector<double>>& covarianceMatrix,
                                             int numPoints = 20) {
    int n = expectedReturns.size();
    vector<vector<double>> frontier;
    
    // Обратная матрица ковариаций (вычисляется один раз)
    vector<vector<double>> inverseCovariance;
    if (!inverseMatrix(covarianceMatrix, inverseCovariance)) {
        return frontier; // Пустая граница в случае ошибки
    }
    
    // Вспомогательные векторы
    vector<double> ones(n, 1.0);
    vector<double> invOnes = multiplyMatrixVector(inverseCovariance, ones);
    vector<double> invReturns = multiplyMatrixVector(inverseCovariance, expectedReturns);
    
    // Константы для расчетов
    double A = dotProduct(ones, invOnes);
    double B = dotProduct(expectedReturns, invOnes);
    double C = dotProduct(expectedReturns, invReturns);
    double discriminant = A * C - B * B;
    
    if (std::abs(discriminant) < 1e-12) {
        return frontier;
    }
    
    // Диапазон доходностей для границы
    double minReturn = *std::min_element(expectedReturns.begin(), expectedReturns.end());
    double maxReturn = *std::max_element(expectedReturns.begin(), expectedReturns.end());
    
    for (int i = 0; i <= numPoints; i++) {
        double targetReturn = minReturn + (maxReturn - minReturn) * i / numPoints;
        
        // Расчет весов для заданной доходности
        double lambda1 = (C - B * targetReturn) / discriminant;
        double lambda2 = (A * targetReturn - B) / discriminant;
        
        vector<double> weights(n);
        for (int j = 0; j < n; j++) {
            weights[j] = lambda1 * invOnes[j] + lambda2 * invReturns[j];
        }
        
        // Расчет риска и доходности
        double risk = calculatePortfolioRisk(weights, covarianceMatrix);
        double returnValue = dotProduct(weights, expectedReturns);
        
        frontier.push_back({returnValue, risk});
    }
    
    return frontier;
}

// === Главная функция ===
int main() {
    cout << "=== СИСТЕМА ОПТИМИЗАЦИИ ПОРТФЕЛЯ ЦЕННЫХ БУМАГ ===\n\n";
    
    // Ввод количества активов
    cout << "Введите количество активов: ";
    int numAssets;
    cin >> numAssets;
    
    if (numAssets <= 0) {
        cout << "Ошибка: количество активов должно быть положительным\n";
        return 1;
    }
    
    // Ввод ожидаемых доходностей
    vector<double> expectedReturns(numAssets);
    cout << "\nВведите ожидаемые доходности (в долях, от -1.0 до 10.0):\n";
    for (int i = 0; i < numAssets; i++) {
        cout << "Актив " << (i + 1) << ": ";
        cin >> expectedReturns[i];
    }
    
    // Ввод ковариационной матрицы
    vector<vector<double>> covarianceMatrix(numAssets, vector<double>(numAssets));
    cout << "\nВведите ковариационную матрицу (" << numAssets << "x" << numAssets << "):\n";
    for (int i = 0; i < numAssets; i++) {
        for (int j = 0; j < numAssets; j++) {
            cout << "Cov[" << (i + 1) << "][" << (j + 1) << "]: ";
            cin >> covarianceMatrix[i][j];
        }
    }
    
    // Ввод безрисковой ставки
    double riskFreeRate;
    cout << "\nБезрисковая ставка (в долях, 0.0-1.0): ";
    cin >> riskFreeRate;
    
    // Ввод ограничений на веса
    vector<double> minWeights(numAssets), maxWeights(numAssets);
    cout << "\nВведите минимальные доли активов (в долях):\n";
    for (int i = 0; i < numAssets; i++) {
        cout << "Мин. доля актива " << (i + 1) << ": ";
        cin >> minWeights[i];
    }
    
    cout << "\nВведите максимальные доли активов (в долях):\n";
    for (int i = 0; i < numAssets; i++) {
        cout << "Макс. доля актива " << (i + 1) << ": ";
        cin >> maxWeights[i];
    }
    
    // Проверка входных данных
    if (!validateInput(expectedReturns, covarianceMatrix, riskFreeRate, minWeights, maxWeights)) {
        cout << "\nОшибка: некорректные входные данные!\n";
        return 1;
    }
    
    try {
        // Расчет оптимального портфеля
        vector<double> optimalWeights = calculateOptimalPortfolio(
            expectedReturns, covarianceMatrix, riskFreeRate, minWeights, maxWeights);
        
        // Расчет характеристик портфеля
        double portfolioReturn = dotProduct(optimalWeights, expectedReturns);
        double portfolioRisk = calculatePortfolioRisk(optimalWeights, covarianceMatrix);
        double sharpeRatio = (portfolioReturn - riskFreeRate) / portfolioRisk;
        
        // Построение эффективной границы
        vector<vector<double>> efficientFrontier = buildEfficientFrontier(expectedReturns, covarianceMatrix);
        
        // Вывод результатов
        cout << std::fixed << std::setprecision(2);
        cout << "\n=== РЕЗУЛЬТАТЫ ОПТИМИЗАЦИИ ===\n";
        
        // Оптимальные доли
        cout << "Оптимальные доли: [";
        for (size_t i = 0; i < optimalWeights.size(); i++) {
            cout << optimalWeights[i] * 100 << "%";
            if (i < optimalWeights.size() - 1) cout << ", ";
        }
        cout << "]\n";
        
        cout << "Ожидаемая доходность: " << portfolioReturn * 100 << "%\n";
        cout << "Риск портфеля: " << portfolioRisk * 100 << "%\n";
        cout << "Коэффициент Шарпа: " << std::setprecision(3) << sharpeRatio << "\n";
        
        // Вывод эффективной границы
        if (!efficientFrontier.empty()) {
            cout << "\nЭффективная граница (первые 5 точек):\n";
            cout << "Доходность%\tРиск%\n";
            for (size_t i = 0; i < std::min(size_t(5), efficientFrontier.size()); i++) {
                cout << efficientFrontier[i][0] * 100 << "\t\t" 
                     << efficientFrontier[i][1] * 100 << "\n";
            }
        }
        
    } catch (const std::exception& e) {
        cout << "\nОшибка при расчетах: " << e.what() << "\n";
        return 1;
    }
    
    return 0;
}
