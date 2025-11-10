#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <iomanip>

using namespace std;

class PortfolioOptimizer {
private:
    vector<double> expectedReturns;
    vector<vector<double>> covarianceMatrix;
    double riskFreeRate;
    vector<vector<double>> constraints;
    int n;

    // Вычисление обратной матрицы 3x3
    bool inverseMatrix3x3(const vector<vector<double>>& A, vector<vector<double>>& invA) {
        double det = A[0][0] * (A[1][1] * A[2][2] - A[1][2] * A[2][1])
                   - A[0][1] * (A[1][0] * A[2][2] - A[1][2] * A[2][0])
                   + A[0][2] * (A[1][0] * A[2][1] - A[1][1] * A[2][0]);
        
        if (fabs(det) < 1e-10) return false;
        
        double invDet = 1.0 / det;
        invA[0][0] = (A[1][1] * A[2][2] - A[1][2] * A[2][1]) * invDet;
        invA[0][1] = (A[0][2] * A[2][1] - A[0][1] * A[2][2]) * invDet;
        invA[0][2] = (A[0][1] * A[1][2] - A[0][2] * A[1][1]) * invDet;
        invA[1][0] = (A[1][2] * A[2][0] - A[1][0] * A[2][2]) * invDet;
        invA[1][1] = (A[0][0] * A[2][2] - A[0][2] * A[2][0]) * invDet;
        invA[1][2] = (A[0][2] * A[1][0] - A[0][0] * A[1][2]) * invDet;
        invA[2][0] = (A[1][0] * A[2][1] - A[1][1] * A[2][0]) * invDet;
        invA[2][1] = (A[0][1] * A[2][0] - A[0][0] * A[2][1]) * invDet;
        invA[2][2] = (A[0][0] * A[1][1] - A[0][1] * A[1][0]) * invDet;
        
        return true;
    }

    // Умножение матрицы на вектор
    vector<double> matrixVectorMultiply(const vector<vector<double>>& A, const vector<double>& v) {
        vector<double> result(n, 0.0);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                result[i] += A[i][j] * v[j];
            }
        }
        return result;
    }

    // Скалярное произведение
    double dotProduct(const vector<double>& a, const vector<double>& b) {
        double result = 0.0;
        for (int i = 0; i < n; i++) {
            result += a[i] * b[i];
        }
        return result;
    }

    // Проверка ограничений
    bool checkConstraints(const vector<double>& weights) {
        for (int i = 0; i < n; i++) {
            if (weights[i] < constraints[i][0] || weights[i] > constraints[i][1]) {
                return false;
            }
        }
        return true;
    }

public:
    PortfolioOptimizer(const vector<double>& returns, const vector<vector<double>>& covMatrix, 
                      double riskFree, const vector<vector<double>>& constr) {
        expectedReturns = returns;
        covarianceMatrix = covMatrix;
        riskFreeRate = riskFree;
        constraints = constr;
        n = returns.size();
    }

    // Расчет оптимальных весов по Марковицу
    vector<double> calculateOptimalWeights() {
        // Вычисляем избыточную доходность
        vector<double> excessReturns(n);
        for (int i = 0; i < n; i++) {
            excessReturns[i] = expectedReturns[i] - riskFreeRate;
        }
        
        // Вычисляем обратную ковариационную матрицу
        vector<vector<double>> invCov(n, vector<double>(n));
        if (!inverseMatrix3x3(covarianceMatrix, invCov)) {
            cout << "Ошибка: матрица вырождена. Используется равномерное распределение." << endl;
            return vector<double>(n, 1.0/n);
        }
        
        // Вычисляем оптимальные веса по формуле Марковица
        vector<double> tempWeights = matrixVectorMultiply(invCov, excessReturns);
        double sumTempWeights = 0.0;
        for (double w : tempWeights) sumTempWeights += w;
        
        vector<double> weights(n);
        for (int i = 0; i < n; i++) {
            weights[i] = tempWeights[i] / sumTempWeights;
        }
        
        // Проверяем и корректируем ограничения
        if (!checkConstraints(weights)) {
            adjustWeightsForConstraints(weights);
        }
        
        return weights;
    }

    // Корректировка весов с учетом ограничений
    void adjustWeightsForConstraints(vector<double>& weights) {
        // Простая корректировка - нормализация с учетом ограничений
        double sum = 0.0;
        for (int i = 0; i < n; i++) {
            if (weights[i] < constraints[i][0]) weights[i] = constraints[i][0];
            if (weights[i] > constraints[i][1]) weights[i] = constraints[i][1];
            sum += weights[i];
        }
        
        // Нормализуем до 100%
        for (int i = 0; i < n; i++) {
            weights[i] /= sum;
        }
    }

    // Расчет риска портфеля
    double calculatePortfolioRisk(const vector<double>& weights) {
        double risk = 0.0;
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                risk += weights[i] * weights[j] * covarianceMatrix[i][j];
            }
        }
        return sqrt(risk);
    }

    // Расчет доходности портфеля
    double calculatePortfolioReturn(const vector<double>& weights) {
        return dotProduct(weights, expectedReturns);
    }

    // Расчет коэффициента Шарпа
    double calculateSharpeRatio(double portfolioReturn, double portfolioRisk) {
        return (portfolioReturn - riskFreeRate) / portfolioRisk;
    }
};

// Функция для ввода данных с проверкой
void inputData(vector<double>& returns, vector<vector<double>>& covMatrix, 
               double& riskFree, vector<vector<double>>& constraints) {
    int n = 3;
    
    cout << "Введите ожидаемые доходности для 3 активов (от -1 до 10):" << endl;
    returns.resize(n);
    for (int i = 0; i < n; i++) {
        cout << "Актив " << i + 1 << ": ";
        cin >> returns[i];
        while (returns[i] < -1.0 || returns[i] > 10.0) {
            cout << "Ошибка! Введите значение от -1 до 10: ";
            cin >> returns[i];
        }
    }
    
    cout << "\nВведите ковариационную матрицу 3x3:" << endl;
    covMatrix.resize(n, vector<double>(n));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            cout << "Элемент [" << i << "][" << j << "]: ";
            cin >> covMatrix[i][j];
        }
    }
    
    cout << "\nВведите безрисковую ставку (от 0 до 1): ";
    cin >> riskFree;
    while (riskFree < 0.0 || riskFree > 1.0) {
        cout << "Ошибка! Введите значение от 0 до 1: ";
        cin >> riskFree;
    }
    
    cout << "\nВведите ограничения для долей активов (мин и макс для каждого):" << endl;
    constraints.resize(n, vector<double>(2));
    for (int i = 0; i < n; i++) {
        cout << "Актив " << i + 1 << " (мин макс): ";
        cin >> constraints[i][0] >> constraints[i][1];
        while (constraints[i][0] < 0.0 || constraints[i][1] > 1.0 || constraints[i][0] > constraints[i][1]) {
            cout << "Ошибка! Мин >= 0, макс <= 1, мин <= макс. Введите снова: ";
            cin >> constraints[i][0] >> constraints[i][1];
        }
    }
}

int main() {
    vector<double> expectedReturns;
    vector<vector<double>> covarianceMatrix;
    double riskFreeRate;
    vector<vector<double>> constraints;
    
    // Ввод данных
    inputData(expectedReturns, covarianceMatrix, riskFreeRate, constraints);
    
    // Создание оптимизатора
    PortfolioOptimizer optimizer(expectedReturns, covarianceMatrix, riskFreeRate, constraints);
    
    // Расчет оптимальных весов
    vector<double> weights = optimizer.calculateOptimalWeights();
    
    // Расчет характеристик портфеля
    double portfolioReturn = optimizer.calculatePortfolioReturn(weights);
    double portfolioRisk = optimizer.calculatePortfolioRisk(weights);
    double sharpeRatio = optimizer.calculateSharpeRatio(portfolioReturn, portfolioRisk);
    
    // Вывод результатов
    cout << fixed << setprecision(1);
    cout << "\n=== РЕЗУЛЬТАТЫ ОПТИМИЗАЦИИ ПОРТФЕЛЯ ===" << endl;
    cout << "Оптимальные доли: [";
    for (size_t i = 0; i < weights.size(); i++) {
        cout << round(weights[i] * 100);
        if (i < weights.size() - 1) cout << "%, ";
        else cout << "%]";
    }
    cout << endl;
    
    cout << "Ожидаемая доходность: " << portfolioReturn * 100 << "%" << endl;
    cout << "Риск портфеля: " << portfolioRisk * 100 << "%" << endl;
    cout << fixed << setprecision(2);
    cout << "Коэффициент Шарпа: " << sharpeRatio << endl;
    
    return 0;
}
