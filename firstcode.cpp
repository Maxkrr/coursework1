#include <iostream>
#include <vector>
#include <iomanip>
#include <cmath>
#include <cstdlib>
#include <ctime>

using namespace std;

// Типы для удобства
typedef vector<double> Weights;
typedef vector<double> Returns;
typedef vector<vector<double>> CovarianceMatrix;

// Функция вычисления доходности портфеля
double calculatePortfolioReturn(const Weights& weights, const Returns& returns) {
    double total_return = 0.0;
    for (size_t i = 0; i < weights.size(); i++) {
        total_return += weights[i] * returns[i];
    }
    return total_return;
}

// Функция вычисления риска портфеля (стандартного отклонения)
double calculatePortfolioRisk(const Weights& weights, const CovarianceMatrix& covariance) {
    double variance = 0.0;
    size_t n = weights.size();
    
    // Вычисляем дисперсию портфеля
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            variance += weights[i] * covariance[i][j] * weights[j];
        }
    }
    
    // Стандартное отклонение = корень из дисперсии
    return sqrt(variance);
}

// Функция проверки ограничений на доли активов
bool isValidPortfolio(const Weights& weights, const Weights& min_limits, const Weights& max_limits) {
    double sum = 0.0;
    
    // Проверяем каждую долю
    for (size_t i = 0; i < weights.size(); i++) {
        // Доля должна быть в пределах [min_limits[i], max_limits[i]]
        if (weights[i] < min_limits[i] || weights[i] > max_limits[i]) {
            return false;
        }
        sum += weights[i];
    }
    
    // Сумма всех долей должна быть равна 1 (100%)
    return abs(sum - 1.0) < 0.000001;
}

// Функция вычисления коэффициента Шарпа
double calculateSharpeRatio(double portfolio_return, double portfolio_risk, double risk_free_rate) {
    return (portfolio_return - risk_free_rate) / portfolio_risk;
}

int main() {
    // Инициализация генератора случайных чисел
    srand(time(0));
    
    cout << "=== ОПТИМИЗАЦИЯ ИНВЕСТИЦИОННОГО ПОРТФЕЛЯ ===" << endl;
    cout << endl;
    
    // Ввод количества активов
    int asset_count;
    cout << "Введите количество активов: ";
    cin >> asset_count;
    
    // Ввод ожидаемых доходностей
    Returns expected_returns(asset_count);
    cout << endl << "Введите ожидаемые доходности активов (от -1 до 10):" << endl;
    for (int i = 0; i < asset_count; i++) {
        do {
            cout << "Доходность актива " << (i + 1) << ": ";
            cin >> expected_returns[i];
            if (expected_returns[i] < -1 || expected_returns[i] > 10) {
                cout << "  Ошибка: доходность должна быть от -1 до 10" << endl;
            }
        } while (expected_returns[i] < -1 || expected_returns[i] > 10);
    }
    
    // Ввод ковариационной матрицы
    CovarianceMatrix covariance(asset_count, vector<double>(asset_count));
    cout << endl << "Введите ковариационную матрицу (построчно):" << endl;
    for (int i = 0; i < asset_count; i++) {
        for (int j = 0; j < asset_count; j++) {
            cout << "Ковариация [" << (i + 1) << "][" << (j + 1) << "]: ";
            cin >> covariance[i][j];
        }
    }
    
    // Ввод безрисковой ставки
    double risk_free_rate;
    cout << endl << "Введите безрисковую ставку (от 0 до 1): ";
    do {
        cin >> risk_free_rate;
        if (risk_free_rate < 0 || risk_free_rate > 1) {
            cout << "  Ошибка: безрисковая ставка должна быть от 0 до 1: ";
        }
    } while (risk_free_rate < 0 || risk_free_rate > 1);
    
    // Ввод ограничений на доли активов
    Weights min_limits(asset_count);
    Weights max_limits(asset_count);
    
    cout << endl << "Введите минимальные доли активов (от 0 до 1):" << endl;
    for (int i = 0; i < asset_count; i++) {
        do {
            cout << "Минимальная доля актива " << (i + 1) << ": ";
            cin >> min_limits[i];
            if (min_limits[i] < 0 || min_limits[i] > 1) {
                cout << "  Ошибка: минимальная доля должна быть от 0 до 1" << endl;
            }
        } while (min_limits[i] < 0 || min_limits[i] > 1);
    }
    
    cout << endl << "Введите максимальные доли активов (от 0 до 1):" << endl;
    for (int i = 0; i < asset_count; i++) {
        do {
            cout << "Максимальная доля актива " << (i + 1) << ": ";
            cin >> max_limits[i];
            if (max_limits[i] < 0 || max_limits[i] > 1 || max_limits[i] < min_limits[i]) {
                cout << "  Ошибка: максимальная доля должна быть от 0 до 1 и не меньше минимальной" << endl;
            }
        } while (max_limits[i] < 0 || max_limits[i] > 1 || max_limits[i] < min_limits[i]);
    }
    
    // Поиск оптимального портфеля
    cout << endl << "Выполняется поиск оптимального портфеля..." << endl;
    
    Weights best_weights(asset_count);
    double best_sharpe_ratio = -1000000; // Очень маленькое начальное значение
    double best_return = 0.0;
    double best_risk = 0.0;
    
    // Количество попыток найти оптимальный портфель
    const int NUMBER_OF_TRIALS = 10000;
    
    for (int trial = 0; trial < NUMBER_OF_TRIALS; trial++) {
        // Генерируем случайные веса в пределах ограничений
        Weights current_weights(asset_count);
        double sum_of_weights = 0.0;
        
        for (int i = 0; i < asset_count; i++) {
            double range = max_limits[i] - min_limits[i];
            current_weights[i] = min_limits[i] + (rand() / (double)RAND_MAX) * range;
            sum_of_weights += current_weights[i];
        }
        
        // Нормализуем веса, чтобы их сумма была равна 1
        for (int i = 0; i < asset_count; i++) {
            current_weights[i] /= sum_of_weights;
        }
        
        // Проверяем, удовлетворяет ли портфель ограничениям
        if (!isValidPortfolio(current_weights, min_limits, max_limits)) {
            continue; // Пропускаем этот портфель
        }
        
        // Вычисляем показатели портфеля
        double current_return = calculatePortfolioReturn(current_weights, expected_returns);
        double current_risk = calculatePortfolioRisk(current_weights, covariance);
        
        // Избегаем деления на ноль
        if (current_risk == 0) {
            continue;
        }
        
        double current_sharpe_ratio = calculateSharpeRatio(current_return, current_risk, risk_free_rate);
        
        // Обновляем лучший портфель, если нашли лучше
        if (current_sharpe_ratio > best_sharpe_ratio) {
            best_sharpe_ratio = current_sharpe_ratio;
            best_weights = current_weights;
            best_return = current_return;
            best_risk = current_risk;
        }
    }
    
    // Вывод результатов
    cout << endl;
    cout << "=== РЕЗУЛЬТАТЫ ОПТИМИЗАЦИИ ===" << endl;
    cout << endl;
    
    cout << "Оптимальные доли: [";
    for (int i = 0; i < asset_count; i++) {
        cout << fixed << setprecision(2) << best_weights[i] * 100 << "%";
        if (i < asset_count - 1) {
            cout << ", ";
        }
    }
    cout << "]" << endl;
    
    cout << "Ожидаемая доходность: " << best_return * 100 << "%" << endl;
    cout << "Риск портфеля: " << best_risk * 100 << "%" << endl;
    cout << "Коэффициент Шарпа: " << fixed << setprecision(3) << best_sharpe_ratio << endl;
    
    return 0;
}
