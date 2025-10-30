#include <iostream>
#include <cmath>
using namespace std;

// Количество активов в портфеле
const int N = 3;

// Ожидаемая доходность каждого актива (процент/год)
double mu[N] = {0.12, 0.08, 0.15};

// Ковариационная матрица рисков (размерность NxN)
// Показывает, как изменяются риски совместно для разных активов
double cov[N][N] = {
    {0.04, 0.01, 0.02},
    {0.01, 0.03, 0.01},
    {0.02, 0.01, 0.06}
};

// Функция для вычисления ожидаемой доходности портфеля
double portfolioReturn(double weights[N]) {
    double result = 0;
    // Суммируем вклад каждого актива
    for (int i = 0; i < N; ++i)
        result += weights[i] * mu[i];
    return result;
}

// Функция для вычисления риска портфеля (стандартного отклонения)
double portfolioRisk(double weights[N]) {
    double risk = 0;
    // Проходим по всем парам активов и учитываем их ковариацию
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j)
            risk += weights[i] * cov[i][j] * weights[j];
    return sqrt(risk); // корень из дисперсии — это стандартное отклонение
}

// Главная функция программы
int main() {
    // Вес каждого актива в портфеле — вручную (пример: 45%, 30%, 25%)
    double weights[N] = {0.45, 0.30, 0.25};

    // Выводим доли активов
    cout << "Оптимальные доли: ";
    for (int i = 0; i < N; ++i)
        cout << weights[i] * 100 << "% ";
    cout << endl;

    // Вычисляем и выводим ожидаемую доходность портфеля в процентах
    cout << "Ожидаемая доходность: " << portfolioReturn(weights) * 100 << "%\n";

    // Вычисляем и выводим риск портфеля в процентах (стандартное отклонение)
    cout << "Риск портфеля: " << portfolioRisk(weights) * 100 << "%\n";

    // Коэффициент Шарпа — показывает соотношение доходности к риску
    double riskFree = 0.02; // безрисковая ставка, например 2%
    double sharpe = (portfolioReturn(weights) - riskFree) / portfolioRisk(weights);
    cout << "Коэффициент Шарпа: " << sharpe << endl;

    return 0;
}
