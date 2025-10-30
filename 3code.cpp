#include <iostream>
#include <vector>
#include <cmath>
using namespace std;

int main() {
    vector<double> r(3); // доходности
    vector<vector<double>> cov(3, vector<double>(3)); // ковариационная матрица
    double rf; // безрисковая ставка
    
    // Ввод данных
    cout << "Введите доходности 3 активов: ";
    cin >> r[0] >> r[1] >> r[2];
    
    cout << "Введите ковариационную матрицу (9 чисел): ";
    for(int i=0; i<3; i++)
        for(int j=0; j<3; j++)
            cin >> cov[i][j];
    
    cout << "Введите безрисковую ставку: ";
    cin >> rf;
    
    // Вычисление обратной матрицы 3x3
    double det = cov[0][0]*(cov[1][1]*cov[2][2]-cov[1][2]*cov[2][1]) 
               - cov[0][1]*(cov[1][0]*cov[2][2]-cov[1][2]*cov[2][0]) 
               + cov[0][2]*(cov[1][0]*cov[2][1]-cov[1][1]*cov[2][0]);
    
    vector<vector<double>> inv(3, vector<double>(3));
    inv[0][0] = (cov[1][1]*cov[2][2]-cov[1][2]*cov[2][1])/det;
    inv[0][1] = (cov[0][2]*cov[2][1]-cov[0][1]*cov[2][2])/det;
    inv[0][2] = (cov[0][1]*cov[1][2]-cov[0][2]*cov[1][1])/det;
    inv[1][0] = (cov[1][2]*cov[2][0]-cov[1][0]*cov[2][2])/det;
    inv[1][1] = (cov[0][0]*cov[2][2]-cov[0][2]*cov[2][0])/det;
    inv[1][2] = (cov[0][2]*cov[1][0]-cov[0][0]*cov[1][2])/det;
    inv[2][0] = (cov[1][0]*cov[2][1]-cov[1][1]*cov[2][0])/det;
    inv[2][1] = (cov[0][1]*cov[2][0]-cov[0][0]*cov[2][1])/det;
    inv[2][2] = (cov[0][0]*cov[1][1]-cov[0][1]*cov[1][0])/det;
    
    // Вычисление оптимальных весов
    vector<double> w(3);
    double sum = 0;
    for(int i=0; i<3; i++) {
        w[i] = inv[i][0]*r[0] + inv[i][1]*r[1] + inv[i][2]*r[2];
        sum += w[i];
    }
    
    for(int i=0; i<3; i++) w[i] /= sum;
    
    // Расчет доходности и риска
    double ret = w[0]*r[0] + w[1]*r[1] + w[2]*r[2];
    
    double risk = 0;
    for(int i=0; i<3; i++)
        for(int j=0; j<3; j++)
            risk += w[i]*w[j]*cov[i][j];
    risk = sqrt(risk);
    
    // Коэффициент Шарпа
    double sharpe = (ret - rf) / risk;
    
    // Вывод результатов
    cout << "Оптимальные доли: [" << w[0]*100 << "%, " << w[1]*100 << "%, " << w[2]*100 << "%]" << endl;
    cout << "Ожидаемая доходность: " << ret*100 << "%" << endl;
    cout << "Риск портфеля: " << risk*100 << "%" << endl;
    cout << "Коэффициент Шарпа: " << sharpe << endl;
    
    return 0;
}
