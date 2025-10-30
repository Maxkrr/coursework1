#include <iostream>
#include <cmath>
using namespace std;

int main() {
    double rf = 0.01; // безрисковая ставка
    double returns[3] = {0.12, 0.08, 0.15};
    double cov[3][3] = {{0.04,0.01,0.02}, {0.01,0.03,0.01}, {0.02,0.01,0.06}};
    
    double best_sharpe = -1e9, best_w[3];

    for(double w1 = 0; w1 <= 1; w1 += 0.1)
    for(double w2 = 0; w2 <= 1 - w1; w2 += 0.1) {
        double w3 = 1 - w1 - w2;
        double er = w1*returns[0] + w2*returns[1] + w3*returns[2];
        
        double risk = 0;
        double w[3] = {w1, w2, w3};
        for(int i = 0; i < 3; i++)
        for(int j = 0; j < 3; j++)
            risk += w[i]*w[j]*cov[i][j];
        risk = sqrt(risk);
        
        double sharpe = (er - rf) / risk;
        if(sharpe > best_sharpe) {
            best_sharpe = sharpe;
            best_w[0] = w1; best_w[1] = w2; best_w[2] = w3;
        }
    }

    cout << "Оптимальные доли: ";
    cout << best_w[0]*100 << "% " << best_w[1]*100 << "% " << best_w[2]*100 << "%\n";
    
    double er = best_w[0]*returns[0] + best_w[1]*returns[1] + best_w[2]*returns[2];
    cout << "Ожидаемая доходность: " << er*100 << "%\n";
    cout << "Коэффициент Шарпа: " << best_sharpe << endl;
    
    return 0;
}
