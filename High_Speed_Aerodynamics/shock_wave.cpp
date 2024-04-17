#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;
int main() {
    ofstream file("data.csv");
    double gamma = 1.4;
    double R = 0.287;
    double Cp = 1.005;
    double M_1 = 1.01;
    double M_2 = sqrt((1 + (( gamma - 1) / 2) * pow(M_1, 2)) / ((gamma * pow(M_1,2)) - ((gamma - 1) / 2)));
    double rho = (gamma + 1)*pow(M_1,2) / (2  + (gamma - 1)* pow(M_1,2));
    double P = (1 + ((2 * gamma) / (gamma + 1)) * (pow(M_1,2) - 1)) ;
    double T = P/rho;
    double S = Cp*log((1+(2*gamma/(gamma+1))*(pow(M_1,2)-1))*((2+(gamma -1 )*pow(M_1,2))/((gamma+1)*pow(M_1,2)))) 
    - R*log(1+(2*gamma/(gamma+1))*(pow(M_1,2)-1));
    double Po = exp(-S/R);
    file<<M_1<<","<<M_2*10<<","<<rho<<","<<P<<","<<T<<","<<Po*10<<endl;
    for (int i =0;i <=20;i++){
        M_1 = 1.1 + i*0.2;
        M_2 = sqrt((1 + (( gamma - 1) / 2) * pow(M_1, 2)) / ((gamma * pow(M_1,2)) - ((gamma - 1) / 2)));
        rho = (gamma + 1)*pow(M_1,2) / (2  + (gamma - 1)* pow(M_1,2));
        P = (1 + ((2 * gamma) / (gamma + 1)) * (pow(M_1,2) - 1)) ;
        T = P/rho;
        S = Cp*log((1+(2*gamma/(gamma+1))*(pow(M_1,2)-1))*((2+(gamma -1 )*pow(M_1,2))/((gamma+1)*pow(M_1,2)))) 
        - R*log(1+(2*gamma/(gamma+1))*(pow(M_1,2)-1));
        Po = exp(-S/R);
        file<<M_1<<","<<M_2*10<<","<<rho<<","<<P<<","<<T<<","<<Po*10<<endl;
    }
    file.close();
    cout<<"Completed"<<endl;

    return 0;
}
