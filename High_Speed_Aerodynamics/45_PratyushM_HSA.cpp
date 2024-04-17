#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#define gamma 1.4
using namespace std;


// Function to compute the desired equations
double computeFunction(double M, double X, int type) {

    if (type == 1){
        return pow(pow(gamma + 1, 2) * M * M / (4 * gamma * M * M - 2 * (gamma - 1)), gamma / (gamma - 1)) * ((1 - gamma + 2 * gamma * M * M) / (gamma + 1)) - X;
    }
    if (type == 2){
        return pow((gamma + 1)/(gamma - 1),0.5) * atan(pow((gamma -1)*(M*M-1)/(gamma + 1),0.5)) - atan(pow(M*M-1,0.5)) - X;
    }

    return 0;
}

// Function to find the root of the equation
double bisection_M( double X, int type) {
    double a = 1.0, b = 5.0;     // Lower & Upper bound of M
    double fa = computeFunction(a, X,type), fb = computeFunction(b, X,type);

    double c = 0.0;
    while ((b - a) / 2 > 1e-6) {
        c = (a + b) / 2;
        double fc = computeFunction(c, X,type);

        if (fc == 0.0)
            break;
        else if (fa * fc < 0)
            b = c;
        else
            a = c;
    }
    return c;
}

double calculate(double P) {
    ofstream datafile;
    datafile.open("HSA.csv", ios::app);
    double p3 = 0.1, theta = 30, beta, M1;
    double pO4_3 = P / p3;

    double M3 = bisection_M(pO4_3, 1);
  
    float v3 = pow((gamma + 1) / (gamma - 1), 0.5) * atan(pow((gamma - 1) * (M3 * M3 - 1) / (gamma + 1), 0.5)) - atan(pow(M3 * M3 - 1, 0.5));

    float v2 = v3 - 30 * M_PI / 180;
  
    double M2 = bisection_M(v2, 2);
  
    // Iterative method to find M1 and beta
    for (double M_i = 1; M_i <= 5; M_i += 0.1) {
        for (double beta_j = 0; beta_j <= 90; beta_j += 0.01) {

            double theta_iterate = (atan(((M_i * M_i * sin(2 * beta_j * M_PI / 180)) - (2 * cos(beta_j * M_PI / 180)) / sin(beta_j * M_PI / 180)) / ((M_i * M_i * (gamma + cos(2 * beta_j * M_PI / 180))) + 2))) * (180 / M_PI);
            
            double diff_theta = abs((theta / 2) - theta_iterate);
            if (diff_theta < 0.01) {
                beta = beta_j;
                break;
            }
        }

        double Mn1 = M_i * sin(beta * M_PI / 180);

        double Mn2_2 = M2 * sin((beta - theta / 2) * M_PI / 180);

        double Mn2_1 = pow((Mn1 * Mn1 + (2 / (gamma - 1))) / (((2 * gamma * Mn1 * Mn1) / (gamma - 1)) - 1), 0.5);

        double diff_Mn2 = abs(Mn2_1 - Mn2_2);
        if (diff_Mn2 < 0.01) {
            M1 = M_i;
            cout<<setw(3) << M1 << endl;
            break;
        }
    }
    datafile << pO4_3 << "," << pO4_3  << "," << M3 << "," << v3 << "," << v2 << "," << M2 << "," << M1 << endl;
    datafile.close();
    return 0;
}

int main() {
    system("cls");
    double sum = 0.1;
    ofstream datafile("HSA.csv");
    datafile << "Po4, Po4/P3, M3, v3, v2, M2, M1" << endl;
    datafile.close();
    cout<<"Pressure(Atm)  Mach Number"<<endl;
    for (double M_i = 2.5; M_i < 3.1; M_i += sum) {
        cout <<setw(3)<< M_i << ",            ";
        calculate (M_i);
    }
    return 0;
}