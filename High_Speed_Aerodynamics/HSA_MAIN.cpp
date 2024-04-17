#include <iostream>
#include <cmath>
#include <fstream>

#define gamma         1.4
using namespace std;

void data(float M, int type){
    float T1 = 1 + ((gamma - 1) / 2) * M * M;
    float P1 = pow(T1, gamma / (gamma - 1));
    float M2 = pow(T1/(gamma*M*M - (gamma - 1)/2),0.5);
    
    float A = 1 / M * pow((2 / (gamma + 1)) * T1, 0.5 * (gamma + 1) / (gamma - 1));
    float v = pow((gamma + 1)/(gamma - 1),0.5) * atan(pow((gamma -1)*(M*M-1)/(gamma + 1),0.5)) - atan(pow(M*M-1,0.5));
    float u = asin(1/M);
    float T2 = 1 + ((gamma - 1) / 2) * M2 * M2;
    float P2 = pow(T2, gamma / (gamma - 1));

    float P2_1 = 1 + (2*gamma/(gamma + 1))*(M*M - 1);
    float rho2_1 = (gamma + 1)*M*M/(2 + (gamma - 1)*M*M);
    float rho1 = pow(T1, 1 / (gamma - 1));
    float T2_1 = P2_1/rho2_1;
    float P0 = P2_1*pow(T2/T1,gamma/(gamma - 1));

    //float T0 = X*X*pow(M2/M,2)*T2/T1;
    //float Po1 = ((1+gamma)/(1+gamma*M*M))*pow((2+(gamma-1)*M*M)/(1+gamma),gamma/(gamma-1));
    //float Po2 = ((1+gamma)/(1+gamma*M2*M2))*pow((2+(gamma-1)*M2*M2)/(1+gamma),gamma/(gamma-1));

    if (type == 1) {
        ofstream myfile;
        myfile.open ("isentropic.csv", ios::app); // open in append mode
        myfile << M << "," << 1/T1 << "," << 1/P1 << "," << 1/rho1 << "," << A << "," << u*180/M_PI << "," << v*180/M_PI << "\n";
        myfile.close();
    }
    
    if (type == 2){
        ofstream myfile;
        myfile.open ("normal.csv", ios::app); // open in append mode
        myfile << M << "," <<M2 << "," << T2_1 << "," << P2_1 << "," << rho2_1 <<","<< P0<<","<<P2_1*P2<< "\n";
        myfile.close();
    }
    
    
}

int main() {
    
    remove("isentropic.csv"); // delete the file if it exists
    remove("normal.csv");
    remove("angle.csv");

    ofstream myfile;
    myfile.open ("isentropic.csv", ios::app);
    myfile << "M,T/To,P/Po,rho/rho0,A/A*,u,v" << "\n";
    myfile.close();

    myfile.open ("normal.csv", ios::app);
    myfile << "M1,M2,T2/T1,P2/P1,rho2/rho1,P02/P01,P0/P0*" << "\n";
    myfile.close();

    double sum = 0.02;
    for (double i = sum; i < 4.8; i+=sum) {
        data(i,1);
    }

    for (double i = 1 + sum; i <= 5.8; i+=sum) {
        data(i,2);
    }


    return 0;
}