#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>

using namespace std;
const double gamma_val = 1.4;

int main() {
    vector<double> beta, a, b;
    for (double angle = M_PI / 2; angle >= 0; angle -= M_PI / 180) {
        beta.push_back(angle);
    }

    ofstream datafile("data.txt");
    for (double M1 = 1; M1 <= 5; M1 += 0.02) {
        double max_theta = 0.0;
        double max_beta = 0.0;

        for (size_t i = 0; i < beta.size(); ++i) {
            double theta = atan(2.0 / tan(beta[i]) * ((M1 * M1 * sin(beta[i]) * sin(beta[i])) - 1.0) /
                                ((gamma_val + cos(2.0 * beta[i])) * M1 * M1 + 2.0));

            datafile << theta * 180 / M_PI << " " << beta[i] * 180 / M_PI << endl;

        if (theta > max_theta) {
                max_theta = theta;
                max_beta = beta[i];
            }
        }

        a.push_back(max_theta);
        b.push_back(max_beta);
    }
    datafile.close();
    ofstream data("max.txt");
    for (size_t i = 0; i < a.size(); ++i) {
        if (i == 0) { 
            data << "0 90" << endl;
        } else {
            data << a[i] * 180 / M_PI << " " << b[i] * 180 / M_PI << endl;
        }
    }
    data.close();
 
    FILE *gnuplotPipe = popen("gnuplot -persist", "w");
    fprintf(gnuplotPipe, "set title 'Theta vs Beta'\n");
    fprintf(gnuplotPipe, "set xlabel 'Theta'\n");
    fprintf(gnuplotPipe, "set ylabel 'Beta'\n");
    fprintf(gnuplotPipe, "set xrange [0:]\n");
    fprintf(gnuplotPipe, "set yrange [0:]\n");
    fprintf(gnuplotPipe, "plot 'data.txt' with lines lc rgb 'blue' notitle, 'max.txt' with lines lt rgb 'red' lw 3 title 'theta max'\n");
    fflush(gnuplotPipe);

    return 0;
}
