#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
const int N = 50, V = 11, T = 2000;
const double dv = 1.0, dx = 10.0, dt = 1.0;

double get_V(int id) { return dv * (id - V / 2) - dv / 2; }

void set_initial_data(double other[N][V]) {

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < V; j++) {
            auto v = get_V(j);
            auto x = (i - N / 2);
            other[i][j] = std::exp(-v * v / 2) * std::exp( - x*x / 2);
        }
    }
}

int main() {

    char filename[100];
    double ar_current[N][V]{0};
    double ar_next[N][V]{0};
    double ar_concentration[N]{0};
    set_initial_data(ar_current);
    std::fstream file_out;
    for (int k = 0; k < T; k++) {
        for (int i = 0; i < N; i++) {

            for (int j = 0; j < V; j++) {
                auto v = get_V(j);
                double gamma = v * dt / dx;

                if (v > 0) {
                    if (i == 0) {
                        ar_next[0][j] = ar_current[N-1][j];
                    } else {
                        ar_next[i][j] = ar_current[i][j] - gamma * (ar_current[i][j] - ar_current[i - 1][j]);
                    }
                } else {
                    if (i == N - 1) {
                        ar_next[N - 1][j] = ar_current[0][j];
                    } else {
                        ar_next[i][j] = ar_current[i][j] - gamma * (ar_current[i + 1][j] - ar_current[i][j]);
                    }
                }
            }
        }
        std::swap(ar_current, ar_next);

        sprintf(filename, "data/out_%03d.txt", k);
        file_out.open(filename, std::fstream::out);
        for (int i = 0; i < N; i++) {
            ar_concentration[i] = 0;
            for (int j = 0; j < V; j++) {
                ar_concentration[i] += ar_current[i][j];
            }
            file_out << i << '\t' << ar_concentration[i] << std::endl;
        }
        file_out.close();
    }
}
