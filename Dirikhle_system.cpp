// Dirikhle_system.cpp : This file contains the 'main' function. Program execution begins and ends there.
//


#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <limits>
#include <iomanip>

const double eps = std::numeric_limits<double>::epsilon();

void writeFinalTable(int n, int m, std::vector<std::vector<double>> v, double a, double b, double c, double d, int S, int S_2, double eps, double eps_max, double eps_max_2, double error_max, double r_norm, double w, double accuracy) {
    std::ofstream outfile("out.txt", std::ofstream::app);
    outfile << "Таблица № " << 1 << std::endl;
    for (size_t j = 0; j < m + 3; j++) {
        if (j == 0) {
            outfile << std::setw(5) << "|" << std::setw(7) << "i" << "|";
            for (int i = 0; i < n + 1; i++) {
                outfile << std::setw(7) << i << "|";
            }

            outfile << std::endl;
        } else if(j == 1) {
            outfile << std::setw(4) << "j" << "|" << std::setw(7) << "Y/X" << "|";
            for (int i = 0; i < n + 1; i++) {
                outfile <<  std::setw(7) << std::setprecision(5) << a + i* (b-a)/n << "|";
            }
            outfile << std::endl;
        } else {
            outfile << std::setw(4) << j - 2 << "|" << std::setw(7) << std::setprecision(4) << c + (j - 2) * (d - c) / m << "|";
            for (int i = 0; i < n + 1; i++) {
                outfile << std::setw(7) << std::setprecision(5) << v[i][j - 2] << "|";
            }
            outfile << std::endl;
        }
    }
    outfile << "\nТаблица № " << 7 << std::endl;
    outfile << std::endl << "Параметр W: " << w << std::endl;
    outfile << std::endl << "Число шагов: " << S << std::endl;
    outfile << std::endl << "Число шагов с двойным шагом: " << S_2 << std::endl;
    outfile << std::endl << "Точность метода: " << eps << std::endl;
    outfile << std::endl << "Достигнутая точность: " << eps_max << std::endl;
    outfile << std::endl << "Достигнутая точность с двойным шагом: " << eps_max_2 << std::endl;
    //outfile << std::endl << "Погрешность: " << error_max << std::endl;
    outfile << std::endl << "Норма невязки: " << r_norm << std::endl;
    outfile << std::endl << "Точность max|v-v_2|: " << accuracy << std::endl;
    outfile.close();
}

void writeTable(int n, int m, std::vector<std::vector<double>> v, double a, double b, double c, double d, int tableNumber) {
    std::ofstream outfile("out.txt", std::ofstream::app);
    outfile << "Таблица № " << tableNumber << std::endl;
    for (size_t j = 0; j < m + 3; j++) {
        if (j == 0) {
            outfile << std::setw(5) << "|" << std::setw(7) << "i" << "|";
            for (int i = 0; i < n + 1; i++) {
                outfile << std::setw(7) << i << "|";
            }

            outfile << std::endl;
        } else if (j == 1) {
            outfile << std::setw(4) << "j" << "|" << std::setw(7) << "Y/X" << "|";
            for (int i = 0; i < n + 1; i++) {
                outfile << std::setw(7) << std::setprecision(5) << a + i * (b - a) / n << "|";
            }
            outfile << std::endl;
        } else {
            outfile << std::setw(4) << j - 2 << "|" << std::setw(7) << std::setprecision(4) << c + (j - 2) * (d - c) / m << "|";
            for (int i = 0; i < n + 1; i++) {
                outfile << std::setw(7) << std::setprecision(5) << v[i][j - 2] << "|";
            }
            outfile << std::endl;
        }
    }
    outfile << std::endl;
    outfile.close();
}

void writeHeader(int n, int m) {
    std::ofstream outfile("out.txt");

    outfile << "\n\nПрименение итерационного метода для решения разностных схем\nна примере задачи Дирихле для уравнения Пуассона\nКоманда №7." << std::endl;
    outfile << "\nu = exp(sin(pi * x * y) * sin(pi * x * y)). Сетка "<<n<<"x"<<m<<". Границы: a=0, b=2, c=0, d=1\n" << std::endl;
}

double equalsZero(double num){
    return (num < eps && num > -eps) ? 0.0 : num;
}

double w_optimal(double a, double b, double c, double d, double h1, double h2) {
    double l1 = b - a;
    double l2 = d - c;

    double lam_min = (2 * h2 * h2 * sin(M_PI * h1 / (2 * l1)) * sin(M_PI * h1 / (2 * l1)) +
        2 * h1 * h1 * sin(M_PI * h2 / (2 * l2)) * sin(M_PI * h2 / (2 * l2))) / (h1 * h1 + h2 * h2);
    return 2 / (1 + sqrt(lam_min * (2 - lam_min)));
}

template <typename Function>
void solve(std::vector<std::vector<double>>& v, Function func, const int n, const int m, double a, double b, double c, double d, int Nmax, int& S, double& eps, double& eps_max, double& error_max) {
    double w = w_optimal(a, b, c, d, (b - a) / n, (d - c) / m);
    int i, j; //индексы
    double a2, k2, h2; // ненулевые элементы матрицы (-A)
    double v_old; // старое значение преобразуемой компоненты вектора v
    double v_new; // новое значение преобразуемой компоненты вектора v
    bool flag = false; // условие остановки
    double eps_cur = 0; // для подсчета текущего значения прироста
    h2 = -(double(n) / (b - a)) * (double(n) / (b - a));
    k2 = -(double(m) / (d - c)) * (double(m) / (d - c));
    a2 = -2 * (h2 + k2);
    for (i = 0; i < n + 1; i++) {
        for (j = 0; j < m + 1; j++) {
            v[i].push_back(0);
        }
    }

    auto m1 = [](double x, double y) {
        return equalsZero(sin(M_PI * y) * sin(M_PI * y));
    };
    auto m2 = [](double x, double y) {
        return equalsZero(sin(M_PI * 2 * y) * sin(M_PI * 2 * y));
    };
    auto m3 = [](double x, double y) {
        return equalsZero(sin(M_PI * x) * sin(M_PI * x));
    };
    auto m4 = [](double x, double y) {
        return equalsZero(sin(M_PI * 2 * x) * sin(M_PI * 2 * x));
    };

    for (j = 0; j < m + 1; j++) {
        v[0][j] = m1(a, c + (d - c) / m * j);
    }
    for (j = 0; j < m + 1; j++) {
        v[n][j] = m2(b, c + (d - c) / m * j);
    }
    for (i = 0; i < n + 1; i++) {
        v[i][0] = m3(a + i * (b - a) / n, c);
    }
    for (i = 0; i < n + 1; i++) {
        v[i][m] = m4(a + i * (b - a) / n, d);
    }

    writeTable(n, m, v, a, b, c, d, 2);

    while (!flag) {
        eps_max = 0;
        for (j = 1; j < m; j++)
            for (i = 1; i < n; i++) {
                v_old = v[i][j];
                v_new = -w * (h2 * (v[i + 1][j] + v[i - 1][j]) + k2 * (v[i][j + 1] + v[i][j - 1]));
                double p = a + i * (b - a) / n;
                double e = c + j * (d - c) / m;
                double func_t = func(p, e);
                v_new = v_new + (1 - w) * a2 * v[i][j] + w * func_t;
                v_new = v_new / a2;
                eps_cur = fabs(v_old - v_new);
                if (eps_cur > eps_max) { eps_max = eps_cur; };
                v[i][j] = v_new;
            }
        if (S < 2) {
            writeTable(n, m, v, a, b, c, d, S * 2 + 4);
        }
        S = S + 1;
        if ((eps_max < eps) or (S >= Nmax)) { flag = true; }
    }
}

int main()
{ 
    auto func = [](double x, double y) {
        return abs(x * x - 2 * y);
    };

    auto u = [](double x, double y) {
        return exp(sin(M_PI * x * y) * sin(M_PI * x * y));
    };

    int Nmax = 10000; // максимальное число итераций (не менее 1)
    int S = 0; // счетчик итераций
    int S_2 = 0; // счетчик итераций
    double eps = 1e-8; // минимально допустимый прирост
    double eps_max = 0; // текущее значение прироста
    double eps_max_2 = 0; // текущее значение прироста
    double eps_cur = 0; // для подсчета текущего значения прироста
    double error_max = 0; // для подсчета текущего значения прироста
    double accuracy = 1000; // для подсчета текущего значения прироста
    double a2, k2, h2; // ненулевые элементы матрицы (-A)


    const int n = 10, m = 10; //размерность сетки
    writeHeader(n,m);

    std::vector<std::vector<double>> v(n + 1); // искомый вектор v
    std::vector<std::vector<double>> v_2(2 * n + 1); // искомый вектор v с половинным шагом
    std::vector<double> r((n - 1) * (m-1)); // невязка
    double a = 0, b = 2, c = 0, d = 1; // границы области определения уравнения
    double w = w_optimal(a, b, c, d, (b - a) / n, (d - c) / m);
    h2 = -(double(n) / (b - a)) * (double(n) / (b - a));
    k2 = -(double(m) / (d - c)) * (double(m) / (d - c));
    a2 = -2 * (h2 + k2);
    bool flag = false;

    int i, j; //индексы
    double v_old; // старое значение преобразуемой компоненты вектора v
    double v_new; // новое значение преобразуемой компоненты вектора v

    solve(v_2, func, 2 * n, 2 * m, a, b, c, d, Nmax, S, eps, eps_max, error_max);
    S_2 = S;
    S = 0;
    eps_max_2 = eps_max;
    eps_max = 0;
    eps_cur = 0;
    error_max = 0;
    solve(v, func, n, m, a, b, c, d, Nmax, S, eps, eps_max, error_max);


    for (j = 1; j < m; j++)
        for (i = 1; i < n; i++) {
            double k_1 = h2 * (v[i + 1][j] * ((i + 1 == n) ? 0 : 1)
                + v[i - 1][j] * ((i - 1 == 0) ? 0 : 1));
            double k_2 = k2 * (v[i][j + 1] * ((j + 1 == m) ? 0 : 1)
                + v[i][j - 1] * ((j - 1 == 0) ? 0 : 1));
            double k_3 = v[i][j] * a2;
            v_new = (k_3 + k_1 + k_2);
            double p = a + i * (b - a) / n;
            double e = c + j * (d - c) / m;
            double func_t = func(p, e);
            double f_t = -(h2 * (v[i + 1][j] * ((i + 1 == n) ? 1 : 0)
                + v[i - 1][j] * ((i - 1 == 0) ? 1 : 0)) +
                k2 * (v[i][j + 1] * ((j + 1 == m) ? 1 : 0)
                    + v[i][j - 1] * ((j - 1 == 0) ? 1 : 0))) + func_t;
            r[(j - 1) * (m - 1) + (i - 1)] = v_new - f_t;
        }

    double r_norm = 0;
    for (j = 0; j < r.size(); j++)
        r_norm += r[j] * r[j];
    r_norm = sqrt(r_norm);
    std::ofstream outfile("test.dat");

    for (j = 1; j < m; j++)
        for (i = 1; i < n; i++) {
            if (abs(v[i][j] - v_2[i][j]) < accuracy) {
                accuracy = abs(v[i][j] - v_2[i][j]);
            }
        }

    for (i = 0; i < n + 1; i++) 
        for (j = 0; j < m+1; j++){        
            v_new = v[i][j];
            double p = a + i * (b - a) / n;
            double e = c + j * (d - c) / m;
            outfile << p << "\t" << e << "\t" << v_new << "\n";
       }
    //writeTable(n, m, vec_u, a, b, c, d, 8);
    writeFinalTable(n, m, v, a, b, c, d, S, S_2, eps, eps_max, eps_max_2, error_max, r_norm, w, accuracy);
    outfile.close();
    system("python show_plot.py");
}

//return -((4 * M_PI * M_PI * y * y * exp(sin(M_PI * x * y) * sin(M_PI * x * y)) * cos(M_PI * x * y) * cos(M_PI * x * y) * sin(M_PI * x * y) * sin(M_PI * x * y) -
//    2 * M_PI * M_PI * y * y * exp(sin(M_PI * x * y) * sin(M_PI * x * y)) * sin(M_PI * x * y) * sin(M_PI * x * y) +
//    2 * M_PI * M_PI * y * y * exp(sin(M_PI * x * y) * sin(M_PI * x * y)) * cos(M_PI * x * y) * cos(M_PI * x * y)) +
//    (4 * M_PI * M_PI * x * x * exp(sin(M_PI * x * y) * sin(M_PI * x * y)) * cos(M_PI * x * y) * cos(M_PI * x * y) * sin(M_PI * x * y) * sin(M_PI * x * y) -
//        2 * M_PI * M_PI * x * x * exp(sin(M_PI * x * y) * sin(M_PI * x * y)) * sin(M_PI * x * y) * sin(M_PI * x * y) +
//        2 * M_PI * M_PI * x * x * exp(sin(M_PI * x * y) * sin(M_PI * x * y)) * cos(M_PI * x * y) * cos(M_PI * x * y)));