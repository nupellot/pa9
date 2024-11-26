#include <cmath>
#include <fstream>
#include <iostream>

using namespace std;
#define Pi 3.141592653589793
#define delta_t_appr 1e-10
#define delta_t_min 1e-12
#define end_time 1e-3
#define epsilon 1e-5 // X
#define E1 1
#define R1 1e3
#define R2 1e3
#define Rb1 20.0
#define Rb2 1e6
#define L1 1e-3
#define Cb1 2e-12
#define C1 1e-6
#define Id1 1e-12
#define MFt 0.026
#define vector_size 70000
#define eps1 1e-15
#define eps2 1e-5
#define appr 0.0
#define C2 1e-6

double U_R1, U_R2, U_Rb1, U_Rb2, U_Id1, I_E1, I_E2, I_L1, I_Cb1, I_C1, I_R1,
I_R2, I_Rb1, I_Rb2, I_Id1, U_E1, U_E2, U_L1, U_Cb1, U_C1, I_L1_pred, U_C2, I_C2,
U_Cb1_pred, U_C1_pred, I_C2_pred, U_C2_pred, U_L1_pred;


double* gauss(double ** a, double * y, int n);
double norm(double * result, int step, double delta_t, double delta_t_pred, int n);
int check_if_vector_less_than_eps(double * x, int n);
void result_to_csv(double * result, int step, int n);
void get_pred_variables(double * result, int n, int step);
void return_past_results(double * result, int n, int step);
void save_step_results(double * result, int n, int step, double time);
void add_phi(double * x);
void make_appr(double * result, int n);
void fill_y(double * y, 
            double U_R1, double U_R2, double U_Rb1, double U_Rb2, double U_Id1, double U_L1, 
            double I_E1, double I_E2, double I_C2, double I_Cb1, double I_C1, 
            double I_R1, double I_R2, double I_Rb1, double I_Rb2, double I_Id1, double I_L1, 
            double U_E1, double U_E2, double U_C2, double U_Cb1, double U_C1, 
            double U_Cb1_pred, double U_C1_pred, double I_C2_pred, double U_C2_pred, double U_L1_pred,
            double delta_t, int n, double time
);
void fill_A(double ** a, double U_Id1, double delta_t);


// Используется итерационный метод для нахождения решения системы дифференциальных уравнений с адаптивным шагом по времени.
int main() {
    int n = 22;
    int n_iter = 0;
    double time = 0.0;
    double ** a = new double * [n];
    double * y = new double[n], * x;
    double * result = new double[vector_size];
    result[0] = time;
    make_appr(result, n);
    double delta_t = delta_t_appr, delta_t_pred = 0.0, eps = 0.0;
    int flag = 1;
    int step = 1;
    int time_step = 0;
    int iter_step = 0;
    while (time < end_time) {
        get_pred_variables(result, n, step);
        time_step++;
        while (flag) {
            for (int i = 0; i < n; i++) {
                a[i] = new double[n];
                y[i] = 0.0;
                for (int j = 0; j < n; j++)
                    a[i][j] = 0.0;
            }
            iter_step++;
            n_iter++;
            fill_A(a, U_Id1, delta_t);
            fill_y(y, 
                        U_R1, U_R2, U_Rb1, U_Rb2, U_Id1, U_L1, 
                        I_E1, I_E2, I_C2, I_Cb1, I_C1, 
                        I_R1, I_R2, I_Rb1, I_Rb2, I_Id1, I_L1, 
                        U_E1, U_E2, U_C2, U_Cb1, U_C1, 
                        U_Cb1_pred, U_C1_pred, I_C2_pred, U_C2_pred, U_L1_pred,
                        delta_t, n, time
            );
            // cout << "time: " << time << endl;
            // cout << "delta_t: " << delta_t << endl;
            x = gauss(a, y, n);
            add_phi(x);
            if (check_if_vector_less_than_eps(x, n) == n) {
                n_iter = 0;
                flag = 0;
            } else {
                if (n_iter > 7) {
                    delta_t = delta_t / 2;
                    n_iter = 0;
                    return_past_results(result, n, step);
                    if (delta_t < delta_t_min) {
                        cout << "Solution doesn't converge" << endl;
                        exit(0);
                    }
                }
            }
        }
        if (step == 1 && time == 0.0) {
            step = 0;
        }
        save_step_results(result, n, step, time);
        delta_t_pred = delta_t;
        // cout << time << endl;
        if (step > 0) {
            eps = norm(result, step, delta_t, delta_t_pred, n);
            // cout<<eps<<endl;
            if (eps < eps1) {
                time += delta_t;
                delta_t_pred = delta_t;
                delta_t = 2 * delta_t;
                for (int i = 1; i < n + 1; i++) {
                    if (step > 1) {
                        result[(step - 2) * (n + 1) + i] = result[(step - 1) * (n + 1) + i];
                        result[(step - 1) * (n + 1) + i] = result[step * (n + 1) + i];
                    }
                }
            }
            if (eps > eps1 && eps < eps2) {
                time += delta_t;
                delta_t_pred = delta_t;
                for (int i = 1; i < n + 1; i++) {
                    if (step > 2) {
                        result[(step - 2) * (n + 1) + i] = result[(step - 1) * (n + 1) + i];
                        result[(step - 1) * (n + 1) + i] = result[step * (n + 1) + i];
                    }
                }
            }
            if (eps > eps2) {
                time += 0;
                delta_t = delta_t / 2;
                for (int i = 1; i < n + 1; i++)
                    if (step >= 1) {
                        result[step * (n + 1) + i] = result[(step - 1) * (n + 1) + i];
                    }
            }
        } else {
            delta_t = delta_t / 2;
            time += delta_t;
        }
        step++;
        flag = 1;
    }
    cout << "success" << endl;
    cout << "Time steps: " << time_step << endl;
    
    result_to_csv(result, step, n);
    delete[] y;
    delete[] result;
    for (int i = 0; i < n; i++)
        delete[] a[i];
    return 0;
}


// Реализует метод Гаусса для решения системы линейных уравнений.
// Принимает на вход матрицу a и вектор y, возвращает вектор x.
double *gauss(double **a, double *y, int n) {
  const double eps = 1e-12;
  double *x, max;
  int k, index;
  x = new double[n];
  k = 0;

  // Прямой ход метода Гаусса
  while (k < n) {
    max = abs(a[k][k]);
    index = k;

    // Находим максимальный элемент в столбце
    for (int i = k + 1; i < n; i++) {
      if (abs(a[i][k]) > max) {
        max = abs(a[i][k]);
        index = i;
      }
    }

    // Проверка на вырожденность матрицы
    if (max < eps) {
      return 0;
    }

    // Перестановка строк для улучшения численной стабильности
    for (int j = 0; j < n; j++) {
      double temp = a[k][j];
      a[k][j] = a[index][j];
      a[index][j] = temp;
    }
    double temp = y[k];
    y[k] = y[index];
    y[index] = temp;

    // Приведение строки к нормальной форме
    for (int i = k; i < n; i++) {
      double temp = a[i][k];
      if (abs(temp) < eps)
        continue;
      for (int j = 0; j < n; j++)
        a[i][j] = a[i][j] / temp;
      y[i] = y[i] / temp;
      if (i == k)
        continue;
      for (int j = 0; j < n; j++)
        a[i][j] = a[i][j] - a[k][j];
      y[i] = y[i] - y[k];
    }
    k++;
  }

  // Обратный ход метода Гаусса
  for (k = n - 1; k >= 0; k--) {
    x[k] = y[k];
    for (int i = 0; i < k; i++)
      y[i] = y[i] - a[i][k] * x[k];
  }
  return x;
}



// Сохраняет результаты в файл CSV.
void result_to_csv(double * result, int step, int n) {
    ofstream fs;
    fs.open("result.csv");
    fs << "time,U_R1,U_R2,U_Rb1,U_Rb2,U_Id1,UL1,I_E1,I_E2,I_C2,I_Cb1,I_C1,I_R1,I_R2,"
    "I_Rb1, I_Rb2, I_Id1, I_L1, U_E1, U_E2, U_C2, U_Cb1, U_C1 " <<
    endl;
    for (int s = 0; s < step; s++) {
        fs << result[s * (n + 1)] << ",";
        fs << result[s * (n + 1) + 1] << ",";
        fs << result[s * (n + 1) + 2] << ",";
        fs << result[s * (n + 1) + 3] << ",";
        fs << result[s * (n + 1) + 4] << ",";
        fs << result[s * (n + 1) + 5] << ",";
        fs << result[s * (n + 1) + 6] << ",";
        fs << result[s * (n + 1) + 7] << ",";
        fs << result[s * (n + 1) + 8] << ",";
        fs << result[s * (n + 1) + 9] << ",";
        fs << result[s * (n + 1) + 10] << ",";
        fs << result[s * (n + 1) + 11] << ",";
        fs << result[s * (n + 1) + 12] << ",";
        fs << result[s * (n + 1) + 13] << ",";
        fs << result[s * (n + 1) + 14] << ",";
        fs << result[s * (n + 1) + 15] << ",";
        fs << result[s * (n + 1) + 16] << ",";
        fs << result[s * (n + 1) + 17] << ",";
        fs << result[s * (n + 1) + 18] << ",";
        fs << result[s * (n + 1) + 19] << ",";
        fs << result[s * (n + 1) + 20] << ",";
        fs << result[s * (n + 1) + 21] << ",";
        fs << result[s * (n + 1) + 22] << endl;
    }
    fs.close();
}



// Проверяет, все ли элементы вектора меньше заданной эпсилон.
int check_if_vector_less_than_eps(double * x, int n) {
    int counter = 0;
    for (int i = 0; i < n; i++) {
        if (x[i] < epsilon)
            counter++;
    }
    return counter;
}

// Вычисляет норму для оценки ошибки.
double norm(double * result, int step, double delta_t, double delta_t_pred,
    int n) {
    double eps = 0.0;
    double dxdt = 0.0;
    for (int i = 1; i < n + 1; i++) {
        dxdt = fabs(result[step * (n + 1) + i] -
            result[(step - 1) * (n + 1) + i] * delta_t / delta_t_pred);
        if (eps < dxdt)
            eps = dxdt;
    }
    return eps * delta_t * delta_t / 2;
}

// Заполняет матрицу A
void fill_A(double ** a, double U_Id1, double delta_t) {
    a[0][0] = 1.0;
    a[1][1] = 1.0;
    a[2][2] = 1.0;
    a[3][3] = 1.0;
    a[4][4] = 1.0;
    a[5][5] = 1.0;
    a[6][6] = 1.0;
    a[7][7] = 1.0;
    a[8][8] = 1.0;
    a[9][9] = 1.0;
    a[10][10] = 1.0;
    a[11][11] = 1.0;
    a[12][12] = 1.0;
    a[13][13] = 1.0;
    a[14][14] = 1.0;
    a[15][15] = 1.0;
    a[16][16] = -L1 / delta_t;
    // a[16][16] = 1;
    a[17][17] = 1.0;
    a[18][18] = 1.0;
    a[19][19] = -C2 / delta_t;
    a[20][20] = -Cb1 / delta_t;
    a[21][21] = -C1 / delta_t;
    
    a[0][18] = - 1;
    a[0][19] = 1;
    a[1][21] = -1;
    a[2][17] = -1;
    a[2][18] = -1;
    a[2][19] = 1;
    a[2][20] = 1;
    a[2][21] = 1;
    a[3][20] = -1;
    a[4][20] = -1;
    a[5][21] = -1;

    a[6][13] = 1;
    a[7][11] = 1;
    a[7][13] = 1;
    a[8][11] = -1;
    a[8][13] = -1;
    a[9][13] = -1;
    a[9][14] = 1;
    a[9][15] = 1;
    a[10][12] = 1;
    a[10][13] = -1;
    a[10][16] = 1;
    
    a[11][0] = -1.0 / R1;
    a[12][1] = -1.0 / R2;
    a[13][2] = -1.0 / Rb1;
    a[14][3] = -1.0 / Rb2;
    a[15][4] = -Id1 * exp(U_Id1 / MFt) / MFt;
    a[16][5] = 1;
    // a[16][5] = -L1 / delta_t;
    
    a[19][8] = 1;
    a[20][9] = 1;
    a[21][10] = 1;
}

// Заполняет вектор правой части системы уравнений.
void fill_y(double * y, 
            double U_R1, double U_R2, double U_Rb1, double U_Rb2, double U_Id1, double U_L1, 
            double I_E1, double I_E2, double I_C2, double I_Cb1, double I_C1, 
            double I_R1, double I_R2, double I_Rb1, double I_Rb2, double I_Id1, double I_L1, 
            double U_E1, double U_E2, double U_C2, double U_Cb1, double U_C1, 
            double U_Cb1_pred, double U_C1_pred, double I_C2_pred, double U_C2_pred, double U_L1_pred,
            double delta_t, int n, double time
) {
    y[0] = U_R1 - U_E2 + U_C2;
    y[1] = U_R2 - U_C1;
    y[2] = U_Rb1 - U_E1 - U_E2 + U_C2 + U_Cb1 + U_C1;
    y[3] = U_Rb2 - U_Cb1;
    y[4] = U_Id1 - U_Cb1;
    y[5] = U_L1 - U_C1;
    
    y[6] = I_E1 + I_Rb1;
    y[7] = I_E2 + I_R1 + I_Rb1;
    y[8] = I_C2 - I_R1 - I_Rb1;
    y[9] = I_Cb1 - I_Rb1 + I_Rb2 + I_Id1;
    y[10] = I_C1 + I_R2 - I_Rb1 + I_L1;
    
    y[11] = I_R1 - U_R1 / R1;
    y[12] = I_R2 - U_R2 / R2;
    y[13] = I_Rb1 - U_Rb1 / Rb1;
    y[14] = I_Rb2 - U_Rb2 / Rb2;
    y[15] = I_Id1 - Id1 * (exp(U_Id1 / MFt) - 1);
    y[16] = I_L1 - L1 / delta_t * (U_L1 - U_L1_pred);
    
    y[17] = U_E1 - 1;
    y[18] = U_E2 - 10 * sin(2 * Pi * time / 0.0001);
    y[18] = I_C2 - C2 * (U_C2 - U_C2_pred) / delta_t;
    y[19] = I_Cb1 - Cb1 * (U_Cb1 - U_Cb1_pred) / delta_t;
    y[21] = I_C1 - C1 * (U_C1 - U_C1_pred) / delta_t;
    
    for (int i = 0; i < n; i++)
        y[i] *= -1.;
}

// Инициализирует начальные условия аппроксимации.
void make_appr(double * result, int n) {
    U_R1 = appr;
    U_R2 = appr;
    U_Rb1 = appr;
    U_Rb2 = appr;
    U_Id1 = appr;
    U_L1 = appr;
    
    I_E1 = appr;
    I_E2 = appr;
    I_C2 = appr;
    I_Cb1 = appr;
    I_C1 = appr;
    
    I_R1 = appr;
    I_R2 = appr;
    I_Rb1 = appr;
    I_Rb2 = appr;
    I_Id1 = appr;
    I_L1 = appr;
    
    U_E1 = appr;
    U_E2 = appr;
    U_C2 = appr;
    U_Cb1 = appr;
    U_C1 = appr;
    
    I_L1_pred = appr;
    U_Cb1_pred = appr;
    U_C1_pred = appr;
    I_C2_pred = appr;
    U_C2_pred = appr;
    U_L1_pred = appr;
    
    for (int i = 1; i < n + 1; i++)
        result[i] = appr;
}

// Корректирует значения переменных согласно решению метода Гаусса.
void add_phi(double * x) {
    U_R1 += x[0];
    U_R2 += x[1];
    U_Rb1 += x[2];
    U_Rb2 += x[3];
    U_Id1 += x[4];
    U_L1 += x[5];
    
    I_E1 += x[6];
    I_E2 += x[7];
    I_C2 += x[8];
    I_Cb1 += x[9];
    I_C1 += x[10];
    
    I_R1 += x[11];
    I_R2 += x[12];
    I_Rb1 += x[13];
    I_Rb2 += x[14];
    I_Id1 += x[15];
    I_L1 += x[16];
    
    U_E1 += x[17];
    U_E2 += x[18];
    U_C2 += x[19];
    U_Cb1 += x[20];
    U_C1 += x[21];
}

// Сохраняет результаты текущего временного шага.
void save_step_results(double * result, int n, int step, double time) {
    result[step * (n + 1)] = time;
    
    result[step * (n + 1) + 1] = U_R1;
    result[step * (n + 1) + 2] = U_R2;
    result[step * (n + 1) + 3] = U_Rb1;
    result[step * (n + 1) + 4] = U_Rb2;
    result[step * (n + 1) + 5] = U_Id1;
    result[step * (n + 1) + 6] = U_L1;
    
    result[step * (n + 1) + 7] = I_E1;
    result[step * (n + 1) + 8] = I_E2;
    result[step * (n + 1) + 9] = I_C2;
    result[step * (n + 1) + 10] = I_Cb1;
    result[step * (n + 1) + 11] = I_C1;
    
    result[step * (n + 1) + 12] = I_R1;
    result[step * (n + 1) + 13] = I_R2;
    result[step * (n + 1) + 14] = I_Rb1;
    result[step * (n + 1) + 15] = I_Rb2;
    result[step * (n + 1) + 16] = I_Id1;
    result[step * (n + 1) + 16] = I_L1;
    
    result[step * (n + 1) + 16] = U_E1;
    result[step * (n + 1) + 17] = U_E2;
    result[step * (n + 1) + 18] = U_C2;
    result[step * (n + 1) + 19] = U_Cb1;
    result[step * (n + 1) + 20] = U_C1;
}

// Возвращает предыдущие результаты для коррекции при неудачной итерации.
void return_past_results(double * result, int n, int step) {
    U_R1 = result[(step - 1) * (n + 1) + 1];
    U_R2 = result[(step - 1) * (n + 1) + 2];
    U_Rb1 = result[(step - 1) * (n + 1) + 3];
    U_Rb2 = result[(step - 1) * (n + 1) + 4];
    U_Id1 = result[(step - 1) * (n + 1) + 5];
    U_L1 = result[(step - 1) * (n + 1) + 6];
    
    I_E1 = result[(step - 1) * (n + 1) + 7];
    I_E2 = result[(step - 1) * (n + 1) + 8];
    I_C2 = result[(step - 1) * (n + 1) + 9];
    I_Cb1 = result[(step - 1) * (n + 1) + 10];
    I_C1 = result[(step - 1) * (n + 1) + 11];
    
    I_R1 = result[(step - 1) * (n + 1) + 12];
    I_R2 = result[(step - 1) * (n + 1) + 13];
    I_Rb1 = result[(step - 1) * (n + 1) + 14];
    I_Rb2 = result[(step - 1) * (n + 1) + 15];
    I_Id1 = result[(step - 1) * (n + 1) + 16];
    I_L1 = result[(step - 1) * (n + 1) + 17];
    
    U_E1 = result[(step - 1) * (n + 1) + 18];
    U_E2 = result[(step - 1) * (n + 1) + 19];
    U_C2 = result[(step - 1) * (n + 1) + 20];
    U_Cb1 = result[(step - 1) * (n + 1) + 21];
    U_C1 = result[(step - 1) * (n + 1) + 22];
}

// Получает предсказанные значения переменных для коррекции.
void get_pred_variables(double * result, int n, int step) {
    I_L1_pred = result[(step - 1) * (n + 1) + 17];
    U_Cb1_pred = result[(step - 1) * (n + 1) + 21];
    U_C1_pred = result[(step - 1) * (n + 1) + 22];
    I_C2_pred = result[(step - 1) * (n + 1) + 9];
    U_C2_pred = result[(step - 1) * (n + 1) + 20];
    U_L1_pred = result[(step - 1) * (n + 1) + 6];
}