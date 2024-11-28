#include <cmath>
#include <fstream>
#include <iostream>

using namespace std;

#define Pi 3.141592653589793  // Значение числа Пи
#define delta_t_appr 1e-10    // Начальный шаг по времени для аппроксимации
#define delta_t_min 1e-14     // Минимальный шаг по времени
#define end_time 1e-3         // Время моделирования (конечное)
#define epsilon 1e-5          // Погрешность для итераций
#define E1 1                   // Напряжение источника E1
#define R1 1e3                 // Сопротивление R1
#define R2 1e3                 // Сопротивление R2
#define Rb1 20.0               // Сопротивление Rb1
#define Rb2 1e6                // Сопротивление Rb2
#define L1 1e-3                // Индуктивность L1
#define Cb1 2e-12              // Ёмкость Cb1
#define C1 1e-6                // Ёмкость C1
#define Id1 1e-12              // Ток Id1
#define MFt 0.026              // Параметр MFt (параметр для экспоненты)
#define vector_size 70000      // Размер вектора для хранения результатов
#define eps1 1e-15             // Минимальная погрешность для оценки точности
#define eps2 1e-5              // Погрешность для принятия шага времени
#define appr 0.0               // Начальная аппроксимация для всех переменных
#define C2 1e-6                // Ёмкость C2

// Объявление всех переменных, которые будут использоваться для хранения значений в процессе моделирования.
double U_R1, U_R2, U_Rb1, U_Rb2, U_Id1, I_E1, I_E2, I_L1, I_Cb1, I_C1, I_R1,
I_R2, I_Rb1, I_Rb2, I_Id1, U_E1, U_E2, U_L1, U_Cb1, U_C1, I_L1_pred, U_C2, I_C2,
U_Cb1_pred, U_C1_pred, I_C2_pred, U_C2_pred, U_L1_pred;

// Прототипы всех функций, которые будут использоваться в программе.
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
            double delta_t, int n, double time);
void fill_A(double ** a, double U_Id1, double delta_t);

// Основная функция для решения системы дифференциальных уравнений с использованием итерационного метода.
int main() {
    int n = 22;               // Размер системы линейных уравнений (22 переменные)
    int n_iter = 0;           // Счётчик числа итераций
    double time = 0.0;        // Текущее время моделирования
    double ** a = new double * [n]; // Матрица коэффициентов для системы уравнений
    double * y = new double[n], * x; // Векторы правой части и решения системы
    double * result = new double[vector_size]; // Массив для хранения результатов
    result[0] = time;          // Инициализация начального времени
    make_appr(result, n);      // Инициализация значений переменных для аппроксимации
    double delta_t = delta_t_appr, delta_t_pred = 0.0, eps = 0.0; // Начальный шаг по времени
    int flag = 1;              // Флаг для проверки сходимости итераций
    int step = 1;              // Счётчик шагов моделирования
    int time_step = 0;         // Счётчик временных шагов
    int iter_step = 0;         // Счётчик числа шагов в каждой итерации
    
    // Основной цикл моделирования
    while (time < end_time) {
        get_pred_variables(result, n, step);  // Получаем предсказанные значения переменных
        time_step++;                          // Увеличиваем счётчик временных шагов
        
        while (flag) {  // Пока не будет достигнута сходимость, продолжаем итерации
            for (int i = 0; i < n; i++) {
                a[i] = new double[n];  // Выделяем память для матрицы коэффициентов
                y[i] = 0.0;            // Инициализируем вектор правой части нулями
                for (int j = 0; j < n; j++)
                    a[i][j] = 0.0;    // Инициализируем элементы матрицы значениями 0
            }
            
            iter_step++;         // Увеличиваем шаг итерации
            n_iter++;            // Увеличиваем общее число итераций
            fill_A(a, U_Id1, delta_t);  // Заполняем матрицу A для системы уравнений
            fill_y(y, 
                        U_R1, U_R2, U_Rb1, U_Rb2, U_Id1, U_L1, 
                        I_E1, I_E2, I_C2, I_Cb1, I_C1, 
                        I_R1, I_R2, I_Rb1, I_Rb2, I_Id1, I_L1, 
                        U_E1, U_E2, U_C2, U_Cb1, U_C1, 
                        U_Cb1_pred, U_C1_pred, I_C2_pred, U_C2_pred, U_L1_pred,
                        delta_t, n, time
            ); // Заполняем вектор правой части системы уравнений
            x = gauss(a, y, n);  // Решаем систему линейных уравнений методом Гаусса
            add_phi(x);          // Корректируем значения переменных на основе решения системы
            
            // Проверка сходимости (если все элементы вектора меньше заданной погрешности)
            if (check_if_vector_less_than_eps(x, n) == n) {
                n_iter = 0;
                flag = 0;  // Если сходимость достигнута, выходим из итераций
            } else {
                // Если не достигнута сходимость за 7 итераций, уменьшаем шаг по времени
                if (n_iter > 7) {
                    delta_t = delta_t / 2;
                    n_iter = 0;
                    return_past_results(result, n, step);  // Возвращаем предыдущие результаты
                    if (delta_t < delta_t_min) {
                        cout << "Solution doesn't converge" << endl;  // Если шаг слишком мал, выводим ошибку
                        exit(0);  // Выход из программы
                    }
                }
            }
        }
        
        // Обработка сохранённых данных и корректировка шага времени в зависимости от погрешности
        if (step == 1 && time == 0.0) {
            step = 0;  // Если это первый шаг моделирования, устанавливаем его как 0
        }
        save_step_results(result, n, step, time);  // Сохраняем результаты текущего шага
        delta_t_pred = delta_t;  // Предсказанный шаг времени для следующего шага
        
        // Оценка погрешности для адаптивного шага времени
        if (step > 0) {
            eps = norm(result, step, delta_t, delta_t_pred, n);  // Вычисляем норму ошибки
            if (eps < eps1) {  // Если погрешность мала, увеличиваем шаг времени
                time += delta_t;  // Переходим к следующему времени
                delta_t_pred = delta_t;
                delta_t = 2 * delta_t;  // Увеличиваем шаг времени
            }
            if (eps > eps1 && eps < eps2) {  // Если погрешность в промежутке, оставляем шаг времени
                time += delta_t;
                delta_t_pred = delta_t;
            }
            if (eps > eps2) {  // Если погрешность велика, уменьшаем шаг времени
                time += 0;  // Время не увеличиваем
                delta_t = delta_t / 2;
            }
        } else {  // Для первого шага времени уменьшаем шаг по времени
            delta_t = delta_t / 2;
            time += delta_t;
        }
        step++;  // Увеличиваем счётчик шагов
        flag = 1;  // Обнуляем флаг для следующей итерации
    }

    // Если моделирование завершено, выводим успех
    cout << "success" << endl;
    cout << "Time steps: " << time_step << endl;  // Выводим количество шагов по времени

    result_to_csv(result, step, n);  // Сохраняем результаты в CSV файл

    // Освобождаем память, выделенную для матриц и векторов
    delete[] y;
    delete[] result;
    for (int i = 0; i < n; i++) {
        delete[] a[i];
    }

    return 0;  // Завершаем выполнение программы
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