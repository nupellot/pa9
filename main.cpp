#include <iostream>
#include <cmath>
#include <fstream>

using namespace std;

// Определение параметров симуляции и констант
#define delta_t_appr 1e-12      // Начальное приближение шага времени
#define delta_t_min 1e-12       // Минимально допустимый шаг времени
#define end_time 1e-3           // Время завершения симуляции
#define epsilon 1e-5            // Порог сходимости для корректировок
#define E1 1                     // Источник напряжения E1
#define R1 1e3                   // Сопротивление R1 в Омах
#define R2 1e3                   // Сопротивление R2 в Омах
#define Rb1 20.0                 // Сопротивление Rb1 в Омах
#define Rb2 1e6                  // Сопротивление Rb2 в Омах
#define L1 1e-3                  // Индуктивность L1 в Генри
#define Cb1 2e-12                // Емкость Cb1 в Фарадах
#define C1 1e-6                  // Емкость C1 в Фарадах
#define Id1 1e-12                // Ток насыщения диода Id1 в Амперах
#define MFt 0.026                // Тепловое напряжение (приблизительно)
#define vector_size 70000        // Размер вектора результатов
#define eps1 1e-15              // Нижний порог нормы для адаптивного шага времени
#define eps2 1e-5               // Верхний порог нормы для адаптивного шага времени
#define appr 0.0                 // Начальное значение приближения

// Объявление глобальных переменных, представляющих напряжения узлов и токи ветвей
double U_R1, U_R2, U_Rb1, U_Rb2, U_Id1, I_E1, I_E2, I_L1, I_Cb1, I_C1,
    I_R1, I_R2, I_Rb1, I_Rb2, I_Id1, U_E1, U_E2, U_L1, U_Cb1, U_C1,
    I_L1_pred, U_Cb1_pred, U_C1_pred;

/**
 * Выполняет метод Гаусса для решения системы линейных уравнений Ax = y.
 * @param a Коэффициентная матрица (изменяется на месте).
 * @param y Вектор правых частей (изменяется на месте).
 * @param n Размер системы (количество уравнений).
 * @return Указатель на вектор решения x, или 0, если систему решить невозможно.
 */
double* gauss(double **a, double *y, int n)
{
    const double eps = 1e-12; // Малая величина для проверки сингулярности
    double *x, max;
    int k, index;
    x = new double[n]; // Выделение памяти под вектор решения
    k = 0;

    // Процесс прямого хода
    while (k < n)
    {
        // Поиск ведущего элемента
        max = abs(a[k][k]);
        index = k;
        for (int i = k + 1; i < n; i++)
        {
            if (abs(a[i][k]) > max)
            {
                max = abs(a[i][k]);
                index = i;
            }
        }

        // Проверка на сингулярную матрицу
        if (max < eps)
        {
            return 0; // Уникальное решение не существует
        }

        // Перестановка текущей строки с ведущей строкой
        for (int j = 0; j < n; j++)
        {
            double temp = a[k][j];
            a[k][j] = a[index][j];
            a[index][j] = temp;
        }
        double temp = y[k];
        y[k] = y[index];
        y[index] = temp;

        // Нормализация ведущей строки и устранение ниже
        for (int i = k; i < n; i++)
        {
            double temp = a[i][k];
            if (abs(temp) < eps)
                continue; // Пропуск, если элемент эффективно равен нулю
            for (int j = 0; j < n; j++)
                a[i][j] = a[i][j] / temp;
            y[i] = y[i] / temp;
            if (i == k)
                continue; // Пропуск ведущей строки
            for (int j = 0; j < n; j++)
                a[i][j] = a[i][j] - a[k][j]; // Устранение текущего столбца
            y[i] = y[i] - y[k];
        }
        k++;
    }

    // Обратный ход для решения системы
    for (k = n - 1; k >= 0; k--)
    {
        x[k] = y[k];
        for (int i = 0; i < k; i++)
            y[i] = y[i] - a[i][k] * x[k];
    }

    return x; // Возврат вектора решения
}

/**
 * Сохраняет результаты симуляции в CSV файл.
 * @param result Указатель на массив результатов.
 * @param step Количество записанных временных шагов.
 * @param n Количество переменных на шаг.
 */
void result_to_csv(double *result, int step, int n)
{
    ofstream fs;
    fs.open("result.csv"); // Открытие (или создание) CSV файла
    // Запись заголовков столбцов
    fs << "time,U_R1,U_R2,U_Rb1,U_Rb2,U_Id1,I_E1,I_E2,I_L1,I_Cb1,I_C1,I_R1,I_R2,I_Rb1,IRb2,I_Id1,U_E1,U_E2,U_L1,U_Cb1,U_C1" << endl;
    // Итерация по каждому временному шагу и запись соответствующих значений
    for (int s = 0; s < step; s++)
    {
        fs << result[s * (n + 1)] << ","; // Время
        fs << result[s * (n + 1) + 1] << ","; // U_R1
        fs << result[s * (n + 1) + 2] << ","; // U_R2
        fs << result[s * (n + 1) + 3] << ","; // U_Rb1
        fs << result[s * (n + 1) + 4] << ","; // U_Rb2
        fs << result[s * (n + 1) + 5] << ","; // U_Id1
        fs << result[s * (n + 1) + 6] << ","; // I_E1
        fs << result[s * (n + 1) + 7] << ","; // I_E2
        fs << result[s * (n + 1) + 8] << ","; // I_L1
        fs << result[s * (n + 1) + 9] << ","; // I_Cb1
        fs << result[s * (n + 1) + 10] << ","; // I_C1
        fs << result[s * (n + 1) + 11] << ","; // I_R1
        fs << result[s * (n + 1) + 12] << ","; // I_R2
        fs << result[s * (n + 1) + 13] << ","; // I_Rb1
        fs << result[s * (n + 1) + 14] << ","; // I_Rb2
        fs << result[s * (n + 1) + 15] << ","; // I_Id1
        fs << result[s * (n + 1) + 16] << ","; // U_E1
        fs << result[s * (n + 1) + 17] << ","; // U_E2
        fs << result[s * (n + 1) + 18] << ","; // U_L1
        fs << result[s * (n + 1) + 19] << ","; // U_Cb1
        fs << result[s * (n + 1) + 20] << endl; // U_C1
    }
    fs.close(); // Закрытие файла
}

/**
 * Проверяет, все ли элементы вектора коррекции меньше порогового значения epsilon.
 * @param x Указатель на вектор коррекции.
 * @param n Количество элементов для проверки.
 * @return Количество элементов, меньших epsilon.
 */
int check_if_vector_less_than_eps(double *x, int n)
{
    int counter = 0;
    for (int i = 0; i < n; i++)
    {
        if (x[i] < epsilon)
            counter++; // Увеличение счетчика, если элемент меньше порога
    }
    return counter;
}

/**
 * Вычисляет норму разности между текущими и предсказанными результатами.
 * @param result Указатель на массив результатов.
 * @param step Текущий временной шаг.
 * @param delta_t Текущий размер шага времени.
 * @param delta_t_pred Предыдущий размер шага времени.
 * @param n Количество переменных.
 * @return Вычисленная норма.
 */
double norm(double *result, int step, double delta_t, double delta_t_pred,
            int n)
{
    double eps = 0.0;
    double dxdt = 0.0;
    // Итерация по каждой переменной для вычисления максимальной разности
    for (int i = 1; i < n + 1; i++)
    {
        dxdt = fabs(result[step * (n + 1) + i] - result[(step - 1) * (n + 1) + i] * delta_t / delta_t_pred);
        if (eps < dxdt)
            eps = dxdt; // Обновление максимальной разности
    }
    return eps * delta_t * delta_t / 2; // Возврат масштабированной нормы
}

/**
 * Заполняет матрицу коэффициентов A на основе текущих переменных состояния и шага времени.
 * @param a Матрица коэффициентов для заполнения.
 * @param U_Id1 Текущее напряжение на диоде Id1.
 * @param delta_t Текущий размер шага времени.
 */
void fill_A(double **a, double U_Id1, double delta_t)
{
    // Инициализация диагональных элементов
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
    a[16][16] = 1.0;
    a[17][17] = 1.0;
    a[18][18] = -Cb1 / delta_t; // Связь с конденсатором Cb1
    a[19][19] = -C1 / delta_t;  // Связь с конденсатором C1

    // Заполнение недиагональных элементов на основе уравнений цепи
    a[0][16] = -1.0; // Связь между U_R1 и U_E2
    a[0][17] = +1.0; // Связь с U_L1

    a[1][19] = -1.0; // Связь с U_C1

    a[2][15] = -1.0; // Связь с U_E1
    a[2][16] = -1.0; // Связь с U_E2
    a[2][17] = 1.0;  // Связь с U_L1
    a[2][18] = 1.0;  // Связь с U_Cb1
    a[2][19] = 1.0;  // Связь с U_C1

    a[3][18] = -1.0; // Связь с U_Cb1

    a[4][18] = -1.0; // Связь с U_Cb1

    a[5][12] = 1.0; // Связь с I_Rb1

    a[6][10] = +1.0; // Связь с I_R1
    a[6][12] = 1.0;  // Связь с I_Rb1

    a[7][10] = -1.0; // Связь с I_R1
    a[7][12] = -1.0; // Связь с I_Rb1

    a[8][12] = -1.0; // Связь с I_Rb1
    a[8][13] = 1.0;  // Связь с I_Rb2
    a[8][14] = 1.0;  // Связь с I_Id1

    a[9][11] = 1.0;  // Связь с I_R2
    a[9][12] = -1.0; // Связь с I_Rb1

    a[10][0] = -1.0 / R1; // Закон Ома для R1

    a[11][1] = -1.0 / R2; // Закон Ома для R2

    a[12][2] = -1.0 / Rb1; // Закон Ома для Rb1

    a[13][3] = -1.0 / Rb2; // Закон Ома для Rb2

    // Линеаризация уравнения диода
    a[14][4] = -Id1 * exp(U_Id1 / MFt) / MFt;

    a[17][7] = -L1 / delta_t; // Связь с индуктивностью

    a[18][8] = 1.0; // Связь с конденсатором Cb1
    a[19][9] = 1.0; // Связь с конденсатором C1
}

/**
 * Заполняет вектор правых частей y на основе текущих переменных состояния и предсказаний.
 * @param y Вектор для заполнения.
 * @param U_R1 Напряжение на резисторе R1.
 * @param U_R2 Напряжение на резисторе R2.
 * @param U_Rb1 Напряжение на резисторе Rb1.
 * @param U_Rb2 Напряжение на резисторе Rb2.
 * @param U_Id1 Напряжение на диоде Id1.
 * @param I_E1 Ток через источник E1.
 * @param I_E2 Ток через источник E2.
 * @param I_L1 Ток через индуктивность L1.
 * @param I_Cb1 Ток через конденсатор Cb1.
 * @param I_C1 Ток через конденсатор C1.
 * @param I_R1 Ток через резистор R1.
 * @param I_R2 Ток через резистор R2.
 * @param I_Rb1 Ток через резистор Rb1.
 * @param I_Rb2 Ток через резистор Rb2.
 * @param I_Id1 Ток через диод Id1.
 * @param U_E1 Напряжение источника E1.
 * @param U_E2 Напряжение источника E2.
 * @param U_L1 Напряжение на индуктивности L1.
 * @param U_Cb1 Напряжение на конденсаторе Cb1.
 * @param U_C1 Напряжение на конденсаторе C1.
 * @param I_L1_pred Предсказанный ток через индуктивность L1 из предыдущего шага.
 * @param U_Cb1_pred Предсказанное напряжение на конденсаторе Cb1 из предыдущего шага.
 * @param U_C1_pred Предсказанное напряжение на конденсаторе C1 из предыдущего шага.
 * @param delta_t Текущий размер шага времени.
 * @param n Количество переменных.
 * @param time Текущее время симуляции.
 */
void fill_y(double *y, double U_R1, double U_R2, double U_Rb1, double U_Rb2,
            double U_Id1, double I_E1, double I_E2, double I_L1, double I_Cb1, double I_C1,
            double I_R1, double I_R2, double I_Rb1, double I_Rb2, double I_Id1, double U_E1, double U_E2, double U_L1, double U_Cb1, double U_C1,
            double I_L1_pred, double U_Cb1_pred, double U_C1_pred, double delta_t, int n, double time)
{
    // Определение системы уравнений на основе законов Кирхгофа и связей элементов цепи
    y[0] = U_R1 - U_E2 + U_L1;
    y[1] = U_R2 - U_C1;
    y[2] = U_Rb1 - U_E1 - U_E2 + U_L1 + U_Cb1 + U_C1;
    y[3] = U_Rb2 - U_Cb1;
    y[4] = U_Id1 - U_Cb1;
    y[5] = I_E1 + I_Rb1;
    y[6] = I_E2 + I_R1 + I_Rb1;
    y[7] = I_L1 - I_R1 - I_Rb1;
    y[8] = I_Cb1 - I_Rb1 + I_Rb2 + I_Id1;
    y[9] = I_C1 + I_R2 - I_Rb1;
    y[10] = I_R1 - U_R1 / R1;
    y[11] = I_R2 - U_R2 / R2;
    y[12] = I_Rb1 - U_Rb1 / Rb1;
    y[13] = I_Rb2 - U_Rb2 / Rb2;
    y[14] = I_Id1 - Id1 * (exp(U_Id1 / MFt) - 1);
    y[15] = U_E1 - 1;
    y[16] = U_E2 - 10 * sin(2 * M_PI * time / 0.0001); // Временно-зависимый источник напряжения E2
    y[17] = U_L1 - L1 * (I_L1 - I_L1_pred) / delta_t; // Связь напряжения на индуктивности
    y[18] = I_Cb1 - Cb1 * (U_Cb1 - U_Cb1_pred) / delta_t; // Связь тока через конденсатор Cb1
    y[19] = I_C1 - C1 * (U_C1 - U_C1_pred) / delta_t; // Связь тока через конденсатор C1

    // Умножение каждого уравнения на -1 для формирования системы Ax = y
    for (int i = 0; i < n; i++)
        y[i] *= -1.;
}

/**
 * Инициализирует все переменные состояния и предсказания значением приближения.
 * @param result Указатель на массив результатов.
 * @param n Количество переменных.
 */
void make_appr(double *result, int n)
{
    // Инициализация всех напряжений и токов значением приближения
    U_R1 = appr;
    U_R2 = appr;
    U_Rb1 = appr;
    U_Rb2 = appr;
    U_Id1 = appr;
    I_E1 = appr;
    I_E2 = appr;
    I_L1 = appr;
    I_Cb1 = appr;
    I_C1 = appr;
    I_R1 = appr;
    I_R2 = appr;
    I_Rb1 = appr;
    I_Rb2 = appr;
    I_Id1 = appr;
    U_E1 = appr;
    U_E2 = appr;
    U_L1 = appr;
    U_Cb1 = appr;
    U_C1 = appr;
    I_L1_pred = appr;
    U_Cb1_pred = appr;
    U_C1_pred = appr;

    // Инициализация массива результатов значениями приближения
    for (int i = 1; i < n + 1; i++)
        result[i] = appr;
}

/**
 * Обновляет переменные состояния, добавляя коррекции из вектора решения.
 * @param x Указатель на вектор коррекции.
 */
void add_phi(double *x)
{
    U_R1 += x[0];
    U_R2 += x[1];
    U_Rb1 += x[2];
    U_Rb2 += x[3];
    U_Id1 += x[4];
    I_E1 += x[5];
    I_E2 += x[6];
    I_L1 += x[7];
    I_Cb1 += x[8];
    I_C1 += x[9];
    I_R1 += x[10];
    I_R2 += x[11];
    I_Rb1 += x[12];
    I_Rb2 += x[13];
    I_Id1 += x[14];
    U_E1 += x[15];
    U_E2 += x[16];
    U_L1 += x[17];
    U_Cb1 += x[18];
    U_C1 += x[19];
}

/**
 * Сохраняет результаты текущего шага в массив результатов.
 * @param result Указатель на массив результатов.
 * @param n Количество переменных.
 * @param step Индекс текущего временного шага.
 * @param time Текущее время симуляции.
 */
void save_step_results(double *result, int n, int step, double time)
{
    result[step * (n + 1)] = time; // Сохранение текущего времени
    result[step * (n + 1) + 1] = U_R1;
    result[step * (n + 1) + 2] = U_R2;
    result[step * (n + 1) + 3] = U_Rb1;
    result[step * (n + 1) + 4] = U_Rb2;
    result[step * (n + 1) + 5] = U_Id1;
    result[step * (n + 1) + 6] = I_E1;
    result[step * (n + 1) + 7] = I_E2;
    result[step * (n + 1) + 8] = I_L1;
    result[step * (n + 1) + 9] = I_Cb1;
    result[step * (n + 1) + 10] = I_C1;
    result[step * (n + 1) + 11] = I_R1;
    result[step * (n + 1) + 12] = I_R2;
    result[step * (n + 1) + 13] = I_Rb1;
    result[step * (n + 1) + 14] = I_Rb2;
    result[step * (n + 1) + 15] = I_Id1;
    result[step * (n + 1) + 16] = U_E1;
    result[step * (n + 1) + 17] = U_E2;
    result[step * (n + 1) + 18] = U_L1;
    result[step * (n + 1) + 19] = U_Cb1;
    result[step * (n + 1) + 20] = U_C1;
}

/**
 * Возвращает переменные состояния к значениям из предыдущего временного шага.
 * Это используется, когда текущая итерация не сходится.
 * @param result Указатель на массив результатов.
 * @param n Количество переменных.
 * @param step Индекс текущего временного шага.
 */
void return_past_results(double *result, int n, int step)
{
    // Восстановление каждой переменной из предыдущего временного шага
    U_R1 = result[(step - 1) * (n + 1) + 1];
    U_R2 = result[(step - 1) * (n + 1) + 2];
    U_Rb1 = result[(step - 1) * (n + 1) + 3];
    U_Rb2 = result[(step - 1) * (n + 1) + 4];
    U_Id1 = result[(step - 1) * (n + 1) + 5];
    I_E1 = result[(step - 1) * (n + 1) + 6];
    I_E2 = result[(step - 1) * (n + 1) + 7];
    I_L1 = result[(step - 1) * (n + 1) + 8];
    I_Cb1 = result[(step - 1) * (n + 1) + 9];
    I_C1 = result[(step - 1) * (n + 1) + 10];
    I_R1 = result[(step - 1) * (n + 1) + 11];
    I_R2 = result[(step - 1) * (n + 1) + 12];
    I_Rb1 = result[(step - 1) * (n + 1) + 13];
    I_Rb2 = result[(step - 1) * (n + 1) + 14];
    I_Id1 = result[(step - 1) * (n + 1) + 15];
    U_E1 = result[(step - 1) * (n + 1) + 16];
    U_E2 = result[(step - 1) * (n + 1) + 17];
    U_L1 = result[(step - 1) * (n + 1) + 18];
    U_Cb1 = result[(step - 1) * (n + 1) + 19];
    U_C1 = result[(step - 1) * (n + 1) + 20];
}

/**
 * Извлекает предсказанные значения некоторых переменных из предыдущего временного шага.
 * @param result Указатель на массив результатов.
 * @param n Количество переменных.
 * @param step Индекс текущего временного шага.
 */
void get_pred_variables(double *result, int n, int step)
{
    I_L1_pred = result[(step - 1) * (n + 1) + 8];
    U_Cb1_pred = result[(step - 1) * (n + 1) + 19];
    U_C1_pred = result[(step - 1) * (n + 1) + 20];
}

/**
 * Главная функция, управляющая симуляцией.
 */
int main()
{
    int n = 20;               // Количество переменных в системе
    int n_iter = 0;           // Счетчик итераций для метода Ньютона-Рафсона
    double time = 0.0;        // Инициализация времени симуляции

    // Выделение памяти под матрицу коэффициентов A и векторы y и result
    double **a = new double *[n];
    double *y = new double[n], *x;
    double *result = new double[vector_size];

    result[0] = time;         // Установка начального времени в массиве результатов
    make_appr(result, n);     // Инициализация всех переменных значениями приближения

    double delta_t = delta_t_appr; // Установка начального размера шага времени
    double delta_t_pred = 0.0;      // Предыдущий размер шага времени
    double eps = 0.0;                // Норма для сходимости

    int flag = 1;        // Флаг для контроля итераций
    int step = 1;        // Индекс текущего временного шага
    int time_step = 0;   // Общее количество выполненных временных шагов
    int iter_step = 0;   // Общее количество итераций

    // Основной цикл симуляции, выполняющийся до достижения end_time
    while (time < end_time)
    {
        get_pred_variables(result, n, step); // Получение предсказанных переменных из предыдущего шага

        time_step++; // Увеличение счетчика временных шагов
        while (flag)
        {
            // Инициализация матрицы коэффициентов A и вектора правых частей y
            for (int i = 0; i < n; i++)
            {
                a[i] = new double[n]; // Выделение памяти под каждую строку A
                y[i] = 0.0;            // Инициализация y нулем
                for (int j = 0; j < n; j++)
                    a[i][j] = 0.0;     // Инициализация A нулями
            }
            iter_step++; // Увеличение счетчика итераций
            n_iter++;    // Увеличение счетчика итераций для Ньютона-Рафсона

            fill_A(a, U_Id1, delta_t); // Заполнение матрицы коэффициентов на основе текущего состояния
            fill_y(y, U_R1, U_R2, U_Rb1, U_Rb2, U_Id1, I_E1, I_E2, I_L1, I_Cb1, I_C1,
                   I_R1, I_R2, I_Rb1, I_Rb2, I_Id1, U_E1, U_E2, U_L1, U_Cb1, U_C1,
                   I_L1_pred, U_Cb1_pred, U_C1_pred, delta_t, n, time); // Заполнение вектора правых частей

            // Раскомментируйте следующие строки для отладки
            // cout << "time: " << time << endl;
            // cout << "delta_t: " << delta_t << endl;

            x = gauss(a, y, n); // Решение линейной системы методом Гаусса
            add_phi(x);          // Обновление переменных состояния с помощью коррекций

            // Проверка, все ли корректировки меньше порога сходимости
            if (check_if_vector_less_than_eps(x, n) == n)
            {
                n_iter = 0; // Сброс счетчика итераций
                flag = 0;    // Выход из цикла итераций
            }
            else
            {
                // Если не сошлось после 7 итераций, уменьшаем шаг времени
                if (n_iter > 7)
                {
                    delta_t = delta_t / 2; // Уменьшение шага времени вдвое
                    n_iter = 0;            // Сброс счетчика итераций
                    return_past_results(result, n, step); // Возврат к предыдущему состоянию

                    // Проверка, не стал ли шаг времени меньше минимального порога
                    if (delta_t < delta_t_min)
                    {
                        cout << "Решение не сходится" << endl;
                        exit(0); // Завершение программы
                    }
                }
            }
        }

        // Специальная обработка для первого шага для корректной инициализации
        if (step == 1 && time == 0.0)
        {
            step = 0;
        }

        save_step_results(result, n, step, time); // Сохранение результатов текущего шага

        delta_t_pred = delta_t; // Обновление предыдущего шага времени
        cout << time << endl;    // Вывод текущего времени симуляции

        if (step > 0)
        {
            eps = norm(result, step, delta_t, delta_t_pred, n); // Вычисление нормы для сходимости

            // Регулировка шага времени на основе значения нормы
            if (eps < eps1)
            {
                time += delta_t;         // Продвижение времени симуляции
                delta_t_pred = delta_t;  // Обновление предыдущего шага времени
                delta_t = 2 * delta_t;   // Удвоение шага времени для повышения эффективности

                // Сдвиг результатов для следующей итерации
                for (int i = 1; i < n + 1; i++)
                {
                    if (step > 1)
                    {
                        result[(step - 2) * (n + 1) + i] = result[(step - 1) * (n + 1) + i];
                        result[(step - 1) * (n + 1) + i] = result[step * (n + 1) + i];
                    }
                }
            }
            if (eps > eps1 && eps < eps2)
            {
                time += delta_t;         // Продвижение времени симуляции
                delta_t_pred = delta_t;  // Обновление предыдущего шага времени

                // Сдвиг результатов для следующей итерации
                for (int i = 1; i < n + 1; i++)
                {
                    if (step > 2)
                    {
                        result[(step - 2) * (n + 1) + i] = result[(step - 1) * (n + 1) + i];
                        result[(step - 1) * (n + 1) + i] = result[step * (n + 1) + i];
                    }
                }
            }
            if (eps > eps2)
            {
                // Если норма слишком большая, уменьшаем шаг времени и возвращаемся к предыдущему состоянию
                time = time; // Время не изменяется
                delta_t = delta_t / 2; // Уменьшение шага времени вдвое
                for (int i = 1; i < n + 1; i++)
                    if (step >= 1)
                    {
                        result[step * (n + 1) + i] = result[(step - 1) * (n + 1) + i];
                    }
            }
        }
        else
        {
            // Если шаг равен нулю, уменьшаем шаг времени и продвигаем время
            delta_t = delta_t / 2;
            time += delta_t;
        }
        step++;      // Переход к следующему временному шагу
        flag = 1;    // Сброс флага итераций для следующего шага
    }

    // Вывод сообщений о завершении симуляции
    cout << "Успех" << endl;
    cout << "Количество временных шагов: " << time_step << endl;
    result_to_csv(result, step, n); // Сохранение всех результатов в CSV

    // Освобождение динамически выделенной памяти
    delete[] y;
    delete[] result;
    for (int i = 0; i < n; i++)
        delete[] a[i];
    delete[] a; // Не забываем удалить массив указателей

    return 0; // Завершение программы
}
