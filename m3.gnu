# Устанавливаем разделитель данных в CSV
set datafile separator ','

# Общие настройки графиков
set xlabel "Time"
set ylabel "Value"
set term pngcairo enhanced font "Arial,10" size 1280,720

# Список имён переменных для заголовков графиков
variables = "U_R1 U_R2 U_Rb1 U_Rb2 U_Id1 I_E1 I_E2 I_C2 I_Cb1 I_C1 I_R1 I_R2 I_Rb1 I_Rb2 I_Id1 U_E1 U_E2 U_C2 U_Cb1 U_C1"

# Создание графиков для каждой переменной
do for [i=2:21] {
    # Извлекаем имя переменной из списка
    variable_name = word(variables, i-1)

    # Настройка заголовка графика
    set title sprintf("Time Dependence of %s", variable_name)

    # Устанавливаем имя выходного файла
    set output sprintf("png/result_%s.png", variable_name)

    # Построение графика для текущей переменной
    plot "result.csv" skip 1 using 1:i with lines title variable_name
}

# Сброс выходного файла
set output
