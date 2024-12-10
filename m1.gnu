set encoding utf8

#Определение стилей линий и маркеров
set style line 1 lt 1 lw 0.5 pt 7 ps 0.5 lc rgb "red"

# Установка заголовка графика
set title "Value of time dependance"

# Установка подписей к осям
set xlabel "Time"
set ylabel "U_C1"

# Определение, что данные в файле разделены запятой
set datafile separator ','

# Отображение легенды
set key top left

# Использование окна x11
set term wxt enhanced

# Построение графика из файла result.csv, начиная со второй строки
plot "result.csv" skip 1 using 1:21 with linespoints linestyle 1 title ""


# Пауза, чтобы удержать окно графика открытым
pause -1