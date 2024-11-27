# Установим стиль графика
set terminal pngcairo size 1920,1080 enhanced font 'Arial,14'
set output 'circuit_simulation_results.png'

# Заголовок графика
set title "Симуляция электрической схемы" font ",18"

# Легенда
set key outside right top box

# Сетки для графика
set grid xtics ytics mxtics mytics linetype -1 linewidth 0.5

# Оси и их формат
set xlabel "Время (с)" font ",14"
set ylabel "Значение (В или А)" font ",14"
set format x "%.1e"
set format y "%.1e"

# Установим линии для разных параметров
set style line 1 linecolor rgb '#FF0000' linewidth 2 pointtype 7 pointsize 1.5
set style line 2 linecolor rgb '#00FF00' linewidth 2 pointtype 7 pointsize 1.5
set style line 3 linecolor rgb '#0000FF' linewidth 2 pointtype 7 pointsize 1.5
set style line 4 linecolor rgb '#FFA500' linewidth 2 pointtype 7 pointsize 1.5
set style line 5 linecolor rgb '#800080' linewidth 2 pointtype 7 pointsize 1.5
set style line 6 linecolor rgb '#00CED1' linewidth 2 pointtype 7 pointsize 1.5
set style line 7 linecolor rgb '#8B4513' linewidth 2 pointtype 7 pointsize 1.5
set style line 8 linecolor rgb '#4682B4' linewidth 2 pointtype 7 pointsize 1.5
set style line 9 linecolor rgb '#FFD700' linewidth 2 pointtype 7 pointsize 1.5

# Читаем данные из файла result.csv
set datafile separator ","

# Параметры схемы
# Колонки: time,U_R1,U_R2,U_Rb1,U_Rb2,U_Id1,UL1,I_E1,I_E2,I_C2,I_Cb1,I_C1,I_R1,I_R2,I_Rb1,I_Rb2,I_Id1,I_L1,U_E1,U_E2,U_C2,U_Cb1,U_C1

# Графики напряжений
set multiplot layout 3,1 title "Графики параметров схемы" font ",16"

set title "Напряжения в цепи"
plot "result.csv" using 1:2 with lines linestyle 1 title "U_R1 (В)", \
     "result.csv" using 1:3 with lines linestyle 2 title "U_R2 (В)", \
     "result.csv" using 1:4 with lines linestyle 3 title "U_Rb1 (В)", \
     "result.csv" using 1:5 with lines linestyle 4 title "U_Rb2 (В)", \
     "result.csv" using 1:6 with lines linestyle 5 title "U_Id1 (В)", \
     "result.csv" using 1:7 with lines linestyle 6 title "UL1 (В)", \
     "result.csv" using 1:20 with lines linestyle 7 title "U_E1 (В)", \
     "result.csv" using 1:21 with lines linestyle 8 title "U_E2 (В)"

# Графики токов
set title "Токи в цепи"
plot "result.csv" using 1:8 with lines linestyle 1 title "I_E1 (А)", \
     "result.csv" using 1:9 with lines linestyle 2 title "I_E2 (А)", \
     "result.csv" using 1:10 with lines linestyle 3 title "I_C2 (А)", \
     "result.csv" using 1:11 with lines linestyle 4 title "I_Cb1 (А)", \
     "result.csv" using 1:12 with lines linestyle 5 title "I_C1 (А)", \
     "result.csv" using 1:13 with lines linestyle 6 title "I_R1 (А)", \
     "result.csv" using 1:14 with lines linestyle 7 title "I_R2 (А)", \
     "result.csv" using 1:15 with lines linestyle 8 title "I_Rb1 (А)", \
     "result.csv" using 1:16 with lines linestyle 9 title "I_Rb2 (А)"

# Сопротивления и индуктивности
set title "Другие параметры цепи"
plot "result.csv" using 1:17 with lines linestyle 1 title "I_Id1 (А)", \
     "result.csv" using 1:18 with lines linestyle 2 title "I_L1 (А)", \
     "result.csv" using 1:19 with lines linestyle 3 title "U_Cb1 (В)"

unset multiplot
