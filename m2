plot: run
	gnuplot m2.gnu

run: compile
	./a.out

compile:
	c++ main.cpp