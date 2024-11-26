plot: run
	gnuplot graph.gnu

run: compile
	./a.out

compile:
	c++ main.cpp