# If you copy and paste this Makefile, change the eight space indent to a tab.
smallworld: smallworld.cpp Makefile
	$(CXX) -fopenmp -O3 -o $@ $<
