CXX = g++
CXXFLAGS = -g -std=c++11 -pthread -march=native

main: main.o MultiDimArray.o Functions.o
	$(CXX) $(CXXFLAGS) -o main main.o MultiDimArray.o Functions.o -lntl -lblitz

main.o: main.cpp MultiDimArray.cpp Functions.h
	$(CXX) $(CXXFLAGS) -c main.cpp

MultiDimArray.o: MultiDimArray.cpp MExponent.o
	$(CXX) $(CXXFLAGS) -c MultiDimArray.cpp

MExponent.o: MExponent.cpp
	$(CXX) $(CXXFLAGS) -c MExponent.cpp

Functions.o: Functions.h