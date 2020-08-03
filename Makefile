CXX = g++
CXXFLAGS = -g -std=c++11 -pthread -march=native

main: main.o MultiDimArray.o
	$(CXX) $(CXXFLAGS) -o main main.o MultiDimArray.o -lntl -lblitz

main.o: main.cpp MultiDimArray.hpp
	$(CXX) $(CXXFLAGS) -c main.cpp

MultiDimArray.o: MultiDimArray.hpp MExponent.o
	$(CXX) $(CXXFLAGS) -c MultiDimArray.hpp

MExponent.o: MExponent.hpp
	$(CXX) $(CXXFLAGS) -c MExponent.hpp
