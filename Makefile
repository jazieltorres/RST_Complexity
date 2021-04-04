CXX = g++
CXXFLAGS = -g -std=c++11 -pthread -march=native

ExpQuadratic_simple: main_ExpQuadratic.o MultiDimArray.o Sequences.o
	$(CXX) $(CXXFLAGS) -o ExpQuadratic_simple main_ExpQuadratic.o MultiDimArray.o Sequences.o -lntl -lblitz -lgmp

main_ExpQuadratic.o: main_ExpQuadratic.cpp MultiDimArray.cpp Sequences.h
	$(CXX) $(CXXFLAGS) -c main_ExpQuadratic.cpp

MultiDimArray.o: MultiDimArray.cpp MultivarPolynomial.o
	$(CXX) $(CXXFLAGS) -c MultiDimArray.cpp

MultivarPolynomial.o: MultivarPolynomial.cpp Monomial.o
	$(CXX) $(CXXFLAGS) -c MultivarPolynomial.cpp

Monomial.o: Monomial.cpp
	$(CXX) $(CXXFLAGS) -c Monomial.cpp

Sequences.o: Sequences.h
