CXX = g++
CXXFLAGS = -g -std=c++11 -pthread -march=native

3D_testBinary: main_MultipleTests.o MultiDimArray.o Sequences.o
	$(CXX) $(CXXFLAGS) -o 3D_testBinary main_MultipleTests.o MultiDimArray.o Sequences.o -lntl -lblitz -lgmp

main_MultipleTests.o: main_MultipleTests.cpp MultiDimArray.cpp Sequences.h
	$(CXX) $(CXXFLAGS) -c main_MultipleTests.cpp

MultiDimArray.o: MultiDimArray.cpp MultivarPolynomial.o
	$(CXX) $(CXXFLAGS) -c MultiDimArray.cpp

MultivarPolynomial.o: MultivarPolynomial.cpp Monomial.o
	$(CXX) $(CXXFLAGS) -c MultivarPolynomial.cpp

Monomial.o: Monomial.cpp
	$(CXX) $(CXXFLAGS) -c Monomial.cpp

Sequences.o: Sequences.h
