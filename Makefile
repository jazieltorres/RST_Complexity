CXX = g++
CXXFLAGS = -g -std=c++11 -pthread -march=native

main_RationalRemoveSimple: main_RationalShift.o MultiDimArray.o Sequences.o
	$(CXX) $(CXXFLAGS) -o main_RationalRemoveSimple main_RationalShift.o MultiDimArray.o Sequences.o -lntl -lblitz -lgmp

main_RationalShift.o: main_RationalShift.cpp MultiDimArray.cpp Sequences.h
	$(CXX) $(CXXFLAGS) -c main_RationalShift.cpp

MultiDimArray.o: MultiDimArray.cpp MultivarPolynomial.o
	$(CXX) $(CXXFLAGS) -c MultiDimArray.cpp

MultivarPolynomial.o: MultivarPolynomial.cpp Monomial.o
	$(CXX) $(CXXFLAGS) -c MultivarPolynomial.cpp

Monomial.o: Monomial.cpp
	$(CXX) $(CXXFLAGS) -c Monomial.cpp

Sequences.o: Sequences.h
