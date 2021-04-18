CXX = g++
CXXFLAGS = -g -std=c++11 -pthread -march=native
INC=-I /usr/include/boost/dynamic_bitset

main_RationalRemoveSimple: main_RationalShift.o MultiDimArray_GF2.o Sequences.o
	$(CXX) $(CXXFLAGS) -o main_RationalRemoveSimple main_RationalShift.o MultiDimArray_GF2.o Sequences.o -lntl -lblitz -lgmp

main_RationalShift.o: main_RationalShift.cpp MultiDimArray_GF2.cpp Sequences.h
	$(CXX) $(CXXFLAGS) -c main_RationalShift.cpp

MultiDimArray_GF2.o: MultiDimArray_GF2.cpp MultivarPolynomial.o
	$(CXX) $(CXXFLAGS) -c MultiDimArray_GF2.cpp

MultivarPolynomial.o: MultivarPolynomial.cpp Monomial.o
	$(CXX) $(CXXFLAGS) -c MultivarPolynomial.cpp

Monomial.o: Monomial.cpp
	$(CXX) $(CXXFLAGS) -c Monomial.cpp

Sequences.o: Sequences.h

clean :
	-rm *.o $(objects)