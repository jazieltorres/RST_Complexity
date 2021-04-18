CXX = g++
CXXFLAGS = -g -std=c++11 -pthread -march=native
INC=-I /usr/include/boost/dynamic_bitset

multiTest_binarySimple: main_MultipleTests.o MultiDimArray_GF2.o Sequences.o
	$(CXX) $(CXXFLAGS) -o multiTest_binarySimple main_MultipleTests.o MultiDimArray_GF2.o Sequences.o -lntl -lblitz -lgmp

main_MultipleTests.o: main_MultipleTests.cpp MultiDimArray_GF2.cpp Sequences.h
	$(CXX) $(CXXFLAGS) -c main_MultipleTests.cpp

MultiDimArray_GF2.o: MultiDimArray_GF2.cpp MultivarPolynomial.o
	$(CXX) $(CXXFLAGS) -c MultiDimArray_GF2.cpp

MultivarPolynomial.o: MultivarPolynomial.cpp Monomial.o
	$(CXX) $(CXXFLAGS) -c MultivarPolynomial.cpp

Monomial.o: Monomial.cpp
	$(CXX) $(CXXFLAGS) -c Monomial.cpp

Sequences.o: Sequences.h

clean :
	-rm *.o $(objects)