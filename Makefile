CXX = g++
CXXFLAGS = -g -std=c++11 -pthread -march=native
INC=-I /usr/include/boost/dynamic_bitset

Test_1: Test_1.o MultiDimArray_GF2.o Sequences.o
	$(CXX) $(CXXFLAGS) -o Test_1 Test_1.o MultiDimArray_GF2.o Sequences.o -lntl -lblitz -lgmp

Test_1.o: Test_1.cpp MultiDimArray_GF2.cpp Sequences.h
	$(CXX) $(CXXFLAGS) -c Test_1.cpp

MultiDimArray_GF2.o: MultiDimArray_GF2.cpp MultivarPolynomial.o
	$(CXX) $(CXXFLAGS) -c MultiDimArray_GF2.cpp

MultivarPolynomial.o: MultivarPolynomial.cpp Monomial.o
	$(CXX) $(CXXFLAGS) -c MultivarPolynomial.cpp

Monomial.o: Monomial.cpp
	$(CXX) $(CXXFLAGS) -c Monomial.cpp

Sequences.o: Sequences.h

clean :
	-rm *.o $(objects)