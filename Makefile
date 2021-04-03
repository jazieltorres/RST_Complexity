CXX = g++
CXXFLAGS = -g -std=c++11 -pthread -march=native

Test_HardCodeSeq: Test_HardCodeSeq.o MultiDimArray.o Sequences.o
	$(CXX) $(CXXFLAGS) -o Test_HardCodeSeq Test_HardCodeSeq.o MultiDimArray.o Sequences.o -lntl -lblitz -lgmp

Test_HardCodeSeq.o: Test_HardCodeSeq.cpp MultiDimArray.cpp Sequences.h
	$(CXX) $(CXXFLAGS) -c Test_HardCodeSeq.cpp

MultiDimArray.o: MultiDimArray.cpp MultivarPolynomial.o
	$(CXX) $(CXXFLAGS) -c MultiDimArray.cpp

MultivarPolynomial.o: MultivarPolynomial.cpp Monomial.o
	$(CXX) $(CXXFLAGS) -c MultivarPolynomial.cpp

Monomial.o: Monomial.cpp
	$(CXX) $(CXXFLAGS) -c Monomial.cpp

Sequences.o: Sequences.h
