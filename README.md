# RST_Complexity
Implementation of the Rubio-Sweedler-Taylor Algorithm for computing the Groebner basis of the ideal of linear recurring polynomials on a multidimensional periodic array with entries in a field.

The code has two dependencies: 
<ul>
  <li> 
    <code>Blitz++</code>. A C++ template class library that provides high-performance multidimensional array containers, found at <a href="https://github.com/blitzpp/blitz">https://github.com/blitzpp/blitz</a>
  </li>
  <li>
    <code>NTL</code> A high-performance, portable C++ library providing data structures and algorithms for manipulating signed, arbitrary length integers, and for vectors, matrices, and polynomials over the integers and over finite fields.
  </li>
</ul>


Most of the action occurs in the file <code>MultiDimArray.cpp</code>. This a class consisting of a multidimensional packed with additional information such as its size, periodicity and complexity.
The <code>RST()</code> method is the member function implementing the Rubio-Sweedler-Taylor algorithms, and in the pocess of obtaining the Groebner basis, it computes the linear complexity.
