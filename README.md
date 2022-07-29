# The Rubio-Sweedler-Taylor Algorithm and Linear Complexity
Implementation of the Rubio-Sweedler-Taylor Algorithm for computing the Groebner basis of the ideal of linear recurring polynomials on a multidimensional periodic array with entries in a field.

The code has two dependencies: 
<ul>
  <li> 
    <code>Blitz++</code>. A C++ template class library that provides high-performance multidimensional array containers, found at <a href="https://github.com/blitzpp/blitz">https://github.com/blitzpp/blitz</a>
  </li>
  <li>
    <code>NTL</code> A high-performance, portable C++ library providing data structures and algorithms for manipulating signed, arbitrary length integers, and for vectors, matrices, and polynomials over the integers and over finite fields, found at 
    <a href="https://libntl.org/">https://libntl.org/</a>
  </li>
</ul>


Most of the action occurs in the file <code>MultiDimArray.cpp</code>. This a class consisting of a multidimensional array packed with additional information such as its size, periodicity, and linear complexity.
The <code>RST()</code> method is the member function implementing the Rubio-Sweedler-Taylor algorithms, and in the process of obtaining the Groebner basis, it computes the linear complexity.


The file <code>MultiDimArray_GF2.cpp</code> is an especialization of <code>MultiDimArray.cpp</code> for arrays with entries in GF(2), i.e., binary arrays. This runs way faster than the general <code>MultiDimArray.cpp</code>.

<code>Monomial.cpp</code> and <code>MultivarPolynomial.cpp</code> are axiliary classes used in <code>MultiDimArray.cpp</code>.
The rest of the <code>cpp</code> files are clients of the <code>MultiDimArray</code> class. 

This software is distributed under the terms of the MIT License.
