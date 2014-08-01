.. 
   py_qcode documentation master file, created by
   sphinx-quickstart on Wed Apr 16 16:30:40 2014.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

py_qcode: decoding topological error-correcting codes in Python
===============================================================

Contents:
---------

.. toctree::
   :maxdepth: 2

   simulation
   lattice
   code
   error
   decoder
   utils

Introduction
============

This library is concerned with simulating the decoding of topological quantum error-correcting codes. This is useful for determining their threshold error rates and logical error rates in the regime where useful error correction can be performed. In order to simulate the act of decoding, the following components are necessary:

 + **Input/Output**: Routines to input necessary values such as the name of the lattice geometry, size, error-correcting code to be used, etc., and save the resulting  output to a file.
 
 + **Lattice**: A set of points on which the qubits and stabilizers of the code are defined. 
 
 + **Error Model**: A set of Gottesman-Knill-compliant operators and probabilities with which they are applied to the qubits in the lattice. 
 
 + **Error-Correcting Code**: A set of stabilizer generators, along with their support, which produces error syndromes given an error.

 + **Decoder**: A rule for inferring the original error, given a syndrome.

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
