py-qcode
========

Decoding topological quantum error-correcting codes in Python

Disclaimer 
---------- 
This code is in constant development as of April
22, 2014. Core functionality may not have yet been written. The
documentation seen here is an attempt at agile development, in which 
the documentation acts as the specification. 

Other Disclaimer 
---------- 
This code is deprecated as of July 22, 2016.
I've managed to publish a paper using it, but it's not good enough. 
The fundamental mistake on which this is all based is that Points can 
be represented as objects.
They can't.
We have to use Lattice objects to write code that's write-able, and 
not have to work around all the time. 
I'm going to try to get a sensible way to simulate topological codes
going with [sparse_pauli](github.com/bcriger/sparse_pauli) and some
(as of yet) private code. 
Feel free to send pull requests or raise issues, I'll even address
them, I just don't think it's worth it. 
