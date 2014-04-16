"""
Code Comparison Project, April 2014
Ben Criger
Shameless plagiarism from Bravyi/Haah
"""

import networkx as nx
from collections import namedtuple

Point = namedtuple('Point', ['x','y','z']) 

dist = lambda d, L: min(abs(d), L-abs(d)) #distance on closed 1d ring

class SquareLattice:
    """
    This class represents a lattice of points with definite 
    co-ordinates on which a stabilizer code can be defined.
    """
def __init__(self):
