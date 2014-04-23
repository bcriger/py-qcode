"""
Code Comparison Project, April 2014
Ben Criger
Shameless plagiarism from Bravyi/Haah
"""

import itertools as it

__all__ = ['Point', 'Lattice', 'SquareLattice', 'SquareOctagonLattice', 'UnionJackLattice']

##constants##
SIDES = ['u', 'd', 'r', 'l', 'f', 'b'] #up, down, left, right, front, back

class Point(object):
    r"""
    Represents a point in two or three dimensions. Normally, I'd use a
    ``namedtuple``, because I'm not a psychopath, but I want to use
    default arguments without getting too fancy. Each point can also 
    contain a value denoting the error which an 
    ``ErrorModel`` has applied to that point, and/or a
    syndrome which results from measurement of an
    ``ErrorCorrectingCode`` on the lattice.
    
    :param coords: co-ordinates of the point in question.

    :type coords: tuple of ints, length 2 or 3
    
    :param error: A value which denotes an error. An ``ErrorCorrectingCode`` must check that this value corresponds to an operator which can be translated into a syndrome.
    
    :type error: any
    
    :param syndrome: A value which denotes an syndrome. A ``Decoder`` must check that this value corresponds to an operator which can be translated into a syndrome.
    
    :type syndrome: any
    """
    def __init__(self, coords, error = None, syndrome = None):
        
        check_int_23_tpl(coords)
        self.coords = coords
        
        self.error = error
        self.syndrome = syndrome 
    
    def __hash__(self):
        """
        A hash function for points is necessary to store
        :class:`py-qcode.Point`s in sets.
        """
        return hash((self.coords, self.error, self.syndrome))

class Lattice:
    """
    A collection of points. Superclass to ``SquareLattice``,
    ``SquareOctagonLattice``, ``UnionJackLattice``, whatever other convenient
    lattices I put in. 

    Represents a 2D/3D lattice of points with integer 
    co-ordinates on which a stabilizer code can be defined.

    :param sz_tpl: linear dimensions of the lattice.
    
    :type sz_tpl: tuple, length 2 or 3, containing integers.
    
    :param is_3D: indicates whether the lattice is to be two- or three- dimensional.
    
    :type is_3D: bool

    :param closed_boundary: Indicates whether to identify the Nth co-ordinate with the zeroth co-ordinate in every dimension.

    :type closed_boundary: bool
    """
    def __init__(self, sz_tpl, is_3D=False, closed_boundary = True):
        check_int_23_tpl(sz_tpl)
        self.size = sz_tpl
        self.is_3D = is_3D


class SquareLattice(Lattice):
    """
    Represents a lattice in which qubits are placed on the edges.

    :param rough_sides: Denotes which, if any, of the sides of the lattice are to have 'rough' boundary conditions. Values in ``rough_sides`` must be drawn from ``['u', 'd', 'r', 'l', 'f', 'b']`` (up, down, left, right, front, back).

    :type rough_sides: tuple of strings
    """
    def __init__(self, sz_tpl, is_3D = False, closed_boundary=True, rough_sides = ('f', 'b')):
        super(SquareLattice, self).__init__(self, sz_tpl, is_3D,
                                            closed_boundary)
        
        if all([side in SIDES for side in rough_sides]):
            self.rough_sides = rough_sides
        else: 
            raise ValueError(("rough_sides must be in the list {0}." +\
                "You entered: {1}").format(SIDES, rough_sides))


class SquareOctagonLattice(Lattice):
    """
    """
    def __init__(self, sz_tpl):
        """
        """
        pass

class UnionJackLattice(Lattice):
    """
    """
    def __init__(self, sz_tpl):
        """
        """
        pass

## Convenience Functions ##
def check_int_23_tpl(coords):
    if len(coords) not in [2, 3]:
        raise ValueError("Input must contain 2 or 3 co-ordinates,"\
            " you entered: {0}".format(coords))
    
    if not all([isinstance(coord,int) for coord in coords]):
        raise ValueError("Input tuple must be nothin' but ints,"\
            " you entered: {0}".format(coords))
    pass