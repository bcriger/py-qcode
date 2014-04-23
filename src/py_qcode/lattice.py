"""
Code Comparison Project, April 2014
Ben Criger
Shameless plagiarism from Bravyi/Haah
"""

__all__ = ['Point', 'Lattice', 'SquareLattice', 'SquareOctagonLattice', 'UnionJackLattice']

class Point(object):
    r"""
    Represents a point in two or three dimensions. Normally, I'd use a
    :type namedtuple:, because I'm not a psychopath, but I want to use
    default arguments without getting too fancy.
    :param coords: co-ordinates of the point in question.
    :type coords: tuple of ints, length 2 or 3
    """
    def __init__(self, coords):

        if len(coords) not in [2, 3]:
            raise ValueError("Point must contain 2 or 3 co-ordinates,"\
                " you entered: {0}".format(coords))
        
        if not all([isinstance(coord,int) for coord in coords]):
            raise ValueError("Input tuple must be nothin' but ints,"\
                " you entered: {0}".format(coords))

        self.coords = coords
    
    def __hash__(self):
        """
        A hash function for points is necessary to store
        :class:`py-qcode.Point`s in sets.
        """
        return hash(self.coords)

class Lattice:
    """
    Wraps a collection of points. Superclass to SquareLattice,
    SquareOctagonLattice, UnionJackLattice, whatever other convenient
    lattices I put in. I have to learn about super/sub-classes before
    I get this going, though.

    Represents a 2D/3D lattice of points with integer 
    co-ordinates on which a stabilizer code can be defined. 
    :param sz_tpl: linear dimensions of the lattice.
    :type sz_tpl: tuple, length 2 or 3, containing integers.
    :param is_3D: indicates whether the lattice is to be two- or three-
    dimensional.
    :type is_3D: bool
    """
    def __init__(self):
        pass

class SquareLattice(Lattice):
    """
    """
    def __init__(self, sz_tpl, is_3D = False):
        """
        """
        pass

class SquareOctagonLattice(Lattice):
    """
    """
    def __init__(self, sz_tpl, is_3D = False):
        """
        """
        pass

class UnionJackLattice(Lattice):
    """
    """
    def __init__(self, sz_tpl, is_3D = False):
        """
        """
        pass

