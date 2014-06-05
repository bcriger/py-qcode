"""
Code Comparison Project, April 2014
Ben Criger
Shameless plagiarism from Bravyi/Haah
"""

import itertools as it

__all__ = ['Point', 'Lattice', 'SquareLattice2D', 'SquareOctagonLattice2D', 'UnionJackLattice2D']

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

    :type coords: tuple of ints
    
    :param error: A value which denotes an error. An ``ErrorCorrectingCode`` must check that this value corresponds to an operator which can be translated into a syndrome.
    
    :type error: any
    
    :param syndrome: A value which denotes an syndrome. A ``Decoder`` must check that this value corresponds to an operator which can be translated into a syndrome.
    
    :type syndrome: any
    """
    def __init__(self, coords, error = None, syndrome = None):
        
        check_int_tpl(coords)
        self.coords = coords
        
        self.error = error
        self.syndrome = syndrome 
    
    def __hash__(self):
        """
        A hash function for points is necessary to store :class:`py_qcode.Point`s in sets or dictionaries.
        """
        return hash((self.coords, self.error, self.syndrome))

    def __repr__(self):
        rtrn_str = "Point at "+str(self.coords)
        if self.error is not None:
            rtrn_str += ", contains error " + str(self.error)
        if self.syndrome is not None:
            rtrn_str += ", contains syndrome " + str(self.syndrome)
        return rtrn_str

    def len(self):
        return len(self.coords)

class Lattice(object):
    """
    A collection of points. Superclass to ``SquareLattice``, ``SquareOctagonLattice``, ``UnionJackLattice``, whatever other convenient lattices I put in. 

    Represents a arbitrary-dimensional lattice of points with integer co-ordinates on which a stabilizer code can be defined. Note that, although the word "lattice" is used to describe these objects, the only requirement is that its constituent points have co-ordinates. No property of the graph structure is assumed, especially planarity. 

    :param points: collection of points on the lattice.
    
    :type points: list of :class:`py_qcode.Point` objects
    
    :param dist: Returns the distance between two points. Note: In order to perform MWPM decoding, only the dual lattice needs a distance function. 

    :type dist: function

    :param is_ft: indicates whether the lattice is to possess an extra dimension for fault-tolerant decoding.
    
    :type is_ft: bool

    :param closed_boundary: Indicates whether to identify the Nth co-ordinate with the zeroth co-ordinate in every dimension.

    :type closed_boundary: bool
    """
    def __init__(self, points, dim, dist=None, is_ft=False):
        
        self.points = points
        self.dim = dim
        self.dist = dist
        self.is_ft = is_ft

class SquareLattice(Lattice):
    """
    Represents a lattice in which qubits are placed on the edges of a grid of squares with size given by `sz_tpl`. 

    :param rough_sides: Denotes which, if any, of the sides of the lattice are to have 'rough' boundary conditions. Values in ``rough_sides`` must be drawn from ``['u', 'd', 'r', 'l', 'f', 'b']`` (up, down, left, right, front, back). Default is `('u','r')`.

    :type rough_sides: tuple of strings
    """
    def __init__(self, sz_tpl, is_dual=False, is_ft = False, closed_boundary=True, rough_sides = ('u', 'r')):
        
        dim = len(sz_tpl)
        x_len, y_len = sz_tpl[:2] 

        #TODO Add convenience functions to make these more legible. 
        if is_dual:
            points_2d = map(Point, list(it.product([2*j for j in range(x_len)], [2*j for j in range(y_len)])) + list(it.product([2*j+1 for j in range(x_len)], [2*j+1 for j in range(y_len)])))
            #TODO define correct distance function given these co-ordinates
            dist=None
        else:
            points_2d = map(Point, list(it.product([2*j for j in range(x_len)], [2*j+1 for j in range(y_len)])) + list(it.product([2*j+1 for j in range(x_len)], [2*j for j in range(y_len)])))
            dist=None

        if is_ft:
            if len(sz_tpl) != 3:
                raise ValueError("Square lattices for fault-tolerant simulations must be 3D.")
            z_len = sz_tpl[2]
            points = layer(points_2d, z_len)
        else:
            if len(sz_tpl) != 2:
                raise ValueError("Square lattices for non-fault-tolerant simulations must be 2D.")
            points = points_2d

        super(SquareLattice, self).__init__(points, dim, dist, is_ft)
        
        if all([side in SIDES for side in rough_sides]):
            self.rough_sides = rough_sides
        else: 
            raise ValueError(("rough_sides must be in the list {0}." +\
                "You entered: {1}").format(SIDES, rough_sides))

class SquareOctagonLattice(Lattice):
    """
    Represents a lattice in which qubits are placed on the corners of squares and octagons. 
    """
    def __init__(self, sz_tpl, is_3D=False, closed_boundary=True):
        """
        """
        pass

class UnionJackLattice(Lattice):
    """
    Represents a lattice in which qubits are placed on the intersections of the diagonals of squares, as well as their corners. 
    """
    def __init__(self, sz_tpl):
        """
        """
        pass

class CubicLattice(Lattice):
    """
    Represents a lattice in which qubits are placed on the intersections of the diagonals of squares, as well as their corners. 
    """
    def __init__(self, sz_tpl):
        """
        """
        pass

## Convenience Functions ##
def check_int_tpl(coords):
    if not all([isinstance(coord,int) for coord in coords]):
        raise ValueError("Input tuple must be nothin' but ints,"\
            " you entered: {0}".format(coords))
    pass

def promote(point, new_coord):
    """
    Adds a new coordinate to a point, preserving the error and syndrome stored therein.
    """
    pt_attrs = point.__dict__.values()
    if type(new_coord) is int:
        pt_attrs[0] += (new_coord,)
    else:
        raise TypeError("New coordinate must be an int.")
    return Point(*pt_attrs)

def layer(points, new_len):
    """
    Naively extends a collection of points into a new dimension by producing `new_len` new copies of the points at integral co-ordinates.
    """
    new_pt_lst = []
    for stratum in range(new_len):
        new_pt_lst.append(map(lambda pt: promote(pt, stratum), points))
    return new_pt_lst