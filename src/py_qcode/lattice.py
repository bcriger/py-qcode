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

    def __len__(self):
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

    def __getitem__(self, key):
        if len(key) != len(self.points[0]):
            raise ValueError("key must be length: " + str(len(self.points[0])))
        for point in self.points:
            if point.coords == tuple(key):
                return point
        raise KeyError("Point not found on lattice; key: "+ str(key))
    
    def __repr__(self):
        pts = map(lambda pt: repr(pt), self.points)
        return '[' + ',\n '.join(pts) + ']'

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
            points_2d = map(Point, sym_coords(x_len, y_len))
            #TODO define correct distance function given these co-ordinates
            dist=None
        else:
            points_2d = map(Point, skew_coords(x_len, y_len))
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
        self.size = sz_tpl
        
        if all([side in SIDES for side in rough_sides]):
            self.rough_sides = rough_sides
        else: 
            raise ValueError(("rough_sides must be in the list {0}." +\
                "You entered: {1}").format(SIDES, rough_sides))
    
    def neighbours(self, location):
        """
        Returns a list of points which are one unit of distance away from a given location.
        Convenience method used to define stars and plaquettes below.

        June 6, 2014: Only supports closed_boundary 

        TODO: Make this depend on the distance function, so that it can be made
        a method of Lattice, instead of SquareLattice.
        """
        x, y = location
        x_sz, y_sz = self.size
        
        left  = (x - 1) % (2 * x_sz)
        right = (x + 1) % (2 * x_sz)
        up    = (y + 1) % (2 * y_sz)
        down  = (y - 1) % (2 * y_sz)
        
        return (self[right,y],self[x,up],self[left,y],self[x,down])

    def stars(self):
        return map(self.neighbours, _even_evens(*self.size))
    
    def plaquettes(self):
        return map(self.neighbours, _odd_odds(*self.size))

class SquareOctagonLattice(Lattice):
    """
    Represents a lattice in which qubits are placed on the corners of squares and octagons. 
    """
    def __init__(self, sz_tpl, is_3D=False, closed_boundary=True):
        """
        """
        pass

    def squares(self):
        pass

    def z_octagons(self):
        pass

    def x_octagons(self):
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

_evens = lambda n: range(0, 2*n, 2)

_odds = lambda n: range(1, 2*n+1, 2)

_even_odds = lambda nx, ny: list(it.product(_evens(nx), _odds(ny)))

_odd_evens = lambda nx, ny: list(it.product(_odds(nx), _evens(ny)))

_even_evens = lambda nx, ny: list(it.product(_evens(nx), _evens(ny)))

_odd_odds = lambda nx, ny: list(it.product(_odds(nx), _odds(ny)))

def sym_coords(nx, ny):
    """
    Convenience function for square lattice definition, returns all pairs of co-ordinates on an n-by-n lattice which are both even or both odd. 
    """
    return  _even_evens(nx, ny) + _odd_odds(nx, ny)


def skew_coords(nx, ny):
    """
    Convenience function for square lattice definition, returns all pairs of co-ordinates on an n-by-n lattice which are "even-odd" or "odd-even". 
    """
    return _even_odds(nx, ny) + _odd_evens(nx, ny)