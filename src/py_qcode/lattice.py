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
    #Magic Methods
    def __init__(self, coords, error = None, syndrome = None, inferred_error=None):
        
        check_int_tpl(coords)
        self.coords = coords
        
        self.error = error
        self.syndrome = syndrome 
        self.inferred_error = inferred_error

    def __hash__(self):
        """
        A hash function for points is necessary to store :class:`py_qcode.Point`s in sets or dictionaries.
        """
        return hash((self.coords, self.error, self.syndrome, self.inferred_error))

    def __repr__(self):
        rtrn_str = "Point at "+str(self.coords)
        
        if self.error is not None:
            rtrn_str += ", contains error " + str(self.error)
        if self.syndrome is not None:
            rtrn_str += ", contains syndrome " + str(self.syndrome)
        if self.inferred_error is not None:
            rtrn_str += ", contains inferred error " + str(self.inferred_error)
        
        return rtrn_str

    def __len__(self):
        return len(self.coords)

    #Pragmatism!
    def clear(self):
        self.syndrome = None
        self.error = None

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
    #Magic Methods
    def __init__(self, points, dim, dist=None, is_ft=False, size = None, is_dual = False):
        
        self.points = points
        self.dim = dim
        self.dist = dist
        self.is_ft = is_ft
        self.size = size
        self.is_dual = is_dual

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

    #Muggle Methods
    def clear(self):
        for point in self.points:
            point.clear()

class SquareLattice(Lattice):
    """
    Represents a lattice in which qubits are placed on the edges of a grid of squares with size given by `sz_tpl`. 

    :param rough_sides: Denotes which, if any, of the sides of the lattice are to have 'rough' boundary conditions. Values in ``rough_sides`` must be drawn from ``['u', 'd', 'r', 'l', 'f', 'b']`` (up, down, left, right, front, back). Default is `('u','r')`.

    :type rough_sides: tuple of strings
    """
    def __init__(self, sz_tpl, is_dual = False, is_ft = False, closed_boundary = True, rough_sides = ('u', 'r')):
        
        dim = len(sz_tpl)
        x_len, y_len = sz_tpl[:2] 

        if is_dual:
            points_2d = map(Point, sym_coords(x_len, y_len))
            dist = lambda coord1, coord2, synd_type: sum([min([abs(a - b) % (2 * sz),
                                                    (2 * sz - abs(a - b)) % (2 * sz)]) 
                                                for a, b, sz in 
                                                zip(coord1, coord2, sz_tpl)])
        else:
            points_2d = map(Point, skew_coords(x_len, y_len))
            dist = None

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
        self.is_dual = is_dual
        
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
        a method of Lattice, overridden by SquareLattice.
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

    def min_distance_path(self, dual_start, dual_end):
        """
        Returns a canonical minimum distance path on the _primal_ lattice 
        between two points of the _dual_ lattice. This is based on cutting
        as few corners as possible on the minimum-size bounding box. The 
        final path is the union of minimum distance paths between start, 
        corners, and end.
        """
        n_d = len(dual_end)
        
        corners = []
        points_on_path = []
        
        #print "dual_start: " + str(dual_start)
        #print "dual_end: " + str(dual_end)

        #make one corner for every co-ordinate that is different:
        corners.append(dual_start)
        for idx in range(n_d - 1):
            if dual_start[idx] != dual_end[idx]:
                corners.append(dual_start[ : idx + 1] + dual_end[ idx + 1 : ])
        corners.append(dual_end)

        #print "corners before removal: " + str(corners)

        #Patches on patches
        for idx, corner in enumerate(corners):
            if corner in corners[idx + 1 : ]:
                corners.remove(corner)
        
        #print "corners after: " + str(corners)

        for idx in range(len(corners) - 1):
            new_points = self._min_between(corners[idx], corners[idx + 1])
            #print new_points
            points_on_path += new_points

        return points_on_path
        
    def _min_between(self, dual_start, dual_end):
        """
        Assuming that two points differ by one co-ordinate only, 
        produces the minimum-length path between the two points.
        This path is straight, iterating over only one co-ordinate.
        """
        for idx, coord in enumerate(dual_start):
            #First index which differs ought to be only
            #index which differs:
            if dual_end[idx] != coord:
                special_idx = idx
                break

        #sort (dual_start, dual_end) using different element as key
        dual_start, dual_end = sorted((dual_start, dual_end),
                                        key=lambda lst: lst[special_idx])

        #coordinates "between" start and end on primal lattice
        betweens_forward = range(dual_start[special_idx] + 1,
                                    dual_end[special_idx], 2)
        
        #coordinates that "loop around" the boundary of the lattice
        betweens_backward = range(dual_start[special_idx]-1, -1, -2) + \
        range(dual_end[special_idx] + 1, 2 * self.size[special_idx], 2)
        
        if len(betweens_forward) < len(betweens_backward):
            betweens = betweens_forward
        else:
            betweens = betweens_backward

        betweens = map(lambda elem: 
                        dual_start[:special_idx] + \
                        (elem,) + \
                        dual_start[special_idx + 1:],
                        betweens)

        return betweens


class SquareOctagonLattice(Lattice):
    """
    Represents a lattice in which qubits are placed on the corners of 
    squares/octagons. This lattice results from allowing each point on 
    the edge of a SquareLattice (TM) to expand into a square. Begins by
    creating a square grid of points according to a size tuple, 
    identically to SquareLattice, with each coordinate passed through 
    an affine map. Each point in this lattice is then 'graduated' to a
    square consiting of nearest neighbours.  
    """
    def __init__(self, sz_tpl, is_dual = False, is_ft = False, closed_boundary = True, rough_sides = ('u', 'r')):
        dim = len(sz_tpl)

        try:
            x_len, y_len = sz_tpl
        except ValueError:
            raise ValueError("Only 2D is supported for now!")

        #We apply this affine map in order to ensure that the leftmost 
        #(bottom) column (row) of points is at co-ordinate 0:
        sq_cntrs = squoct_square_centers(x_len, y_len)

        squoct_coords = []
        for coord_pair in sq_cntrs:
            x, y = coord_pair
            squoct_coords.extend(_square_neighbourhood(x, y))

        points = map(Point, squoct_coords)
        super(SquareOctagonLattice, self).__init__(points, dim, dist, is_ft)
        
        #max coordinate value is derived by a change of co-ordinates, 
        #adding 1 to account for neighbourhoods
        
        self.size = x_len, y_len
        max_x, max_y = 2 * x_len - 1, 2 * y_len - 1
        
        #total_size is the values to mod by for the boundary conditions
        #Largest center co-ordinate + 1 (for the neighbour) + 1
        #(for the boundary)
        
        self.total_size = sq2oct(max_x) + 2, sq2oct(max_x) + 2
    
    def squares(self):
        nx, ny = self.size
        square_centers = _squoct_affine_map(skew_coords(nx, ny))
        point_list = []
        
        s_x, s_y = self.total_size
        for pt in square_centers:
            x, y = pt
            
            left, right = (x - 1) % s_x, (x + 1) % s_x
            down, up    = (y - 1) % s_y, (y + 1) % s_y
            
            point_list.append([self[(left, down)],  self[(left, up)],
                    self[(right, down)], self[(right, up)]])
        
        return point_list

    def z_octagons(self):
        nx, ny = self.size
        squoct_z_oct_centers = _squoct_affine_map(_odd_odds(nx, ny))
        s_x, s_y = self.total_size
        point_list = []
        for pt in squoct_z_oct_centers:
            x, y = pt
            
            xm2, xm1, xp1, xp2 = map(lambda x: x % s_x,
                                [(x - 2), (x - 1), (x + 1), (x + 2)])
            
            ym2, ym1, yp1, yp2 = map(lambda y: y % s_y,
                                [(y - 2), (y - 1), (y + 1), (y + 2)])
            
            point_list.append([(xm2, ym1), (xm2, yp1), (xm1, ym2),
                                (xp1, ym2), (xp1, yp2), (xm1, yp2),
                                (xp2, ym1), (xp2, yp1)])
        return point_list

    def x_octagons(self):
        #TODO: Un-copy this code, invoke a private function.
        nx, ny = self.size
        squoct_x_oct_centers = _squoct_affine_map(_even_evens(nx, ny))
        s_x, s_y = self.total_size
        point_list = []
        for pt in squoct_x_oct_centers:
            x, y = pt
            
            xm2, xm1, xp1, xp2 = map(lambda x: x % s_x,
                                [(x - 2), (x - 1), (x + 1), (x + 2)])
            
            ym2, ym1, yp1, yp2 = map(lambda y: y % s_y,
                                [(y - 2), (y - 1), (y + 1), (y + 2)])
            
            point_list.append([(xm2, ym1), (xm2, yp1), (xm1, ym2),
                                (xp1, ym2), (xp1, yp2), (xm1, yp2),
                                (xp2, ym1), (xp2, yp1)])
        return point_list

class UnionJackLattice(Lattice):
    """
    Gives the dual lattice to the SquareOctagonLattice above. 
    """
    def __init__(self, sz_tpl, is_dual = False, is_ft = False, closed_boundary = True, rough_sides = ('u', 'r')):
        dim = len(sz_tpl)
        
        try:
            x_len, y_len = sz_tpl
        except ValueError:
            raise ValueError("Only 2D is supported for now!")

        point_list = _squoct_affine_map(all_coords(x_len, y_len))
        points = map(Point, point_list)

        def dist(pt1, pt2, synd_type):
            """
            This function is complicated because the number of errors 
            required to produce a chain from one point to another is 
            not proportional to any actual distance; it has to ensure
            that octagons on the path are satisfied. To determine the 
            weight of a prospective error, we begin by isolating an
            easy case, then we relate the more difficult cases back to
            it.
            """
            
            #Return appropriate distance given co-ordinate types
            if is_sq_cent(pt1):
                if is_sq_cent(pt2):
                    return square_square_dist(pt1, pt2, synd_type,
                                                sz_tpl)
                else:
                    return square_octagon_dist(pt1, pt2, synd_type,
                                                sz_tpl)
            else:
                if is_sq_cent(pt2):
                    return square_octagon_dist(pt2, pt1, synd_type,
                                                sz_tpl)
                else:
                    return octagon_octagon_dist(pt1, pt2, synd_type,
                                                sz_tpl)
            
            raise ValueError("Co-ordinates {0} could not be"+\
                " identified as square or octagon centers."\
                .format([pt1, pt2]))


        super(UnionJackLattice, self).__init__(points, dim, dist, is_ft)
        
        #max coordinate value is derived by a change of co-ordinates, 
        #adding 1 to account for neighbourhoods
        self.size = x_len, y_len

        max_x, max_y = 2 * x_len - 1, 2 * y_len - 1
        
        #total_size is the values to mod by for the boundary conditions
        #Largest center co-ordinate + 1 (for the boundary)
        
        self.total_size = sq2oct(max_x) + 1, sq2oct(max_x) + 1

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

_evens = lambda n: range(0, 2 * n, 2)

_odds = lambda n: range(1, 2 * n + 1, 2)

_even_odds = lambda nx, ny: list(it.product(_evens(nx), _odds(ny)))

_odd_evens = lambda nx, ny: list(it.product(_odds(nx), _evens(ny)))

_even_evens = lambda nx, ny: list(it.product(_evens(nx), _evens(ny)))

_odd_odds = lambda nx, ny: list(it.product(_odds(nx), _odds(ny)))

def sym_coords(nx, ny):
    """
    Convenience function for square lattice definition, returns all pairs of co-ordinates on an n-by-n lattice which are both even or both odd. 
    """
    return _even_evens(nx, ny) + _odd_odds(nx, ny)


def skew_coords(nx, ny):
    """
    Convenience function for square lattice definition, returns all pairs of co-ordinates on an n-by-n lattice which are "even-odd" or "odd-even". 
    """
    return _even_odds(nx, ny) + _odd_evens(nx, ny)

def all_coords(nx, ny):
    """
    Returns all points on a square lattice, used to determine the dual lattice of the SquareOctagonLattice.
    """
    return sym_coords(nx, ny) + skew_coords(nx, ny)

def sq2oct(coord):
    """
    Takes a coordinate from the initial square lattice and produces a 
    new coordinate such that square neighbourhoods don't collide.
    """
    return 3 * coord + 1

def oct2sq(coord):
    """
    Is the inverse of sq2oct. This function can only be called on 
    points in the dual of the SquareOctagonLattice, so it throws an 
    error if it is fed other co-ordinates.
    """
    if (coord - 1) % 3 == 0:
        return (coord - 1)/3
    else:
        raise ValueError('Invalid co-ordinate: {0}'.format(coord))

def is_sq_cent(coord):
    """
    Squares are defined on the skew_coords of the virtual square 
    lattice. With this in mind, we invert the affine map from earlier,
    and sum the coordinates to determine if the result is odd.
    """
    return bool(sum(oct2sq(coord))%2)

#If it's not an square center, it's an octagon center.
is_oct_cent = lambda coord: not(is_sq_cent(coord))

def _squoct_affine_map(tpl_lst):
    """
    Maps sq2oct onto all elements in a list of tuples. Used in the big
    coordinate change. 
    """
    return map(lambda tpl: map(sq2oct, tpl), tpl_lst)

def straight_octagon_dist(x1, x2, sz):
    """
    Input two numbers which are raw one-d coordinates of octagons 
    whose other coordinates are identical. The function outputs the 
    weight of an error chain joining the two octagons.
    """
    diff = abs(x2 - x1)
    log_op_weight = 2 * sz #virtual lattice size
    return min([diff / 3, (log_op_weight - diff) / 3]) #Correcting for affine map

def octagon_octagon_dist(coord1, coord2, sz_tpl):
    """
    maps straight_octagon_dist to pairs of n-D coordinates, providing
    a consistent distance between any two octagons.
    """
    return sum(map(lambda tpl: straight_octagon_dist(*tpl),
                    zip(coord1, coord2, sz_tpl)))

def square_octagon_dist(sq_coord, oct_coord, synd_type, sz_tpl):
    """
    This function finds the number of errors necessary to form a chain
    between a square and an octagon, given the syndrome type (which
    determines the appropriate neighbouring octagons) and the 
    coordinates of the two points in question.

    It does this by minimization over the two cases, each corresponding
    to a neighbouring octagon.
    """

    test_octagons = appropriate_neighbours(sq_coord)
    return min(map(lambda o_c: 
        octagon_octagon_dist(o_c, oct_coord, sz_tpl), 
        test_octagons))

def square_square_dist(coord1, coord2, synd_type, sz_tpl):
    """
    This function finds the number of errors necessary to form a chain
    between two squares, given the syndrome type (which
    determines the appropriate neighbouring octagons) and the 
    coordinates of the two points in question.

    It does this by minimization over four cases, each corresponding
    to a pair of neighbouring octagons.
    """
    #Non-descriptive variable names for the octagons we're going to use
    ao1, ao2 = test_octagons(coord1)
    bo1, bo2 = test_octagons(coord2)
    return min(map(lambda se: octagon_octagon_dist(*se, sz_tpl), 
        [(ao1, bo1), (ao1, bo2), (ao2, bo1), (ao2, bo2)]))

def appropriate_neighbours(sq_coord, synd_type, sz_tpl):
    """
    If the square is in an even row on the virtual square lattice, 
    then the 'X'-stabilizer neighbours are to the left and right, with
    the 'Z'-stabilizer neighbours up and down. For squares on odd rows
    of the virtual square lattice, the opposite is true.
    """
    #Make sure the syndrome makes sense
    if synd_type not in 'xXzZ':
        raise ValueError('Weird syndrome type on '\
                        +'SquareOctagonLattice: {0}'.format(synd_type))

    synd_type = synd_type.upper()
    
    x, y = oct2sq(sq_coord)
    sz_x, sz_y = sz_tpl
    sz_x *= 2; sz_y *= 2 #Value to mod by on virtual lattice.
    
    #use Y-value and synd_type to determine neighbours
    y_even = bool(y % 2)

    if y_even:
        if synd_type = 'X':
            virtual_neighbourhood = [(x, (y - 1) % sz_y),(x, (y + 1) % sz_y)]
        else:
            virtual_neighbourhood = [((x + 1) % sz_x, y),((x - 1) % sz_x, y)]
    else:
        if synd_type = 'X':
            virtual_neighbourhood = [((x + 1) % sz_x, y),((x - 1) % sz_x, y)]
        else:
            virtual_neighbourhood = [(x, (y - 1) % sz_y),(x, (y + 1) % sz_y)]

    return _squoct_affine_map(virtual_neighbourhood)