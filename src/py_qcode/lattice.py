"""
Code Comparison Project, April 2014
Ben Criger
"""

from itertools import product
from math import floor
from ctypes import c_ushort, c_char, cdll
import os.path
# FIXME: Detect whether lattice is being imported or py_qcode is being
# installed:
path_to_lib = os.path.join(os.path.dirname(
    __file__), 'libqcode_dist.so')
# path_to_lib = os.path.abspath('libqcode_dist.so')
# libqcode_dist = cdll.LoadLibrary(os.path.join(me, 'libqcode_dist.so'))
libqcode_dist = cdll.LoadLibrary(path_to_lib)

__all__ = ['Point', 'Lattice', 'SquareLattice', 'SquareOctagonLattice',
           'UnionJackLattice']

#"""
__all__.extend(
    ['skew_coords', '_squoct_affine_map', 'straight_octagon_dist',
     'straight_octagon_path', 'octagon_octagon_path',
     'octagon_octagon_dist', 'square_octagon_dist',
     'square_octagon_path', 'square_square_dist',
     'square_square_path', 'appropriate_neighbours', 'nwes_pairs'])
#"""

# constants##
SIDES = ['u', 'd', 'r', 'l', 'f', 'b']
    # up, down, left, right, front, back


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

    :param error: A value which denotes an error. An 
    ``ErrorCorrectingCode`` must check that this value corresponds to 
    an operator which can be translated into a syndrome.

    :type error: any

    :param syndrome: A value which denotes an syndrome. A ``Decoder``
     must check that this value corresponds to an operator which can 
     be translated into a syndrome.

    :type syndrome: any
    """
    # Magic Methods
    def __init__(self, coords, error=None,
                 syndrome=None, inferred_error=None):

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
        rtrn_str = "Point at " + str(self.coords)

        if self.error is not None:
            rtrn_str += ", contains error " + str(self.error)
        if self.syndrome is not None:
            rtrn_str += ", contains syndrome " + str(self.syndrome)
        if self.inferred_error is not None:
            rtrn_str += ", contains inferred error " + \
                str(self.inferred_error)

        return rtrn_str

    def __len__(self):
        return len(self.coords)

    # Pragmatism!
    def clear(self):
        self.syndrome = None
        self.error = None


class Lattice(object):

    """
    A collection of points. Superclass to ``SquareLattice``,
    ``SquareOctagonLattice``, ``UnionJackLattice``, whatever other
    convenient lattices I put in.

    Represents a arbitrary-dimensional lattice of points with integer
    co-ordinates on which a stabilizer code can be defined. Note that,
    although the word "lattice" is used to describe these objects, the
    only requirement is that its constituent points have co-ordinates.
    No property of the graph structure is assumed, especially planarity.

    :param points: collection of points on the lattice.

    :type points: list of :class:`py_qcode.Point` objects

    :param dist: Returns the distance between two points. Note: In
    order to perform MWPM decoding, only the dual lattice needs a
    distance function.

    :type dist: function

    :param closed_boundary: Indicates whether to identify the Nth co-
    ordinate with the zeroth co-ordinate in every dimension.

    :type closed_boundary: bool
    """
    # Magic Methods
    def __init__(self, points, dim, dist=None, size=None,
                 is_dual=False):

        self.points = points
        self.dim = dim
        self.dist = dist
        self.size = size
        self.is_dual = is_dual

    def __getitem__(self, key):
        if len(key) != len(self.points[0]):
            raise ValueError("key must be length: " +
                             str(len(self.points[0])))
        for point in self.points:
            if point.coords == tuple(key):
                return point
        raise KeyError("Point not found on lattice; key: " + str(key))

    def __repr__(self):
        pts = map(lambda pt: repr(pt), self.points)
        return '[' + ',\n '.join(pts) + ']'

    # Muggle Methods
    def clear(self):
        for point in self.points:
            point.clear()


class SquareLattice(Lattice):

    """
    Represents a lattice in which qubits are placed on the edges of a
    grid of squares with size given by `sz_tpl`.

    :param rough_sides: Denotes which, if any, of the sides of the
    lattice are to have 'rough' boundary conditions. Values in
    ``rough_sides`` must be drawn from
    ``['u', 'd', 'r', 'l', 'f', 'b']`` (up, down, left, right, front,
    back). Default is `('u','r')`.

    :type rough_sides: tuple of strings
    """
    def __init__(self, sz_tpl, is_dual=False, closed_boundary=True,
                 rough_sides=('u', 'r')):

        dim = len(sz_tpl)
        x_len, y_len = sz_tpl[:2]

        if is_dual:
            points_2d = map(Point, sym_coords(x_len, y_len))
            """
            dist = lambda coord1, coord2, synd_type: \
                    sum([
                        min([abs(a - b) % (2 * sz),
                            (2 * sz - abs(a - b)) % (2 * sz)])
                                        for a, b, sz in
                                        zip(coord1, coord2, sz_tpl)])
            """
            def dist(coord1, coord2, synd_type):
                return sum([min([abs(a - b) % (2 * sz),
                            (2 * sz - abs(a - b)) % (2 * sz)])
                            for a, b, sz in
                            zip(coord1, coord2, sz_tpl)]) / 2
            '''
            def dist(coord1, coord2, synd_type):
                x1, y1 = map(c_ushort, coord1)
                x2, y2 = map(c_ushort, coord2)
                sz_x, sz_y = map(c_ushort, [x_len, y_len])
                return libqcode_dist.toric_dist(x1, y1, x2, y2,
                    sz_x, sz_y, c_char(synd_type))
            '''
        else:
            points_2d = map(Point, skew_coords(x_len, y_len))
            dist = None

        if len(sz_tpl) != 2:
            raise ValueError("Square lattices for " +
                             "non-fault-tolerant simulations must be 2D.")
        points = points_2d

        super(SquareLattice, self).__init__(points, dim, dist)
        self.size = sz_tpl
        self.is_dual = is_dual

        if all([side in SIDES for side in rough_sides]):
            self.rough_sides = rough_sides
        else:
            raise ValueError(("rough_sides must be in the list {0}." +
                              "You entered: {1}").format(SIDES, rough_sides))

    # Overwriting __getitem__ for increased speed:
    def __getitem__(self, coord_pair):
        x, y = coord_pair
        sz_x = self.size[0]
        shift = -(x % 2) if self.is_dual else x % 2 - 1
        return self.points[x * sz_x + (y + shift) / 2]

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

        left = (x - 1) % (2 * x_sz)
        right = (x + 1) % (2 * x_sz)
        up = (y + 1) % (2 * y_sz)
        down = (y - 1) % (2 * y_sz)

        #This determines gate order during noisy decoding
        return (self[x, up], self[left, y], self[right, y], self[x, down])

    def stars(self):
        return map(self.neighbours, _even_evens(*self.size))

    def plaquettes(self):
        return map(self.neighbours, _odd_odds(*self.size))

    def min_distance_path(self, dual_start, dual_end, synd_type=None):
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

        # print "dual_start: " + str(dual_start)
        # print "dual_end: " + str(dual_end)

        # make one corner for every co-ordinate that is different:
        corners.append(dual_start)
        for idx in range(n_d - 1):
            if dual_start[idx] != dual_end[idx]:
                corners.append(dual_start[
                               : idx + 1] + dual_end[idx + 1:])
        corners.append(dual_end)

        # print "corners before removal: " + str(corners)

        # Patches on patches
        for idx, corner in enumerate(corners):
            if corner in corners[idx + 1:]:
                corners.remove(corner)

        # print "corners after: " + str(corners)

        for idx in range(len(corners) - 1):
            new_points = self._min_between(
                corners[idx], corners[idx + 1])
            # print new_points
            points_on_path += new_points

        return points_on_path

    def _min_between(self, dual_start, dual_end):
        """
        Assuming that two points differ by one co-ordinate only,
        produces the minimum-length path between the two points.
        This path is straight, iterating over only one co-ordinate.
        """
        for idx, coord in enumerate(dual_start):
            # First index which differs ought to be only
            # index which differs:
            if dual_end[idx] != coord:
                special_idx = idx
                break

        # sort (dual_start, dual_end) using different element as key
        dual_start, dual_end = sorted((dual_start, dual_end),
                                      key=lambda lst: lst[special_idx])

        # coordinates "between" start and end on primal lattice
        betweens_forward = range(dual_start[special_idx] + 1,
                                 dual_end[special_idx], 2)

        # coordinates that "loop around" the boundary of the lattice
        betweens_backward = range(dual_start[special_idx] - 1, -1, -2) + \
            range(dual_end[special_idx] + 1, 2 * self.size[
                  special_idx], 2)

        if len(betweens_forward) < len(betweens_backward):
            betweens = betweens_forward
        else:
            betweens = betweens_backward

        betweens = map(lambda elem:
                       dual_start[:special_idx] +
                      (elem,) +
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
    def __init__(self, sz_tpl, is_dual=False,
                 closed_boundary=True, rough_sides=('u', 'r')):

        dim = len(sz_tpl)

        try:
            x_len, y_len = sz_tpl
        except ValueError:
            raise ValueError("Only 2D is supported for now!")

        dist = None  # Primal lattices don't need distance functions for now

        # Co-ordinates are placed in a sorted list, beginning with the
        # x-coordinates:
        total_x = 6 * x_len
        total_y = 6 * y_len
        coords = []
        x_cs = [x for x in range(total_x) if (x - 1) % 3 != 0]
        # primitive y-list requires x-dependent shift
        y_list = [y for y in range(2, total_y, 2) if y % 3 != 0]

        for x in x_cs:
            shift = -2 if x % 2 else 1
            coords.extend([(x, y + shift) for y in y_list])

        points = map(Point, coords)
        super(SquareOctagonLattice, self).__init__(points, dim, dist)

        self.size = x_len, y_len

        self.total_size = (total_x, total_y)

    def __getitem__(self, coord_pair):
        x, y = coord_pair
        sz_y = self.size[1] * 2
        num_blocks_ahead = (x + int(floor(x / 3))) / 2
        shift = -2 if x % 2 else 1
        num_elems_ahead = y - shift  # 2, 4, 8, 10, ...
        num_elems_ahead -= 2 * (1 + int(floor(y / 6)))  # 0, 2, 4, 6 ...
        num_elems_ahead /= 2  # 0, 1, 2, 3
        return self.points[num_blocks_ahead * sz_y + num_elems_ahead]

    def squares(self):
        nx, ny = self.size
        square_centers = _squoct_affine_map(skew_coords(nx, ny))
        point_list = []

        s_x, s_y = self.total_size
        for pt in square_centers:
            x, y = pt

            left, right = (x - 1) % s_x, (x + 1) % s_x
            down, up = (y - 1) % s_y, (y + 1) % s_y

            point_list.append([self[(left, down)], self[(left, up)],
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

            point_list.append(
                map(self.__getitem__, [
                    (xm2, ym1), (xm2, yp1), (xm1, ym2),
                    (xp1, ym2), (xp1, yp2), (xm1, yp2),
                    (xp2, ym1), (xp2, yp1)]))
        return point_list

    def x_octagons(self):
        # TODO: Un-copy this code, invoke a private function.
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

            point_list.append(
                map(self.__getitem__, [
                    (xm2, ym1), (xm2, yp1), (xm1, ym2),
                    (xp1, ym2), (xp1, yp2), (xm1, yp2),
                    (xp2, ym1), (xp2, yp1)]))
        return point_list

    def min_distance_path(self, dual_start, dual_end, synd_type):
        """
        Returns a canonical minimum distance path on the _primal_
        lattice between two points of the _dual_ lattice. This is based
        on finding the path between neighbouring octagons of whatever
        dual coordinates are fed in, and testing to see which sets of
        neighbours are optimal.
        """
        total_size = self.total_size
        if is_sq_cent(dual_start):
            if is_sq_cent(dual_end):
                path = square_square_path(
                    dual_start, dual_end, synd_type,
                    total_size)
            else:
                path = square_octagon_path(
                    dual_start, dual_end, synd_type,
                    total_size)
        else:
            if is_sq_cent(dual_end):
                path = square_octagon_path(
                    dual_end, dual_start, synd_type,
                    total_size)
            else:
                path = octagon_octagon_path(
                    dual_start, dual_end, total_size)

        return path


class UnionJackLattice(Lattice):

    """
    Gives the dual lattice to the SquareOctagonLattice above.
    """
    def __init__(self, sz_tpl, is_dual=True,
                 closed_boundary=True, rough_sides=('u', 'r')):
        dim = len(sz_tpl)

        try:
            x_len, y_len = sz_tpl
        except ValueError:
            raise ValueError("Only 2D is supported for now!")

        point_list = _squoct_affine_map(all_coords(x_len, y_len))
        points = map(Point, point_list)

        # max coordinate value is derived by a change of co-ordinates,
        # adding 1 to account for neighbourhoods
        max_x, max_y = 2 * x_len - 1, 2 * y_len - 1
        total_size = sq2oct(max_x) + 2, sq2oct(max_y) + 2

        def __getitem__(self, coord_pair):
            sz_y = self.size[1] * 2
            x, y = (x - 1) / 3, (y - 1) / 3
            return self.points[x * sz_y + y]
        '''
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
            if pt1 == pt2:
                return 0
            #Return appropriate distance given co-ordinate types
            if is_sq_cent(pt1):
                if is_sq_cent(pt2):
                    return square_square_dist(pt1, pt2, synd_type,
                                                total_size)
                else:
                    return square_octagon_dist(pt1, pt2, synd_type,
                                                total_size)
            else:
                if is_sq_cent(pt2):
                    return square_octagon_dist(pt2, pt1, synd_type,
                                                total_size)
                else:
                    return octagon_octagon_dist(pt1, pt2, total_size)

            raise ValueError("Co-ordinates {0} could not be"+\
                " identified as square or octagon centers."\
                .format([pt1, pt2]))
        '''
        def dist(pt1, pt2, synd_type):
            x1, y1 = map(c_ushort, pt1)
            x2, y2 = map(c_ushort, pt2)
            sz_x, sz_y = map(c_ushort, total_size)
            return libqcode_dist.squoct_dist(x1, y1, x2, y2, sz_x, sz_y,
                                             c_char(synd_type))
        #'''
        super(UnionJackLattice, self).__init__(points, dim, dist)

        self.size = x_len, y_len
        self.total_size = total_size


class CubicLattice(Lattice):

    """
    Represents a lattice in which qubits are placed on the intersections of the diagonals of squares, as well as their corners.
    """
    def __init__(self, sz_tpl):
        """
        """
        pass


# Convenience Functions ##
def check_int_tpl(coords):
    if not all([isinstance(coord, int) for coord in coords]):
        raise ValueError("Input tuple must be nothin' but ints,"
                         " you entered: {0}".format(coords))
    pass

_evens = lambda n: range(0, 2 * n, 2)

_odds = lambda n: range(1, 2 * n + 1, 2)

_even_odds = lambda nx, ny: map(tuple, list(
    product(_evens(nx), _odds(ny))))

_odd_evens = lambda nx, ny: map(
    tuple, list(product(_odds(nx), _evens(ny))))

_even_evens = lambda nx, ny: map(tuple, list(
    product(_evens(nx), _evens(ny))))

_odd_odds = lambda nx, ny: map(
    tuple, list(product(_odds(nx), _odds(ny))))


def sym_coords(nx, ny):
    """
    Convenience function for square lattice definition, returns all
    pairs of co-ordinates on an n-by-n lattice which are both even or
    both odd. Note that it iterates over all of the points in the 2D
    grid which is nx-by-ny large, this is so the list of returned
    coordinates is sorted.
    """
    symmetric_coordinates = []
    for x in range(2 * nx):
        if x % 2 == 0:
            for y in range(2 * ny):
                if y % 2 == 0:
                    symmetric_coordinates.append((x, y))
        else:
            for y in range(2 * ny):
                if y % 2 == 1:
                    symmetric_coordinates.append((x, y))
    return symmetric_coordinates


def skew_coords(nx, ny):
    """
    Convenience function for square lattice definition, returns all
    pairs of co-ordinates on an n-by-n lattice which are "even-odd" or
    "odd-even".Note that it iterates over all of the points in the 2D
    grid which is nx-by-ny large, this is so the list of returned
    coordinates is sorted.
    """
    skewed_coordinates = []
    for x in range(2 * nx):
        if x % 2 == 0:
            for y in range(2 * ny):
                if y % 2 == 1:
                    skewed_coordinates.append((x, y))
        else:
            for y in range(2 * ny):
                if y % 2 == 0:
                    skewed_coordinates.append((x, y))
    return skewed_coordinates


def all_coords(nx, ny):
    """
    Returns all points on a square lattice, used to determine the dual
    lattice of the SquareOctagonLattice.
    """
    all_coordinates = []
    for x in range(2 * nx):
        for y in range(2 * ny):
            all_coordinates.append((x, y))
    return all_coordinates


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
        return (coord - 1) / 3
    else:
        raise ValueError('Invalid co-ordinate: {0}'.format(coord))


def is_sq_cent(coord):
    """
    Squares are defined on the skew_coords of the virtual square
    lattice. With this in mind, we invert the affine map from earlier,
    and sum the coordinates to determine if the result is odd.
    """
    return bool(sum(map(oct2sq, coord)) % 2)

# If it's not a square center, it's an octagon center.
is_oct_cent = lambda coord: not(is_sq_cent(coord))


def _squoct_affine_map(tpl_lst):
    """
    Maps sq2oct onto all elements in a list of tuples. Used in the big
    coordinate change.
    """
    return map(lambda tpl: tuple(map(sq2oct, tpl)), tpl_lst)


def straight_octagon_dist(x1, x2, sz):
    """
    Input two numbers which are raw one-d coordinates of octagons
    whose other coordinates are identical. The function outputs the
    weight of an error chain joining the two octagons.
    """
    diff = abs(x2 - x1)
    diff /= 3
    log_op_weight = sz / 3  # virtual lattice size
    return min([diff, abs(log_op_weight - diff)])


def octagon_octagon_dist(coord1, coord2, total_size):
    """
    maps straight_octagon_dist to pairs of n-D coordinates, providing
    a consistent distance between any two octagons.
    """
    return sum(map(lambda tpl: straight_octagon_dist(*tpl),
                   zip(coord1, coord2, total_size)))


def square_octagon_dist(sq_coord, oct_coord, synd_type, total_size):
    """
    This function finds the number of errors necessary to form a chain
    between a square and an octagon, given the syndrome type (which
    determines the appropriate neighbouring octagons) and the
    coordinates of the two points in question.

    It does this by minimization over the two cases, each corresponding
    to a neighbouring octagon.
    """

    test_octagons = appropriate_neighbours(
        sq_coord, synd_type, total_size)
    neighbour_dist = min(map(lambda o_c:
                             octagon_octagon_dist(
                                 o_c, oct_coord, total_size),
                             test_octagons))

    return neighbour_dist + 1  # 1-qubit op to step from oct to square.


def square_square_dist(coord1, coord2, synd_type, total_size):
    """
    This function finds the number of errors necessary to form a chain
    between two squares, given the syndrome type (which
    determines the appropriate neighbouring octagons) and the
    coordinates of the two points in question.

    It does this by minimization over four cases, each corresponding
    to a pair of neighbouring octagons.
    """

    # print synd_type

    # Non-descriptive variable names for the octagons we're going to use
    ao1, ao2 = appropriate_neighbours(coord1, synd_type, total_size)
    bo1, bo2 = appropriate_neighbours(coord2, synd_type, total_size)
    neighbour_dist = min(
        map(lambda se: octagon_octagon_dist(se[0], se[1], total_size),
            [(ao1, bo1), (ao1, bo2), (ao2, bo1), (ao2, bo2)]))
    return neighbour_dist + 2  # 2 oct-to-square steps necessary.


def appropriate_neighbours(sq_coord, synd_type, total_size):
    """
    If the square is in an even row on the virtual square lattice,
    then the 'X'-stabilizer neighbours are to the left and right, with
    the 'Z'-stabilizer neighbours up and down. For squares on odd rows
    of the virtual square lattice, the opposite is true.
    """
    # Make sure the syndrome makes sense
    if not (synd_type in 'xXzZ'):
        raise ValueError('Weird syndrome type on '
                         + 'SquareOctagonLattice: {0}'.format(synd_type))

    synd_type = synd_type.upper()

    x, y = map(oct2sq, sq_coord)
    sz_x, sz_y = total_size[0] / 3, total_size[
        1] / 3  # size of virtual lattice.

    # use Y-value and synd_type to determine neighbours
    y_even = not(bool(y % 2))

    if y_even:
        if synd_type == 'X':
            virtual_neighbourhood = [(x, (
                y - 1) % sz_y), (x, (y + 1) % sz_y)]
        else:
            virtual_neighbourhood = [((
                x + 1) % sz_x, y), ((x - 1) % sz_x, y)]
    else:
        if synd_type == 'X':
            virtual_neighbourhood = [((
                x + 1) % sz_x, y), ((x - 1) % sz_x, y)]
        else:
            virtual_neighbourhood = [(x, (
                y - 1) % sz_y), (x, (y + 1) % sz_y)]

    return _squoct_affine_map(virtual_neighbourhood)


def straight_octagon_path(c_1, c_2, sz):
    """
    This gives the 1-D list of coordinates on which to place an error
    in order to traverse the space between two octagonal checks of the
    same type.
    """
    # July 22, 2014: I'm on a sorting kick today.
    c_1, c_2 = sorted([c_1, c_2])

    frwrd = sorted(range(c_1 + 2, c_2, 6) + range(c_1 + 4, c_2, 6))

    # first qubit in row/col might not be @ 0:
    start = 0 if oct2sq(c_1) % 2 else 3
    rvrs = sorted(range(start, c_1, 6) + range(start + 2, c_1, 6) +
                  range(c_2 + 2, sz, 6) + range(c_2 + 4, sz, 6))

    path_list = [frwrd, rvrs]
    path_list.sort(key=len)
    return path_list[0]


def octagon_octagon_path(oct_1, oct_2, total_size):
    """
    Finds a path between octagons which are not necessarily on a
    gridline.
    """
    # Find the corner octagon
    corner_oct = oct_2[0], oct_1[1]
    # Paths in 1D
    x_path = straight_octagon_path(oct_1[
                                   0], corner_oct[0], total_size[0])
    y_path = straight_octagon_path(
        corner_oct[1], oct_2[1], total_size[1])
    # Promote paths by appending corner co-ords to all points
    x_path = [(num, oct_1[1] + 1) for num in x_path]
    y_path = [(oct_2[0] + 1, num) for num in y_path]

    return x_path + y_path


def sq_oct_shift(sq_c, oct_c, total_size):
    """
    Takes the centers of a square and an adjacent octagon as arguments.

    Returns a single point on a square which shifts a syndrome from the
    octagon to the square. Will always return the 'north-east' or
    'south-west' point, since they are equivalent to their counterparts
    up to gauge operators (lazy).
    """
    # Collect sign to add to every coordinate in sq_c
    for idx in range(len(sq_c)):

        if sq_c[idx] != oct_c[idx]:
            # print abs(sq_c[idx] - oct_c[idx])
            if abs(sq_c[idx] - oct_c[idx]) == 3:
                delta = cmp(oct_c[idx] - sq_c[idx], 0)
            elif abs(sq_c[idx] - oct_c[idx]) == total_size[idx] - 3:
                delta = -cmp(oct_c[idx] - sq_c[idx], 0)
            else:
                raise ValueError("Input square and octagon to " +
                                 "sq_oct_shift are not adjacent.")

    return tuple(sq_c[j] + delta for j in range(len(sq_c)))


def square_octagon_path(sq_1, oct_2, synd_type, total_size):
    """
    Finds the appropriate nearest-neighbour path using the
    appropriate_neighbours function, returning the path instead of the
    distance.
    """
    square_neighbours = appropriate_neighbours(
        sq_1, synd_type, total_size)

    # Determine which neighbour to use by whichever is closer to oct_2
    oct_1 = min(square_neighbours,
                key=lambda o_1: octagon_octagon_dist(o_1, oct_2, total_size))

    intersection = sq_oct_shift(sq_1, oct_1, total_size)

    return octagon_octagon_path(oct_1, oct_2, total_size) + [intersection]


def square_square_path(sq_1, sq_2, synd_type, total_size):
    """
    Same as the square_square_dist above, but returns the path.
    """
    # Non-descriptive variable names for the octagons we're going to use
    neighb_1 = appropriate_neighbours(sq_1, synd_type, total_size)
    neighb_2 = appropriate_neighbours(sq_2, synd_type, total_size)

    oct_1, oct_2 = min(product(neighb_1, neighb_2), key=lambda o_tpl:
                       octagon_octagon_dist(o_tpl[0], o_tpl[1], total_size))

    intersection_1 = sq_oct_shift(sq_1, oct_1, total_size)
    intersection_2 = sq_oct_shift(sq_2, oct_2, total_size)

    # short circuit to avoid length 2 paths which hit the same point
    # twice:
    if intersection_1 == intersection_2:
        return []

    oct_path = octagon_octagon_path(oct_1, oct_2, total_size)

    return [intersection_1] + oct_path + [intersection_2]


def hoelzer_dist(x1, x2, y1, y2, sz_x, sz_y):
    """
    Here, I test the square-octagon distance metric given by Tobias
    Hoelzer in his Bachelor's Thesis. It is meant to be numerically
    identical to the current UnionJackLattice.dist.
    """
    # First, we map the co-ordinates back to the virtual lattice:
    x1, x2, y1, y2 = map(oct2sq, [x1, x2, y1, y2])
    sz_x, sz_y = sz_x / 3, sz_y / 3

    dist = min([abs(x1 - x2), sz_x - abs(x1 - x2)]) + \
        min([abs(y1 - y2), sz_y - abs(y1 - y2)])

    if x1 != x2:
        dist += 2 * ((x1 + 1) % 2)
    if y1 != y2:
        dist += 2 * ((y1 + 1) % 2)

    return dist

def nwes_pairs(lattice, dual_lattice, label, odd_even=None):
    """
    Lists the pairs corresponding to compass directions on a primal 
    and dual square lattice, corresponding to gate/twirl locations for
    the toric code. 
    """

    if odd_even:
        if odd_even not in ['odd', 'even']:
            raise ValueError("odd_even must be odd, even or None."
                "{} entered.".format(odd_even))
        elif odd_even == 'odd':
            dual_coord_set = _odd_odds(*dual_lattice.size)
        elif odd_even == 'even':
            dual_coord_set = _even_evens(*dual_lattice.size)
    else:
        dual_coord_set = sym_coords(*dual_lattice.size)

    if label not in 'nwes':
        raise ValueError('label must be n, w, e, or s.'
            ' {} entered.'.format(label))
    
    shift = {'n': [0, 1], 'w': [-1, 0],
             'e': [1, 0], 's': [0, -1]}[label]
    
    output = []
    for crds in dual_coord_set:
        nb_crds = [(crds[k] + shift[k]) % (2 * lattice.size[k])
                    for k in range(2)]
        #output order: [primal, dual]
        output.append([lattice[nb_crds], dual_lattice[crds]])

    return output