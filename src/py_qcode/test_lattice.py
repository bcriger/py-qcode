import lattice as lt
import code as cd
import utils as ut
from qecc import I, X, Z

#Global lattices to use in various tests
g_sq_lat = lt.SquareLattice((8,8))
g_d_sq_lat = lt.SquareLattice((8,8), is_dual=True)
g_squoct_lat = lt.SquareOctagonLattice((4,4))
g_uj_lat = lt.UnionJackLattice((4,4))

#Codes to be measured as part of min_distance_path tests:
toric_code = cd.toric_code(g_sq_lat, g_d_sq_lat)
squoct_code = cd.square_octagon_code(g_squoct_lat, g_uj_lat)

crds = lambda a: a.coords

def square_dist_path_pred(pt, o_pt, synd_type):
    """
    Predicate which tests equality between distance and path length on
    the square lattice.
    """
    return g_d_sq_lat.dist(pt, o_pt, synd_type) == \
                len(g_sq_lat.min_distance_path(pt, o_pt, synd_type))

def squoct_dist_path_pred(pt, o_pt, synd_type):
    """
    Predicate which tests equality between distance and path length on
    the square lattice.
    """
    return g_uj_lat.dist(pt, o_pt, synd_type) == \
            len(g_squoct_lat.min_distance_path(pt, o_pt, synd_type))

def pt_pred_loop(pt_lst, o_pt_lst, synd_type, pred_fun):
    """
    Template for the first few unit tests, asserts a predicate for a 
    pair of point lists. 
    """
    for pt in pt_lst:
        for o_pt in o_pt_lst:
            if pred_fun(pt, o_pt, synd_type):
                pass
            else:
                print pt, o_pt
                assert False
    return None

def plaquette_test():
    """
    For all pairs of plaquettes on the dual lattice of the toric code, 
    tests that the length of the minimum-length path as returned by the
    class method is the same as the distance between the points. 
    """
    return pt_pred_loop(lt._odd_odds(*g_d_sq_lat.size), 
                        lt._odd_odds(*g_d_sq_lat.size),
                        None, square_dist_path_pred)
            
def star_test():
    """
    For all pairs of stars on the dual lattice of the toric code, 
    tests that the length of the minimum-length path as returned by the
    class method is the same as the distance between the points. 
    """
    return pt_pred_loop(lt._even_evens(*g_d_sq_lat.size), 
                        lt._even_evens(*g_d_sq_lat.size),
                        None, square_dist_path_pred)

#In this section, we break square-octagon distance testing into a 
#pedantic number of cases:

def oct_oct_z_test():
    z_oct_lst = lt._squoct_affine_map(lt._odd_odds(*g_squoct_lat.size))
    return pt_pred_loop(z_oct_lst, z_oct_lst, 'X',
                                squoct_dist_path_pred)

def oct_oct_x_test():
    x_oct_lst = lt._squoct_affine_map(lt._even_evens(*g_squoct_lat.size))
    return pt_pred_loop(x_oct_lst, x_oct_lst, 'Z',
                                squoct_dist_path_pred)

def squ_oct_z_test():
    z_oct_lst = lt._squoct_affine_map(lt._odd_odds(*g_squoct_lat.size))
    sq_lst = lt._squoct_affine_map(lt.skew_coords(*g_squoct_lat.size))
    return pt_pred_loop(z_oct_lst, sq_lst, 'X',
                                squoct_dist_path_pred)

def squ_oct_x_test():
    x_oct_lst = lt._squoct_affine_map(lt._even_evens(*g_squoct_lat.size))
    sq_lst = lt._squoct_affine_map(lt.skew_coords(*g_squoct_lat.size))
    return pt_pred_loop(x_oct_lst, sq_lst, 'Z',
                                squoct_dist_path_pred)

def squ_squ_x_test():
    sq_lst = lt._squoct_affine_map(lt.skew_coords(*g_squoct_lat.size))
    return pt_pred_loop(sq_lst, sq_lst, 'X',
                                squoct_dist_path_pred)

def squ_squ_z_test():
    sq_lst = lt._squoct_affine_map(lt.skew_coords(*g_squoct_lat.size))
    return pt_pred_loop(sq_lst, sq_lst, 'X',
                                squoct_dist_path_pred)

#__getitem__ tests, we confirm that, for all points in a lattice, 
#lattice[(x, y)].coords == (x, y)

def toric_items_test():
    for pt in g_sq_lat.points:
        assert g_sq_lat[pt.coords] == pt

def dual_toric_items_test():
    for pt in g_d_sq_lat.points:
        assert g_d_sq_lat[pt.coords] == pt

def squoct_items_test():
    for pt in g_squoct_lat.points:
            assert g_squoct_lat[pt.coords] == pt

def uj_items_test():
    for pt in g_uj_lat.points:
        assert g_uj_lat[pt.coords] == pt

#Path commutation/anticommutation tests: 
#Check that min_distance_path produces Paulis which anticommute with 
#stabilizers at their endpoints and not elsewhere:

def path_com_pred(lat, d_lat, code, crd_lst, o_crd_lst, pauli):
    """
    Given a lattice, dual lattice, code, pair of coordinate lists, and
    a pauli, returns true if laying down that Pauli on paths between 
    every pair of points across the lists results in syndromes of the 
    matching type on those points, and no syndromes elsewhere.
    """
    #Put the identity down on the lattice
    for prml_pt in lat.points:
        prml_pt.error = I

    for crd in crd_lst:
        for o_crd in o_crd_lst:
            path = lat.min_distance_path(crd, o_crd, pauli.op)
            if len(path) > 0:
                for coord in path:
                    lat[coord].error *= pauli
                code.measure()
                if d_lat[crd].syndrome != pauli.op:
                    return False
                if d_lat[o_crd].syndrome != pauli.op:
                    return False
                d_lat[crd].syndrome = None
                d_lat[o_crd].syndrome = None
                for dual_pt in d_lat.points:
                    if dual_pt.syndrome:
                        return False
            for coord in path:
                lat[coord].error = I

    #If False hasn't yet been returned, return True
    return True

def toric_path_code_test():
    for pauli, crd_lst in zip([X, Z], [lt._odd_odds(*g_sq_lat.size),
                                    lt._even_evens(*g_sq_lat.size)]):
        
        assert path_com_pred(g_sq_lat, g_d_sq_lat, toric_code, 
                        crd_lst, crd_lst, pauli)
    return None
    
def squoct_path_code_test():
    x_oct_lst = lt._squoct_affine_map(lt._even_evens(*g_squoct_lat.size))
    z_oct_lst = lt._squoct_affine_map(lt._odd_odds(*g_squoct_lat.size))
    sq_lst = lt._squoct_affine_map(lt.skew_coords(*g_squoct_lat.size))
    crd_1 = [sq_lst, x_oct_lst, z_oct_lst, sq_lst, 
                    sq_lst, x_oct_lst, z_oct_lst, sq_lst]
    crd_2 = [sq_lst, sq_lst, sq_lst, sq_lst, x_oct_lst,
                    x_oct_lst, z_oct_lst, z_oct_lst]
    p_lst = [X, Z, X, Z, Z, Z, X, X]
    count = 0
    for crd_lst, o_crd_lst, pauli in zip(crd_1, crd_2, p_lst):
        print count
        count += 1
        assert path_com_pred(g_squoct_lat, g_uj_lat, squoct_code, 
                                        crd_lst, o_crd_lst, pauli)