from qecc import Pauli, commutes_with
from lattice import _even_evens, _odd_odds 
from types import FunctionType

__all__ = ['ErrorCorrectingCode', 'ErrorCheck', 'StabilizerCheck', 'toric_code']

class ErrorCheck(object):
    """
    This is the primitive operation of measurement for error-correcting codes; it takes a list of errors on a subset of the primal lattice of the code and translates it into syndromes on a subset of the dual lattice.

    :param primal_sets: co-ordinates from which the error check will collect syndromes. I'll add input checking so that tuples of co-ordinates can be entered on their own instead of the objects which wrap them.

    :type primal_sets: collection of whatever maps to :class:`py_qcode.Point` objects

    :param dual_points: points on the dual lattice to which the syndrome for this error check will be written

    :type dual_points: set/list of tuples or :class:`py_qcode.Point` objects.

    :param rule: lookup table or other mechanism that maps errors to syndromes.

    :type rule: function or dict
    """
    def __init__(self, primal_sets, dual_points, rule):

        self.primal_sets = primal_sets
        self.dual_points = dual_points
        self.rule = rule

    def evaluate(self):
        for idx, point in enumerate(self.dual_points):
            error_str = ''.join([pt.error for pt in self.primal_sets[idx]])
            if isinstance(self.rule, dict):
                try:
                    point.syndrome = self.rule[error_str]
                except KeyError:
                    raise KeyError("There is no entry in the lookup table for the error " + error_str)
                finally:
                    pass
            elif isinstance(self.rule, FunctionType):
                point.syndrome = self.rule(error_str)
            else:
                raise TypeError("Rule used by error check must be a function or dict, you entered a value of type: " + type(self.rule))

class StabilizerCheck(ErrorCheck):
    """
    subclass of :class:`py_qcode.ErrorCheck`, takes anything that can be cast to a :class:`qecc.Pauli` instead of a rule, and uses commutation to determine the syndrome. 
    """
    def __init__(self, primal_sets, dual_points, stabilizer):
        
        if type(stabilizer) is str:
            stabilizer = Pauli(stabilizer)
        
        #returns 0 if error commutes with stabilizer, 1 if it anti-commutes
        stab_rule = lambda err_str: 1 - int(commutes_with(stabilizer)(Pauli(err_str)))
        
        super(StabilizerCheck, self).__init__(primal_sets, dual_points, stab_rule)

        self.stabilizer = stabilizer
        

class ErrorCorrectingCode():
    """
    Wraps a bunch of parity checks. 

    :param parity_check_list: A list of :class:`py_qcode.ErrorCheck` objects, which can be a mix of any subclass of :class:`py_qcode.ErrorCheck`.

    :type parity_check_list: list  
    """
    def __init__(self, parity_check_list):
        
        self.parity_check_list = parity_check_list

    def measure(self):
        """
        Evaluates all the parity checks.
        """
        for check in self.parity_check_list:
            check.evaluate()
        
#UTILITY FUNCTIONS
def toric_code(primal_grid, dual_grid):
    """
    Uses a few convenience functions to produce the toric code on a set of square lattices.
    """
    star_coords = _even_evens(*dual_grid.size)    
    star_duals = [dual_grid[coord] for coord in star_coords]
    star_primal = [primal_grid.neighbours(coord) for coord in star_coords]
    star_check = StabilizerCheck(star_primal, star_duals, 'XXXX')
    
    plaq_coords = _odd_odds(*dual_grid.size)    
    plaq_duals = [dual_grid[coord] for coord in plaq_coords]
    plaq_primal = [primal_grid.neighbours(coord) for coord in plaq_coords]
    plaq_check = StabilizerCheck(plaq_primal, plaq_duals, 'ZZZZ')

    return ErrorCorrectingCode([star_check, plaq_check])