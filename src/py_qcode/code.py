from qecc import Pauli

__all__ = ['ErrorCorrectingCode', 'StabilizerCode', 'ErrorCheck', 'StabilizerCheck']

class ErrorCheck():
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
        self.primal_set = primal_set
        self.dual_point = dual_point
        self.rule = rule

    def evaluate(self):
        pass
class StabilizerCheck(ErrorCheck):
    """
    subclass of :class:`py_qcode.ErrorCheck`, takes anything that can be cast to a :class:`qecc.Pauli` instead of a rule, and uses commutation to determine the syndrome. 
    """
    def __init__(self, arg):
        super(StabilizerCheck).__init__()
        self.arg = arg


class ErrorCorrectingCode():
    """
    An error-correcting code, for the purpose of this module, is a rule
    for taking errors on sets of points, and turning them into
    discrete-valued syndromes which are interpreted by the decoder.
    Normally, we would immediately make the restriction to stabilizer
    codes which return binary-valued syndromes, but we want to make
    room for codes which correct non-Pauli errors, and return fuzzy
    syndromes.

    :param primal_lattice: The lattice on which the qubits live.

    :type primal_lattice: :class:`py_qcode.Lattice` 
    
    :param dual_lattice: The lattice on which the parity checks live.

    :type dual_lattice: :class:`py_qcode.Lattice`

    :param parity_check: A rule for mapping errors on the primal lattice to measurement results on the dual lattice.

    :type parity_check: function 
    """
    def __init__(self, primal_lattice, dual_lattice, parity_check):
        
        self.primal_lattice = primal_lattice
        self.dual_lattice = dual_lattice
        self.parity_check = parity_check
        self.parity_check_list

class StabilizerCode(ErrorCorrectingCode):
    """ 
    subclass of
    ErrorCorrectingCode for which syndromes are determined by commutation
    /anti-commutation with a Pauli stabilizer.
    """
    def __init__(self, arg):
        super(StabilizerCode, self).__init__()
        self.arg = arg
        
#UTILITY FUNCTIONS
def commutation_rule(point_set, pauli):
    """
    This function determines whether an error on the set of points `point_set` commutes with an input Pauli. It's used to produce anonymous functions for input to :class:`py_qcode.ErrorCheck`'s `.__init__` method when initializing a :class:`py_qcode.StabilizerCheck` instance.
    """
    pass