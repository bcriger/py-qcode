from qecc import Pauli

__all__ = ['ErrorCorrectingCode', 'StabilizerCode', 'ErrorCheck', 'StabilizerCheck']

class ErrorCheck():
    """
    This is the primitive operation of measurement for error-correcting codes; it takes a list of errors on a subset of the primal lattice of the code and translates it into a syndrome on the dual lattice.

    :param primal_set: co-ordinates from which the error check will collect syndromes. I'll add input checking so that tuples of co-ordinates can be entered on their own instead of the objects which wrap them.

    :type primal_set: collection of whatever maps to :class:`py_qcode.Point` objects

    :param dual_point: point on the dual lattice to which the syndrome will be written

    :type dual_point: tuple or :class:`py_qcode.Point`

    :param rule: lookup table or other mechanism that maps errors to syndromes.

    :type rule: function
    """
    def __init__(self, primal_set, dual_point, rule):
        self.primal_set = primal_set
        self.dual_point = dual_point
        self.rule = rule

class StabilizerCheck(ErrorCheck):
    """

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

class StabilizerCode(ErrorCorrectingCode):
    """ 
    subclass of
    ErrorCorrectingCode for which syndromes are determined by commutation
    /anti-commutation with a Pauli stabilizer.
    """
    def __init__(self, arg):
        super(StabilizerCode, self).__init__()
        self.arg = arg
        
