__all__ = ['ErrorCorrectingCode', 'StabilizerCode'] 

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
    def __init__(self, arg):
        super(ErrorCorrectingCode, self).__init__()
        self.arg = arg

class StabilizerCode(ErrorCorrectingCode):
    """ 
    subclass of
    ErrorCorrectingCode for which syndromes are determined by commutation
    /anti-commutation with a Pauli stabilizer.
    """
    def __init__(self, arg):
        super(StabilizerCode, self).__init__()
        self.arg = arg
        
