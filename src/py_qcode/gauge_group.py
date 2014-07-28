from qecc import Pauli, PauliList

__all__ = ['GaugeGroup']

class GaugeGroup(object):
    """
    This class serves two purposes. Firstly, it enables us to reduce
    the weight of a post-inference error 
    """
    def __init__(self, lattice, paulis):
        self.lattice = lattice
        self.paulis = paulis
    
    def apply(self):