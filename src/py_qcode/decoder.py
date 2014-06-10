__all__ = ['Decoder', 'MWPMDecoder', 'RGBPDecoder', 'BHRGDecoder']

class Decoder():
    """
    The role of a decoder is to infer an error given a syndrome. The 
    overall pattern is  
    """
    def __init__(self, primal_lattice, dual_lattice, algorithm):
        self.primal_lattice = primal_lattice
        self.dual_lattice = dual_lattice
        self.algorithm = algorithm

    def infer(self):
        """
        Uses `self.algorithm` to update the error on the primal_lattice,
        given the syndromes on the dual lattice.
        """
        #Pseudocode, requires refactor on lattice.py
        self.primal_lattice.update(self.algorithm(self.dual_lattice),'error')

class MWPMDecoder(Decoder):
    """
    Decoder based on minimum-weight perfect matching using the blossom algorithm. 
    """
    def __init__(self, arg):
        self.arg = arg

class RGBPDecoder(Decoder):
    """
    A renormalization group / belief propagation decoder based on Duclos-Cianci and Poulin.
    """
    def __init__(self, arg):
        self.arg = arg

class BHRGDecoder(Decoder):
    """
    Bravyi/Haah renormalization group decoder.
    """
    def __init__(self, arg):
        self.arg = arg
