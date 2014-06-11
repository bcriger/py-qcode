__all__ = ['Decoder', 'MWPMDecoder', 'RGBPDecoder', 'BHRGDecoder']

class Decoder():
    """
    The role of a decoder is to infer an error given a syndrome. This 
    requires the presence of the dual lattice (where the syndromes are
    stored), the primal lattice (where the real error is stored)
    """
    def __init__(self, algorithm, primal_lattice, dual_lattice):
        self.algorithm = algorithm
        self.update_rule = update_rule
        self.primal_lattice = primal_lattice
        self.dual_lattice = dual_lattice

    def infer(self):
        """
        Uses `self.algorithm` to update the error on the primal_lattice,
        given the syndromes on the dual lattice.
        """
        self.algorithm(self.update_rule, self.primal_lattice, self.dual_lattice)

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
