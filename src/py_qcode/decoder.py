import networkx as nx

__all__ = ['Decoder', 'MWPMDecoder', 'RGBPDecoder', 'BHRGDecoder']

class Decoder():
    """
    The role of a decoder is to infer an error given a syndrome. This 
    requires the presence of the dual lattice (where the syndromes are
    stored), the primal lattice (where the real error is stored)
    """
    def __init__(self, algorithm, primal_lattice, dual_lattice):
        self.algorithm = algorithm
        self.primal_lattice = primal_lattice
        self.dual_lattice = dual_lattice

    def infer(self):
        """
        Uses `self.algorithm` to update the error on the primal_lattice,
        given the syndromes on the dual lattice.
        """
        self.algorithm(self.primal_lattice, self.dual_lattice)

class MWPMDecoder(Decoder):
    """
    Decoder based on minimum-weight perfect matching using the blossom algorithm,
    implemented in networkx.  
    """
    def __init__(self, primal_lattice, dual_lattice):
        
        def matching_alg(primal_lattice, dual_lattice):
            """
            There are two steps to this algorithm. First, we solve the
            matching problem on the dual lattice, identifying pairs of
            points with minimum-weight error chains between them. Then,
            we use a simple rule to produce an error chain from the
            matching.

            This decoder is only implemented with toric codes in mind.
            It treats X and Z errors as completely independent.

            TODO: Put some kind of flag to denote X and Z errors. 
            Change syndromes from {0,1} to a larger set. 
            """
            
            #First, construct a pair of graphs given syndrome data:
            x_graph = nx.Graph(); z_graph = nx.Graph()
            


        super(MWPMDecoder, self).__init__(matching_alg, primal_lattice, dual_lattice)

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
