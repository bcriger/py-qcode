import networkx as nx

__all__ = ['Decoder', 'mwpm_decoder', 'RGBPDecoder', 'BHRGDecoder']

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

def mwpm_decoder(primal_lattice, dual_lattice):
    """
    Decoder based on minimum-weight perfect matching using the blossom algorithm,
    implemented in networkx.  
    """
    def matching_alg(primal_lattice, dual_lattice):
        """
        There are two steps to this algorithm. First, we solve the
        matching problem on the dual lattice, identifying pairs of
        points with minimum-weight error chains between them. Then,
        we use a simple rule to produce an error chain from the
        matching.

        This decoder is only implemented with toric codes in mind.
        It treats X and Z errors as completely independent.
        """
        
        #First, construct a pair of graphs given syndrome data:
        x_graph = nx.Graph(); z_graph = nx.Graph()
        #For all points on the dual lattice, add a vertex to the
        #appropriate graph and weighted edges connecting it to 
        #every prior vertex.
        
        for point in dual_lattice.points:
            if any([ltr in point.syndrome for ltr in 'xX']):
                x_graph.add_node(point.coords)
            if any([ltr in point.syndrome for ltr in 'zZ']):
                z_graph.add_node(point.coords)
        
        for g in [x_graph, z_graph]:
            for node in g.nodes():
                other_nodes = g.nodes()
                other_nodes.remove(node)
                for other_node in other_nodes:
                    edge_tuple = (node, other_node,
                                    -dual_lattice.dist(node,other_node))
                    g.add_weighted_edges_from([edge_tuple])

        x_mate_dict, z_mate_dict = \
        map(nx.max_weighted_matching, (x_graph, z_graph))

        #Produce error chains according to min-length path between
        #mated points
        for pair in x_mate_dict.items():
            coord_set = min_length_path(*pair, dual_lattice.dist)
            pass #what happens if error chains cross? check if threshold is low.
        pass #This function is secretly a subroutine

    return Decoder(matching_alg, primal_lattice, dual_lattice)

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

#Convenience Fuctions

def min_length_path(start, end, sz_tpl, closed_boundary=True):
    """
    Finds the shortest path between two points on a union-jack lattice
    ('squares' inscribed with 45-degree-rotated 'squares'), given the
    start point, the end point and a flag to indicate whether the 
    underlying lattice is closed. 

    The basic idea is to use the union of disjoint straight paths 
    between 'corner' sites on the dual lattice. This minimizes
    total path length while minimizing the number of diagonal moves, 
    and (hopefully) the complexity of the code.
    """
    #Determine corner sites:
    

    pass

