import networkx as nx
from qecc import X, Z

__all__ = ['Decoder', 'mwpm_decoder', 'RGBPDecoder', 'BHRGDecoder']

class Decoder():
    """
    The role of a decoder is to infer an error given a syndrome. This 
    requires the presence of the dual lattice (where the syndromes are
    stored), the primal lattice (where the real error is stored)
    """
    def __init__(self, algorithm, primal_lattice, dual_lattice, name='Un-named'):
        self.algorithm = algorithm
        self.primal_lattice = primal_lattice
        self.dual_lattice = dual_lattice
        self.name = name

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
            if point.syndrome: #exists
                if any([ltr in point.syndrome for ltr in 'xX']):
                    x_graph.add_node(point.coords)
                if any([ltr in point.syndrome for ltr in 'zZ']):
                    z_graph.add_node(point.coords)

        #set an additive constant large enough for all weights to be positive:
        size_constant = 2 * len(primal_lattice.size) * max(primal_lattice.size)

        for g in [x_graph, z_graph]:
            for node in g.nodes():
                other_nodes = g.nodes()
                other_nodes.remove(node)
                for other_node in other_nodes:
                    #Negative weights are no good for networkx
                    edge_tuple = (node, other_node,
                        size_constant - dual_lattice.dist(node,other_node))
                    g.add_weighted_edges_from([edge_tuple])

        x_mate_dict, z_mate_dict = \
        map(nx.max_weight_matching, (x_graph, z_graph))
        x_mate_tuples = x_mate_dict.items()
        z_mate_tuples = z_mate_dict.items()

        #NetworkX assumes directional graph, includes reversed edges.
        #This will "correct errors twice", leaving stray errors on the
        #lattice.
        for tpl_lst in [x_mate_tuples, z_mate_tuples]:
                for tpl in tpl_lst:
                    rvrs = tuple(reversed(tpl))
                    if rvrs in tpl_lst:
                        tpl_lst.remove(rvrs)

        #Produce error chains according to min-length path between
        #mated points
        for pauli, tpl_lst in zip([X,Z],[x_mate_tuples, z_mate_tuples]):
            for pair in tpl_lst:    
                coord_set = primal_lattice.min_distance_path(*pair)
                for coord in coord_set:
                    primal_lattice[coord].error *= pauli   

        pass #This function is secretly a subroutine

    return Decoder(matching_alg, primal_lattice, dual_lattice, name='Minimum-Weight Matching')

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