import networkx as nx
from qecc import X, Z
import pdb
from scipy import weave
from os import getcwd
from os.path import abspath
from numpy import zeros, int16
from py_qcode import __path__ as install_path
install_path = abspath(install_path[0]) # str instead of str list.

__all__ = ['Decoder', 'mwpm_decoder', 'ft_mwpm_decoder']

#"""
__all__.extend(['matching_alg'])
#"""


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

def mwpm_decoder(primal_lattice, dual_lattice, blossom=True):
    """
    Decoder based on minimum-weight perfect matching using the blossom algorithm,
    implemented in networkx.  
    """
    if blossom:
        return Decoder(blossom_matching_alg, primal_lattice,
                        dual_lattice, name='Minimum-Weight Matching')
    else:
        return Decoder(matching_alg, primal_lattice, dual_lattice, 
                                        name='Minimum-Weight Matching')

def ft_mwpm_decoder(primal_lattice, dual_lattice_list):
    """
    Fault-tolerant decoder based on minimum-weight perfect matching 
    using the blossom algorithm, implemented in networkx. This decoder 
    follows the scheme in Dennis/Kitaev/Landahl/Preskill.
    Key Points:
    -----------
     + We treat X and Z syndromes as being completely independent.
     + We produce a :math:`d+1`-dimensional lattice as an intermediate
       object, recording differences between the syndromes as vertices 
       on the graph.
    """
    return Decoder(hi_d_matching_alg, primal_lattice, dual_lattice_list, 
                    name = '(n+1)-dimensional Minimum-Weight Matching')

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

    for g, synd_type in zip([x_graph, z_graph], ['X','Z']):
        for node in g.nodes():
            other_nodes = g.nodes()
            other_nodes.remove(node)
            for other_node in other_nodes:
                #Negative weights are no good for networkx
                edge_tuple = (node, other_node,
                    size_constant - dual_lattice.dist(node, other_node, synd_type))
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
        #pdb.set_trace()
        for pair in tpl_lst:
            coord_set = primal_lattice.min_distance_path(*pair, synd_type=str(pauli.op))
            for coord in coord_set:
                #print coord
                try:
                    primal_lattice[coord].error *= pauli
                except KeyError: 
                    print pair
                    print coord_set

    pass #This function is secretly a subroutine

def blossom_matching_alg(primal_lattice, dual_lattice):
    """
    This algorithm uses the C++ library 'Blossom V' to perform the 
    perfect matching algorithm, hopefully saving a bit of time on the 
    most expensive part of decoding.
    """
    x_verts = []; z_verts = []
    for point in dual_lattice.points:
        if point.syndrome: #exists
            if any([ltr in point.syndrome for ltr in 'xX']):
                x_verts.append(point.coords)
            if any([ltr in point.syndrome for ltr in 'zZ']):
                z_verts.append(point.coords)

    num_x_verts, num_z_verts = len(x_verts), len(z_verts) 
    num_x_edges = num_x_verts * (num_x_verts - 1) / 2
    num_z_edges = num_z_verts * (num_z_verts - 1) / 2

    x_edges = zeros((num_x_edges, 3), dtype = int16)
    z_edges = zeros((num_z_edges, 3), dtype = int16)
    x_partners = zeros((num_x_verts,), dtype = int16)
    z_partners = zeros((num_z_verts,), dtype = int16)
    
    dist = dual_lattice.dist
    for verts, edges, synd_type in zip([x_verts, z_verts],
                                         [x_edges, z_edges], 'XZ'):
        edge_count = 0
        for vert_idx, vert in enumerate(verts):
            for other_idx, o_vert in enumerate(verts[vert_idx + 1:]):
                edges[edge_count,:] = vert_idx,\
                                        other_idx + vert_idx + 1, \
                                        dist(vert, o_vert, synd_type)
                edge_count += 1

    #Bring in code
    c_code = '''
    int edge_idx, vert_idx;
    int return_val[num_verts];
    PerfectMatching *pm = new PerfectMatching(num_verts, num_edges);

    struct PerfectMatching::Options options;
    options.verbose = false; //suppress printing from c++

    pm->options = options;

    for ( edge_idx = 0; edge_idx < num_edges; edge_idx++ )
    {
        pm->AddEdge(edges(edge_idx,0), edges(edge_idx,1), edges(edge_idx,2));
    }
    pm->Solve();
    for (vert_idx = 0; vert_idx < num_verts; ++vert_idx)
        {
            int partner = pm->GetMatch(vert_idx);
            partners(vert_idx) = partner;
        }
    delete pm;
    '''
    
    #Auxiliary arguments to scipy.weave.inline:
    arg_names = ['num_verts', 'num_edges', 'edges', 'partners']
    headers = ['<PerfectMatching.h>']
    libraries = ["rt"]
    #print install_path
    include_dirs = [install_path + '/blossom5-v2.04.src/']
    extra_objects = [include_dirs[0] + 'blossom.o']

    #The heavy lifting:
    """
    print 'x_edges\n======='
    print x_edges
    print 'z_edges\n======='
    print z_edges
    """
    for num_verts, num_edges, edges, partners in \
    zip([num_x_verts, num_z_verts],[num_x_edges, num_z_edges],
        [x_edges, z_edges], [x_partners, z_partners]):
        
        weave.inline(c_code, arg_names = arg_names, 
            headers = headers, include_dirs = include_dirs, 
            type_converters = weave.converters.blitz, 
            extra_objects = extra_objects, 
            compiler='gcc', libraries=libraries)

    #Post-process 1D partner lists to avoid duplicate paths:

    x_mate_tuples, z_mate_tuples = [], []
    for verts, partners, mate_tuples in zip([x_verts, z_verts],
                                        [x_partners, z_partners],
                                        [x_mate_tuples, z_mate_tuples]):
        partnered_verts = []
        for vert, partner in enumerate(partners):
            if vert in partnered_verts:
                continue
            else:
                mate_tuples.append((verts[vert], verts[partners[vert]]))
                partnered_verts.append(partners[vert])

    #Produce error chains according to min-length path between
    #mated points
    for pauli, tpl_lst in zip([X,Z],[x_mate_tuples, z_mate_tuples]):
        #pdb.set_trace()
        for pair in tpl_lst:
            coord_set = primal_lattice.min_distance_path(*pair, synd_type=str(pauli.op))
            for coord in coord_set:
                #print coord
                try:
                    primal_lattice[coord].error *= pauli
                except KeyError: 
                    print pair
                    print coord_set

    pass #This function is secretly a subroutine

def hi_d_matching_alg(primal_lattice, dual_lattice_list):
    #First, construct a pair of graphs given syndrome data:
    x_graph = nx.Graph(); z_graph = nx.Graph()
    #For all points on each dual_lattice, compare with the point at
    #the previous time step, and add a node to the appropriate 
    #graph if they differ:
    for point in dual_lattice_list[0].points:
        if point.syndrome: #exists
            if any([ltr in point.syndrome for ltr in 'xX']):
                x_graph.add_node(point.coords + (0, ))
            if any([ltr in point.syndrome for ltr in 'zZ']):
                z_graph.add_node(point.coords + (0, ))

    for idx in range(1, len(dual_lattice_list)):
        curr_lattice = dual_lattice_list[idx]
        prev_lattice = dual_lattice_list[idx - 1]
        for point in curr_lattice.points:
            crds = point.coords
            #print "co-ordinates: {0}".format(str(crds))
            #print "syndrome comparison: {0} vs. {1}".format(point.syndrome, prev_lattice[crds].syndrome)
            if point.syndrome != prev_lattice[crds].syndrome:
                if any([ltr in point.syndrome for ltr in 'xX']):
                    x_graph.add_node(crds + (idx, ))
                if any([ltr in point.syndrome for ltr in 'zZ']):
                    z_graph.add_node(crds + (idx, ))                    
    print 'x_graph = ' + str(x_graph.edges())    
    #set an additive constant large enough for all weights to be 
    #positive:
    size_constant = 2 * len(primal_lattice.size) * \
                max(primal_lattice.size) + len(dual_lattice_list)

    for g, synd_type in zip([x_graph, z_graph], ['X','Z']):
        for node in g.nodes():
            other_nodes = g.nodes()
            other_nodes.remove(node)
            for other_node in other_nodes:
                #Negative weights are no good for networkx
                edge_tuple = (node, other_node,
                    size_constant - dual_lattice_list[0].dist(node, other_node, synd_type)
                    + abs(node[-1]-other_node[-1]))
                g.add_weighted_edges_from([edge_tuple])

    x_mate_dict, z_mate_dict = \
    map(nx.max_weight_matching, (x_graph, z_graph))
    x_mate_temps = x_mate_dict.items()
    z_mate_temps = z_mate_dict.items()

    #NetworkX assumes directional graph, includes reversed edges.
    #This will "correct errors twice", leaving stray errors on the
    #lattice.
    for tpl_lst in [x_mate_temps, z_mate_temps]:
            for tpl in tpl_lst:
                rvrs = tuple(reversed(tpl))
                if rvrs in tpl_lst:
                    tpl_lst.remove(rvrs)
    
    x_mate_tuples = []
    z_mate_tuples = []

    #Eliminate vertical paths
    for lst_in, lst_out in zip([x_mate_temps, z_mate_temps],
                                [x_mate_tuples, z_mate_tuples]):
        for mate_tuple in lst_in:
            if mate_tuple[0][:-1] != mate_tuple[1][:-1]:
                lst_out.append(mate_tuple)

    #Project remaining paths onto n-dimensions:
    for lst in ([x_mate_tuples, z_mate_tuples]):
        for item in lst:
            item[0], item[1] = item[0][:-1], item[1][:-1]

    #Produce error chains according to min-length path between
    #mated points
    for pauli, tpl_lst in zip([X,Z],[x_mate_tuples, z_mate_tuples]):
        #pdb.set_trace()
        for pair in tpl_lst:
            coord_set = primal_lattice.min_distance_path(*pair, synd_type=str(pauli.op))
            for coord in coord_set:
                #print coord
                try:
                    primal_lattice[coord].error *= pauli
                except KeyError: 
                    print pair
                    print coord_set

    pass #Subroutine