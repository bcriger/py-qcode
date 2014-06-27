import py_qcode as pq
import networkx as nx
from qecc import X, Z, I

test_lattice = pq.SquareLattice((4,4))
test_dual_lattice = pq.SquareLattice((4,4), is_dual=True)
test_model = pq.depolarizing_model(0.05)
test_code = pq.toric_code(test_lattice, test_dual_lattice)
test_decoder = pq.mwpm_decoder(test_lattice, test_dual_lattice)

test_model.act_on(test_lattice)

print "Real Errors: \n============"
for point in test_lattice.points:
    if not(point.error == I):
        print point

test_code.measure()

print "Syndromes: \n=========="
for point in test_dual_lattice.points:
    if point.syndrome:
        print point

x_graph = nx.Graph(); z_graph = nx.Graph()
#For all points on the dual lattice, add a vertex to the
#appropriate graph and weighted edges connecting it to 
#every prior vertex.

for point in test_dual_lattice.points:
    if point.syndrome: #exists
        if any([ltr in point.syndrome for ltr in 'xX']):
            x_graph.add_node(point.coords)
        if any([ltr in point.syndrome for ltr in 'zZ']):
            z_graph.add_node(point.coords)

print "Syndrome Graphs:\n================"
print "X Graph Nodes:\n--------------"
print x_graph.nodes()

print "Z Graph Nodes:\n--------------"
print z_graph.nodes()

size_constant = 2 * len(test_lattice.size) * max(test_lattice.size)
print "\nsize_constant = " + str(size_constant) + '\n'

for g in [x_graph, z_graph]:
    for node in g.nodes():
        other_nodes = g.nodes()
        other_nodes.remove(node)
        for other_node in other_nodes:
            #Negative weights are no good for networkx
            edge_tuple = (node, other_node,
                size_constant - test_dual_lattice.dist(node,other_node))
            g.add_weighted_edges_from([edge_tuple])

print "X Graph:\n--------"
for item in x_graph.adjacency_iter():
    print item

print "Z Graph:\n--------"
for item in z_graph.adjacency_iter():
    print item

x_mate_dict, z_mate_dict = \
map(nx.max_weight_matching, (x_graph, z_graph))

print "Matched Pairs:\n=============="

x_mate_tuples = x_mate_dict.items()
z_mate_tuples = z_mate_dict.items()

for tpl_lst in [x_mate_tuples, z_mate_tuples]:
    for tpl in tpl_lst:
        if tuple(reversed(tpl)) in tpl_lst:
            print 'removing'
            tpl_lst.remove(tpl)

print "X Errors:\n========="
print x_mate_tuples
print "Z Errors:\n========="
print z_mate_tuples

#Produce error chains according to min-length path between
#mated points
print "X Error Paths:\n=============="
for pair in x_mate_tuples:    
    coord_set = test_lattice.min_distance_path(*pair)
    for coord in coord_set:
        print coord
        test_lattice[coord].error *= X   

print "Z Error Paths:\n=============="
for pair in z_mate_tuples:
    coord_set = test_lattice.min_distance_path(*pair)
    for coord in coord_set:
        print coord
        test_lattice[coord].error *= Z
