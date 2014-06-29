import py_qcode as pq
import networkx as nx
from qecc import X, Z, I

test_lattice = pq.SquareLattice((4,4))
test_dual_lattice = pq.SquareLattice((4,4), is_dual=True)
test_model = pq.depolarizing_model(0.15)
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

test_decoder.infer()

test_code.measure()

print "Syndromes After Decoding: \n=========="
for point in test_dual_lattice.points:
    if point.syndrome:
        print point