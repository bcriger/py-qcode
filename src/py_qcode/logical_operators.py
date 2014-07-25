from qecc import Pauli, com
from lattice import _evens, _odds

__all__ = ['LogicalOperator', 'toric_log_ops', 'squoct_log_ops']

class LogicalOperator():
    """
    This class wraps a function which tests anticommutation on a 
    lattice.
    """
    def __init__(self, name, pauli, coord_list):
        self.name = name
        self.pauli = pauli
        self.coord_list = coord_list
    
    def test(self, lattice):
        test_pauli = reduce(lambda a, b: a.tens(b), 
                        [lattice[coord].error for coord in self.coord_list])
        return self.name, com(test_pauli, self.pauli)

def toric_log_ops(sz_tpl):
    """
    This is really poor design, but I'm kind of tired of this code. 
    Maybe that's causative. Anyway, this returns a list of 
    LogicalOperator objects that matches the toric code defined on a 
    rectangular lattice.
    """
    if len(sz_tpl) != 2:
        raise ValueError('Only 2D codes are supported in this'+\
         'function for now, maybe define your operators manually')
    x_1 = LogicalOperator("X_1", Pauli('X' * sz_tpl[0]),
                            [(x, 1) for x in _evens(sz_tpl[0])])
    x_2 = LogicalOperator("X_2", Pauli('X' * sz_tpl[1]),
                            [(1, x) for x in _evens(sz_tpl[1])])
    z_1 = LogicalOperator("Z_1", Pauli('Z' * sz_tpl[1]),
                            [(0, x) for x in _odds(sz_tpl[1])])
    z_2 = LogicalOperator("Z_2", Pauli('Z' * sz_tpl[0]),
                            [(x, 0) for x in _odds(sz_tpl[0])])
    
    return [x_1, x_2, z_1, z_2]

def squoct_log_ops(total_size):
    """
    Returns a list of LogicalOperator objects matching the non-gauge
    qubits from a concatenated [[4,2,2]]/toric code.
    """
    if len(total_size) != 2:
        raise ValueError('Only 2D codes are supported in this'+\
         'function for now, maybe define your operators manually')
    
    edge_0 = range(0, total_size[0], 6) + range(2, total_size[0], 6)
    edge_1 = range(0, total_size[1], 6) + range(2, total_size[1], 6)
    shift_edge_0 = range(3, total_size[0], 6) + range(5, total_size[0], 6)
    shift_edge_1 = range(3, total_size[1], 6) + range(5, total_size[1], 6)
    
    x_1 = LogicalOperator("X_1", Pauli('X' * (total_size[1] / 3)),
                            [(0, x) for x in shift_edge_1])
    x_2 = LogicalOperator("X_2", Pauli('X' * (total_size[0] / 3)),
                            [(x, 3) for x in edge_0])
    z_1 = LogicalOperator("Z_1", Pauli('Z' * (total_size[0] / 3)),
                            [(x, 3) for x in edge_0])
    z_2 = LogicalOperator("Z_2", Pauli('Z' * (total_size[1] / 3)),
                            [(0, x) for x in shift_edge_1])
    
    return [x_1, x_2, z_1, z_2]