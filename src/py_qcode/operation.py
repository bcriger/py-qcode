from collections import Iterable
from qecc import Pauli, com
from qecc import Clifford as clfrd #collision avoidance
#circular dep below?
from lattice import Lattice
from error import depolarizing_model, two_bit_twirl

__all__ = ['Clifford', 'Measurement', 'Timestep']

class Clifford(object):
    """
    a `py_qcode.Clifford` is a wrapper for a `qecc.Clifford` of one of
    the three fundamental types: CNot, Hadamard and Phase. This 
    Clifford is defined over many locations on a lattice, and acts 
    simultaneously on all of them.  
    """
    def __init__(self, gate, point_sets):
        #sanity checks

        if not(isinstance(gate, clfrd)):
            raise ValueError("Input 'gate' must be a "+\
                "`qecc.Clifford`, you entered:\n {}".format(gate))
        
        if not(gate.is_valid()):
            raise ValueError("Input 'gate' must be a valid "+\
                "Clifford, input performs the transformation:\n"+\
                "{}".format(gate.str_sparse()))
        
        _check_lengths(point_sets, gate.nq)
        
        #end sanity checks
        self.gate = gate
        self.point_sets = point_sets
        self._support = None

    @property
    def support(self):
        if self._support == None:
            self._support = set().union(*map(set, self.point_sets))
        return self._support
    
    def apply(self, length=None):
        """
        Transforms a Pauli error on a lattice or pair of lattices. 
        This is made possible by the fact that a `py_qcode.Point` can 
        store an error, even if it is on a dual lattice.

        Note: This is field expedient, we either iterate over points or
        sets of two points, this should be refactored if three-qubit 
        Cliffords become elements of the fundamental gate set. 
        """
        if not(length):
            length = self.gate.nq

        if length == 1:
            for point in self.point_sets:
                point.error = self.gate(point.error)
        elif length == 2:
            for lst in self.point_sets:
                lst[0].error, lst[1].error = \
                    self.gate(lst[0].error & lst[1].error)
        else:
            raise ValueError("Input length is forbidden: "
                "{}".format(length))

class Measurement():
    """
    Pauli measurements to be done on `py_qcode.Point`s. Each object 
    consists of a Pauli type, which it will measure, and a set of 
    `py_qcode.Point`s on which the measurement will take place.  
    """
    #TODO: Change this to have a callback function with Pauli
    #measurement as a subclass
    def __init__(self, pauli, outputs, point_set):
        
        #sanity chex
        if not(isinstance(pauli, Pauli)):
            raise ValueError("input pauli must be a `qecc.Pauli`, "+\
                "{} entered.".format(pauli))
        
        self.pauli = pauli
        self.outputs = outputs
        self.point_set = point_set
    
    def apply(self):
        for pt in self.point_set:
            symbol = self.outputs[com(self.pauli, pt.error)]
            if pt.syndrome:
                pt.syndrome += symbol
            else:
                pt.syndrome = symbol

class Timestep():
    """
    Wraps a bunch of Cliffords for noisy application in circuit models.
    This way, I can calculate which qubits to depolarize at 
    initialization, and running the circuit is just looping over a list
    of Timesteps.
    """
    def __init__(self, cliff_list):
        #SANITIZE
        bad_list_err = ValueError("Timestep wraps an iterable of"
            " py_qcode.Clifford's, you entered: {}".format(cliff_list))
        
        if not isinstance(cliff_list, Iterable):
            raise bad_list_err
        for cliff in cliff_list:
            if not isinstance(cliff, Clifford):
                raise bad_list_err

        self.cliff_list = cliff_list
        #FIXME: Only works for two-qubit gates
        self.twirl_support = set().union(*[x.support for x in cliff_list])
    
    def noisy_apply(self, dat_lat, anc_lat, p):
        """
        Applies the gates in the time step to the data/ancilla qubits, 
        applying a depolarizing map after one-qubit gates, a two-qubit 
        twirl after two-qubit gates, and a depolarizing map after all 
        memory locations.
        """
        dep_model = depolarizing_model(p)
        twirl = two_bit_twirl(p)        

        dep_support = set(dat_lat.points + anc_lat.points)
        dep_support -= self.twirl_support

        for cliff in self.cliff_list:
            cliff.apply()
            if cliff.gate.nq == 2:
                twirl.act_on(cliff.point_sets)

        if dep_support != set():
            dep_model.act_on(list(dep_support))

        pass        

#Convenience Functions
def _check_lengths(coord_sets, length):
    for coord_set in coord_sets:
        if len(coord_set) != length:
            raise ValueError(("length of coord_set must match gate "+\
                "length. CNot has length {}, 'coord_sets' contains "+\
                "{}. ").format(length, coord_set))

def _check_lat(lat, dual):
    
    lat_type = "Dual" if dual else "Primal"
    err_str = "{} lat must be a `py_qcode.Lat`, you "+\
        "entered: {}".format(lat_type, lat)
    
    if lat and not(isinstance(lat, Lattice)):
        raise ValueError(err_str)
    
    if lat.is_dual != dual:
        input_name = '`dual_lat`' if dual else '`lat`'
        lat_status = 'dual' if lat.is_dual else 'primal'
        raise ValueError("The `is_dual` of the input lat must "+\
            "match its assignment to 'lat' or 'dual_lat'. A "+\
            "{} lattice was input to {}.".format(lat_status, 
                input_name))
