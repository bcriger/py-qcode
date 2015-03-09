from qecc import Clifford, Pauli, com
#circular dep below?
from lattice import Lattice

class Clifford(object):
    """
    a `py_qcode.Clifford` is a wrapper for a `qecc.Clifford` of one of
    the three fundamental types: CNot, Hadamard and Phase. This 
    Clifford is defined over many locations on a lattice, and acts 
    simultaneously on all of them.  
    """
    def __init__(self, gate, coord_sets):
        #sanity checks
        
        if not(isinstance(gate, Clifford)):
            raise ValueError("Input 'gate' must be a "+\
                "`qecc.Clifford`, you entered: {}".format(gate))
        
        if not(gate.is_valid()):
            raise ValueError("Input 'gate' must be a valid "+\
                "Clifford, input performs the transformation:\n"+\
                "{}".format(gate.str_sparse()))
        
        _check_lengths(coord_sets, gate.nq)
        
        #end sanity checks
        self.gate = gate
        self.coord_sets = coord_sets

    def apply(self, lat=None, dual_lat=None):
        """
        Transforms a Pauli error on a lattice or pair of lattices. 
        This is made possible by the fact that a `py_qcode.Point` can 
        store an error, even if it is on a dual lattice. 
        """
        #sanity checks
        if !(lat or dual_lat):
            raise ValueError("At least one of `lat` and "+\
                "`dual_lat` must be set.")

        _check_lat(lat, dual=False)
        _check_lat(dual_lat, dual=True)
        
        if !(lat):
            for crd in coord_sets:
                dual_lat[crd].error = self.gate(dual_lat[crd].error)
        elif !(dual_lat):
            for crd in coord_sets:
                lat[crd].error = self.gate(lat[crd].error)
        else:
            for crd in coord_sets:
                lat[crd].error, dual_lat[crd].error = \
                    self.gate(lat[crd].error & dual_lat[crd].error)

class Measurement():
    """
    Pauli measurements to be done on `py_qcode.Point`s. Each object 
    consists of a Pauli type, which it will measure, and a set of 
    coordinates on which the measurement will take place.  
    """
    def __init__(self, pauli, coord_set):
        
        #sanity chex
        if !isinstance(pauli, Pauli):
            raise ValueError("input pauli must be a `qecc.Pauli`, "+\
                "{} entered.".format(pauli))
        
        self.pauli = pauli
        self.coord_set = coord_set
    
    def apply(self, dual_lat):
        _check_lat(dual_lat, dual=True)
        
        for x in self.coord_set:
            dual_lat[x].syndrome = com(self.pauli, dual_lat[x].error)

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
        "entered: {}".format(err_str, lat)
    
    if lat and not(isinstance(lat, Lattice)):
        raise ValueError(err_str)
    
    if lat.is_dual != dual:
        input_name = '`dual_lat`' if dual else '`lat`'
        lat_status = 'dual' if lat.is_dual else 'primal'
        raise ValueError("The `is_dual` of the input lat must "+\
            "match its assignment to 'lat' or 'dual_lat'. A "+\
            "{} lattice was input to {}.".format(lat_status, 
                input_name))
