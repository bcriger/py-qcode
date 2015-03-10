import py_qcode as pq, qecc as q

class HardCodeToricSim():
    """
    Sandbox for trying to simulate the toric code directly using 
    Cliffords and small error models. This should really be a subclass
    of Simulation, but I'll figure that out later.  
    """
    def __init__(self, size, p, n_trials):
        self.size = size
        self.p = p
        self.n_trials = n_trials
        self.logical_errors = []
    
    def run(self):
        primal_lattice = pq.SquareLattice((size,size))
        dual_lattice = pq.SquareLattice((size,size), is_dual=True)
        for trial in self.n_trials:
            primal_lattice.clear()
            dual_lattice.clear()
        pass

    def save(self, filename):
        pass