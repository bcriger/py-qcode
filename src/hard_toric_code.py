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
        self.logical_error = []
    
    def run(self):
        sz = self.size

        lat = pq.SquareLattice((sz, sz))
        d_lat = pq.SquareLattice((sz, sz), is_dual=True)

        d_lat_lst = [pq.SquareLattice((sz, sz), is_dual=True)
                        for _ in range(sz)]

        decoder = pq.ft_mwpm_decoder(lat, d_lat_lst)

        log_ops = pq.toric_log_ops((sz, sz))
        
        x_flip = pq.PauliErrorModel({q.I : 1. - p, q.X : p})
        z_flip = pq.PauliErrorModel({q.I : 1. - p, q.Z : p})
        twirl = pq.two_bit_twirl(p)

        odd_prs = {drctn : pq.nwes_pairs(lat, d_lat, drctn, 'odd') 
                    for drctn in 'nwes'}
        even_prs = {drctn : pq.nwes_pairs(lat, d_lat, drctn, 'odd') 
                    for drctn in 'nwes'}
        cnot = {drctn : pq.Clifford(q.cnot(2,0,1), odd_prs[drctn])
                for drctn in 'nwes'}
        notc = {drctn : pq.Clifford(q.cnot(2,1,0), even_prs[drctn])
                for drctn in 'nwes'}
        
        x_meas = pq.Measurement(q.X, ['I', 'Z'], d_lat.star_centers())
        z_meas = pq.Measurement(q.X, ['I', 'X'], d_lat.plaq_centers())

        for _ in self.n_trials:
            #clear last sim
            pq.error_fill(lat, q.I)
            pq.error_fill(d_lat, q.I)
            pq.syndrome_fill(d_lat, 0)
            for ltc in d_lat_lst:
                pq.syndrome_fill(ltc, 0) #may break
            
            #fill d_lat_lst with syndromes by copying
            for idx in range(sz):
                meas_cycle(lat, d_lat, x_flip, z_flip, twirl, odd_prs,
                             even_prs, cnot, notc, x_meas, z_meas)
                pq.syndrome_copy(d_lat, d_lat_lst[idx])

            #run decoder, with no final lattice check (laaaaater)
            decoder.infer()

            com_relation_list = []
            for operator in log_ops:
                com_relation_list.append(operator.test(lat))
            self.logical_error.append(com_relation_list)

        pass

    def save(self, filename):
        big_dict = {}
        big_dict['lattice_class'] = 'SquareLattice'
        big_dict['lattice_size'] = self.size
        big_dict['dual_lattice_class'] = 'Dual SquareLattice'
        big_dict['dual_lattice_size'] = self.size
        big_dict['error_model'] = 'custom hard-coded'
        big_dict['code'] = 'Toric Code'
        big_dict['decoder'] = 'FT MWPM'
        big_dict['n_trials'] = self.n_trials
        big_dict['logical_errors'] = self.logical_error

        with open(filename, 'w') as phil:
            pkl.dump(big_dict, phil)

def meas_cycle(lat, d_lat, x_flip, z_flip, twirl, odd_prs, even_prs, 
                cnot, notc, x_meas, z_meas):
    
    """
    Does one cycle of measurement for the toric code on a SquareLattice
    + Flip syndrome qubits with appropriate error type
    + Apply north-facing CNots
    + Two-qubit twirl on north links
    + Apply west-facing CNots
    + Two-qubit twirl on west links
    + Apply east-facing CNots
    + Two-qubit twirl on east links
    + Apply south-facing CNots
    + Two-qubit twirl on south links
    + Flip syndrome qubits with appropriate error type
    + Measure syndrome qubits.
    """

    x_flip.act_on(d_lat.plaq_centers())
    z_flip.act_on(d_lat.star_centers())

    for drctn in 'nwes':
        cnot[drctn].apply()
        notc[drctn].apply()
        twirl.act_on(odd_prs[drctn])
        twirl.act_on(even_prs[drctn])
    
    x_flip.act_on(d_lat.plaq_centers())
    z_flip.act_on(d_lat.star_centers())

    x_meas.apply()
    z_meas.apply()

