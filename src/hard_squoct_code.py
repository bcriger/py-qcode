import py_qcode as pq, qecc as q, cPickle as pkl
from hard_toric_code import SIM_TYPES, _sanitize_sim_type

#clockwise around center, starting from upper left
oct_drctns = [(-1, 2), (1, 2), (2, 1), (2, -1),
                    (1, -2), (-1, -2), (-2, -2), (-2, 1)]

sq_drctns = [(-1, 1), (1, 1), (1, -1), (-1, -1)]

class HardCodeSquoctSim():
    """
    Sandbox for trying to simulate the square-octagon code directly using 
    Cliffords and small error models. Copypasta of HardToricCodeSim.  
    """
    def __init__(self, size, p, n_trials):
        self.size = size
        self.p = p
        self.n_trials = n_trials
        self.logical_error = []
    
    def run(self, sim_type='cb'):
        #sanitize input
        _sanitize_sim_type(sim_type)

        sz = self.size

        lat = pq.SquareOctagonLattice((sz, sz))
        d_lat = pq.UnionJackLattice((sz, sz), is_dual=True)

        d_lat_lst = [pq.UnionJackLattice((sz, sz), is_dual=True)
                        for _ in range(sz + 1)]

        decoder = pq.ft_mwpm_decoder(lat, d_lat_lst)

        log_ops = pq.squoct_log_ops(lat.total_size)
        
        x_flip = pq.PauliErrorModel({q.I : 1. - self.p, q.X : self.p})
        z_flip = pq.PauliErrorModel({q.I : 1. - self.p, q.Z : self.p})
        dep = pq.depolarizing_model(self.p)
        twirl = pq.two_bit_twirl(self.p)

        z_prs = {shft : pq.oct_pairs(lat, d_lat, shft, oct_type='z')
                    for shft in oct_drctns}
        x_prs = {shft : pq.oct_pairs(lat, d_lat, shft, oct_type='x')
                    for shft in oct_drctns}
        sq_prs = {shft : pq.sq_pairs(lat, d_lat, shft)
                    for shft in sq_drctns}

        z_oct_cx = {drctn : pq.Clifford(q.cnot(2, 0, 1), z_prs[drctn])
                for drctn in oct_drctns}
        x_oct_xc = {drctn : pq.Clifford(q.cnot(2, 1, 0), x_prs[drctn])
                for drctn in oct_drctns}
        
        x_sq_xc = {drctn : pq.Clifford(q.cnot(2, 1, 0), sq_prs[drctn])
                for drctn in sq_drctns}

        z_sq_cx = {drctn : pq.Clifford(q.cnot(2, 0, 1), sq_prs[drctn])
                for drctn in sq_drctns}

        x_oct_meas = pq.Measurement(q.X, ['', 'Z'], d_lat.octagon_centers(oct_type='X'))
        z_oct_meas = pq.Measurement(q.Z, ['', 'X'], d_lat.octagon_centers(oct_type='Z'))
        x_sq_meas = pq.Measurement(q.X, ['', 'Z'], d_lat.square_centers())
        z_sq_meas = pq.Measurement(q.Z, ['', 'X'], d_lat.square_centers())
        
        noiseless_code = pq.squoct_code(lat, d_lat_lst[-1])
        
        for _ in range(self.n_trials):
            #clear last sim
            pq.error_fill(lat, q.I)
            d_lat.clear()
            for ltc in d_lat_lst:
                ltc.clear() #may break
            pq.error_fill(d_lat, q.I)
            
            #fill d_lat_lst with syndromes by copying
            for idx in range(sz - 1):
                meas_cycle(lat, d_lat, x_flip, z_flip, dep, twirl, 
                            z_prs, x_prs, sq_prs, z_oct_cx, x_oct_xc,
                            z_sq_cx, x_sq_xc, z_oct_meas, x_oct_meas,
                            z_sq_meas, x_sq_meas, sim_type=sim_type)
                pq.syndrome_copy(d_lat, d_lat_lst[idx])

            #noise is now already on the syndrome qubits (including meas. noise)
            noiseless_code.measure()
            #run decoder, with no final lattice check (laaaaater)
            decoder.infer()

            # Error checking, if the resulting Pauli is not in the
            # normalizer, chuck an error:

            d_lat_lst[-1].clear()
            noiseless_code.measure()
            for point in d_lat_lst[-1].points:
                if point.syndrome:
                    raise ValueError('Product of "inferred error"' 
                                     ' with actual error anticommutes'
                                     ' with some stabilizers.')
            
            com_relation_list = []
            for operator in log_ops:
                com_relation_list.append(operator.test(lat))
            self.logical_error.append(com_relation_list)

        pass

    def save(self, filename):
        big_dict = {}
        big_dict['lattice_class'] = 'SquareOctagonLattice'
        big_dict['lattice_size'] = self.size
        big_dict['dual_lattice_class'] = 'UnionJackLattice'
        big_dict['dual_lattice_size'] = self.size
        big_dict['error_model'] = 'custom hard-coded'
        big_dict['code'] = 'Square-Octagon Code'
        big_dict['decoder'] = 'FT MWPM'
        big_dict['n_trials'] = self.n_trials
        big_dict['logical_errors'] = self.logical_error

        with open(filename, 'w') as phil:
            pkl.dump(big_dict, phil)

def meas_cycle(lat, d_lat, x_flip, z_flip, dep, twirl, 
                z_prs, x_prs, sq_prs, z_oct_cx, x_oct_xc,
                z_sq_cx, x_sq_xc, z_oct_meas, x_oct_meas,
                z_sq_meas, x_sq_meas, sim_type):
    
    """
    Does one cycle of measurement for the square-octagon code on a
    SquareOctagonLattice. There are three possible error models, these
    are selected using the input variable `sim_type`. One in
    which only the data qubits get errors (`sim_type='p'`), one in 
    which the data and syndrome qubits are subject to independent flip
    errors (`sim_type='pq'`), and an all-serial circuit-based model
    (`sim_type='cb'`). The all-serial measurement schedule is as 
    follows:
    
    + Flip syndrome qubits with appropriate error type
    + Apply Z-octagon CNots, alternated with two-qubit twirl and
      depolarization on remaining qubits 
    + Apply X-octagon CNots, alternated with two-qubit-twirl and
      depolarization on remaining qubits
    + Apply Z-square CNots, alternated with two-qubit-twirl and
      depolarization on remaining qubits
    + Apply X-square CNots, alternated with two-qubit-twirl and
      depolarization on remaining qubits
    + Flip syndrome qubits with appropriate error type
    + Measure syndrome qubits.

    If the model is an RPGM, that's:
    + Flip data qubits with IIDXZ errors
    + Apply clean gates in the same sequence as above
    + Flip syndrome qubits with appropriate error type
    + Measure syndrome qubits. 
    """
    d_lat.clear()
    pq.error_fill(d_lat, q.I)
    
    if sim_type == 'cb':
        x_flip.act_on(d_lat.oct_centers('Z'))
        x_flip.act_on(d_lat.sq_centers())
        z_flip.act_on(d_lat.oct_centers('X'))
        z_flip.act_on(d_lat.sq_centers())
        
    if sim_type in ['pq', 'p']:
        x_flip.act_on(lat)
        z_flip.act_on(lat)

    for gate_set, drctns in zip([z_oct_cx, x_oct_xc, z_sq_cx, x_sq_xc],
                                [oct_drctns, oct_drctns, sq_drctns, sq_drctns])
        #ALSO ZIP PAIRS
        for drctn in drctns:
            gate_set[drctn].apply()
            if sim_type == 'cb':
                twirl.act_on(odd_prs[drctn])
            #DEPOLARIZE COMPLEMENT OF CNOT
        if sim_type in ['cb', 'pq']:
            x_flip.act_on(d_lat.plaq_centers())
            z_flip.act_on(d_lat.star_centers())
        elif sim_type == 'p':
            pass

    x_meas.apply()
    z_meas.apply()