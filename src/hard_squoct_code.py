import py_qcode as pq, qecc as q, cPickle as pkl
from hard_toric_code import SIM_TYPES, _sanitize_sim_type
from itertools import chain

#clockwise around center, starting from upper left
oct_clck_dirs = [(-1, 2), (1, 2), (2, 1), (2, -1),
                    (1, -2), (-1, -2), (-2, -1), (-2, 1)]

sq_clck_dirs = [(-1, 1), (1, 1), (1, -1), (-1, -1)]

#directions for interleaved simulation
# z_oct_il_dirs = [(2, 1), (2, -1), (1, -2), (-1, -2), 
#                     (-2, -1), (-2, 1), (-1, 2), (1, 2)]

# x_oct_il_dirs = [(2, 1), (2, -1), (1, -2), (1, 2),
#                     (-1, -2), (-2, 1), (-2, -1), (-1, 2)]
z_oct_il_dirs = [(1, 2), (-1, 2), (-2, 1), (-2, -1),
                    (2, 1), (1, -2), (2, -1), (-1, -2)]

x_oct_il_dirs = z_oct_il_dirs

#gates are arranged so that square CNots are always parallel.
#square measurement is broken up into two rounds
sq_il_dirs = [[(1, 1), (1, -1), (-1, 1), (-1, -1)][x] 
                for x in [0, 1, 2, 3, 3, 1, 2, 0]]

def pair_complements(lat, pair_list):
    """
    In order to figure out which qubits to depolarize while twirling a
    bunch of others, I present a function to determine which points 
    from a lattice are not in the support of a list of pairs.
    """
    pair_support = set(chain.from_iterable(pair_list))
    return filter(lambda pt: pt not in pair_support, lat.points)

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

        sz = int(self.size / 2.)

        lat = pq.SquareOctagonLattice((sz, sz))
        d_lat = pq.UnionJackLattice((sz, sz), is_dual=True)
        # We won't be able to measure multiple syndromes on the same
        # point, so I make a spare dual lattice for the X square 
        # syndromes:
        d_lat_x_sq = pq.UnionJackLattice((sz,sz), is_dual=True)

        d_lat_lst = [pq.UnionJackLattice((sz, sz), is_dual=True)
                        for _ in range(self.size + 1)]

        decoder = pq.ft_mwpm_decoder(lat, d_lat_lst, blossom=False)

        log_ops = pq.squoct_log_ops(lat.total_size)
        
        x_flip = pq.PauliErrorModel({q.I : 1. - self.p, q.X : self.p})
        z_flip = pq.PauliErrorModel({q.I : 1. - self.p, q.Z : self.p})
        dep = pq.depolarizing_model(self.p)
        twirl = pq.two_bit_twirl(self.p)

        z_prs = {shft : pq.oct_pairs(lat, d_lat, shft, oct_type='z')
                    for shft in oct_clck_dirs}
        z_deps = {shft : pair_complements(lat, z_prs[shft])
                    for shft in oct_clck_dirs}
        
        x_prs = {shft : pq.oct_pairs(lat, d_lat, shft, oct_type='x')
                    for shft in oct_clck_dirs}
        x_deps = {shft : pair_complements(lat, x_prs[shft])
                    for shft in oct_clck_dirs}
        
        sq_prs = {shft : pq.sq_pairs(lat, d_lat, shft)
                    for shft in sq_clck_dirs}

        x_sq_prs = {shft : pq.sq_pairs(lat, d_lat_x_sq, shft)
                    for shft in sq_clck_dirs}
        
        sq_deps = {shft : pair_complements(lat, sq_prs[shft])
                    for shft in sq_clck_dirs}
        
        z_oct_cx = {drctn : pq.Clifford(q.cnot(2, 0, 1), z_prs[drctn])
                for drctn in oct_clck_dirs}
        
        x_oct_xc = {drctn : pq.Clifford(q.cnot(2, 1, 0), x_prs[drctn])
                for drctn in oct_clck_dirs}
        
        x_sq_xc = {drctn : pq.Clifford(q.cnot(2, 1, 0), x_sq_prs[drctn])
                for drctn in sq_clck_dirs}

        z_sq_cx = {drctn : pq.Clifford(q.cnot(2, 0, 1), sq_prs[drctn])
                for drctn in sq_clck_dirs}

        x_oct_meas = pq.Measurement(q.X, ['', 'Z'], d_lat.octagon_centers(oct_type='X'))
        z_oct_meas = pq.Measurement(q.Z, ['', 'X'], d_lat.octagon_centers(oct_type='Z'))
        x_sq_meas = pq.Measurement(q.X, ['', 'Z'], d_lat_x_sq.square_centers())
        z_sq_meas = pq.Measurement(q.Z, ['', 'X'], d_lat.square_centers())
        
        noiseless_code = pq.square_octagon_code(lat, d_lat_lst[-1])
        
        for _ in range(self.n_trials):
            #clear last sim
            pq.error_fill(lat, q.I)
            d_lat.clear()
            d_lat_x_sq.clear()
            
            for ltc in d_lat_lst:
                ltc.clear() #may break
                
            pq.error_fill(d_lat, q.I)
            pq.error_fill(d_lat_x_sq, q.I)
            
            #fill d_lat_lst with syndromes by copying
            for idx in range(len(d_lat_lst) - 1):
                
                meas_cycle(lat, d_lat, d_lat_x_sq, x_flip, z_flip, dep, twirl, 
                            z_prs, x_prs, sq_prs, z_deps, x_deps, 
                            sq_deps, z_oct_cx, x_oct_xc,
                            z_sq_cx, x_sq_xc, z_oct_meas, x_oct_meas,
                            z_sq_meas, x_sq_meas, sim_type=sim_type)
                    
                pq.syndrome_copy(d_lat, d_lat_lst[idx])
                pq.syndrome_copy(d_lat_x_sq, d_lat_lst[idx], append=True)

            #noise is now already on the syndrome qubits (including meas. noise)
            noiseless_code.measure()
            #run decoder, with no final lattice check (laaaaater)
            decoder.infer()

            # Error checking, if the resulting Pauli is not in the
            # normalizer, chuck an error:

            d_lat_lst[-1].clear()
            pq.syndrome_fill(d_lat_lst[-1], '')
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

class InterleavedSquoctSim(HardCodeSquoctSim):
    """
    Marginally smarter than copypasta, I rewrite the run method for 
    HardCodeSquoctSim. 
    """
    def __init__(self, size, p, n_trials, vert_dist=None):
        HardCodeSquoctSim.__init__(self, size, p, n_trials)
        self.vert_dist = vert_dist
    
    def run(self, sim_type='cb'):
        #sanitize input
        _sanitize_sim_type(sim_type)
        if sim_type != 'cb':
            raise ValueError("InterleavedSquoctSim only supports "
                            "sim_type=cb, use HardCodeSquoctSim.")
        
        sz = int(self.size / 2.)

        lat = pq.SquareOctagonLattice((sz, sz))
        d_lat = pq.UnionJackLattice((sz, sz), is_dual=True)
        # We won't be able to measure multiple syndromes on the same
        # point, so I make a spare dual lattice for the X square 
        # syndromes:
        d_lat_x_sq = pq.UnionJackLattice((sz,sz), is_dual=True)

        d_lat_lst = [pq.UnionJackLattice((sz, sz), is_dual=True)
                        for _ in range(self.size + 1)]

        decoder = pq.ft_mwpm_decoder(lat, d_lat_lst, blossom=False)
        noiseless_code = pq.square_octagon_code(lat, d_lat_lst[-1])

        log_ops = pq.squoct_log_ops(lat.total_size)
        
        dep = pq.depolarizing_model(self.p)
        x_flip = pq.PauliErrorModel({q.I : 1. - self.p, q.X : self.p})
        z_flip = pq.PauliErrorModel({q.I : 1. - self.p, q.Z : self.p})
        synd_flip = {q.X : z_flip, q.Z : x_flip}
        
        v_x_prs = [pq.sq_pairs(lat, d_lat_x_sq, d, 'v') for d in sq_il_dirs[:4]]
        v_z_prs = [pq.sq_pairs(lat, d_lat, d, 'v') for d in sq_il_dirs[4:]]
        h_x_prs = [pq.sq_pairs(lat, d_lat_x_sq, d, 'h') for d in sq_il_dirs[4:]]
        h_z_prs = [pq.sq_pairs(lat, d_lat, d, 'h') for d in sq_il_dirs[:4]]
        o_x_prs = [pq.oct_pairs(lat, d_lat, d, 'x') for d in x_oct_il_dirs]
        o_z_prs = [pq.oct_pairs(lat, d_lat, d, 'z') for d in z_oct_il_dirs]
        
        v_x_cnots = [pq.Clifford(q.cnot(2, 1, 0), pr) for pr in v_x_prs]
        v_z_cnots = [pq.Clifford(q.cnot(2, 0, 1), pr) for pr in v_z_prs]
        h_x_cnots = [pq.Clifford(q.cnot(2, 1, 0), pr) for pr in h_x_prs]
        h_z_cnots = [pq.Clifford(q.cnot(2, 0, 1), pr) for pr in h_z_prs]
        o_x_cnots = [pq.Clifford(q.cnot(2, 1, 0), pr) for pr in o_x_prs]
        o_z_cnots = [pq.Clifford(q.cnot(2, 0, 1), pr) for pr in o_z_prs]

        x_o_meas = pq.Measurement(q.X, ['', 'Z'], d_lat.octagon_centers(oct_type='X'))
        z_o_meas = pq.Measurement(q.Z, ['', 'X'], d_lat.octagon_centers(oct_type='Z'))
        x_v_meas = pq.Measurement(q.X, ['', 'Z'], d_lat_x_sq.square_centers('v'))
        z_v_meas = pq.Measurement(q.Z, ['', 'X'], d_lat.square_centers('v'))
        x_h_meas = pq.Measurement(q.X, ['', 'Z'], d_lat_x_sq.square_centers('h'))
        z_h_meas = pq.Measurement(q.Z, ['', 'X'], d_lat.square_centers('h'))

        cycle = map(pq.Timestep, zip(v_x_cnots + v_z_cnots,
                                        h_z_cnots + h_x_cnots,
                                        o_x_cnots, o_z_cnots))

        for _ in xrange(self.n_trials):
            #clear last sim
            pq.error_fill(lat, q.I)
            d_lat.clear()
            d_lat_x_sq.clear()
            
            for ltc in d_lat_lst:
                ltc.clear() #may break
                
            pq.error_fill(d_lat, q.I)
            pq.error_fill(d_lat_x_sq, q.I)
            
            #fill d_lat_lst with syndromes by copying
            for idx in range(len(d_lat_lst) - 1):
                d_lat.clear()
                d_lat_x_sq.clear()
                pq.error_fill(d_lat, q.I)
                pq.error_fill(d_lat_x_sq, q.I)
                pq.syndrome_fill(d_lat, '')
                pq.syndrome_fill(d_lat_x_sq, '')
                #flip first round of ancillas
                for meas in [x_o_meas, z_o_meas, x_v_meas, z_h_meas]:
                    synd_flip[meas.pauli].act_on(meas.point_set)
                #first four noisy gates
                for tdx in range(4):
                    cycle[tdx].noisy_apply(lat, None, self.p, False)
                    for pt in lat.points:
                        if pt not in cycle[tdx].twirl_support:
                            dep.act_on(pt)
                    for pt in d_lat.points + d_lat_x_sq.points:
                        if not any([pt in sprt 
                                    for sprt in map(lambda a: a.support,
                                        [v_x_cnots[tdx], h_z_cnots[tdx],
                                         o_x_cnots[tdx], o_z_cnots[tdx]] )]):
                            dep.act_on(pt)

                #first measurement round
                for meas in [x_v_meas, z_h_meas]:
                    synd_flip[meas.pauli].act_on(meas.point_set)
                    meas.apply()
                
                pq.syndrome_copy(d_lat, d_lat_lst[idx])
                pq.syndrome_copy(d_lat_x_sq, d_lat_lst[idx], append=True)
                
                for pt in x_v_meas.point_set + z_h_meas.point_set:
                    pt.error = q.I
                    pt.syndrome = ''

                #depolarisation during measurement
                dep.act_on(lat)
                dep.act_on(d_lat.octagon_centers())

                #prep new square measurements
                for meas in [x_h_meas, z_v_meas]:
                    synd_flip[meas.pauli].act_on(meas.point_set)

                #next 4 noisy gates
                for tdx in range(4, 8):
                    cycle[tdx].noisy_apply(lat, None, self.p, False)
                    for pt in lat.points:
                        if pt not in cycle[tdx].twirl_support:
                            dep.act_on(pt)
                    for pt in d_lat.points + d_lat_x_sq.points:
                        if not any([pt in sprt 
                                    for sprt in map(lambda a: a.support,
                                        [v_z_cnots[tdx - 4], h_x_cnots[tdx - 4],
                                         o_x_cnots[tdx], o_z_cnots[tdx]] )]):
                            dep.act_on(pt)

                #flip and measure remaining ancillas
                for meas in [x_h_meas, z_v_meas, x_o_meas, z_o_meas]:
                    synd_flip[meas.pauli].act_on(meas.point_set)
                    meas.apply()

                #copy syndromes onto 3D lattice.
                pq.syndrome_copy(d_lat, d_lat_lst[idx], append=True)
                pq.syndrome_copy(d_lat_x_sq, d_lat_lst[idx], append=True)

            #noise is now already on the syndrome qubits (including meas. noise)
            noiseless_code.measure()
            #run decoder, with no final lattice check (laaaaater)
            decoder.infer()

            # Error checking, if the resulting Pauli is not in the
            # normalizer, chuck an error:

            d_lat_lst[-1].clear()
            pq.syndrome_fill(d_lat_lst[-1], '')
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


def meas_cycle(lat, d_lat, d_lat_x_sq, x_flip, z_flip, dep, twirl, 
                z_prs, x_prs, sq_prs, z_deps, x_deps, 
                sq_deps, z_oct_cx, x_oct_xc,
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
    d_lat_x_sq.clear()
    pq.error_fill(d_lat, q.I)
    pq.error_fill(d_lat_x_sq, q.I)
    pq.syndrome_fill(d_lat, '')
    pq.syndrome_fill(d_lat_x_sq, '')
            
    if sim_type in ['pq', 'p']:
        x_flip.act_on(lat)
        z_flip.act_on(lat)

    synd_flip = {q.X : z_flip, q.Z : x_flip}

    for gate_set, clck_dirs, meas, deps in zip([z_oct_cx, x_oct_xc, z_sq_cx, x_sq_xc],
                                            [oct_clck_dirs, oct_clck_dirs, sq_clck_dirs, sq_clck_dirs],
                                            [z_oct_meas, x_oct_meas, z_sq_meas, x_sq_meas],
                                            [z_deps, x_deps, sq_deps, sq_deps]):
        for pt in meas.point_set:
            pt.syndrome = ''
    
        if sim_type == 'cb':
            synd_flip[meas.pauli].act_on(meas.point_set)
    
        for drctn in clck_dirs:
            gate_set[drctn].apply()
            if sim_type == 'cb':
                twirl.act_on(gate_set[drctn].point_sets)
                dep.act_on(deps[drctn])
    
        if sim_type in ['pq', 'cb']:
            synd_flip[meas.pauli].act_on(meas.point_set)
        
        meas.apply()
