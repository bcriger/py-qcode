import py_qcode as pq, qecc as q, cPickle as pkl
from hard_toric_code import SIM_TYPES, _sanitize_sim_type
from itertools import chain

#clockwise around center, starting from upper left
oct_clck_dirs = [(-1, 2), (1, 2), (2, 1), (2, -1),
                    (1, -2), (-1, -2), (-2, -1), (-2, 1)]

sq_clck_dirs = [(-1, 1), (1, 1), (1, -1), (-1, -1)]

#directions for interleaved simulation
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

        if hasattr(p, 'items'):
            self.p = p
        else:
            self.p = {err_type: p for err_type in 
                        ['dep', 'twirl', 'prep', 'meas']}
        
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
        
        x_flip = {key : pq.PauliErrorModel({q.I : 1. - self.p[key],
                                            q.X : self.p[key]})
                                            for key in ['prep', 'meas']
                                            }
        z_flip = {key : pq.PauliErrorModel({q.I : 1. - self.p[key],
                                            q.Z : self.p[key]})
                                            for key in ['prep', 'meas']
                                            }
        dep = pq.depolarizing_model(self.p['dep'])
        twirl = pq.two_bit_twirl(self.p['twirl'])

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
        big_dict = _save_dict(self)
        with open(filename, 'w') as phil:
            pkl.dump(big_dict, phil)

class InterleavedSquoctSim(HardCodeSquoctSim):
    """
    Marginally smarter than copypasta, I rewrite the run method for 
    HardCodeSquoctSim. 
    """
    def __init__(self, size, p, n_trials, vert_dist=None, oct_factor=1.):
        HardCodeSquoctSim.__init__(self, size, p, n_trials)
        self.vert_dist = vert_dist
        self.data_errors = {'X': 0, 'Y': 0, 'Z': 0}
        self.syndrome_errors = {'xv': 0, 'zv': 0, 'xh': 0, 'zh': 0,
                                'xo': 0, 'zo': 0}
        self.sim_type = None
        self.oct_factor = oct_factor
    
    def run(self, sim_type='cb'):
        #sanitize input
        _sanitize_sim_type(sim_type)
        if sim_type not in ['cb', 'stats']:
            raise ValueError("InterleavedSquoctSim only supports "
                            "sim_type='cb' or 'stats', use"
                            " HardCodeSquoctSim.")
        self.sim_type = sim_type
        sz = int(self.size / 2.)

        lat = pq.SquareOctagonLattice((sz, sz))
        d_lat = pq.UnionJackLattice((sz, sz), is_dual=True)
        # We won't be able to measure multiple syndromes on the same
        # point, so I make a spare dual lattice for the X square 
        # syndromes:
        d_lat_x_sq = pq.UnionJackLattice((sz,sz), is_dual=True)

        if sim_type == 'cb':
            d_lat_lst = [pq.UnionJackLattice((sz, sz), is_dual=True)
                        for _ in range(self.size + 1)]
        elif sim_type == 'stats':
            d_lat_lst = [pq.UnionJackLattice((sz, sz), is_dual=True)
                        for _ in range(2)]

        decoder = pq.ft_mwpm_decoder(lat, d_lat_lst, blossom=False, 
                                        vert_dist=self.vert_dist)
        noiseless_code = pq.square_octagon_code(lat, d_lat_lst[-1])

        log_ops = pq.squoct_log_ops(lat.total_size)
        
        dep = pq.depolarizing_model(self.p['dep'])
        x_flip = {key : pq.PauliErrorModel({q.I : 1. - self.p[key],
                                            q.X : self.p[key]})
                                            for key in ['prep', 'meas']
                                            }
        z_flip = {key : pq.PauliErrorModel({q.I : 1. - self.p[key],
                                            q.Z : self.p[key]})
                                            for key in ['prep', 'meas']
                                            }
        synd_flip = {}
        synd_flip['prep'] = {q.X : z_flip['prep'], q.Z : x_flip['prep']}
        synd_flip['meas'] = {q.X : z_flip['meas'], q.Z : x_flip['meas']}
        
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

        sq_cycle = map(pq.Timestep, zip(v_x_cnots + v_z_cnots,
                                        h_z_cnots + h_x_cnots))
        oct_cycle = map(pq.Timestep, zip(o_x_cnots, o_z_cnots))

        if sim_type == 'stats':
            synd_keys = ['xv', 'zv', 'xh', 'zh', 'xo', 'zo']
            synd_types = ['Z', 'X', 'Z', 'X', 'Z', 'X']
            crd_sets = [
                        pq._square_centers((sz, sz), 'v'),
                        pq._square_centers((sz, sz), 'v'),
                        pq._square_centers((sz, sz), 'h'),
                        pq._square_centers((sz, sz), 'h'),
                        pq._octagon_centers((sz, sz), 'x'),
                        pq._octagon_centers((sz, sz), 'z')
                        ]

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
                    synd_flip['prep'][meas.pauli].act_on(meas.point_set)
                #first four noisy gates
                for tdx in range(4):
                    sq_cycle[tdx].noisy_apply(lat, None, self.p['twirl'], False)
                    oct_cycle[tdx].noisy_apply(lat, None, 
                            self.oct_factor * self.p['twirl'], False)
                    for pt in lat.points:
                        if not any([pt in seq.twirl_support for seq in sq_cycle[tdx], oct_cycle[tdx]]):
                            dep.act_on(pt)
                    for pt in d_lat.points + d_lat_x_sq.points:
                        if not any([pt in sprt 
                                    for sprt in map(lambda a: a.support,
                                        [v_x_cnots[tdx], h_z_cnots[tdx],
                                         o_x_cnots[tdx], o_z_cnots[tdx]] )]):
                            dep.act_on(pt)

                #first measurement round
                for meas in [x_v_meas, z_h_meas]:
                    synd_flip['meas'][meas.pauli].act_on(meas.point_set)
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
                    synd_flip['prep'][meas.pauli].act_on(meas.point_set)

                #next 4 noisy gates
                for tdx in range(4, 8):
                    sq_cycle[tdx].noisy_apply(lat, None, self.p['twirl'], False)
                    oct_cycle[tdx].noisy_apply(lat, None, 
                            self.oct_factor * self.p['twirl'], False)
                    for pt in lat.points:
                        if not any([pt in seq.twirl_support for seq in sq_cycle[tdx], oct_cycle[tdx]]):
                            dep.act_on(pt)
                    for pt in d_lat.points + d_lat_x_sq.points:
                        if not any([pt in sprt 
                                    for sprt in map(lambda a: a.support,
                                        [v_z_cnots[tdx - 4], h_x_cnots[tdx - 4],
                                         o_x_cnots[tdx], o_z_cnots[tdx]] )]):
                            dep.act_on(pt)

                #flip and measure remaining ancillas
                for meas in [x_h_meas, z_v_meas, x_o_meas, z_o_meas]:
                    synd_flip['meas'][meas.pauli].act_on(meas.point_set)
                    meas.apply()

                #copy syndromes onto 3D lattice.
                pq.syndrome_copy(d_lat, d_lat_lst[idx], append=True)
                pq.syndrome_copy(d_lat_x_sq, d_lat_lst[idx], append=True)

            #noise is now already on the syndrome qubits (including meas. noise)
            noiseless_code.measure()
            if sim_type == 'cb':
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
            elif sim_type == 'stats':
                for key in self.data_errors.keys():
                    for point in lat.points:
                        if point.error.op == key:
                            self.data_errors[key] += 1
                #check syndromes
                for synd_key, synd_type, crd_set in zip(synd_keys,
                                                        synd_types,
                                                        crd_sets):
                    for crd in crd_set:
                        if (synd_type in d_lat_lst[0][crd].syndrome) != \
                            (synd_type in d_lat_lst[1][crd].syndrome):
                            self.syndrome_errors[synd_key] += 1
                        
    def save(self, filename):
        if self.sim_type == 'cb':
            big_dict = _save_dict(self)
            big_dict['oct_factor'] = self.oct_factor
        elif self.sim_type == 'stats':
            big_dict = {}
            big_dict['lattice_class'] = 'SquareOctagonLattice'
            big_dict['lattice_size'] = self.size
            big_dict['dual_lattice_class'] = 'UnionJackLattice'
            big_dict['dual_lattice_size'] = self.size
            big_dict['error_model'] = 'custom hard-coded'
            big_dict['code'] = 'Square-Octagon Code'
            big_dict['n_trials'] = self.n_trials
            big_dict['data_errors'] = self.data_errors
            big_dict['syndrome_errors'] = self.syndrome_errors
            big_dict['oct_factor'] = self.oct_factor

            with open(filename, 'w') as phil:
                pkl.dump(big_dict, phil)


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
        x_flip['prep'].act_on(lat)
        z_flip['prep'].act_on(lat)

    synd_flip = {}
    synd_flip['prep'] = {q.X : z_flip['prep'], q.Z : x_flip['prep']}
    synd_flip['meas'] = {q.X : z_flip['meas'], q.Z : x_flip['meas']}

    for gate_set, clck_dirs, meas, deps in zip([z_oct_cx, x_oct_xc, z_sq_cx, x_sq_xc],
                                            [oct_clck_dirs, oct_clck_dirs, sq_clck_dirs, sq_clck_dirs],
                                            [z_oct_meas, x_oct_meas, z_sq_meas, x_sq_meas],
                                            [z_deps, x_deps, sq_deps, sq_deps]):
        for pt in meas.point_set:
            pt.syndrome = ''
    
        if sim_type == 'cb':
            synd_flip['prep'][meas.pauli].act_on(meas.point_set)
    
        for drctn in clck_dirs:
            gate_set[drctn].apply()
            if sim_type == 'cb':
                twirl.act_on(gate_set[drctn].point_sets)
                dep.act_on(deps[drctn])
    
        if sim_type in ['pq', 'cb']:
            synd_flip['meas'][meas.pauli].act_on(meas.point_set)
        
        meas.apply()

def _save_dict(sim):
    big_dict = {}
    big_dict['lattice_class'] = 'SquareOctagonLattice'
    big_dict['lattice_size'] = sim.size
    big_dict['dual_lattice_class'] = 'UnionJackLattice'
    big_dict['dual_lattice_size'] = sim.size
    big_dict['error_model'] = 'custom hard-coded'
    big_dict['code'] = 'Square-Octagon Code'
    big_dict['decoder'] = 'FT MWPM'
    big_dict['n_trials'] = sim.n_trials
    big_dict['logical_errors'] = sim.logical_error
    big_dict['oct_factor'] = sim.oct_factor
    return big_dict