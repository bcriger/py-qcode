import py_qcode as pq, qecc as q, cPickle as pkl

DRCTNS = 'nwes'
SIM_TYPES = ['cb', 'pq', 'p', 'stats']

class HardCodeToricSim():
    """
    Sandbox for trying to simulate the toric code directly using 
    Cliffords and small error models. This should really be a subclass
    of Simulation, but I'll figure that out later.  
    """
    def __init__(self, size, p, n_trials, vert_dist=None):
        self.size = size
        #setup dict of probabilities
        
        if hasattr(p, 'items'):
            self.p = p
        else:
            self.p = {err_type: p for err_type in 
                        ['dep', 'twirl', 'prep', 'meas']}
        
        self.n_trials = n_trials
        self.logical_error = []
        self.vert_dist = vert_dist
        #for stats
        self.data_errors = {'X': 0, 'Y': 0, 'Z': 0}
        self.syndrome_errors = {'x': 0, 'z': 0}
        self.sim_type = None
    
    def run(self, sim_type='cb'):
        #sanitize input
        _sanitize_sim_type(sim_type)
        self.sim_type = sim_type

        sz = self.size
        lat = pq.SquareLattice((sz, sz))
        d_lat = pq.SquareLattice((sz, sz), is_dual=True)

        if sim_type in ['cb', 'pq']:
            d_lat_lst = [pq.SquareLattice((sz, sz), is_dual=True)
                        for _ in range(sz)]
        elif sim_type == 'stats':
            d_lat_lst = [pq.SquareLattice((sz, sz), is_dual=True)
                        for _ in range(2)]

        
        decoder = pq.ft_mwpm_decoder(lat, d_lat_lst, blossom=False, 
                                        vert_dist=self.vert_dist)

        log_ops = pq.toric_log_ops((sz, sz))
        
        x_flip = {key : pq.PauliErrorModel({q.I : 1. - self.p[key],
                                            q.X : self.p[key]})
                                            for key in ['prep', 'meas']
                                            }
        z_flip = {key : pq.PauliErrorModel({q.I : 1. - self.p[key],
                                            q.Z : self.p[key]})
                                            for key in ['prep', 'meas']
                                            }
        twirl = pq.two_bit_twirl(self.p['twirl'])

        odd_prs = {drctn : pq.nwes_pairs(lat, d_lat, drctn, 'odd') 
                    for drctn in DRCTNS}
        even_prs = {drctn : pq.nwes_pairs(lat, d_lat, drctn, 'even') 
                    for drctn in DRCTNS}
        cx = {drctn : pq.Clifford(q.cnot(2,0,1), odd_prs[drctn])
                for drctn in DRCTNS}
        xc = {drctn : pq.Clifford(q.cnot(2, 1, 0), even_prs[drctn])
                for drctn in DRCTNS}
        
        x_meas = pq.Measurement(q.X, ['', 'Z'], d_lat.star_centers())
        z_meas = pq.Measurement(q.Z, ['', 'X'], d_lat.plaq_centers())
        
        noiseless_code = pq.toric_code(lat, d_lat_lst[-1])
        
        if sim_type == 'stats':
            synd_keys = ['x', 'z']
            synd_types = ['Z', 'X']
            crd_sets = [pq._even_evens(sz, sz), 
                        pq._odd_odds(sz, sz)]

        for _ in range(self.n_trials):
            #clear last sim
            pq.error_fill(lat, q.I)
            d_lat.clear()
            for ltc in d_lat_lst:
                ltc.clear() #may break
            pq.error_fill(d_lat, q.I)
            
            #fill d_lat_lst with syndromes by copying
            for idx in range(len(d_lat_lst) - 1):
                meas_cycle(lat, d_lat, x_flip, z_flip, twirl, odd_prs,
                            even_prs, cx, xc, x_meas, z_meas, 
                            sim_type=sim_type)
                pq.syndrome_copy(d_lat, d_lat_lst[idx])

            #print d_lat 
            noiseless_code.measure()
            if sim_type in ['cb', 'pq']:
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
            elif sim_type == 'stats':
                for key in self.data_errors.keys():
                    # raise Exception
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
        pass

    def save(self, filename):
        big_dict = {}
        big_dict['lattice_class'] = 'SquareLattice'
        big_dict['lattice_size'] = self.size
        big_dict['dual_lattice_class'] = 'Dual SquareLattice'
        big_dict['dual_lattice_size'] = self.size
        big_dict['error_model'] = 'custom hard-coded'
        big_dict['code'] = 'Toric Code'
        big_dict['n_trials'] = self.n_trials
        
        if self.sim_type == 'cb':
            big_dict['decoder'] = 'FT MWPM'
            big_dict['logical_errors'] = self.logical_error    
        elif self.sim_type == 'stats':
            big_dict['data_errors'] = self.data_errors
            big_dict['syndrome_errors'] = self.syndrome_errors
        
        with open(filename, 'w') as phil:
            pkl.dump(big_dict, phil)

class BellStateToricSim(HardCodeToricSim):
    """
    Testing out a hunch from Landahl et al, you're okay to serialize 
    the syndrome extraction circuit as long as you measure one syndrome
    type at a time; X or Z.  
    """
    def run(self, sim_type='cb'):
        #sanitize input
        _sanitize_sim_type(sim_type)
        self.sim_type = sim_type

        sz = self.size
        lat = pq.SquareLattice((sz, sz))
        d_lat = pq.SquareLattice((sz, sz), is_dual=True)

        if sim_type in ['cb', 'pq']:
            d_lat_lst = [pq.SquareLattice((sz, sz), is_dual=True)
                        for _ in range(sz)]
        elif sim_type == 'stats':
            d_lat_lst = [pq.SquareLattice((sz, sz), is_dual=True)
                        for _ in range(2)]
            #get correct syndrome at z measurement time
            early_d_lat = pq.SquareLattice((sz, sz), is_dual=True)

        
        decoder = pq.ft_mwpm_decoder(lat, d_lat_lst, blossom=False, 
                                        vert_dist=self.vert_dist)

        log_ops = pq.toric_log_ops((sz, sz))
        
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

        z_prs = {d : pq.nwes_pairs(lat, d_lat, d) for d in 'en'}
        x_prs = {d : pq.nwes_pairs(lat, d_lat, d) for d in 'sw'}
        
        p_prs = zip(d_lat.star_centers(), d_lat.plaq_centers())

        z_cycle = map(pq.Timestep, [[pq.Clifford(q.cnot(2, 0, 1), z_prs[d])]
                for d in 'en'])
        x_cycle = map(pq.Timestep, [[pq.Clifford(q.cnot(2, 1, 0), x_prs[d])]
                for d in 'sw'])
        
        prep_step_x = pq.Timestep([pq.Clifford(q.cnot(2, 0, 1), p_prs)])
        prep_step_z = pq.Timestep([pq.Clifford(q.cnot(2, 1, 0), p_prs)])

        x_meas = pq.Measurement(q.X, ['', 'Z'], d_lat.points)
        z_meas = pq.Measurement(q.Z, ['', 'X'], d_lat.points)
        
        noiseless_code = pq.toric_code(lat, d_lat_lst[-1])
        
        if sim_type == 'stats':
            synd_keys = ['x', 'z']
            synd_types = ['Z', 'X']
            crd_sets = [pq._even_evens(sz, sz), 
                        pq._odd_odds(sz, sz)]

        for _ in range(self.n_trials):
            #clear last sim
            pq.error_fill(lat, q.I)
            for ltc in d_lat_lst:
                ltc.clear() #may break
            
            #fill d_lat_lst with syndromes by copying
            for idx in range(len(d_lat_lst) - 1):
                d_lat.clear()
                pq.error_fill(d_lat, q.I)
                #z stabilisers
                
                #flip bell state ancillas
                z_flip['prep'].act_on(d_lat.star_centers())
                x_flip['prep'].act_on(d_lat.plaq_centers())
                prep_step_z.noisy_apply(None, None, self.p['twirl'], 0., False)
                dep.act_on(lat)
                for step in z_cycle:
                    step.noisy_apply(None, None, self.p['twirl'], 0., False)
                z_flip['meas'].act_on(d_lat)
                z_meas.apply()
                for star, plaq, o_plaq in zip(d_lat.star_centers(),
                                                d_lat.plaq_centers(),
                                                d_lat_lst[idx].plaq_centers()):
                    o_plaq.syndrome = 'X' if star.syndrome != plaq.syndrome else ''
                if sim_type == 'stats':
                    noiseless_code.measure()
                    pq.syndrome_copy(d_lat_lst[-1], early_d_lat)        
                    d_lat_lst[-1].clear()        

                d_lat.clear()
                pq.error_fill(d_lat, q.I)
                #x stabilisers
                
                #flip bell state ancillas
                x_flip['prep'].act_on(d_lat.star_centers())
                z_flip['prep'].act_on(d_lat.plaq_centers())
                prep_step_x.noisy_apply(None, None, self.p['twirl'], 0., False)
                dep.act_on(lat)
                for step in x_cycle:
                    step.noisy_apply(None, None, self.p['twirl'], 0., False)
                x_flip['meas'].act_on(d_lat)
                x_meas.apply()
                for star, plaq, o_star in zip(d_lat.star_centers(),
                                                d_lat.plaq_centers(),
                                                d_lat_lst[idx].star_centers()):
                    o_star.syndrome = 'Z' if star.syndrome != plaq.syndrome else ''
                    
                
            #print d_lat 
            noiseless_code.measure()
            if sim_type in ['cb', 'pq']:
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
            elif sim_type == 'stats':
                for key in self.data_errors.keys():
                    # raise Exception
                    for point in lat.points:
                        if point.error.op == key:
                            self.data_errors[key] += 1
                #check syndromes
                for synd_key, synd_type, crd_set in zip(synd_keys,
                                                        synd_types,
                                                        crd_sets):
                    if synd_key == 'x':
                        for crd in crd_set:
                            if (synd_type in d_lat_lst[0][crd].syndrome) != \
                                (synd_type in d_lat_lst[1][crd].syndrome):
                                self.syndrome_errors[synd_key] += 1
                    elif synd_key == 'z':
                        for crd in crd_set:
                            if (synd_type in d_lat_lst[0][crd].syndrome) != \
                                (synd_type in early_d_lat[crd].syndrome):
                                self.syndrome_errors[synd_key] += 1
        pass

def meas_cycle(lat, d_lat, x_flip, z_flip, twirl, odd_prs, even_prs, 
                cx, xc, x_meas, z_meas, sim_type):
    
    """
    Does one cycle of measurement for the toric code on a
    SquareLattice. If the error model is circuit-based, that's: 
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

    If the model is an RPGM, that's:
    + Apply north-facing CNots
    + Apply west-facing CNots
    + Apply east-facing CNots
    + Apply south-facing CNots
    + Flip syndrome qubits with appropriate error type
    + Measure syndrome qubits. 
    """
    d_lat.clear()
    pq.error_fill(d_lat, q.I)
    
    if sim_type in ['cb', 'stats']:
        x_flip['prep'].act_on(d_lat.plaq_centers())
        z_flip['prep'].act_on(d_lat.star_centers())
    elif sim_type in ['pq', 'p']:
        x_flip['prep'].act_on(lat)
        z_flip['prep'].act_on(lat)

    for drctn in DRCTNS:
        cx[drctn].apply()
        if sim_type in ['cb', 'stats']:
            twirl.act_on(odd_prs[drctn])
        
        xc[drctn].apply()
        if sim_type in ['cb', 'stats']:
            twirl.act_on(even_prs[drctn])
    
    if sim_type in ['cb', 'pq', 'stats']:
        x_flip['meas'].act_on(d_lat.plaq_centers())
        z_flip['meas'].act_on(d_lat.star_centers())
    elif sim_type == 'p':
        pass

    x_meas.apply()
    z_meas.apply()

def _sanitize_sim_type(sim_type):
    if sim_type not in SIM_TYPES:
        raise ValueError("Simulation type (sim_type) must be "
            "either 'cb' (circuit-based), 'pq' (3D Z2 RPGM), or 'p' "
            "(no syndrome errors). {} entered.".format(sim_type))
