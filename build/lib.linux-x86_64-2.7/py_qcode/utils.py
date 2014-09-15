import cPickle as pkl

from simulation import Simulation, FTSimulation
from lattice import SquareLattice, SquareOctagonLattice, UnionJackLattice
from error import depolarizing_model
from decoder import mwpm_decoder, ft_mwpm_decoder
from code import toric_code, noisy_toric_code, square_octagon_code, noisy_squoct_code
from logical_operators import toric_log_ops, squoct_log_ops

import cPickle as pkl

__all__ = ['sim_from_file', 'square_toric_code_sim', 'error_print',
            'syndrome_print', 'squoct_sim', 'noisy_toric_code_sim',
            'noisy_squoct_sim']

def error_print(lattice):
    printnone = True
    for point in lattice.points:
        if (point.error) and (not(point.error.op == 'I')):
            print point
            printnone = False

    if printnone:
        print "No Errors"

    pass

def syndrome_print(lattice):
    printnone = True
    for point in lattice.points:
        if point.syndrome:
            print point
            printnone = False

    if printnone:
        print "No Syndromes"
    pass

def sim_from_file(filename):
    """
    The purpose of this function is to:

    + open a file containing a pickled dictionary of input values to a simulation,

    + initialize the objects which the corresponding `py_qcode.Simulation` takes as input,
    
    + run the simulation, and 

    + save the results to a file of the same name as the input, with a different extension.  
    """
    #Obsolete, scavenging code for pickle-independent implementation
    #(you can't pickle functions).
    with open(filename,'r') as phil:
        sim_dict = pkl.load(phil)
    sim = Simulation(**sim_dict)
    sim.run()

    split_name = filename.split('.')
    try:
        file_prefix, file_ext = split_name
    except ValueError:
        raise ValueError('Filenames are assumed to be of the form'+\
        ' "prefix.ext".')

    output_name = '.'.join([file_prefix, 'out'])

    sim.save(output_name)

def square_toric_code_sim(size, error_rate, n_trials, filename):
    """
    This function is square in more than one sense; it does everything
    the most vanilla way possible, and it uses a square grid to define 
    the torus. You put in an integer size, an error rate and a number
    of trials to execute, and it produces a pickled dict containing 
    the input to a simulation object in a file.
    """
    
    sim_lattice = SquareLattice((size,size))
    sim_dual_lattice = SquareLattice((size,size), is_dual=True)
    sim_model = depolarizing_model(error_rate)
    sim_code = toric_code(sim_lattice, sim_dual_lattice)
    sim_decoder = mwpm_decoder(sim_lattice, sim_dual_lattice)
    sim_log_ops = toric_log_ops(sim_lattice.total_size)

    sim_keys = ['lattice', 'dual_lattice', 'error_model', 'code', 
                            'decoder', 'logical_operators', 'n_trials']

    sim_values = [sim_lattice, sim_dual_lattice, sim_model, sim_code, 
                                sim_decoder, sim_log_ops, n_trials]
    
    sim_dict = dict(zip(sim_keys, sim_values))

    sim = Simulation(**sim_dict)
    sim.run()
    sim.save(filename + '.sim')

def noisy_toric_code_sim(size, error_rate, n_trials, filename):
    """
    Slightly less vanilla than `square_toric_code_sim`, this function 
    assigns `error_rate` to both the `error_model` and `code` 
    properties of the simulation, such that the wrong syndrome is 
    returned from the code at that error rate. 
    """
    
    sim_lattice = SquareLattice((size,size))
    sim_dual_lattice_list = []
    for idx in xrange(size):
        sim_dual_lattice_list.append(SquareLattice((size,size), is_dual=True))
    sim_model = depolarizing_model(error_rate)
    sim_code_func = noisy_toric_code
    sim_decoder = ft_mwpm_decoder(sim_lattice, sim_dual_lattice_list)
    sim_log_ops = toric_log_ops((size,size))

    sim_keys = ['lattice', 'dual_lattice_list', 'error_model', 
                'synd_noise', 'code_func', 'decoder', 
                'logical_operators', 'n_trials']

    sim_values = [sim_lattice, sim_dual_lattice_list, sim_model, 
                    error_rate, sim_code_func, 
                    sim_decoder, sim_log_ops, n_trials]
    
    sim_dict = dict(zip(sim_keys, sim_values))

    sim = FTSimulation(**sim_dict)
    sim.run()
    sim.save(filename + '.sim')


def squoct_sim(size, error_rate, n_trials, filename):
    """
    Copypasta from square_toric_code_sim, maybe I'll un-pasta this 
    later. 
    """
    
    sim_lattice = SquareOctagonLattice((size,size))
    sim_dual_lattice = UnionJackLattice((size,size), is_dual=True)
    sim_model = depolarizing_model(error_rate)
    sim_code = square_octagon_code(sim_lattice, sim_dual_lattice)
    sim_decoder = mwpm_decoder(sim_lattice, sim_dual_lattice)
    sim_log_ops = squoct_log_ops(sim_lattice.total_size)

    sim_keys = ['lattice', 'dual_lattice', 'error_model', 'code', 
                            'decoder', 'logical_operators', 'n_trials']

    sim_values = [sim_lattice, sim_dual_lattice, sim_model, sim_code, 
                                sim_decoder, sim_log_ops, n_trials]
    
    sim_dict = dict(zip(sim_keys, sim_values))

    sim = Simulation(**sim_dict)
    sim.run()
    sim.save(filename + '.sim')

def noisy_squoct_sim(size, error_rate, n_trials, filename):
    """
    This function provides a convenient means of doing a typical 
    fault-tolerant simulation using the square-octagon code with noisy
    syndrome measurements. 
    """
    
    sim_lattice = SquareOctagonLattice((size,size))
    sim_dual_lattice_list = []
    for idx in xrange(size):
        sim_dual_lattice_list.append(UnionJackLattice((size,size), is_dual=True))
    sim_model = depolarizing_model(error_rate)
    sim_code_func = noisy_squoct_code
    sim_decoder = ft_mwpm_decoder(sim_lattice, sim_dual_lattice_list)
    sim_log_ops = squoct_log_ops(sim_lattice.total_size)

    sim_keys = ['lattice', 'dual_lattice_list', 'error_model', 
                'synd_noise', 'code_func', 'decoder', 
                'logical_operators', 'n_trials']

    sim_values = [sim_lattice, sim_dual_lattice_list, sim_model, 
                    error_rate, sim_code_func, 
                    sim_decoder, sim_log_ops, n_trials]
    
    sim_dict = dict(zip(sim_keys, sim_values))

    sim = FTSimulation(**sim_dict)
    sim.run()
    sim.save(filename + '.sim')

