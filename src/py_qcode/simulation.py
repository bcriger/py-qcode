import cPickle as pkl

__all__ = ['Simulation']

class Simulation():     
    """
    `Simulation` is the top-level class
    for py_qcode, the user is meant to set up, execute and save results
    using these objects.

    :param lattice: A grid of points to contain errors and syndromes. 

    :type lattice: :class:`py_qcode.Lattice`, input

    :param error_model: A description of the errors to be applied independently to the qubits of `lattice`.

    :type error_model: :class:`py_qcode:ErrorModel`, input 

    :param code: An error-correcting code which translates errors on a lattice to syndromes on the dual lattice.

    :type code: :class:`py_qcode.ErrorCorrectingCode`, input

    :param decoder: A protocol for inferring errors given syndromes.

    :type decoder: :class:`py_qcode.Decoder`, input

    :param n_trials: a number of simulations to be performed in series. This can be used to organize batch jobs so that one can submit more than one simulation per job.

    :type n_trials: integer, input

    :param true_coset: The actual coset to which the random error applied in the simulation belongs.

    :type true_coset: str, output

    :param inferred_coset: The coset assigned by the error-correcting code during the simulation. 
    """
    def __init__(self, lattice, error_model, code, decoder, n_trials):
        
        #Initial Values
        self.lattice = lattice
        self.error_model = error_model
        self.code = code
        self.decoder = decoder
        self.n_trials = n_trials
        
        #Final Values
        self.true_coset = None
        self.inferred_coset = None

    def run(self):
        """
        The main routine in this library, follows the recipe `n_trials` times in series:

        + Apply the error model to the lattice, assigning values to the `error` attributes of the :class:`py_qcode.Point` objects within. 

        + Obtain the true coset of the error with respect to the :class:`py_qcode.ErrorCorrectingCode` being used.

        + Perform a measurement, using the attributes of the error-correcting code to generate syndromes on the dual lattice.

        + Infer the coset of the error by acting the decoder on the dual lattice. 

        + Record both cosets
        """
        for idx in range(n_trials):
            continue
        pass
    
    def save(self):
        pass
#Convenience Functions
def sim_from_file(filename):
    """
    The purpose of this function is to:

    + open a file containing a pickled dictionary of input values to a simulation,

    + initialize the objects which the corresponding `py_qcode.Simulation` takes as input,
    
    + run the simulation, and 

    + save the results to a file of the same name as the input, with a different extension.  
    """
    pass