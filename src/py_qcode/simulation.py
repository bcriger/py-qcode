import cPickle as pkl

__all__ = ['Simulation']

class Simulation():     
    """
    `Simulation` is the top-level class
    for py_qcode, the user is meant to set up, execute and save results
    using these objects.

    :param lattice: A grid of points to contain errors and syndromes. 

    :type lattice: :class:`py_qcode.Lattice`

    :param error_model: A description of the errors to be applied independently to the qubits of `lattice`.

    :type error_model: :class:`py_qcode:ErrorModel` 

    :param code: An error-correcting code which translates errors on a lattice to syndromes on the dual lattice.

    :type code: :class:`py_qcode.ErrorCorrectingCode`

    :param decoder: A protocol for inferring errors given syndromes.

    :type decoder: :class:`py_qcode.Decoder`

    :param n_trials: a number of simulations to be performed in series. This can be used to organize batch jobs so that one can submit more than one simulation per job.

    :type n_trials: integer
    """
    def __init__(self, n_trials):
        self.arg = arg
    def run(self):
        """

        """

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