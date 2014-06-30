import cPickle as pkl
from collections import Iterable

__all__ = ['Simulation']

class Simulation():     
    """
    `Simulation` is the top-level class
    for py_qcode, the user is meant to set up, execute and save results
    using these objects.

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
    #Magic Methods
    def __init__(self, lattice, dual_lattice, error_model, code, 
                                decoder, logical_operators, n_trials):

        #Defined objects
        self.lattice = lattice
        self.dual_lattice = dual_lattice
        self.error_model = error_model
        self.code = code
        self.decoder = decoder
        
        if not isinstance(logical_operators, Iterable):
            self.logical_operators = [logical_operators]
        else:
            self.logical_operators = logical_operators

        #Integer
        self.n_trials = n_trials
        
        #Final Values
        self.logical_error = None

    def run(self):
        """
        The main routine in this library, follows the recipe `n_trials` times in series:

        + Apply the error model to the primary lattice, assigning values to the `error` attributes of the :class:`py_qcode.Point` objects within. 

        + Obtain the true coset of the error with respect to the :class:`py_qcode.ErrorCorrectingCode` being used.

        + Perform a measurement, using the attributes of the error-correcting code to generate syndromes on the dual lattice.

        + Infer the coset of the error by acting the decoder on the dual lattice. 

        + Record both cosets. 
        """
        self.logical_error = []
        for idx in range(self.n_trials):
            self.error_model.act_on(self.lattice)
            self.code.measure()
            self.decoder.infer()

            #Error checking, if the resulting Pauli is not in the 
            #normalizer, chuck an error:
            self.code.measure()
            for point in self.dual_lattice.points:
                if point.syndrome is not None:
                    raise ValueError('Product of "inferred error"'+\
                        ' with actual error anticommutes with some'+\
                        ' stabilizers.')

            for operator in self.logical_operators:
                com_relation_list = []
                com_relation_list.append(operator.test(self.lattice))
            self.logical_error.append(com_relation_list)
    
    def save(self, filename):
        
        big_dict =  {}
        big_dict['lattice_class'] = \
            str(self.lattice.__class__).split('.')[-1][:-2]
        big_dict['lattice_size'] = self.lattice.size
        big_dict['dual_lattice_class'] = \
            str(self.dual_lattice.__class__).split('.')[-1][:-2]
        big_dict['dual_lattice_size'] = self.dual_lattice.size
        big_dict['error_model'] = repr(self.error_model)
        big_dict['code'] = self.code.name
        big_dict['decoder'] = self.decoder.name
        big_dict['n_trials'] = self.n_trials
        big_dict['logical_errors'] = self.logical_error
        
        with open(filename,'w') as phil:
            pkl.dump(big_dict, phil)