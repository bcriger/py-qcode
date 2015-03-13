import cPickle as pkl
from collections import Iterable
from qecc import I
# from utils import syndrome_print, error_print

__all__ = ['Simulation', 'FTSimulation']


class Simulation():

    """
    `Simulation` is the top-level class
    for py_qcode, the user is meant to set up, execute and save results
    using these objects.

    :param lattice: a set of points, identified with 
    :math:`n`-dimensional co-ordinates, each corresponding to a 
    physical qubit.

    :type lattice: :class:`py_qcode.Lattice`

    :param dual_lattice: a second set of :math:`n`-dimensional points,
    corresponding to check operators.

    :type dual_lattice: :class:`py_qcode.Lattice`

    :param error_model: A description of the errors to be applied 
    independently to the qubits of `lattice`.

    :type error_model: :class:`py_qcode:ErrorModel`

    :param code: An error-correcting code which translates errors on a 
    lattice to syndromes on the dual lattice.

    :type code: :class:`py_qcode.ErrorCorrectingCode`

    :param decoder: A protocol for inferring errors given syndromes.

    :type decoder: :class:`py_qcode.Decoder`

    :param n_trials: a number of simulations to be performed in series.
    This can be used to organize batch jobs so that one can submit more
    than one simulation per job.

    :type n_trials: integer

    :param logical_error: Commutation relations between the logical 
    operators and the product of the guessed error and the actual 
    error.

    :type logical_error: list
    """
    # Magic Methods
    def __init__(self, lattice, dual_lattice, error_model, code,
                 decoder, logical_operators, n_trials):

        # Defined objects
        self.lattice = lattice
        self.dual_lattice = dual_lattice
        self.error_model = error_model
        self.code = code
        self.decoder = decoder

        if not isinstance(logical_operators, Iterable):
            self.logical_operators = [logical_operators]
        else:
            self.logical_operators = logical_operators

        # Integer
        self.n_trials = n_trials

        # Final Values
        self.logical_error = None

    def run(self):
        """
        The main routine in this library, follows the recipe `n_trials` 
        times in series:

        + Apply the error model to the primary lattice, assigning 
          values to the `error` attributes of the 
          :class:`py_qcode.Point` objects within.

        + Obtain the true coset of the error with respect to the 
          :class:`py_qcode.ErrorCorrectingCode` being used.

        + Perform a measurement, using the attributes of the 
          error-correcting code to generate syndromes on the dual 
          lattice.

        + Infer the error by acting the decoder on the dual lattice, 
          applying the resulting operator to the primary lattice.

        + Record the commutation relations of the resulting operator
          with the logical operators.
        """
        self.logical_error = []
        for idx in range(self.n_trials):
            # Clean up results from previous simulation
            self.lattice.clear()
            self.dual_lattice.clear()

            # The bulk of the work
            self.error_model.act_on(self.lattice)
            self.code.measure()
            self.decoder.infer()

            # Error checking, if the resulting Pauli is not in the
            # normalizer, chuck an error:

            self.dual_lattice.clear()
            # syndrome_print(self.dual_lattice)
            self.code.measure()
            # syndrome_print(self.dual_lattice)

            for point in self.dual_lattice.points:
                if point.syndrome:
                    raise ValueError('Product of "inferred error"'
                                     ' with actual error anticommutes'
                                     'with some stabilizers.')

            com_relation_list = []
            for operator in self.logical_operators:
                # print operator
                com_relation_list.append(operator.test(self.lattice))
            self.logical_error.append(com_relation_list)

        # Clean up results from final simulation
        self.lattice.clear()
        self.dual_lattice.clear()

    def save(self, filename):
        big_dict = {}
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

        with open(filename, 'w') as phil:
            pkl.dump(big_dict, phil)


class FTSimulation():

    """
    `FTSimulation` is the class for representing simulations of 
    fault-tolerant error-correction protocols.
    """
    #Magic Methods
    def __init__(self, lattice, dual_lattice_list, error_model, 
                 synd_noise, code_func, last_code_func, decoder,
                 logical_operators, n_trials):

        # Defined objects
        self.lattice = lattice
        self.dual_lattice_list = dual_lattice_list
        self.error_model = error_model
        self.code_func = code_func
        self.last_code_func = last_code_func
        self.decoder = decoder
        self.synd_noise = synd_noise

        if not isinstance(logical_operators, Iterable):
            self.logical_operators = [logical_operators]
        else:
            self.logical_operators = logical_operators

        # Integer
        self.n_trials = n_trials

        # Final Values
        self.logical_error = None

    def run(self):
        """
        The main routine in this library, follows the recipe 
        `n_trials` times in series:

        + For every dual lattice provided:
          - Apply the error map to the primary lattice
          - perform a measurement onto the specified dual lattice.
        + Perform a fault-tolerant decoding, applying the resulting
          operator to the primary lattice.
        + Measure commutation relations with the logical operators
        + Record the results.
        """
        self.logical_error = []
        for idx in range(self.n_trials):
            # Clean up results from previous simulation
            self.lattice.clear()
            for dual_lattice in self.dual_lattice_list:
                dual_lattice.clear()

            # The bulk of the work
            #self.error_model.act_on(self.lattice)
            #Begin in code state
            for point in self.lattice.points:
                point.error = I
            
            for dual_lattice in self.dual_lattice_list[:-1]:
                #No new memory errors, just errors from code 
                #back-action
                #self.error_model.act_on(self.lattice)
                
                # New code object created for every iteration:
                current_code = self.code_func(dual_lattice)
                current_code.measure()
            
            #In order to guarantee that the resulting lattice operator 
            #is in the normalizer, we assume that the final round of 
            #error correction is perfect. Theoretically, errors hidden
            #by this round survive to the next. In practice, thus far,
            #These errors are cleared.  
            noiseless_code = self.last_code_func(
                                            self.dual_lattice_list[-1])

            noiseless_code.measure()

            self.decoder.infer()

            # Error checking, if the resulting Pauli is not in the
            # normalizer, chuck an error:

            self.dual_lattice_list[-1].clear()
            # syndrome_print(self.dual_lattice_list[-1])
            noiseless_code.measure()
            for point in self.dual_lattice_list[-1].points:
                if point.syndrome:
                    raise ValueError('Product of "inferred error"' 
                                     ' with actual error anticommutes'
                                     'with some stabilizers.')

            com_relation_list = []
            for operator in self.logical_operators:
                # print operator
                com_relation_list.append(operator.test(self.lattice))
            self.logical_error.append(com_relation_list)

        # Clean up results from final simulation
        self.lattice.clear()
        for dual_lattice in self.dual_lattice_list:
            dual_lattice.clear()

    def save(self, filename):
        big_dict = {}
        big_dict['lattice_class'] = \
            str(self.lattice.__class__).split('.')[-1][:-2]
        big_dict['lattice_size'] = self.lattice.size
        big_dict['dual_lattice_class'] = \
            str(self.dual_lattice_list[0].__class__).split('.')[-1][:-2]
        big_dict['dual_lattice_size'] = self.dual_lattice_list[0].size
        big_dict['n_measurements'] = len(self.dual_lattice_list)
        big_dict['error_model'] = repr(self.error_model)
        big_dict['code'] = self.code_func.func_name
        big_dict['decoder'] = self.decoder.name
        big_dict['n_trials'] = self.n_trials
        big_dict['logical_errors'] = self.logical_error

        with open(filename, 'w') as phil:
            pkl.dump(big_dict, phil)
