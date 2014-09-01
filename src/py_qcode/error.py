from numpy.random import rand
from qecc import Pauli

## ALL ##
__all__ = ['ErrorModel', 'PauliErrorModel', 'depolarizing_model', 'iidxz_model']

## CONSTANTS ##
PAULIS = ['I', 'X', 'Y', 'Z']
ACCEPTABLE_OPERATORS = PAULIS + ['H', 'P']

class ErrorModel(object):
    """
    Wraps a list of tuples corresponding to a discrete probability and 
    an operator. This assumes independent identically-distributed 
    noise, though not necessarily Pauli.

    :param prob_op_list: A list of probabilites and associated 
    operators. The probabilites are floats which must sum to 1 to 
    within :math:`10^{-12}`. The operators are represented by strings
    which must be drawn from the list of acceptable operators 
    `['I','X','Y','Z','H','P']` Each pair (a probability and its 
    associated operator) is stored in a tuple.

    :type prob_op_list: list
    """
    #Magic Methods
    def __init__(self, prob_op_list):
        #Sanitize input
        probs, ops = zip(*prob_op_list)
        
        if abs(1. - sum(probs)) > 10. ** -12:
            raise ValueError('The probabilites of different errors'+\
                ' must sum to 1, the probabilites entered here sum'+\
                 ' to {0}'.format(sum(probs)))
        if any([op not in ACCEPTABLE_OPERATORS for op in ops]):
            raise ValueError('Received operator outside set of acceptable operators {0}. You entered: {1}'\
                .format(ACCEPTABLE_OPERATORS, ops))

        self.prob_op_list = prob_op_list
    
    def __repr__(self):
        return ', '.join(['{1} with probability {0}'.format(
                *prob_op)
                for prob_op in self.prob_op_list])

    #Other Methods
    def act_on(self, lattice):
        for point in lattice.points:
            point.error = _action(self.prob_op_list, rand())

class PauliErrorModel(ErrorModel):
    def __init__(self, prob_op_list):

        for prob_op in prob_op_list:
            if prob_op[1] not in PAULIS:
                raise ValueError('Received operator {0} outside of Pauli set {1}'\
                    .format(prob_op[1], PAULIS))

        super(PauliErrorModel, self).__init__(prob_op_list)

    def act_on(self, lattice):
        for point in lattice.points:
            new_pauli = Pauli(_action(self.prob_op_list, rand()))
            if point.error is None:
                point.error = new_pauli
            else:
                point.error *= new_pauli

#Convenience functions

def depolarizing_model(p):
    """
    The depolarizing model applies the identity with probability 
    :math:`1-p`, and each of the single qubit Pauli operators 
    :math:`X`, :math:`Y`, and :math:`Z` with probability 
    :math:`\dfrac{p}{3}`. 
    """
    return PauliErrorModel([(1. - p, 'I'), (p / 3., 'X'), (p / 3., 'Y'), (p / 3., 'Z')])

def iidxz_model(px, pz=None):
    """
    The independent identically-distributed X/Z model applies a bit 
    and phase flip to each site, resulting in a reduced probability of
    Y errors.
    """
    if pz is None:
        pz = px

    return PauliErrorModel([((1. - px) * (1. - pz), 'I'), (px * (1. - pz), 'X'),
                        ((1. - px) * pz, 'Z'), (px * pz, 'Y')])

rolling_sum = lambda lst: [sum(lst[:idx+1]) for idx in range(len(lst))]

def _action(prob_op_list, sample):
    """
    returns the operator from a prob_op_list corresponding to a number between 0 and 1. 
    """
    if (sample < 0.) or (sample > 1.):
        raise ValueError("`sample` must be between 0 and 1, preferably uniform.")

    probs, ops = zip(*prob_op_list)
    cum_probs = rolling_sum(probs)
    for idx in range(len(cum_probs)):
        if sample < cum_probs[idx]:
            return ops[idx]
