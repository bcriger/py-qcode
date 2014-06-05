from numpy.random import rand

## ALL ##
__all__ = ['ErrorModel', 'DepolarizingModel']

## CONSTANTS ##
ACCEPTABLE_OPERATORS = ['I','X','Y','Z','H','P']

class ErrorModel(object):
    """
    Wraps a list of tuples corresponding to a discrete probability and an operator. This assumes independent identically-distributed noise, though not necessarily Pauli.

    :param prob_op_list: A list of probabilites and associated operators. The probabilites are floats which must sum to 1 to within :math:`10^{-12}`. The operators are represented by strings which must be drawn from the list of acceptable operators `['I','X','Y','Z','H','P']` Each pair (a probability and its associated operator) is stored in a tuple.

    :type prob_op_list: list
    """

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
    
    def act_on(self, lattice):
        for point in lattice.points:
            point.error = action(self.prob_op_list, rand())

class DepolarizingModel(ErrorModel):
    """
    The depolarizing model applies the identity with probability :math:`1-p`, and each of the single qubit Pauli operators :math:`X`, :math:`Y`, and :math:`Z` with probability :math:`\dfrac{p}{3}`. 

    :param p: A probability, between 0 and 1.

    :type p: float
    """
    def __init__(self, p):
        super(DepolarizingModel, self).__init__([(1. - p, 'I'),
                                                 (p / 3., 'X'),
                                                 (p / 3., 'Y'),
                                                 (p / 3., 'Z')])

#Convenience functions

rolling_sum = lambda lst: [sum(lst[:idx+1]) for idx in range(len(lst))]

def action(prob_op_list, sample):
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
