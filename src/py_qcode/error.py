## ALL ##
__all__ = ['ErrorModel', 'BoundingErrorModel']

## CONSTANTS ##
ACCEPTABLE_OPERATORS = ['I','X','Y','Z','H','P']

class ErrorModel():
    """
    Wraps a list of tuples corresponding to a discrete probability and 
    an operator. This assumes independent identically-distributed 
    noise, though not necessarily Pauli.
    """

    def __init__(self, prob_op_list):
        #Sanitize input
        self.prob_op_list = prob_op_list
        
class BoundingErrorModel():
    """
    Wraps a pair of ErrorModel objects, one of which minimally
    overestimates the fidelity of a 'weird' map (think non-unital,
    non-Pauli, etc.), while the other minimally underestimates the
    fidelity.
    """

    def __init__(self, upper, lower):
        #Sanitize input
        self.upper = upper
        self.lower = lower
