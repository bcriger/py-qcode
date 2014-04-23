
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
		