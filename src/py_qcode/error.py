from numpy.random import rand
from qecc import Pauli, pauli_group
import qecc as q
#from lattice import Lattice, Point #?
from collections import Iterable

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
        
        test_len = len(ops[0])
        for string in ops:
            _whitelist_string(ACCEPTABLE_OPERATORS, string)
            if len(string) != test_len:
                raise ValueError("All errors must act on registers "+\
                    ("of the same size, operators {0} and {1} act on"+\
                        " registers of different size")\
                            .format(ops[0], string))

        self.prob_op_list = prob_op_list
    
    def __repr__(self):
        return ', '.join(['{1} with probability {0}'.format(*prob_op)
                            for prob_op in self.prob_op_list])

    #Other Methods
    def act_on(self, register):
        """
        Acts an error model on a register. First detects whether the 
        register in question is a py_qcode.Lattice and the model is 
        single-qubit, or the model has operations the same size as the
        register. 
        """
        if isinstance(register, Lattice):

            #Check that model has single-qubit errors
            for op in zip(*self.prob_op_list)[1]:
                if len(op) != 1:
                    raise ValueError("Length of operators applied "+\
                        "to lattice must be 1, operator "+\
                        "{0} has length {1}.".format(op, len(op)))
            
            for point in register.points:
                #FIXME: Needs to update the error, not assign a new value
                # see PauliErrorModel.act_on
                point.error = _action(self.prob_op_list, rand())
        
        elif isinstance(register, Iterable):
            #TODO: provide support for multi-qubit non-pauli errors
            #TODO: check that all elements of iterable are Points.
            #TODO: check weight of errors == length of register
            pass

        else:
            raise ValueError("ErrorModel objects must act on a "+\
                "Lattice object or an iterable of Point objects. "+\
                "You entered: {}".format(register))

class PauliErrorModel(ErrorModel):
    def __init__(self, prob_op_list):

        for string in zip(*prob_op_list)[1]:
            _whitelist_string(PAULIS, string)

        super(PauliErrorModel, self).__init__(prob_op_list)

    def act_on(self, register):
        
        ops = zip(*self.prob_op_list)[1]
        
        if isinstance(register, Lattice):
            
            for op in ops:
                if len(op) != 1:
                    raise ValueError("Only weight-1 Paulis may be "+\
                                        "used on whole Lattices")
            
            for point in register.points:
                new_pauli = Pauli(_action(self.prob_op_list, rand()))
                if point.error is None:
                    point.error = new_pauli
                else:
                    point.error *= new_pauli
        
        elif isinstance(register, Iterable):
            #Test register to see that it contains points
            for point in register:
                if not isinstance(point, Point):
                    raise ValueError("Registers which are not "+\
                        "entire lattices must be iterables of "+\
                        "points, this input contains non-points.")

            #Test model to see if ops are same len as register
            for op in ops:
                if len(op) != len(register):
                    raise ValueError("Pauli must be same length "+\
                        "as register")

            for pt in register:
                if pt.error is None:
                    pt.error = I
            
            error = Pauli("".join(pt.error for pt in register))
            error = Pauli(_action(self.prob_op_list , rand())) * error

            for idx, pt in enumerate(register):
                pt.error = error[idx]
        
        else:
            raise ValueError("register must be either Lattice or "+\
                                "iterable of points")

        pass

    def compress(self):
        """
        If `self` has non-unique operators, this method sums the 
        probabilities of those operators and assigns the sum to the 
        operator.  
        """
        probs, ops = zip(*self.prob_op_list)
        unique_ops = list(set(ops))
        unique_probs = [0 for elem in unique_ops]
        for idx, unique_op in enumerate(unique_ops):
            unique_probs[idx] = sum(probs[jdx] 
                                    for jdx in range(len(ops))
                                        if ops[jdx] == unique_op) 

        self.prob_op_list = zip(unique_probs, unique_ops)
    
    def __mul__(self, other):
        
        if not isinstance(other, PauliErrorModel):
            raise ValueError("PauliErrorModel instances can only be "+\
                "multiplied with each other.")
        
        new_prob_ops=[]
        for self_p, self_o in self.prob_op_list:
            for oth_p, oth_o in other.prob_op_list:
                new_prob_ops.append((self_p * oth_p, 
                                (q.Pauli(self_o) * q.Pauli(oth_o)).op))

        new_model = PauliErrorModel(new_prob_ops)
        new_model.compress()
        return new_model 

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

def two_bit_twirl(p):
    """
    With probability p, selects a Pauli at random from the 
    non-identity Paulis on two qubits.
    """
    prob_op_list = [(1. - p, 'II')]
    for op in map(lambda a: a.op, list(pauli_group(2)))[1:]:
        prob_op_list.append((p / 15., op))
    
    return PauliErrorModel(prob_op_list)

def fowler_meas_model(stab_type, nq, p):
    """
    Given a description of a circuit which measures a stabilizer of a 
    CSS code and an error probability, returns the error model on 
    `nq + 1` bits that results from conjugating all the wait location 
    faults and CNOT faults as well as a measurement fault and a state 
    preparation fault, which are determined by the stabilizer type.  
    """
    ident = q.Pauli('I' * (nq + 1))
    meas_circ = css_meas_circuit(stab_type, nq)
    err_dict = {}

    #State preparation flip
    flip_type = 'X' if stab_type == 'Z' else 'Z'
    
    prep_fault = q.propagate_fault(meas_circ, 
                    q.Pauli.from_sparse({nq: flip_type}, nq = nq + 1))
    
    err_dict_update(err_dict, prep_fault, p)
    #err_dict_update(err_dict, ident, (1. - p))
    
    #Faults resulting from wait locations and CNOTs, provided by QuaEC:
    for idx, sub_circuit in enumerate(meas_circ):
        faults = q.possible_faults(sub_circuit)
        for fault in faults:
            if fault.wt == 0:
                #identity
                pass
                #err_dict_update(err_dict, fault, (1. - p))
            elif fault.wt == 1:
                end_fault = q.propagate_fault(meas_circ[idx + 1:], fault)
                err_dict_update(err_dict, end_fault, p / (3.))
            elif fault.wt == 2:
                end_fault = q.propagate_fault(meas_circ[idx + 1:], fault)
                err_dict_update(err_dict, end_fault, p / (15.))
            else:
                raise ValueError("Something weird has happened.")

    #measurement fault
    meas_fault = q.Pauli.from_sparse({nq: flip_type}, nq = nq + 1)
    err_dict_update(err_dict, meas_fault, p)
    #err_dict_update(err_dict, ident, (1. - p) )

    return [(tpl[1], tpl[0]) for tpl in err_dict.items()]

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

def _whitelist_string(whitelist, string):
    """
    Common pattern throughout this file, test a string to see if all 
    its letters are members of a whitelist. 
    """
    if not all(letter in whitelist for letter in string):
        raise ValueError(("Input string {0} contains letter {1} "+\
                            "not found in whitelist {2}")\
                                .format(string, letter, whitelist))
    pass

def css_meas_circuit(stab_type, nq):
    """
    Produces a circuit which measures a CSS stabilizer of the required 
    type, acting on a given number of bits. 
    """
    
    if not isinstance(nq, int):
        raise ValueError("nq must be integer")
    elif nq < 0:
        raise ValueError("nq must be positive")
    
    if stab_type == 'Z':
        cnot_pairs = [(idx, nq) for idx in range(nq)]
    elif stab_type == 'X':
        cnot_pairs = [(nq, idx) for idx in range(nq)]
    else:
        raise ValueError("stab_type must be X or Z")

    pr_to_cnot = lambda pr: ('CNOT', pr[0], pr[1])

    circ = q.Circuit(*map(pr_to_cnot, cnot_pairs))
    
    return list(circ.group_by_time(pad_with_waits=True))

def err_dict_update(err_dict, pauli, prob):
    """
    Takes a dictionary of the form {qecc.Pauli: float}, and performs 
    an update. To do this, it checks whether the input pauli is in the
    dictionary already. If so, it adds `prob` to the value of `pauli`,
    if not it adds `pauli` to the dictionary with value `prob`.
    """

    if pauli.op in err_dict.keys():
        err_dict[pauli.op] += prob
    else:
        err_dict[pauli.op] = prob
    
    return err_dict