from numpy.random import rand
import numpy as np
from qecc import Pauli, pauli_group, I, X, Y, Z
from lattice import Lattice, Point #?
from collections import Iterable
from math import fsum, log
from scipy.weave import inline

## ALL ##
__all__ = ['ErrorModel', 'PauliErrorModel', 'depolarizing_model',
           'iidxz_model', 'DensePauliErrorModel']

## TEMPORARY ADDITIONS TO ALL ##
__all__.extend(['ensure_probabilities', 'mask_from_bits', 'num_ys', 
                'hamming_weight', 'weight_from_idx', 'two_bit_twirl',
                '_action'])

## CONSTANTS ##
PAULIS = ['I', 'X', 'Y', 'Z']
ACCEPTABLE_OPERATORS = PAULIS + ['H', 'P']
float_type = np.float64

class ErrorModel(dict):

    """
    Dictionary with operator objects as keys and floating- point
    probabilites as values.

    :param op_pr_dict: A dictionary. The probabilites are floats which
    must sum to 1 to within :math:`10^{-12}`. The operators are
    represented by arbitrary objects (so that either strings or
    :class:`qecc.Pauli` objects can be used).

    :type op_pr_dict: dict
    """
    # Magic Methods
    def __init__(self, op_pr_dict={}):
        # Sanitize optional input
        if op_pr_dict:
            ops, probs = zip(*op_pr_dict.items())

            if not all([0. <= p <= 1. for p in probs]):
                raise ValueError("All input probabilites must be " +
                                 "between 0 and 1, you entered: {}".format(probs))

            if abs(1. - fsum(probs)) > 10. ** -12:
                raise ValueError('The probabilites of different errors' +
                                 ' must sum to 1, the probabilites entered here sum' +
                                 ' to {0}'.format(sum(probs)))

            test_len = len(ops[0])
            for operator in ops:
                # FIXME: Find out what methods an object must support in
                # order to be a valid key in an error model dict and
                # sanitize accordingly:
                #_whitelist_string(ACCEPTABLE_OPERATORS, string)
                if len(operator) != test_len:
                    raise ValueError("All errors must act on registers " +
                                    ("of the same size, operators {0} and {1} act on" +
                                     " registers of different size")
                                     .format(ops[0], operator))

        # Initialize a dict
        super(ErrorModel, self).__init__()
        self.update(op_pr_dict)

    # Convenience Properties
    @property
    def ops(self):
        return self.keys()

    @property
    def probs(self):
        return self.values()

    @property
    def o_p_pairs(self):
        return self.items()

    # Other Methods
    def act_on(self, register):
        """
        Acts an error model on a register. First detects whether the
        register in question is a py_qcode.Lattice and the model is
        single-qubit, or the model has operations the same size as the
        register.
        """
        if hasattr(register, 'points'):

            # Check that model has single-qubit errors
            for op in self.ops:
                if len(op) != 1:
                    raise ValueError("Length of operators applied " +
                                     "to lattice must be 1, operator " +
                                     "{} has length {}.".format(op, len(op)))

            for point in register.points:
                # FIXME: Needs to update the error, not assign a new value
                # see PauliErrorModel.act_on
                point.error = _action(self, rand())

        elif isinstance(register, Iterable):
            # TODO: provide support for multi-qubit non-pauli errors
            # TODO: check that all elements of iterable are Points.
            # TODO: check weight of errors == length of register
            _point_iter_apply(self, register)

        else:
            raise ValueError("ErrorModel objects must act on a "
                             "Lattice object or an iterable of Point"
                             " objects. You entered: "
                             "{}".format(register))


class PauliErrorModel(ErrorModel):

    """
    subclass of ErrorModel, restricts the keys of the optional input
    to be :class:`qecc.Pauli` objects
    """
    def __init__(self, op_pr_dict={}):
        """
        Sanitizes the optional input dictionary, and initializes the
        superclass object.
        """
        if op_pr_dict:
            for pauli in op_pr_dict.keys():
                if not isinstance(pauli, Pauli):
                    raise ValueError("Input keys to PauliErrorModel " +
                                     "must be Paulis.")

        super(PauliErrorModel, self).__init__(op_pr_dict)

    def __len__(self):
        return len(self.ops[0])

    def act_on(self, register):
        """
        Ungodly version of act_on, distinguishes between three cases.
        Either:
         + We apply a single-qubit map to the entire lattice, which is
           passed in as the register.
         + We apply a single-qubit map to a subset of the lattice, 
           which is passed in as a list or other iterable.
         + We apply a multi-qubit map to subsets of the lattice; the 
           register is an iterable of iterables.
        """
        if hasattr(register, 'points'):
            _full_lattice_apply(self, register)
        
        elif isinstance(register, Iterable):
            if isinstance(register[0], Iterable):    
                for pt_set in register:
                    mul_error = _action(self, rand())
                    for pdx, pt in enumerate(pt_set):
                        if pt.error is None:
                            pt.error = I
                        pt.error *= mul_error[pdx]
            else:
                for pt in register:
                    if pt.error is None:
                        pt.error = Pauli('I')
                    pt.error *= _action(self, rand())    
        
        elif isinstance(register, Point):
            if pt.error is None:
                pt.error = Pauli('I')
            pt.error *= _action(self, rand())
        
        else:
            raise ValueError("Could not determine how to act error "
                "model {} on register {}.".format(self, register))


class DensePauliErrorModel(object):

    """
    When an error model acts on a large number of qubits, and has a
    non-negligible probability of producing any Pauli in the group,
    dictionaries are an extremely slow way to store and manipulate the
    error model. Here, I store an array, called `vec`, whose indices
    correspond to Paulis by the binary symplectic mapping
    (where [001101] corresponds to ZIY, etc.), and each value is a
    probability.
    """
    def __init__(self, optional_vec=None, nq=1):

        if optional_vec is not None:
            # ensure 1d array of probabilities
            optional_vec = sanitize_vector(optional_vec)
            ensure_probabilities(optional_vec)
            self.vec = optional_vec
            self.nq = int(log(len(optional_vec), 4))
        else:
            self.vec = np.zeros((4 ** nq,), dtype=float_type)
            self.vec[0] = 1.
            self.nq = nq

    def sample(self):
        """
        Lazily selects a Pauli at random from the vector.
        """
        prob = rand()
        idx = 0
        for elem in np.nditer(self.vec):
            if prob < elem:
                return pauli_from_int(idx, self.nq)
            else:
                prob -= elem
                idx += 1

    def __mul__(self, other):

        if not isinstance(other, DensePauliErrorModel):
            try:
                other = DensePauliErrorModel(other)
            except Exception, e:
                raise TypeError("You can only multiply DensePauliErrorModels" +
                                " by other DensePauliErrorModels")
        if self.nq != other.nq:
            raise ValueError("self and other must have the same " +
                             "number of qubits")

        new_vec = np.empty((4 ** self.nq,), dtype=float_type)
        c_code = """
        int s_idx, out_idx;
        double temp_p;

        for (int out_idx = 0; out_idx < vec_len; out_idx++)
        {
            temp_p = 0;
            for (int s_idx = 0; s_idx < vec_len; s_idx++)
            {
                temp_p += s_vec[s_idx] * o_vec[s_idx ^ out_idx];
            }
            new_vec[out_idx] = temp_p;
        }

        """
        vec_len = 4 ** self.nq
        s_vec = self.vec
        o_vec = other.vec
        arg_names = ['vec_len', 's_vec', 'o_vec', 'new_vec']

        inline(c_code, arg_names=arg_names, compiler='gcc', 
            extra_compile_args=['-O3'], libraries=['rt'])

        return DensePauliErrorModel(new_vec)

    def __pow__(self, exponent):
        """
        Exponentiates an error model. model**n is equivalent to acting 
        a model n times in a row.
        """
        if not isinstance(exponent, int):
            raise TypeError("Exponent must be int")

        if exponent == 0:
            return DensePauliErrorModel(nq=self.nq)
        elif exponent == 1:
            return self
        else:
            model_copy = self
            for idx in range(exponent - 1):
                model_copy *= self

            return model_copy

    def cnot(self, ctrl, targ):
        """
        The action of any Clifford on a Pauli probability vector can be
        expressed as a permutation. This function permutes the elements
        of the vector according to a CNOT acting on two qubits.
        """

        nq = self.nq

        if any([ctrl > nq - 1, targ > nq - 1, ctrl < 0,
                targ < 0, targ == ctrl]):
            raise ValueError("ctrl and/or targ unsuitable: " +
                             "{} {}".format(ctrl, targ))

        for idx in xrange(4 ** nq):
            # split index into x and z portions
            x_int, z_int = xz_split(idx, nq)

            # act cnot on each portion
            new_x_int = c_xor_int(x_int, ctrl, targ, nq)
            new_z_int = c_xor_int(z_int, targ, ctrl, nq)

            # rejoin and assign
            if (new_x_int != x_int) or (new_z_int != z_int):
                new_idx = xz_join(new_x_int, new_z_int, nq)
                self.vec[new_idx], self.vec[idx] = \
                    self.vec[idx], self.vec[new_idx]

        pass

    def depolarize(self, q, p):
        """
        For something a little faster than __mul__, I hard-code
        multiplication by a depolarizing model on a select qubit.
        """
        nq = self.nq
        new_vec = np.empty(len(self.vec), dtype=float_type)
        c_code = '''
        int mask_01 = 1 << (nq - q - 1);
        int mask_10 = 1 << (2 * nq - q - 1);
        int mask_11 = mask_10 + mask_01;
        for (int idx = 0; idx < vec_len; ++idx)
        {
            new_vec[idx] = (1. - double(p)) * s_vec[idx] + \
                           double(p) / 3. * (s_vec[idx ^ mask_01] +
                                     s_vec[idx ^ mask_10] +
                                     s_vec[idx ^ mask_11]);
        }
        '''
        s_vec = self.vec
        vec_len = len(self.vec)
        arg_names = ['s_vec', 'new_vec', 'vec_len', 'nq', 'p', 'q']

        inline(c_code, arg_names=arg_names, compiler='gcc', 
            extra_compile_args=[''], libraries=['rt'])
        #self = DensePauliErrorModel(new_vec) 
        self.vec = new_vec

        pass

    def twirl(self, q1, q2, p):
        """
        For something a little faster than __mul__, I hard-code
        multiplication by a full twirl on two qubits.
        """
        nq = self.nq
        new_vec = np.zeros((len(self.vec),), dtype=float_type)
        c_code = '''
        int mask_0001 = 1 << (nq - q2 - 1);
        int mask_0010 = 1 << (nq - q1 - 1);
        int mask_0100 = 1 << (2 * nq - q2 - 1);
        int mask_1000 = 1 << (2 * nq - q1 - 1);

        int masks[15] = {mask_0001, mask_0010, mask_0010 + mask_0001,
                         mask_0100, mask_0100 + mask_0001,
                         mask_0100 + mask_0010,
                         mask_0100 + mask_0010 + mask_0001, mask_1000,
                         mask_1000 + mask_0001,
                         mask_1000 + mask_0010,
                         mask_1000 + mask_0010 + mask_0001,
                         mask_1000 + mask_0100,
                         mask_1000 + mask_0100 + mask_0001,
                         mask_1000 + mask_0100 + mask_0010,
                         mask_1000 + mask_0100 + mask_0010 + mask_0001};

        for (int idx = 0; idx < vec_len; ++idx)
        {
            new_vec[idx] = (1. - double(p)) * s_vec[idx];
            for (int mask_dx = 0; mask_dx < 15; ++mask_dx)
            {
                new_vec[idx] += double(p) / 15. * s_vec[idx ^ (mask_dx + 1)];
            }
        }
        '''

        s_vec = self.vec
        vec_len = len(self.vec)
        arg_names = ['s_vec', 'new_vec', 'vec_len',
                                 'nq', 'p', 'q1', 'q2']
        
        inline(c_code, arg_names=arg_names, compiler='gcc', 
            extra_compile_args=['-O3'], libraries=['rt'])
        
        #self = DensePauliErrorModel(new_vec)
        self.vec = new_vec

        pass

    def flip(self, q, p, flip_type):
        """
        Applies a bit or phase flip (depending on 'flip_type') to qubit
        q with probability p.
        """
        nq = self.nq

        if flip_type == 'X':
            mask = 1 << 2 * nq - q - 1
        else:
            mask = 1 << nq - q - 1
        pass
        new_vec = np.empty((len(self.vec)), dtype=float_type)
        for idx in xrange(len(new_vec)):
            new_vec[idx] = (1. - p) * self.vec[idx] +\
                p * self.vec[idx ^ mask]

        #self = DensePauliErrorModel(new_vec)
        self.vec = new_vec

    def net_prob(self, qubit=None, err_type=None):
        """
        Returns the total probability that an error of a given type 
        occurs on a given bit. Default qubit is the last bit in the
        error model (for syndrome flips). Default error type is 'any',
        note that this does not make sense for syndrome flips.
        """
        
        if qubit is None:
            qubit = self.nq - 1
        
        if not(err_type is None):
            if err_type not in 'xyzXYZ':
                raise ValueError(("Error type '{}' should be "+\
                    "X, Y, Z, or None").format(err_type))
            err_type = err_type.lower()

        qubit_mask = (1 << qubit) + (1 << (self.nq + qubit))

        if err_type == 'x':
            pauli_val = 1 << (self.nq + qubit)
        if err_type == 'y':
            pauli_val = (1 << (self.nq + qubit)) + (1 << qubit)
        if err_type == 'z':
            pauli_val = 1 << qubit
        
        output_prob = 0.
        if err_type is None:
            for idx, prob in enumerate(self.vec):
                if idx & qubit_mask != 0:
                    output_prob += prob
        else:
            for idx, prob in enumerate(self.vec):
                if idx & qubit_mask == pauli_val:
                    output_prob += prob

        return output_prob
    
    def av_weight(self,subset=None):
        """
        Returns the average weight of a sampled error on a subset of 
        the qubits. 
        """
        if subset is None:
            subset = range(self.nq) #All bits
        mean = fsum(self.vec[j] * weight_from_idx(j, self.nq, subset) 
                                            for j in bits(2 * self.nq))
        return mean

    @staticmethod
    def x_flip(p):
        vec = np.zeros((4,), dtype=float_type)
        vec[0] = 1. - p
        vec[2] = p
        return DensePauliErrorModel(vec)

    @staticmethod
    def z_flip(p):
        vec = np.zeros((4,), dtype=float_type)
        vec[0] = 1. - p
        vec[1] = p
        return DensePauliErrorModel(vec)

    @staticmethod
    def dense_dep_model(p):
        vec = np.empty((4,), dtype=float_type)
        vec[0] = 1. - p
        vec[1:] = p / 3.
        return DensePauliErrorModel(vec)

    @staticmethod
    def two_bit_twirl(p):
        vec = np.empty((16,), dtype=float_type)
        vec[0] = 1. - p
        vec[1:] = p / 15.
        return DensePauliErrorModel(vec)

    @staticmethod
    def fowler_meas_model(p, nq, stab_type):
        """
        produces an `nq + 1` - qubit error model representing the 
        output from independent 1- and 2-qubit pauli noise after each
        gate in the measurement of a CSS stabilizer on `nq` bits.
        """
        if stab_type not in 'xzXZ':
            raise ValueError(("stab_type '{}' " + \
                                "is not in 'xzXZ'").format(stab_type))
        stab_type = stab_type.upper()

        flip_type = 'Z' if stab_type == 'X' else 'X'

        output_model = DensePauliErrorModel(nq=nq + 1)
        # state prep error
        output_model.flip(nq, p, flip_type)
        for cnot_idx in reversed(range(nq)):

            # All wait locations prior to CNOT
            dep_p_before = (1. -
                           (DensePauliErrorModel.dense_dep_model(p) **
                            cnot_idx).vec[0])
            output_model.depolarize(cnot_idx, dep_p_before)

            # act CNOT
            if stab_type == 'X':
                output_model.cnot(nq, cnot_idx)
            elif stab_type == 'Z':
                output_model.cnot(cnot_idx, nq)

            # Two-qubit noise after CNOT
            output_model.twirl(nq, cnot_idx, p)

            # All wait locations after CNOT
            dep_p_after = (1. -
                          (DensePauliErrorModel.dense_dep_model(p) **
                          (nq - 1 - cnot_idx)).vec[0])
            output_model.depolarize(cnot_idx, dep_p_after)

        # measurement error
        output_model.flip(nq, p, flip_type)
        return output_model

# ---------------------------Bit Manipulation--------------------------#

def xz_split(num, nb):
    """
    Splits an integer into two, each on `nb` bits.
    """
    x_int = num >> nb
    z_int = num - (x_int << nb)
    return x_int, z_int


def xz_join(x_num, z_num, nb):
    """
    Joins two `nb`-bit bitstrings.
    """
    return (x_num << nb) + z_num


def bit_of(num, k, nb):
    """
    Returns the `k`th bit of some `nb`-bit integer, treating the
    leftmost bit as the 0'th.
    """
    return (num >> (nb - k - 1)) & 1


def bits(n):
    """
    Iterator of length 2**n. If a negative n is given, returns an empty
    iterator.
    """
    return xrange(int(2 ** n))


def c_xor_int(num, i, o, nb):
    """
    Performs an XOR on the `o`'th bit of an n-bit integer if the `i`'th
    bit is equal to 1.
    """
    ith_bit = bit_of(num, i, nb)
    if ith_bit:
        return num ^ (1 << (nb - 1 - o))
    else:
        return num
    pass


def pad_int(num, locs, old_n, new_n):
    """
    Places the bits of an integer into specific locations in a larger
    register.
    """
    return sum(((num >> (old_n - k - 1)) & 1) << (new_n - locs[k] - 1)
               for k in range(len(locs)))


def mask_from_bits(bit_tpl, nb):
    """
    Uses the indices in bit_tpl to generate a mask integer on nb bits.
    Accepts the values on the bits included in bit_tpl under bitwise 
    AND.
    """
    return pad_int(2 ** len(bit_tpl) - 1, bit_tpl, len(bit_tpl), nb)


def pauli_from_int(p_int, nq):
    """
    Field-expedient function, splits an integer into two n-bit halves,
    then produces a phase-free Pauli from them.
    """
    x_half, z_half = xz_split(p_int, nq)
    x_half = int_to_str(x_half, nq)
    z_half = int_to_str(z_half, nq)
    output_pauli = Pauli.from_string(x_half, 'X') *\
        Pauli.from_string(z_half, 'Z')
    output_pauli.ph = 0
    return output_pauli

int_to_str = lambda n_int, nb: bin(n_int).lstrip('0b').zfill(nb)

def weight_from_idx(idx, nq, subset):
    """
    Given an integer index, returns the weight of a Pauli assuming that
    the index, read as a bitstring, is a binary symplectic vector for 
    that Pauli.
    """
    #total hamming weight - number of pauli-Ys (double counting)
    return hamming_weight(idx, nq, subset) - num_ys(idx, nq, subset)

def hamming_weight(index, nq, subset):
    """
    Given an integer and a list containing bit indices, returns the 
    Hamming weight of the index restricted to the subset. 
    """
    mask = mask_from_bits(subset, nq)
    mask += mask << nq
    return bin(index & mask).count('1')

def num_ys(idx, nq, subset):
    mask = mask_from_bits(subset, nq)
    x_half, z_half = xz_split(idx, nq)
    return hamming_weight((x_half & mask) & (z_half & mask), nq, subset)

# ---------------------------------------------------------------------#

# Convenience functions

def _full_lattice_apply(err_mod, register):
    for op in err_mod.ops:
        if len(op) != 1:
            raise ValueError("Only weight-1 Paulis may be " +
                             "used on whole Lattices")

    for point in register.points:
        nu_pauli = _action(err_mod, rand())
        if point.error is None:
            point.error = nu_pauli
        else:
            try:
                if isinstance(point.error, Pauli):
                    point.error = point.error * nu_pauli
                else:
                    point.error = Pauli(point.error) * nu_pauli
            except ValueError as e:
                e.args += (" (pauli: {})".format(point.error),)
                raise e
    pass

def _point_iter_apply(err_mod, register):

    # Test register to see that it contains points
    for point in register:
        if not isinstance(point, Point):
            raise ValueError("Registers which are not "
                             "entire lattices must be iterables of "
                             "points, this input contains non-points.")

    # Test model to see if ops are same len as register
    for op in err_mod.ops:
        if len(op) != len(register):
            raise ValueError("Pauli must be same length " +
                             "as register")

    for pt in register:
        if pt.error is None:
            pt.error = Pauli('I')
        pt.error *= _action(err_mod, rand())


def _point_set_iter_apply(err_mod, register):
    
    for pt_set in register:
        for pt in pt_set:
            if pt.error is None:
                pt.error = I

    for pt_set in register:
        
        error = reduce(lambda a, b: a.tens(b),
                   [pt.error for pt in pt_set])

        error = _action(err_mod, rand()) * error
        
        for idx, pt in enumerate(pt_set):
            pt.error = error[idx]
    pass

def depolarizing_model(p):
    """
    The depolarizing model applies the identity with probability
    :math:`1-p`, and each of the single qubit Pauli operators
    :math:`X`, :math:`Y`, and :math:`Z` with probability
    :math:`\dfrac{p}{3}`.
    """
    return PauliErrorModel({I: 1. - p, X: p / 3.,
                            Y: p / 3., Z: p / 3.})


def iidxz_model(px, pz=None):
    """
    The independent identically-distributed X/Z model applies a bit
    and phase flip to each site, resulting in a reduced probability of
    Y errors.
    """
    if pz is None:
        pz = px

    return PauliErrorModel({I: (1. - px) * (1. - pz), X: px * (1. - pz),
                            Z: (1. - px) * pz, Y: px * pz})


def two_bit_twirl(p):
    """
    With probability p, selects a Pauli at random from the
    non-identity Paulis on two qubits.
    """
    err_dict = {Pauli('II'): 1. - p}
    err_dict.update({pauli: p / 15. for pauli in list(
        pauli_group(2))[1:]})

    return PauliErrorModel(err_dict)


def _action(err_mod, sample):
    """
    returns the operator from an :class:`ErrorModel` corresponding to a
    number between 0 and 1.
    """
    if not (0. <= sample <= 1.):
        raise ValueError("`sample` must be between 0 and 1, " +
                         "preferably uniform.")

    probs, ops = err_mod.probs, err_mod.ops
    for idx, prob in enumerate(probs):
        if sample > prob:
            sample -= prob
        else:
            return ops[idx]
    raise Exception("No operator was selected during sampling.")


def _whitelist_string(whitelist, string):
    """
    Common pattern throughout this file, test a string to see if all
    its letters are members of a whitelist.
    """
    if not all(letter in whitelist for letter in string):
        raise ValueError(("Input string {0} contains letter {1} " +
                          "not found in whitelist {2}")
                         .format(string, letter, whitelist))
    pass

# -------------------------Exception Handling--------------------------#


def sanitize_vector(vector):

    if not isinstance(vector, np.ndarray):
        try:
            vector = np.array(vector, dtype=float_type)
        except Exception, e:
            raise TypeError("Input vector must cast to 1d array, you" +
                            " entered {}".format(type(vector)))

    if len(vector.shape) is not 1:
        raise ValueError("Vector must be 1d, vector has shape" +
                         " {}".format(vector.shape))

    if not np.issubdtype(vector.dtype, np.float):
        raise TypeError("Only floating-point data accepted, you " +
                        "entered {}".format(vector.dtype))

    return vector


def ensure_probabilities(vec, tol=10 ** (-12)):
    """
    Assuming that the input is a 1d array, checks that the elements
    are probabilities; their minimum is greater than 0, their maximum
    is less than 1, and their sum is close to 1.
    """
    minn, maxx, total = np.amin(vec), np.amax(vec), fsum(vec)

    if minn < 0.:
        raise ValueError("Minimum {} should be > 0".format(minn))

    if maxx > 1.:
        raise ValueError("Maximum {} should be < 1".format(maxx))

    if np.abs(total - 1.) > tol:
        raise ValueError("Total {} should == 1".format(total))

    pass
# ---------------------------------------------------------------------#
