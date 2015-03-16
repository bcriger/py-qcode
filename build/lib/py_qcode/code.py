from qecc import Pauli, commutes_with, eye_p
from lattice import _even_evens, _odd_odds, _squoct_affine_map, skew_coords
from types import FunctionType
from numpy.random import rand

__all__ = ['ErrorCorrectingCode', 'ErrorCheck', 'StabilizerCheck',
           'toric_code', 'square_octagon_code', 'noisy_squoct_code']

class ErrorCheck(object):

    """
    This is the primitive operation of measurement for
    error-correcting codes; it takes a list of errors on a subset of
    the primal lattice of the code and translates it into syndromes on
    a subset of the dual lattice.

    :param primal_sets: co-ordinates from which the error check will
    collect syndromes. I'll add input checking so that tuples of co-
    ordinates can be entered on their own instead of the objects which
    wrap them.

    :type primal_sets: collection of whatever maps to
    :class:`py_qcode.Point` objects

    :param dual_points: points on the dual lattice to which the
    syndrome for this error check will be written

    :type dual_points: set/list of tuples or :class:`py_qcode.Point`
    objects.

    :param rule: lookup table or other mechanism that maps errors to
    syndromes.

    :type rule: function or dict

    :param noise_model: a probability/function pair which is applied to
    the correct syndrome to yield a noisy syndrome.

    :type noise_model: tuple
    """
    def __init__(self, primal_sets, dual_points, rule,
                    noise_model=None, fault_model=None):

        self.primal_sets = primal_sets
        self.dual_points = dual_points
        self.rule = rule

        if noise_model is None:
            noise_model = (0., lambda a: a)

        def noise_func(syndrome):
            prob, func = noise_model
            val = rand()
            if val < prob:
                syndrome = func(syndrome)
            return syndrome

        self.noise_func = noise_func
        self.fault_model = fault_model

    def evaluate(self):
        for idx, point in enumerate(self.dual_points):
            error_str = _sum([pt.error for pt in self.primal_sets[idx]])
            if isinstance(self.rule, dict):
                try:
                    if point.syndrome is None:
                        point.syndrome = \
                            self.noise_func(self.rule[error_str])
                    else:
                        point.syndrome += \
                            self.noise_func(self.rule[error_str])

                except KeyError:
                    raise KeyError("There is no entry in the lookup "
                        "table for the error " + error_str)
                finally:
                    pass

            elif isinstance(self.rule, FunctionType):
                if point.syndrome is None:
                    point.syndrome = \
                        self.noise_func(self.rule(error_str))
                else:
                    point.syndrome += \
                        self.noise_func(self.rule(error_str))
            else:
                raise TypeError("Rule used by error check must be a " +
                                "function or dict, you entered a value of type: "
                                + type(self.rule))


class StabilizerCheck(ErrorCheck):

    """
    subclass of :class:`py_qcode.ErrorCheck`, takes anything that can
    be cast to a :class:`qecc.Pauli` instead of a rule, and uses
    commutation to determine the syndrome.
    """
    def __init__(self, primal_sets, dual_points, stabilizer,
                 noise_model=None, fault_model=None,
                 indy_css=False):

        if (noise_model and fault_model):
            raise SyntaxWarning("Presence of argument 'fault_model' " +
                                "over-rides the argument " +
                                "'noise_model'.")

        if type(stabilizer) is str:
            stabilizer = Pauli(stabilizer)

        # returns 0 if error commutes with stabilizer, 1 if it
        # anti-commutes
        if indy_css == False:
            stab_rule = lambda err_obj: 1 - int(
                commutes_with(stabilizer)(Pauli(err_obj)))
        else:
            # Returns the appropriate letter, X or Z
            def stab_rule(err_obj):
                if type(err_obj) is str:
                    err_pauli = Pauli(err_obj)
                elif type(err_obj) is Pauli:
                    err_pauli = err_obj
                else:
                    raise TypeError(
                        "Input type to stabilizer rule not understood.")

                if all([ltr in 'xX' for ltr in stabilizer.op]):
                    syn_str = 'Z'
                elif all([ltr in 'zZ' for ltr in stabilizer.op]):
                    syn_str = 'X'
                else:
                    raise ValueError("CSS Stabilizers must be all-X " +
                                     "or all-Z; you entered: {0}".format(stabilizer))

                if not(commutes_with(stabilizer)(err_pauli)):
                    return syn_str
                else:
                    return ''

        super(StabilizerCheck, self).__init__(primal_sets, dual_points,
                                              stab_rule, noise_model,
                                              fault_model)

        self.stabilizer = stabilizer
        self.indy_css = indy_css

    def evaluate(self):
        # Use error on first point to typecheck
        test_error = self.primal_sets[0][0].error
        if type(test_error) is str:
            super(StabilizerCheck, self).evaluate()
        elif type(test_error) is Pauli:
            for idx, point in enumerate(self.dual_points):

                multi_bit_error = reduce(lambda p1, p2: p1.tens(p2),
                                         [pt.error for pt in self.primal_sets[idx]])

                if self.fault_model:
                    # Generate Sample, switch syndrome if sample causes
                    # syndrome error
                    if not self.indy_css:
                        raise NotImplementedError("Only CSS " +
                                                  "stabilizers are permitted.")
                    try:
                        big_pauli = self.fault_model.sample()
                    except:
                        raise NotImplementedError(
                            "Use DensePauliErrorModel")

                    nq = self.stabilizer.nq
                    synd_err = eye_p(nq - 1).tens(big_pauli[-1])
                    if not(commutes_with(self.stabilizer)(synd_err)):
                        stab_type = self.stabilizer.op[0]
                        flip_type = 'X' if stab_type == 'Z' else 'Z'
                        self.noise_func = lambda synd: letter_flip(
                            synd, flip_type)
                    else:
                        self.noise_func = lambda synd: synd

                # Collect syndrome with updated noise function
                if point.syndrome == None:
                    point.syndrome = self.noise_func(
                        self.rule(multi_bit_error))
                else:
                    point.syndrome += self.noise_func(
                        self.rule(multi_bit_error))

                # Record remainder of error onto lattice, where it will
                # be picked up in the next time step.
                if self.fault_model:
                    multi_bit_error *= big_pauli[:-1]
                    for jdx, pt in enumerate(self.primal_sets[idx]):
                        pt.error = multi_bit_error[jdx]


class ErrorCorrectingCode():

    """
    Wraps a bunch of parity checks.

    :param parity_check_list: A list of :class:`py_qcode.ErrorCheck`
    objects, which can be a mix of any subclass of
    :class:`py_qcode.ErrorCheck`.

    :type parity_check_list: list

    :param fault_model: a map from :math:`n+1` - bit Pauli operators to
    probabilities, sampled in order to determine whether the syndrome
    is to be flipped, and what the post-measurement error is.

    :type fault_model: :class:`py_qcode.DensePauliErrorModel`

    :param name: a title for the error-correcting code, saved as part
    of simulation logging.

    :type name: str
    """
    def __init__(self, parity_check_list, name='Un-named'):

        self.parity_check_list = parity_check_list
        self.name = name

    def measure(self):
        """
        Evaluates all the parity checks.
        """
        for check in self.parity_check_list:
            check.evaluate()

# UTILITY FUNCTIONS


def toric_code(primal_grid, dual_grid, error_rate=None,
               star_fault_mod=None, plaq_fault_mod=None):
    """
    Uses a few convenience functions to produce the toric code on a set
    of square lattices.

    Necessary arguments are primal_grid and dual_grid, errors from
    primal_grid will be read onto dual_grid.

    If none of error_rate, star_fault_mod, or plaq_fault_mod are
    specified, this function returns a 'clean' toric code, with perfect
    certainty in the results of the measurement and no back-action on
    the lattice.

    Of the optional arguments, two configurations are permitted. Either
     + error_rate is specified, and neither of the fault models are, or
     + *both* fault models are specified and no error_rate is.

    This is to avoid confusion between situations in which the syndrome
    can be corrupted without back-action, and the action of a fault
    model, which affects both the state of the grid and the syndrome.
    """
    # Check for permissible input
    if error_rate and (star_fault_mod or plaq_fault_mod):
        raise NotImplementedError("If error_rate is specified, "
                                    "no fault model can be.")
    elif not(bool(star_fault_mod) == bool(plaq_fault_mod)):
        raise NotImplementedError("If one fault model is specified, "
                                  "they must both be.")

    # The basics: Where to locate the star/plaquette checks
    star_coords = _even_evens(*dual_grid.size)
    star_duals = [dual_grid[coord] for coord in star_coords]
    star_primal = [primal_grid.neighbours(coord)
                   for coord in star_coords]

    plaq_coords = _odd_odds(*dual_grid.size)
    plaq_duals = [dual_grid[coord] for coord in plaq_coords]
    plaq_primal = [primal_grid.neighbours(coord)
                   for coord in plaq_coords]

    # noise models
    star_noise_mod = (error_rate, z_flip) if error_rate else None
    plaq_noise_mod = (error_rate, x_flip) if error_rate else None

    # sanitize inputs that are false-like but not None
    star_fault_mod = star_fault_mod if star_fault_mod else None
    plaq_fault_mod = plaq_fault_mod if plaq_fault_mod else None

    star_check = StabilizerCheck(star_primal, star_duals, 'XXXX',
                                 star_noise_mod, star_fault_mod,
                                 indy_css=True)

    plaq_check = StabilizerCheck(plaq_primal, plaq_duals, 'ZZZZ',
                                 star_noise_mod, star_fault_mod,
                                 indy_css=True)

    name = ""
    name += "Faulty " if star_fault_mod else ""
    name += "Noisy " if error_rate else ""
    name += "Toric Code"

    return ErrorCorrectingCode([star_check, plaq_check], name=name)

# TODO: Refactor Square-Octagon Codes to match single Toric Code
# function above.


def square_octagon_code(primal_grid, dual_grid):
    """
    Uses a few convenience functions to produce the concatenated
    toric/[[4,2,2]] code on a set of square lattices.
    """
    nx, ny = primal_grid.size

    sq_coords = _squoct_affine_map(skew_coords(nx, ny))
    sq_duals = [dual_grid[coord] for coord in sq_coords]
    sq_primal = primal_grid.squares()
    sq_check_X = StabilizerCheck(sq_primal, sq_duals, 'XXXX',
                                 indy_css=True)

    sq_check_Z = StabilizerCheck(sq_primal, sq_duals, 'ZZZZ',
                                 indy_css=True)

    x_oct_coords = _squoct_affine_map(_even_evens(nx, ny))
    x_oct_duals = [dual_grid[coord] for coord in x_oct_coords]
    x_oct_primal = primal_grid.x_octagons()
    x_oct_check = StabilizerCheck(x_oct_primal, x_oct_duals,
                                  'XXXXXXXX', indy_css=True)

    z_oct_coords = _squoct_affine_map(_odd_odds(nx, ny))
    z_oct_duals = [dual_grid[coord] for coord in z_oct_coords]
    z_oct_primal = primal_grid.z_octagons()
    z_oct_check = StabilizerCheck(z_oct_primal, z_oct_duals,
                                  'ZZZZZZZZ', indy_css=True)

    return ErrorCorrectingCode([sq_check_Z, sq_check_X,
                                x_oct_check, z_oct_check],
                               name="Square-Octagon Code")


def noisy_squoct_code(primal_grid, dual_grid, error_rate):
    """
    Uses a few convenience functions to produce the concatenated
    toric/[[4,2,2]] code on a set of square lattices.
    """
    nx, ny = primal_grid.size

    sq_coords = _squoct_affine_map(skew_coords(nx, ny))
    sq_duals = [dual_grid[coord] for coord in sq_coords]
    sq_primal = primal_grid.squares()
    sq_check_X = StabilizerCheck(sq_primal, sq_duals, 'XXXX',
                                (error_rate, z_flip), indy_css=True)

    sq_check_Z = StabilizerCheck(sq_primal, sq_duals, 'ZZZZ',
                                (error_rate, x_flip), indy_css=True)

    x_oct_coords = _squoct_affine_map(_even_evens(nx, ny))
    x_oct_duals = [dual_grid[coord] for coord in x_oct_coords]
    x_oct_primal = primal_grid.x_octagons()
    x_oct_check = StabilizerCheck(x_oct_primal, x_oct_duals,
                                  'XXXXXXXX', (error_rate, z_flip), indy_css=True)

    z_oct_coords = _squoct_affine_map(_odd_odds(nx, ny))
    z_oct_duals = [dual_grid[coord] for coord in z_oct_coords]
    z_oct_primal = primal_grid.z_octagons()
    z_oct_check = StabilizerCheck(z_oct_primal, z_oct_duals,
                                  'ZZZZZZZZ', (error_rate, x_flip), indy_css=True)

    return ErrorCorrectingCode([sq_check_Z, sq_check_X,
                                x_oct_check, z_oct_check],
                               name="Noisy Square-Octagon Code")


def letter_flip(synd, letter):
    """
    This is a convenience function used to 'noise up' input syndromes.
    """
    # print synd
    if synd == letter:
        return ''
    elif synd == '':
        return letter
    else:
        raise ValueError("Unknown syndrome: {0}".format(synd))

z_flip = lambda synd: letter_flip(synd, 'Z')
x_flip = lambda synd: letter_flip(synd, 'X')

_sum = lambda iterable: reduce(lambda a, b: a + b, iterable)
