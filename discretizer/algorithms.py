import itertools
import sympy
import numpy as np
from collections import defaultdict


# ************************** Some globals *********************************
momentum_operators = sympy.symbols('k_x k_y k_z', commutative=False)
coord = sympy.symbols('x y z', commutative=False)
wf = sympy.Symbol('Psi')(*coord)

lattice_constants = sympy.symbols('a_x a_y a_z')
a = sympy.Symbol('a')


# **************** Operation on sympy expressions **************************
def substitute_functions(expression, space_dependent={}):
    """ Substitute `AppliedUndef` functions into expression.

    Symbols defined in space_dependent will be substitute with a
    function of specified coordinates.

    Parameters:
    -----------
    expression : sympy expression
    space_dependent : dict
        Dictionary list which keys are symbols standing for a space
        dependent function and keys are coordinates on which it depeneds

    Returns:
    --------
    expression

    Note:
    -----
    This function may became part of preprocessing function.
    """
    subs = {}
    for s, v in space_dependent.items():
        if isinstance(v, (tuple, list)):
            for i in v:
                assert i in coord, \
                    "Argument '{}' should be symbol from discretizer.coord.".format(i)
            subs[s] = s(*v)
        else:
            assert v in coord, \
                "Argument '{}' should be symbol from discretizer.coord.".format(v)
            subs[s] = s(v)

    return expression.subs(subs)


def derivate(expression, operator):
    """ Calculate derivate of expression for given momentum operator:

    Parameters:
    -----------
    expression : sympy expression
        Valid sympy expression containing functions to to be derivated
    operator : sympy symbol
        Sympy symbol representing momentum operator

    Returns:
    --------
    expression : derivate of input expression

    Examples:
    ---------
    >>> A = sympy.Symbol('A')
    >>> x = discretizer.coord[0]
    >>> kx = discretizer.momentum_operators[0]
    >>> derivate(A(x), kx)
    -I*(-A(-a_x + x)/(2*a_x) + A(a_x + x)/(2*a_x))
    """
    assert operator in momentum_operators, \
            "Operator '{}' does not belong to [kx, ky, kz].".format(operator)

    if isinstance(expression, (int, float)):
        return 0
    else:
        ind = momentum_operators.index(operator)
        expr1 = expression.subs(coord[ind], coord[ind] + lattice_constants[ind])
        expr2 = expression.subs(coord[ind], coord[ind] - lattice_constants[ind])
        output = (expr1 - expr2) / 2 / lattice_constants[ind]
        return -sympy.I * sympy.expand(output)


def split_factors(expression):
    """ Split symbolic `expression` for a discretization step.

    Parameters:
    -----------
    expression : sympy.Expr instance
        The expression to be split. It should represents single summand.

    Output:
    -------
    lhs : sympy.Expr instance
        Part of expression standing to the left from operators
        that acts in current discretization step.

    operators: sympy.Expr instance
        Operator that perform discretization in current step.

    rhs : sympy.Expr instance
        Part of expression that is derivated in current step.

    Raises:
    -------
    AssertionError
        if input `expression` is of type ``sympy.Add``
    """
    assert not isinstance(expression, sympy.Add), \
        'Input expression must not be sympy.Add. It should be a single summand.'
    output = {'rhs': [1], 'operator': [1], 'lhs': [1]}

    if isinstance(expression, sympy.Pow):
        output['operator'].append(expression.args[0])
        output['lhs'].append(sympy.Pow(expression.args[0], expression.args[1]-1))

    elif isinstance(expression, (int, float, sympy.Integer, sympy.Float)):
        output['rhs'].append(expression)

    elif isinstance(expression, (sympy.Symbol, sympy.Function)):
        if expression in momentum_operators:
            output['operator'].append(expression)
        else:
            output['rhs'].append(expression)

    elif isinstance(expression, sympy.Mul):
        iterator = iter(expression.args[::-1])
        for factor in iterator:
            if factor in momentum_operators:
                output['operator'].append(factor)
                break
            elif factor.func == sympy.Pow and factor.args[0] in momentum_operators:
                operator = factor.args[0]
                power = factor.args[1]
                output['operator'].append(operator)
                output['lhs'].append(sympy.Pow(operator, power-1))
                break
            else:
                output['rhs'].append(factor)

        for factor in iterator:
            output['lhs'].append(factor)

    output = tuple(sympy.Mul(*output[key][::-1])
                   for key in ['lhs', 'operator', 'rhs'])
    return output


def _discretize_summand(summand):
    """ Discretize one summand. """
    assert not isinstance(summand, sympy.Add), "Input should be one summand."

    def do_stuff(expr):
        """ Derivate expr recursively. """
        expr = sympy.expand(expr)

        if isinstance(expr, sympy.Add):
            return do_stuff(expr.args[-1]) + do_stuff(sympy.Add(*expr.args[:-1]))

        lhs, operator, rhs = split_factors(expr)
        if rhs == 1 and operator != 1:
            return 0
        elif operator == 1:
            return lhs*rhs
        elif lhs == 1:
            return derivate(rhs, operator)
        else:
            return do_stuff(lhs*derivate(rhs, operator))

    return do_stuff(summand)


def _discretize_expression(expression):
    """ Discretize continous `expression` into discrete tb representation.

    Parameters:
    -----------
    expression : sympy.Expr instance
        The expression to be discretized.

    Returns:
    --------
    discrete_expression: dict
        dict in which key is offset of hopping ((0, 0, 0) for onsite)
        and value is corresponding hopping (onsite) value.

    Note:
    -----
    Recursive derivation implemented in _discretize_summand is applied
    on every summand. Shortening is applied before return on output.
    """

    if not isinstance(expression, sympy.Expr):
        raise TypeError('Input hamiltonian should be a valid sympy expression.')

    if wf in expression.atoms(sympy.Function):
        raise ValueError("Hamiltonian must not contain {}.".format(wf))

    expression = sympy.expand(expression * wf)

    if expression.func == sympy.Add:
        summands = expression.args
    else:
        summands = [expression]

    outputs = []
    for summand in summands:
        outputs.append(_discretize_summand(summand))

    outputs = [extract_hoppings(summand) for summand in outputs]
    outputs = [shortening(summand) for summand in outputs]

    discrete_expression = defaultdict(int)
    for summand in outputs:
        for k, v in summand.items():
                discrete_expression[k] += v

    return dict(discrete_expression)


def discretize(hamiltonian):
    """ Discretize continous `expression` into discrete tb representation.

    Parameters:
    -----------
    hamiltonian : sympy.Expr instance
        The expression for the Hamiltonian.

    Returns:
    --------
    discrete_hamiltonian: dict
        dict in which key is offset of hopping ((0, 0, 0) for onsite)
        and value is corresponding hopping (onsite) value.

    Note:
    -----
    Recursive derivation implemented in _discretize_summand is applied
    on every summand. Shortening is applied before return on output.
    """

    if not isinstance(hamiltonian, sympy.Matrix):
        return _discretize_expression(hamiltonian)

    shape = hamiltonian.shape

    discrete_hamiltonian = defaultdict(lambda: sympy.zeros(*shape))
    for i,j in itertools.product(range(shape[0]), repeat=2):
        hoppings = _discretize_expression(hamiltonian[i, j])

        for offset, hop in hoppings.items():
            discrete_hamiltonian[offset][i,j] += hop
    return discrete_hamiltonian


# ****** extracting hoppings ***********
def read_hopping_from_wf(inp_psi):
    """ Read hopping from an input wave function.

    Parameters:
    -----------
    inpt_psi : ``~discretizer.algorithms.wf`` like object
        example could be: wf(x + ax, y + 2 * ay, z + 3 * nz)

    Returns:
    --------
    offset: tuple
        offset of a wave function in respect to (x, y, z)

    Examples:
    ---------
    >>> import discretizer
    >>> wf = discretizer.algorithms.wf
    >>> ax, ay, az = discretizer.algorithms.lattice_constants
    >>> x, y, z = discretizer.algorithms.coord
    >>> subs = {x: x+ax, y: y + 2 * ay, z: z + 3 * az}
    >>> expr = wf.subs(subs)
    >>> read_hopping_from_wf(expr)
    (1, 2, 3)
    """
    assert inp_psi.func == wf.func, 'Input should be correct wf used in module.'
    offset = []
    for argument in inp_psi.args:
        temp = sympy.expand(argument)
        if temp in sympy.symbols('x y z', commutative = False):
            offset.append(0)
        elif temp.func == sympy.Add:
            for arg_summands in temp.args:
                if arg_summands.func == sympy.Mul:
                    if len(arg_summands.args) > 2:
                        print('More than two factors in an argument of inp_psi')
                    if not arg_summands.args[0] in sympy.symbols('a_x a_y a_z'):
                        offset.append(arg_summands.args[0])
                    else:
                        offset.append(arg_summands.args[1])
                elif arg_summands in sympy.symbols('a_x a_y a_z'):
                    offset.append(1)
        else:
            print('Argument of \inp_psi is neither a sum nor a single space variable.')
    return tuple(offset)


# This function will not have the shortening included
def extract_hoppings(expr):
    # this line should be unneccessary in the long run. Now I want to avoid errors due to wrong formats of the input.
    expr = sympy.expand(expr)

    # The output will be stored in a dictionary, with terms for each hopping kind
    # This format is probably not good in the long run
    hoppings = {}

    if expr.func == sympy.Add:
        for summand in expr.args:
            #find a way to make it readable
            if not summand.func == wf.func:
                for i in range(len(summand.args)):
                    if summand.args[i].func == wf.func:
                        index = i
                if index < len(summand.args) - 1:
                    print('Psi is not in the very end of the term. Output will be wrong!')

                try:
                    hoppings[read_hopping_from_wf(summand.args[-1])] += sympy.Mul(*summand.args[:-1])
                except:
                    hoppings[read_hopping_from_wf(summand.args[-1])] = sympy.Mul(*summand.args[:-1])
            else:
                try:
                    hoppings[read_hopping_from_wf(summand)] += 1
                except:
                    hoppings[read_hopping_from_wf(summand)] = 1

    else:
        if not expr.func == wf.func:
            for i in range(len(expr.args)):
                if expr.args[i].func == wf.func:
                    index = i
            if index < len(expr.args) - 1:
                print('Psi is not in the very end of the term. Output will be wrong!')

            try:
                hoppings[read_hopping_from_wf(expr.args[-1])] += sympy.Mul(*expr.args[:-1])
            except:
                hoppings[read_hopping_from_wf(expr.args[-1])] = sympy.Mul(*expr.args[:-1])
        else:
            try:
                hoppings[read_hopping_from_wf(expr)] += 1
            except:
                hoppings[read_hopping_from_wf(expr)] = 1
    return hoppings


def shortening(hoppings):
    # make a list of all hopping kinds we have to consider during the shortening
    hops_kinds = np.array(list(hoppings))
    # find the longest hopping range in each direction
    longest_ranges = [np.max(hops_kinds[:,i]) for i in range(len(hops_kinds[0,:]))]
    # define an array in which we are going to store by which factor we
    # can shorten the hoppings in each direction
    shortening_factors = np.ones_like(longest_ranges)
    # Loop over the direction and each potential shortening factor.
    # Inside the loop test whether the hopping distances are actually
    # multiples of the potential shortening factor.
    for dim in np.arange(len(longest_ranges)):
        for factor in np.arange(longest_ranges[dim])+1:
            modulos = np.mod(hops_kinds[:, dim], factor)
            if np.sum(modulos) < 0.1:
                shortening_factors[dim] = factor
    # Apply the shortening factors on the hopping.
    short_hopping = {}
    for hopping_kind in hoppings.keys():
        short_hopping_kind = tuple(np.array(hopping_kind) / shortening_factors)

        for i in short_hopping_kind:
            if isinstance(i, float):
                assert i.is_integer()
        short_hopping_kind = tuple(int(i) for i in short_hopping_kind)

        short_hopping[short_hopping_kind] = hoppings[hopping_kind]
        for lat_const, factor in zip(lattice_constants, shortening_factors):
            factor = int(factor)
            subs = {lat_const: lat_const/factor}
            short_hopping[short_hopping_kind] = short_hopping[short_hopping_kind].subs(subs)

    # We don't need separate a_x, a_y and a_z anymore.
    for key, val in short_hopping.items():
        short_hopping[key] = val.subs({i: a for i in lattice_constants})

    return short_hopping
