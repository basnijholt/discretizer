import itertools
import sympy
import numpy as np
from collections import defaultdict

from sympy.utilities.lambdify import lambdastr
from sympy.printing.lambdarepr import LambdaPrinter
from sympy.core.function import AppliedUndef

class NumericPrinter(LambdaPrinter):
    def _print_ImaginaryUnit(self, expr):
        return "1.j"


# ************************** Some globals *********************************
wavefunction_name = 'Psi'


# **************** Operation on sympy expressions **************************
def substitute_functions(expression, space_dependent={}):
    """ Substitute `AppliedUndef` functions into expression.

    Symbols defined in space_dependent will be substitute with a
    function of specified coordinates.

    Parameters:
    -----------
    expression : sympy.Expr instance
    space_dependent : dict
        Dictionary list which keys are symbols standing for a space dependent
        parameters and values are list or tuple with strings representing names
        of space dependent coordinates.

    Returns:
    --------
    expression : sympy.Expr instance
        Expression containing space dependent functions.

    """
    subs = {}
    for s, v in space_dependent.items():
        if isinstance(v, (tuple, list)):
            for i in v:
                assert i in ['x', 'y', 'z'], \
                    "Argument '{}' should be list of ['x', 'y', 'z'].".format(i)
            subs[s] = s(*sympy.symbols(v, commutative=False))
        else:
            assert v in ['x', 'y', 'z'], \
                "Argument '{}' should be list of ['x', 'y', 'z'].".format(v)
            subs[s] = s(sympy.Symbol(v, commutative=False))

    return expression.subs(subs)


def split_factors(expression, discrete_coordinates=('x', 'y', 'z')):
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

    momentum_names = ['k_{}'.format(s) for s in discrete_coordinates]
    momentum_operators = sympy.symbols(momentum_names, commutative=False)

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


def derivate(expression, operator):
    """ Calculate derivate of expression for given momentum operator:

    Parameters:
    -----------
    expression : sympy.Expr instance
        Sympy expression containing functions to to be derivated.
    operator : sympy.Symbol
        Sympy symbol representing momentum operator.

    Returns:
    --------
    output : sympy.Expr instance
        Derivated input expression.

    Examples:
    ---------
    >>> from discretizer.algorithms import derivate
    >>> import sympy
    >>> A = sympy.Function('A')
    >>> x = sympy.Symbol('x')
    >>> kx = sympy.Symbol('k_x')
    >>> derivate(A(x), kx)
    -I*(-A(-a_x + x)/(2*a_x) + A(a_x + x)/(2*a_x))
    """
    if not isinstance(operator, sympy.Symbol):
        raise TypeError("Input operator '{}' is not type sympy.Symbol.")

    if operator.name not in ['k_x', 'k_y', 'k_z']:
        raise ValueError("Input operator '{}' unkown.".format(operator))

    if isinstance(expression, (int, float, sympy.Symbol)):
        return 0
    else:
        coordinate_name = operator.name.split('_')[1]
        ct = sympy.Symbol(coordinate_name, commutative=True)
        cf = sympy.Symbol(coordinate_name, commutative=False)
        h = sympy.Symbol('a_'+coordinate_name)

        expr1 = expression.subs({ct: ct + h, cf: cf + h})
        expr2 = expression.subs({ct: ct - h, cf: cf - h})
        output = (expr1 - expr2) / 2 / h
        return -sympy.I * sympy.expand(output)


def _discretize_summand(summand, discrete_coordinates=('x', 'y', 'z')):
    """ Discretize one summand. """
    assert not isinstance(summand, sympy.Add), "Input should be one summand."

    def do_stuff(expr):
        """ Derivate expr recursively. """
        expr = sympy.expand(expr)

        if isinstance(expr, sympy.Add):
            return do_stuff(expr.args[-1]) + do_stuff(sympy.Add(*expr.args[:-1]))

        lhs, operator, rhs = split_factors(expr, discrete_coordinates)
        if rhs == 1 and operator != 1:
            return 0
        elif operator == 1:
            return lhs*rhs
        elif lhs == 1:
            return derivate(rhs, operator)
        else:
            return do_stuff(lhs*derivate(rhs, operator))

    return do_stuff(summand)


def _discretize_expression(expression, discrete_coordinates=('x', 'y', 'z')):
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
        raise TypeError('Input expression should be a valid sympy expression.')

    expression = sympy.expand(expression)

    if expression.func == sympy.Add:
        summands = expression.args
    else:
        summands = [expression]

    outputs = []
    for summand in summands:
        outputs.append(_discretize_summand(summand, discrete_coordinates))

    outputs = [extract_hoppings(summand) for summand in outputs]
    outputs = [shortening(summand, discrete_coordinates) for summand in outputs]

    discrete_expression = defaultdict(int)
    for summand in outputs:
        for k, v in summand.items():
                discrete_expression[k] += v

    return dict(discrete_expression)


def discretize(hamiltonian, discrete_coordinates=('x', 'y', 'z')):
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
    coordinates = sympy.symbols(discrete_coordinates, commutative=False)
    wf = sympy.Function(wavefunction_name)(*coordinates)

    if wf in hamiltonian.atoms(sympy.Function):
        raise ValueError("Hamiltonian must not contain {}.".format(wf))

    if not isinstance(hamiltonian, sympy.Matrix):
        return _discretize_expression(hamiltonian * wf, discrete_coordinates)

    shape = hamiltonian.shape

    discrete_hamiltonian = defaultdict(lambda: sympy.zeros(*shape))
    for i,j in itertools.product(range(shape[0]), repeat=2):
        expression = hamiltonian[i, j] * wf
        hoppings = _discretize_expression(expression, discrete_coordinates)

        for offset, hop in hoppings.items():
            discrete_hamiltonian[offset][i,j] += hop
    return discrete_hamiltonian


# ************ Making kwant functions ***********
def make_return_string(expr):
    """Process a sympy expression into an evaluatable Python return statement.

    Parameters:
    -----------
    expr : sympy.Expr instance

    Returns:
    --------
    output : string
        A return string that can be used to assemble a Kwant value function.
    func_symbols : set of sympy.Symbol instances
        All space dependent functions that appear in the expression.
    const_symbols : set of sympy.Symbol instances
        All constants that appear in the expression.
    """
    func_symbols = {sympy.Symbol(i.func.__name__) for i in
                    expr.atoms(AppliedUndef)}

    free_symbols = {i for i in expr.free_symbols if i not in coordinates}
    const_symbols = free_symbols - func_symbols

    output = lambdastr((), expr, printer=NumericPrinter)[len('lambda : '):]
    output = output.replace('MutableDenseMatrix', 'np.array')

    return 'return {}'.format(output), func_symbols, const_symbols


def assign_symbols(func_symbols, const_symbols, onsite=True,
                   discrete_coordinates=('x', 'y', 'z')):
    """Generate a series of assingments defining a set of symbols.

    Parameters:
    -----------
    func_symbols : set of sympy.Symbol instances
        All space dependent functions that appear in the expression.
    const_symbols : set of sympy.Symbol instances
        All constants that appear in the expression.

    Returns:
    --------
    assignments : list of strings
        List of lines used for including in a function.

    Notes:
    where A, B, C are all the free symbols plus the symbols that appear on the
    ------
    The resulting lines begin with a coordinates assignment of form
    `x,y,z = site.pos` when onsite=True, or
    `x,y,z = site2.pos` when onsite=False

    followed by two lines of form
    `A, B, C = p.A, p.B, p.C`
    `f, g, h = p.f, p.g, p.h`
    where A, B, C are symbols representing constants and f, g, h are symbols
    representing functions. Separation of constant and func symbols is probably
    not necessary but I leave it for now, just in case.
    """
    lines = []
    func_names = [i.name for i in func_symbols]
    const_names = [i.name for i in const_symbols]

    lines.insert(0, ', '.join(func_names) + ' = p.' +
                 ', p.'.join(func_names))

    lines.insert(0, ', '.join(const_names) + ' = p.' +
                 ', p.'.join(const_names))

    if onsite:
        site = 'site'
    else:
        site = 'site2'

    lines.insert(0, '{} = {}.pos'.format(', '.join(discrete_coordinates), site))

    return lines


def value_function(content, name='_anonymous_func', onsite=True, verbose=False):
    """Generate a Kwant value function from a list of lines containing its body.

    Parameters:
    -----------
    content : list of lines
        Lines forming the body of the function.
    name : string
        Function name (not important).
    onsite : bool
        If True, the function call signature will be `f(site, p)`, otherwise
        `f(site1, site2, p)`.
    verbose : bool
        Whether the function bodies should be printed.

    Returns:
    --------
    f : function
        The function defined in a namespace containing only the cached values
        of Pauli matrices.
    """
    if not content[-1].startswith('return'):
        raise ValueError('The function does not end with a return statement')

    separator = '\n' + 4 * ' '
    site_string = 'site' if onsite else 'site1, site2'
    header = 'def {0}({1}, p):'.format(name, site_string)
    func_code = separator.join([header] + list(content))

    namespace = {}
    if verbose:
        print(func_code)
    exec(func_code, namespace)
    return namespace[name]


def make_kwant_functions(discrete_hamiltonian, verbose=False):
    """Transform discrete hamiltonian into valid kwant functions.

    Parameters:
    -----------


    """
    pass


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
    >>> x, y, z = discretizer.algorithms.coordinates
    >>> subs = {x: x+ax, y: y + 2 * ay, z: z + 3 * az}
    >>> expr = wf.subs(subs)
    >>> read_hopping_from_wf(expr)
    (1, 2, 3)
    """
    assert inp_psi.func.__name__ == wavefunction_name, \
        'Input should be function that represents wavefunction in module.'
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
    hoppings = defaultdict(int)

    if expr.func == sympy.Add:
        for summand in expr.args:
            #find a way to make it readable
            if not summand.func.__name__ == wavefunction_name:
                for i in range(len(summand.args)):
                    if summand.args[i].func.__name__ == wavefunction_name:
                        index = i
                if index < len(summand.args) - 1:
                    print('Psi is not in the very end of the term. Output will be wrong!')
                hoppings[read_hopping_from_wf(summand.args[-1])] += sympy.Mul(*summand.args[:-1])
            else:
                hoppings[read_hopping_from_wf(summand)] += 1

    else:
        if not expr.func.__name__ == wavefunction_name:
            for i in range(len(expr.args)):
                if expr.args[i].func.__name__ == wavefunction_name:
                    index = i
            if index < len(expr.args) - 1:
                print('Psi is not in the very end of the term. Output will be wrong!')

            hoppings[read_hopping_from_wf(expr.args[-1])] += sympy.Mul(*expr.args[:-1])
        else:
            hoppings[read_hopping_from_wf(expr)] += 1
    return hoppings


def shortening(hoppings, discrete_coordinates=('x', 'y', 'z')):
    """ Perform shortening of hoppings."""
    tmps = ['a_{}'.format(s) for s in discrete_coordinates]
    lattice_constants = sympy.symbols(tmps)
    a = sympy.Symbol('a')

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
        short_hopping[key] = val.subs({i: a for i in sympy.symbols('a_x a_y a_z')})

    return short_hopping
