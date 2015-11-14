import itertools
import sympy
import numpy as np
from math import factorial
from functools import reduce
from operator import mul


# ************* modified Anton's functions ********************
def make_tb_coefs():
    cache = {}

    def tb_coefs_1d(n):
        """Return tight-binding coefficients of k**n."""
        if n in cache:
            return cache[n]

        d = 1 + (n+1) // 2
        a = (1j * np.arange(1-d, d))**np.arange(2*d - 1).reshape(-1, 1)
        a = sympy.Matrix(a)
        a = a.inv().T
        result = list(factorial(n) * (a.row(n)) / sympy.Symbol('a')**n)

        cache[n] = result
        return result
    return tb_coefs_1d

tb_coefs_1d = make_tb_coefs()
del make_tb_coefs


def tb_coefs(k_powers):
    """Calculate tight-binding coefficients of a monomial of momenta.

    Parameters:
    -----------
    k_powers : list of integers
       List with powers of momenta in each direction, e.g. `k_x**2 * k_y` would
       correspond to `(2, 1)`.

    Returns:
    --------
    coefficients : iterator
        An iterator of (vector, coefficient) tuples, where vector is an integer
        tuple corresponding to the hopping, and coefficient is a sympy
        constant expression (complex rational).
    """
    dir_coefs = [tb_coefs_1d(n)[::-1] for n in k_powers]
    deltas = [range(-len(i)//2 + 1, len(i)//2 + 1) for i in dir_coefs]
    hopping_vecs = (i for i in itertools.product(*deltas))
    hopping_coefs = (i for i in itertools.product(*dir_coefs))
    prod = lambda x, y=1: reduce(mul, x, y)
    result = ((i, prod(j)) for i, j in zip(hopping_vecs, hopping_coefs))
    return ((i, j) for i, j in result if j != 0)


# **************** Operation on sympy expressions **************************
# Some globals
coord = sympy.symbols('x y z', commutative=False)
momentum_operators = sympy.symbols('k_x k_y k_z', commutative=False)
a = sympy.Symbol('a')


def substitute_functions(expr, space_dependent=[]):
    """ Substitute space_dependent symbols with function of (x, y, z) """
    symbols = [s for s in expr.atoms(sympy.Symbol) if s.name in space_dependent]
    subs = {s: sympy.Function(s.name)(*coord) for s in symbols}

    return expr.subs(subs)


def derivate(expression, k_powers):
    """ Calculate derivate of expression for given momentum powers:

    Parameters:
    -----------
    expression : sympy expression
        Valid sympy expression containing functions to to be derivated
    k_powers : list of integers
       List with powers of momenta in each direction, e.g. `k_x**2 * k_y` would
       correspond to `(2, 1, 0)`.

    Returns:
    --------
    expression : derivate of input expression
    """
    coefs = tb_coefs(k_powers)
    output = []
    for v, coef in coefs:
        subs = {c: c+v*a for c, v in zip(coord, v)}
        output.append(coef * expression.subs(subs))
    return sympy.expand(sympy.Add(*output))


def split_factors(expr):
    """ Split symbolic expression for a discretization step.

    Parameters:
    -----------
    expr : sympy expression
        symbolic expression to be split

    Output:
    -------
    lhs : sympy expression
        part of expression standing to the left from operators
        that acts in current discretization step

    operators: sympy expression
        operators that perform discretization in current step

    rhs : sympy expression
        part of expression that is derivated in current step

    """
    assert isinstance(expr, sympy.Mul), 'input expr is not sympy.Mul'
    output = {'rhs': [], 'operators': [], 'lhs': []}
    normal = True

    iterator = iter(expr.args[::-1])
    for factor in iterator:
        if factor in momentum_operators:
            output['operators'].append(factor)
            break
        elif factor.func == sympy.Pow and factor.args[0] in momentum_operators:
            output['operators'].append(factor)
            break
        else:
            output['rhs'].append(factor)

    for factor in iterator:
        if factor in momentum_operators:
            output['operators'].append(factor)
        elif factor.func == sympy.Pow and factor.args[0] in momentum_operators:
            output['operators'].append(factor)
        else:
            output['lhs'].append(factor)
            break

    for factor in iterator:
        output['lhs'].append(factor)

    output = tuple(sympy.Mul(*output[key][::-1])
                   for key in ['lhs', 'operators', 'rhs'])
    return output


def get_powers(expr):
    """ Read powers of momentum operators.

    Parameters:
    -----------
    expr : sympy expression
       Expression that is built of momentum operators, e.g. `k_x**2 * k_y`.

    Returns:
    --------
    powers : tuple
        Tuple of momentum powers, e.g. `k_x**2 * k_y` would return `(2, 1, 0)`.
    """

    gens = [sympy.Symbol(s.name) for s in momentum_operators]
    subs = {c: n for c,n in zip(momentum_operators, gens)}

    expr = sympy.poly(expr.subs(subs), gens)
    monomials = [(power, coef)
        for power, coef in zip(expr.monoms(), expr.coeffs())]

    if len(monomials) != 1:
        raise ValueError('Momentum operator is not built from simple Mul-s.')

    powers, coef = monomials[0]
    if coef != 1:
        raise ValueError('Momentum operator contains sth except momentum operators.')

    return powers


def recursive(expr):
    """ Recursive derivation.

    Function moved out of notebook because Sebastian wants to work on it.
    """
    expr = sympy.expand(expr)

    def do_stuff(expr):
        lhs, operators, rhs = split_factors(expr)
        if lhs == 1:
            return derivate(rhs, get_powers(operators))
        else:
            return recursive(lhs*derivate(rhs, get_powers(operators)))

    if expr.func == sympy.Mul:
        return do_stuff(expr)

    elif expr.func == sympy.Add:
        return do_stuff(expr.args[-1]) + recursive(sympy.Add(*expr.args[:-1]))

    else:
        raise ValueError('Incorrect input', expr)
