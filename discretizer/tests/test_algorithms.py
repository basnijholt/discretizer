import sympy
import discretizer
from discretizer.algorithms import split_factors
from discretizer.algorithms import derivate

from nose.tools import raises
import numpy as np

kx, ky, kz = discretizer.algorithms.momentum_operators
x, y, z = discretizer.algorithms.coord
Psi = discretizer.algorithms.wf
A, B, C = sympy.symbols('A B C', commutative=False)


def test_split_factors_1():
    test = {
        kz * Psi                      : (1, kz, Psi),
        C * kx**2 * Psi               : (C * kx, kx, Psi),
        C * kx**2 * ky * Psi          : (C * kx**2, ky, Psi),
        ky * A * kx * B * Psi         : (ky * A, kx, B * Psi),
        kx                            : (1, kx, 1),
        kx**2                         : (kx, kx, 1),
        A                             : (1, 1, A),
        A(x, y, z)                    : (1, 1, A(x, y, z)),
        Psi                           : (1, 1, Psi),
        np.int(5)                     : (1, 1, np.int(5)),
        np.float(5)                   : (1, 1, np.float(5)),
        sympy.Integer(5)              : (1, 1, sympy.Integer(5)),
        sympy.Float(5)                : (1, 1, sympy.Float(5)),
        1                             : (1, 1, 1),
        1.0                           : (1, 1, 1.0),
        5                             : (1, 1, 5),
        5.0                           : (1, 1, 5.0),
    }

    for inp, out in test.items():
        got = split_factors(inp)
        assert  got == out,\
            "Should be: split_factors({})=={}. Not {}".format(inp, out, got)


@raises(AssertionError)
def test_split_factors_2():
    split_factors(A+B)


def test_derivate_1():
    test = {
        (A(x), kx)       : "-I*(-A(-a_x + x)/(2*a_x) + A(a_x + x)/(2*a_x))",
        (A(x), ky)       : "0",
        (A(x) * B, kx)   : "-I*(-A(-a_x + x)*B/(2*a_x) + A(a_x + x)*B/(2*a_x))"

    }

    for inp, out in test.items():
        got = sympy.srepr(sympy.expand(sympy.sympify(str(derivate(*inp)))))
        out = sympy.srepr(sympy.expand(sympy.sympify(out)))
        assert  got == out,\
            "Should be: derivate({})=={}. Not {}".format(inp, out, got)
