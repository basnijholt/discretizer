import sympy
import discretizer
from discretizer.algorithms import split_factors
from discretizer.algorithms import wavefunction_name
from discretizer.algorithms import derivate
from discretizer.algorithms import _discretize_summand
from discretizer.algorithms import read_hopping_from_wf

from nose.tools import raises
import numpy as np

kx, ky, kz = sympy.symbols('k_x k_y k_z', commutative=False)
x, y, z = sympy.symbols('x y z', commutative=False)
ax, ay, az = sympy.symbols('a_x a_y a_z')

wf =  sympy.Function(wavefunction_name)
Psi = sympy.Function(wavefunction_name)(x, y, z)
A, B = sympy.symbols('A B', commutative=False)

ns = {'A': A, 'B': B, 'a_x': ax, 'a_y': ay, 'az': az, 'x': x, 'y': y, 'z': z}

def test_split_factors_1():
    test = {
        kz * Psi                      : (1, kz, Psi),
        A * kx**2 * Psi               : (A * kx, kx, Psi),
        A * kx**2 * ky * Psi          : (A * kx**2, ky, Psi),
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
        got = split_factors(inp, discrete_coordinates={'x', 'y', 'z'})
        assert  got == out,\
            "Should be: split_factors({})=={}. Not {}".format(inp, out, got)


@raises(AssertionError)
def test_split_factors_2():
    split_factors(A+B, discrete_coordinates={'x', 'y', 'z'})


def test_derivate_1():
    test = {
        (A(x), kx): '-I*(-A(-a_x + x)/(2*a_x) + A(a_x + x)/(2*a_x))',
        (A(x), ky): '0',
        (A(x)*B, kx): '-I*(-A(-a_x + x)*B/(2*a_x) + A(a_x + x)*B/(2*a_x))',
        (A(x) + B(x), kx): '-I*(-A(-a_x + x)/(2*a_x) + A(a_x + x)/(2*a_x) - B(-a_x + x)/(2*a_x) + B(a_x + x)/(2*a_x))',
        (A, kx): '0',
        (5, kx): '0',
        (A(x) * B(x), kx): '-I*(-A(-a_x + x)*B(-a_x + x)/(2*a_x) + A(a_x + x)*B(a_x + x)/(2*a_x))',
        (Psi, ky): '-I*(-Psi(x, -a_y + y, z)/(2*a_y) + Psi(x, a_y + y, z)/(2*a_y))',
    }

    for inp, out in test.items():
        got = (derivate(*inp))
        out = sympy.sympify(out, locals=ns)
        assert  sympy.simplify(sympy.expand(got - out)) == 0,\
            "Should be: derivate({0[0]}, {0[1]})=={1}. Not {2}".format(inp, out, got)


@raises(TypeError)
def test_derivate_2():
    derivate(A(x), kx**2)


@raises(ValueError)
def test_derivate_3():
    derivate(A(x), sympy.Symbol('A'))


def test_discretize_summand_1():
    test = {
        kx * A(x): '-I*(-A(-a_x + x)/(2*a_x) + A(a_x + x)/(2*a_x))',
        kx * Psi: '-I*(-Psi(-a_x + x, y, z)/(2*a_x) + Psi(a_x + x, y, z)/(2*a_x))',
        kx**2 * Psi: 'Psi(x, y, z)/(2*a_x**2) - Psi(-2*a_x + x, y, z)/(4*a_x**2) - Psi(2*a_x + x, y, z)/(4*a_x**2)',
        kx * A(x) * kx * Psi: 'A(-a_x + x)*Psi(x, y, z)/(4*a_x**2) - A(-a_x + x)*Psi(-2*a_x + x, y, z)/(4*a_x**2) + A(a_x + x)*Psi(x, y, z)/(4*a_x**2) - A(a_x + x)*Psi(2*a_x + x, y, z)/(4*a_x**2)',
    }

    for inp, out in test.items():
        got = _discretize_summand(inp, discrete_coordinates={'x', 'y', 'z'})
        out = sympy.sympify(out, locals=ns)
        assert  sympy.simplify(sympy.expand(got - out)) == 0,\
            "Should be: _discretize_summand({})=={}. Not {}".format(inp, out, got)

@raises(AssertionError)
def test_discretize_summand_2():
    _discretize_summand(kx*A(x)+ B(x), discrete_coordinates={'x', 'y', 'z'})


def test_read_hoppings_from_wf_1():
    offsets = [(0,0,0), (1,0,0), (0,1,0), (0,0,1), (1,1,1), (1,2,3)]
    test = {}

    for offset in offsets:
        nx, ny, nz = offset
        key = Psi.subs({x: x + nx*ax, y: y + ny * ay, z: z + nz * az})
        test[key] = offset

    for inp, out in test.items():
        got = read_hopping_from_wf(inp)
        assert got == out,\
            "Should be: read_hopping_from_wf({}) == {}. Not {}".format(inp, out, got)


def test_read_hoppings_from_wf_2():
    test = {
        wf(x, y, z): (0,0,0),
        wf(x, y): (0, 0),
        wf(x, z): (0, 0),
        wf(x+ax, y-2*ay): (1, -2),
        wf(x, az+3*az): (0, 3)
    }

    for inp, out in test.items():
        got = read_hopping_from_wf(inp)
        assert got == out,\
            "Should be: read_hopping_from_wf({}) == {}. Not {}".format(inp, out, got)

@raises(AssertionError)
def test_read_hoppings_from_wf_2():
    read_hopping_from_wf(A)
    read_hopping_from_wf(5*Psi)
    read_hopping_from_wf(Psi+2)
    read_hopping_from_wf(Psi + A)
