import sympy
import discretizer
from discretizer.algorithms import wavefunction_name
from discretizer.algorithms import read_hopping_from_wf

from nose.tools import raises
from nose.tools import assert_raises
import numpy as np


kx, ky, kz = sympy.symbols('k_x k_y k_z', commutative=False)
x, y, z = sympy.symbols('x y z', commutative=False)
ax, ay, az = sympy.symbols('a_x a_y a_z')

wf =  sympy.Function(wavefunction_name)
Psi = sympy.Function(wavefunction_name)(x, y, z)
A, B = sympy.symbols('A B', commutative=False)

ns = {'A': A, 'B': B, 'a_x': ax, 'a_y': ay, 'az': az, 'x': x, 'y': y, 'z': z}


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
        wf(x, z+3*az): (0, 3),
        wf(y, z): (0, 0),
        wf(y, z+az): (0, 1),
    }

    for inp, out in test.items():
        got = read_hopping_from_wf(inp)
        assert got == out,\
            "Should be: read_hopping_from_wf({}) == {}. Not {}".format(inp, out, got)


def test_test_read_hoppings_from_wf_ValueError():
    tests = {
        wf(x+ay, ay),
        wf(y, x),
        wf(x, x),
        wf(x, ax),
        wf(y, A),
    }
    for inp in tests:
        assert_raises(ValueError, read_hopping_from_wf, inp)


def test_test_read_hoppings_from_wf_TypeError():
    tests = {
        wf(x,y,z) + A,
        A(x,y),
        5*Psi,
        Psi+2,
        A,
    }
    for inp in tests:
        assert_raises(TypeError, read_hopping_from_wf, inp)
