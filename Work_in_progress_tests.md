{
 "cells": [],
 "metadata": {},
 "nbformat": 4,
 "nbformat_minor": 0
}

```python
>>> import sympy
>>> import discretizer
>>> from discretizer.algorithms import split_factors
>>> from discretizer.algorithms import wavefunction_name
>>> from discretizer.algorithms import derivate
>>> from discretizer.algorithms import _discretize_summand
>>> from discretizer.algorithms import read_hopping_from_wf
>>> from discretizer.algorithms import extract_hoppings
...
>>> from nose.tools import raises
>>> import numpy as np
...
>>> kx, ky, kz = sympy.symbols('k_x k_y k_z', commutative=False)
>>> x, y, z = sympy.symbols('x y z', commutative=False)
>>> ax, ay, az = sympy.symbols('a_x a_y a_z')
...
>>> wf =  sympy.Function(wavefunction_name)
>>> Psi = sympy.Function(wavefunction_name)(x, y, z)
>>> A, B = sympy.symbols('A B', commutative=False)
...
>>> ns = {'A': A, 'B': B, 'a_x': ax, 'a_y': ay, 'az': az, 'x': x, 'y': y, 'z': z}
```

```python
>>> temp = _discretize_summand((kx*A(x)*kx+kx)*Psi, ('x y')).subs(sympy.Symbol('a_x'), sympy.Symbol('ax')).subs(sympy.Symbol('a_y'), sympy.Symbol('ay')).subs(sympy.Symbol('az'), sympy.Symbol('a'))
>>> temp
-I*(-Psi(-ax + x, y, z)/(2*ax) + Psi(ax + x, y, z)/(2*ax)) + A(-ax + x)*Psi(x, y, z)/(4*ax**2) - A(-ax + x)*Psi(-2*ax + x, y, z)/(4*ax**2) + A(ax + x)*Psi(x, y, z)/(4*ax**2) - A(ax + x)*Psi(2*ax + x, y, z)/(4*ax**2)
```

```python
>>> out = dict(extract_hoppings(temp))
...
>>> out
{(-2, 0, 0): -A(-ax + x)/(4*ax**2),
 (-1, 0, 0): I/(2*ax),
 (0, 0): -I/(2*ax),
 (0, 0, 0): A(-ax + x)/(4*ax**2) + A(ax + x)/(4*ax**2),
 (2, 0, 0): -A(ax + x)/(4*ax**2)}
```

```python

```

```python
>>> def test_extract_hoppings():
...     test = {1j*(Psi.subs(x, x + ax) - Psi.subs(x, x - ax)):
...             {(-1, 0, 0): -1.0*1j, (1, 0, 0): 1.0*1j},
...             A(-ax + x)*Psi/(4*ax**2) - A(-ax + x)*Psi.subs(x, -2*ax + x)/(4*ax**2)
...             + A(ax + x)*Psi/(4*ax**2) - A(ax + x)*Psi.subs(x, 2*ax + x)/(4*ax**2):
...             {(-2, 0, 0): -A(-ax + x)/(4*ax**2), (0, 0, 0): A(-ax + x)/(4*ax**2)
...             + A(ax + x)/(4*ax**2), (2, 0, 0): -A(ax + x)/(4*ax**2)}}
...
...     for inp, out in test.items():
...         got = dict(extract_hoppings(inp))
...         assert  got == out,\
...             "Should be: extract_hoppings({})=={}. Not {}".format(inp, out, got)
```

```python
>>> test_extract_hoppings()
```

```python

```
