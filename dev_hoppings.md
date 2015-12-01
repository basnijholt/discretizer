---
dictitems:
  kernelspec:
    display_name: dev2
    language: python2
    name: dev2
  language_info:
    codemirror_mode:
      name: ipython
      version: 2
    file_extension: .py
    mimetype: text/x-python
    name: python
    nbconvert_exporter: python
    pygments_lexer: ipython2
    version: 2.7.6
state:
  _allownew: true
---

# Possible inputs to hopping functions as outputs from discretization

```python
>>> import sympy
>>> sympy.init_printing(use_latex='mathjax')
>>> from discretizer import discretize
```

```python
>>> from discretizer.algorithms import wavefunction_name
>>> from discretizer.algorithms import _discretize_summand
```

```python
>>> from discretizer.algorithms import read_hopping_from_wf
>>> from discretizer.algorithms import extract_hoppings
```

```python
>>> def get_wf(discrete_coordinates):
...     """Rturn Psi depending only on discrete_coordinates. Help function."""
...     coordinates_names = sorted(list(discrete_coordinates))
...     coordinates = [sympy.Symbol(c, commutative=False) for c in coordinates_names]
...     Psi = sympy.Function(wavefunction_name)(*coordinates)
...
...     return Psi
```

# Support stuff

```python
>>> kx, ky, kz = sympy.symbols('k_x k_y k_z', commutative=False)
>>> x, y, z = sympy.symbols('x y z', commutative=False)
>>> A, B = sympy.symbols('A B', commutative=False)
```

```python
>>> discrete_coordinates = {'x', 'y'}
>>> Psi = get_wf(discrete_coordinates)
>>> summand = kx*A(x)*kx * Psi
```

```python
>>> inp = _discretize_summand(summand, discrete_coordinates)
>>> summand
kₓ⋅A(x)⋅kₓ⋅Ψ(x, y)
```

```python
>>> inp
A(-aₓ + x)⋅Ψ(x, y)   A(-aₓ + x)⋅Ψ(-2⋅aₓ + x, y)   A(aₓ + x)⋅Ψ(x, y)   A(aₓ + x
────────────────── - ────────────────────────── + ───────────────── - ────────
          2                        2                        2                 
      4⋅aₓ                     4⋅aₓ                     4⋅aₓ                  

)⋅Ψ(2⋅aₓ + x, y)
────────────────
     2          
 4⋅aₓ
```

```python
>>> print(inp)
A(-a_x + x)*Psi(x, y)/(4*a_x**2) - A(-a_x + x)*Psi(-2*a_x + x, y)/(4*a_x**2) + A(a_x + x)*Psi(x, y)/(4*a_x**2) - A(a_x + x)*Psi(2*a_x + x, y)/(4*a_x**2)
```

```python
>>> hoppings = extract_hoppings(inp, discrete_coordinates); hoppings
⎧           ⎛  a    ⎞            ⎛  a    ⎞    ⎛a    ⎞            ⎛a    ⎞ ⎫
⎪         -A⎜- ─ + x⎟           A⎜- ─ + x⎟   A⎜─ + x⎟          -A⎜─ + x⎟ ⎪
⎪           ⎝  2    ⎠            ⎝  2    ⎠    ⎝2    ⎠            ⎝2    ⎠ ⎪
⎨(-1, 0): ────────────, (0, 0): ────────── + ────────, (1, 0): ──────────⎬
⎪               2                    2           2                  2    ⎪
⎪              a                    a           a                  a     ⎪
⎩                                                                        ⎭
```

```python
>>> for key, val in hoppings.items():
...     print("{}: '{}',".format(key, val))
(1, 0): '-A(a/2 + x)/a**2',
(0, 0): 'A(-a/2 + x)/a**2 + A(a/2 + x)/a**2',
(-1, 0): '-A(-a/2 + x)/a**2',
```
