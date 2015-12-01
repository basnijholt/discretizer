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
>>> discrete_coordinates = {'x', 'y', 'z'}
>>> Psi = get_wf(discrete_coordinates)
>>> summand = kx**2*A(x, y)*ky * Psi
```

```python
>>> inp = _discretize_summand(summand, discrete_coordinates)
>>> summand
  2                       
kₓ ⋅A(x, y)⋅k_y⋅Ψ(x, y, z)
```

```python
>>> inp
ⅈ⋅A(x, y)⋅Ψ(x, -a_y + y, z)   ⅈ⋅A(x, y)⋅Ψ(x, a_y + y, z)   ⅈ⋅A(-2⋅aₓ + x, y)⋅Ψ
─────────────────────────── - ────────────────────────── - ───────────────────
             2                            2                                   
         4⋅aₓ ⋅a_y                    4⋅aₓ ⋅a_y                             8⋅

(-2⋅aₓ + x, -a_y + y, z)   ⅈ⋅A(-2⋅aₓ + x, y)⋅Ψ(-2⋅aₓ + x, a_y + y, z)   ⅈ⋅A(2⋅
──────────────────────── + ────────────────────────────────────────── - ──────
  2                                            2                              
aₓ ⋅a_y                                    8⋅aₓ ⋅a_y                          

aₓ + x, y)⋅Ψ(2⋅aₓ + x, -a_y + y, z)   ⅈ⋅A(2⋅aₓ + x, y)⋅Ψ(2⋅aₓ + x, a_y + y, z)
─────────────────────────────────── + ────────────────────────────────────────
              2                                          2                    
          8⋅aₓ ⋅a_y                                  8⋅aₓ ⋅a_y
```

```python
>>> print(inp)
I*A(x, y)*Psi(x, -a_y + y, z)/(4*a_x**2*a_y) - I*A(x, y)*Psi(x, a_y + y, z)/(4*a_x**2*a_y) - I*A(-2*a_x + x, y)*Psi(-2*a_x + x, -a_y + y, z)/(8*a_x**2*a_y) + I*A(-2*a_x + x, y)*Psi(-2*a_x + x, a_y + y, z)/(8*a_x**2*a_y) - I*A(2*a_x + x, y)*Psi(2*a_x + x, -a_y + y, z)/(8*a_x**2*a_y) + I*A(2*a_x + x, y)*Psi(2*a_x + x, a_y + y, z)/(8*a_x**2*a_y)
```

```python
>>> hoppings = extract_hoppings(inp, discrete_coordinates); hoppings
⎧             -ⅈ⋅A(-a + x, y)               ⅈ⋅A(-a + x, y)              ⅈ⋅A(x,
⎪(-1, -1, 0): ────────────────, (-1, 1, 0): ──────────────, (0, -1, 0): ──────
⎨                      3                            3                        3
⎪                   2⋅a                          2⋅a                        a 
⎩                                                                             

 y)             -ⅈ⋅A(x, y)               -ⅈ⋅A(a + x, y)              ⅈ⋅A(a + x
───, (0, 1, 0): ───────────, (1, -1, 0): ───────────────, (1, 1, 0): ─────────
                      3                           3                          3
                     a                         2⋅a                        2⋅a 
                                                                              

, y)⎫
────⎪
    ⎬
    ⎪
    ⎭
```

```python
>>> for key, val in hoppings.items():
...     print("{}: '{}',".format(key, val))
(-1, 1, 0): 'I*A(-a + x, y)/(2*a**3)',
(1, 1, 0): 'I*A(a + x, y)/(2*a**3)',
(0, -1, 0): 'I*A(x, y)/a**3',
(0, 1, 0): '-I*A(x, y)/a**3',
(1, -1, 0): '-I*A(a + x, y)/(2*a**3)',
(-1, -1, 0): '-I*A(-a + x, y)/(2*a**3)',
```

```python

```
