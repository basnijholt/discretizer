# Showcase

This will be notebook showing how our stuff works

## Imports

```python
>>> import sympy
>>> from sympy.interactive import printing
>>> printing.init_printing(use_latex='mathjax')
...
>>> import discretizer
>>> from discretizer import discretize
```

```python
>>> from discretizer.algorithms import substitute_functions
>>> from discretizer.postprocessing import make_kwant_functions
```

```python
>>> from discretizer.algorithms import magic
```

## Defining sample expression

```python
>>> kx, ky, kz = discretizer.momentum_operators
>>> x, y, z = discretizer.coordinates
...
>>> A, B, C = sympy.symbols('A B C', commutative=False)
>>> # H = sympy.Matrix([[kx*A*kx, B*kx], [kx*B, ky*C*ky+kx*sympy.sin(x)]])
... H = kx*A*kx + kx*B + C
```

```python
>>> H
C + kₓ⋅A⋅kₓ + kₓ⋅B
```

```python
>>> hamiltonian = H
>>> space_dependent = {'A', 'B', 'C'}
>>> discrete_coordinates = {'x'}
```

```python
>>> magic(hamiltonian, space_dependent, discrete_coordinates, verbose=True, symbolic_output=True)
Discrete coordinates set to:  ['x']
⎧                      ⎛  a    ⎞                ⎛  a    ⎞    ⎛a    ⎞          
⎪                     A⎜- ─ + x⎟               A⎜- ─ + x⎟   A⎜─ + x⎟          
⎪       ⅈ⋅B(-a + x)    ⎝  2    ⎠                ⎝  2    ⎠    ⎝2    ⎠          
⎨(-1,): ─────────── - ──────────, (0,): C(x) + ────────── + ────────, (1,): - 
⎪           2⋅a            2                        2           2             
⎪                         a                        a           a              
⎩                                                                             

              ⎛a    ⎞⎫
             A⎜─ + x⎟⎪
ⅈ⋅B(a + x)    ⎝2    ⎠⎪
────────── - ────────⎬
   2⋅a           2   ⎪
                a    ⎪
                     ⎭
```

```python
>>> magic(hamiltonian, space_dependent, discrete_coordinates, verbose=True, symbolic_output=False)
Discrete coordinates set to:  ['x']

Function generated for (0,):
def _anonymous_func(site, p):
    x = site.pos
    a = p.a
    C, A = p.C, p.A
    return (C(x) + A(-a/2 + x)/a**2 + A(a/2 + x)/a**2)

Function generated for (-1,):
def _anonymous_func(site1, site2, p):
    x = site2.pos
    a = p.a
    B, A = p.B, p.A
    return (1.j*B(-a + x)/(2*a) - A(-a/2 + x)/a**2)

Function generated for (1,):
def _anonymous_func(site1, site2, p):
    x = site2.pos
    a = p.a
    B, A = p.B, p.A
    return (-1.j*B(a + x)/(2*a) - A(a/2 + x)/a**2)
{(-1,): <function _anonymous_func>,
 (0,): <function _anonymous_func>,
 (1,): <function _anonymous_func>}
```

```python

```
