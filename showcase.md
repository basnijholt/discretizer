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
>>> from discretizer.algorithms import make_kwant_functions
```

## Defining sample expression

```python
>>> # def magic(hamiltonian, space_dependent=None, discrete_coordinates=None,
... #           verbose=False, symbolic_output=False):
...
... #     tmp = substitute_functions(hamiltonian, discrete_coordinates, space_dependent)
... #     hamiltonian, discrete_coordinates = tmp
...
... #     if verbose:
... #         print('Discrete coordinates set to: ', sorted(discrete_coordinates))
...
... #     discrete_hamiltonian = discretize(hamiltonian, discrete_coordinates)
...
... #     if symbolic_output:
... #         return discrete_hamiltonian
...
... #     tb = make_kwant_functions(discrete_hamiltonian, discrete_coordinates, verbose)
... #     return tb
```

```python
>>> kx, ky, kz = discretizer.momentum_operators
>>> x, y, z = discretizer.coordinates
...
>>> A, B, C = sympy.symbols('A B C', commutative=False)
>>> H = sympy.Matrix([[kx*A*kx, B*kx], [kx*B, ky*C*ky+sympy.sin(x)]])
```

```python
>>> H
⎡kₓ⋅A⋅kₓ         B⋅kₓ       ⎤
⎢                           ⎥
⎣ kₓ⋅B    k_y⋅C⋅k_y + sin(x)⎦
```

```python
>>> from discretizer.algorithms import magic
```

```python
>>> magic(H, space_dependent=[B, C], discrete_coordinates={'x', 'y'}, verbose=True,
...             symbolic_output=True)
Discrete coordinates set to:  ['x', 'y']
defaultdict(<function discretizer.algorithms.discretize.<locals>.<lambda>>,
            ⎧         ⎡     -A         ⅈ⋅B(x, y)⎤                                         
⎪(-1, 0): ⎢     ───        ─────────⎥, (0, -1): ⎡0         0       ⎤, (0, 0): 
⎪         ⎢       2           2⋅a   ⎥           ⎢                  ⎥          
⎪         ⎢      a                  ⎥           ⎢     ⎛     a    ⎞ ⎥          
⎪         ⎢                         ⎥           ⎢   -C⎜x, - ─ + y⎟ ⎥          
⎨         ⎢ⅈ⋅B(-a + x, y)           ⎥           ⎢     ⎝     2    ⎠ ⎥          
⎪         ⎢──────────────      0    ⎥           ⎢0  ───────────────⎥          
⎪         ⎣     2⋅a                 ⎦           ⎢           2      ⎥          
⎪                                               ⎣          a       ⎦          
⎪                                                                             
⎩                                                                             

⎡2⋅A                                      ⎤                                   
⎢───                   0                  ⎥, (0, 1): ⎡0        0      ⎤, (1, 0
⎢  2                                      ⎥          ⎢                ⎥       
⎢ a                                       ⎥          ⎢     ⎛   a    ⎞ ⎥       
⎢                                         ⎥          ⎢   -C⎜x, ─ + y⎟ ⎥       
⎢               ⎛     a    ⎞    ⎛   a    ⎞⎥          ⎢     ⎝   2    ⎠ ⎥       
⎢              C⎜x, - ─ + y⎟   C⎜x, ─ + y⎟⎥          ⎢0  ─────────────⎥       
⎢               ⎝     2    ⎠    ⎝   2    ⎠⎥          ⎢          2     ⎥       
⎢ 0   sin(x) + ───────────── + ───────────⎥          ⎣         a      ⎦       
⎢                     2              2    ⎥                                   
⎣                    a              a     ⎦                                   

   ⎡      -A         -ⅈ⋅B(x, y) ⎤⎫
): ⎢      ───        ───────────⎥⎪
   ⎢        2            2⋅a    ⎥⎪
   ⎢       a                    ⎥⎪
   ⎢                            ⎥⎪
   ⎢-ⅈ⋅B(a + x, y)              ⎥⎬
   ⎢───────────────       0     ⎥⎪
   ⎣      2⋅a                   ⎦⎪
                                 ⎪
                                 ⎪
                                 ⎭)
```

```python
>>> tb = magic(H, space_dependent=[A, B, C], discrete_coordinates={'x'}, verbose=True)
Discrete coordinates set to:  ['x']
Function generated for (0,):
def _anonymous_func(site, p):
    x = site.pos
    a, k_y = p.a, p.k_y
    C, A = p.C, p.A
    return (np.array([[A(-a/2 + x)/a**2 + A(a/2 + x)/a**2, 0], [0, k_y**2*C(x) + sin(x)]]))

Function generated for (-1,):
def _anonymous_func(site1, site2, p):
    x = site2.pos
    a = p.a
    B, A = p.B, p.A
    return (np.array([[-A(-a/2 + x)/a**2, I*B(x)/(2*a)], [I*B(-a + x)/(2*a), 0]]))

Function generated for (1,):
def _anonymous_func(site1, site2, p):
    x = site2.pos
    a = p.a
    B, A = p.B, p.A
    return (np.array([[-A(a/2 + x)/a**2, -I*B(x)/(2*a)], [-I*B(a + x)/(2*a), 0]]))
```
