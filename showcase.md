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

## Defining sample expression

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
>>> magic(H, space_dependent=[A, B, C], discrete_coordinates={'x', 'y'}, verbose=True,
...             symbolic_output=True)
Discrete coordinates set to:  ['x', 'y']
defaultdict(<function discretizer.algorithms.discretize.<locals>.<lambda>>,
            ⎧         ⎡  ⎛  a       ⎞            ⎤                                        
⎪         ⎢-A⎜- ─ + x, y⎟            ⎥                                        
⎪         ⎢  ⎝  2       ⎠   ⅈ⋅B(x, y)⎥                                        
⎪(-1, 0): ⎢───────────────  ─────────⎥, (0, -1): ⎡0         0       ⎤, (0, 0):
⎪         ⎢        2           2⋅a   ⎥           ⎢                  ⎥         
⎪         ⎢       a                  ⎥           ⎢     ⎛     a    ⎞ ⎥         
⎨         ⎢                          ⎥           ⎢   -C⎜x, - ─ + y⎟ ⎥         
⎪         ⎢ⅈ⋅B(-a + x, y)            ⎥           ⎢     ⎝     2    ⎠ ⎥         
⎪         ⎢──────────────       0    ⎥           ⎢0  ───────────────⎥         
⎪         ⎣     2⋅a                  ⎦           ⎢           2      ⎥         
⎪                                                ⎣          a       ⎦         
⎪                                                                             
⎩                                                                             

 ⎡ ⎛  a       ⎞    ⎛a       ⎞                                      ⎤          
 ⎢A⎜- ─ + x, y⎟   A⎜─ + x, y⎟                                      ⎥          
 ⎢ ⎝  2       ⎠    ⎝2       ⎠                                      ⎥          
 ⎢───────────── + ───────────                   0                  ⎥, (0, 1): 
 ⎢       2              2                                          ⎥          
 ⎢      a              a                                           ⎥          
 ⎢                                                                 ⎥          
 ⎢                                       ⎛     a    ⎞    ⎛   a    ⎞⎥          
 ⎢                                      C⎜x, - ─ + y⎟   C⎜x, ─ + y⎟⎥          
 ⎢                                       ⎝     2    ⎠    ⎝   2    ⎠⎥          
 ⎢             0               sin(x) + ───────────── + ───────────⎥          
 ⎢                                             2              2    ⎥          
 ⎣                                            a              a     ⎦          

                            ⎡   ⎛a       ⎞               ⎤⎫
                            ⎢ -A⎜─ + x, y⎟               ⎥⎪
                            ⎢   ⎝2       ⎠    -ⅈ⋅B(x, y) ⎥⎪
⎡0        0      ⎤, (1, 0): ⎢ ─────────────   ───────────⎥⎪
⎢                ⎥          ⎢        2            2⋅a    ⎥⎪
⎢     ⎛   a    ⎞ ⎥          ⎢       a                    ⎥⎪
⎢   -C⎜x, ─ + y⎟ ⎥          ⎢                            ⎥⎬
⎢     ⎝   2    ⎠ ⎥          ⎢-ⅈ⋅B(a + x, y)              ⎥⎪
⎢0  ─────────────⎥          ⎢───────────────       0     ⎥⎪
⎢          2     ⎥          ⎣      2⋅a                   ⎦⎪
⎣         a      ⎦                                        ⎪
                                                          ⎪
                                                          ⎭)
```

```python
>>> magic(H, space_dependent=[A, B, C], discrete_coordinates={'x', 'y'}, verbose=True);
Discrete coordinates set to:  ['x', 'y']

Function generated for (0, 1):
def _anonymous_func(site1, site2, p):
    y, x = site2.pos
    a = p.a
    C = p.C
    return (np.array([[0, 0], [0, -C(x, a/2 + y)/a**2]]))

Function generated for (0, -1):
def _anonymous_func(site1, site2, p):
    y, x = site2.pos
    a = p.a
    C = p.C
    return (np.array([[0, 0], [0, -C(x, -a/2 + y)/a**2]]))

Function generated for (1, 0):
def _anonymous_func(site1, site2, p):
    y, x = site2.pos
    a = p.a
    A, B = p.A, p.B
    return (np.array([[-A(a/2 + x, y)/a**2, -1.j*B(x, y)/(2*a)], [-1.j*B(a + x, y)/(2*a), 0]]))

Function generated for (0, 0):
def _anonymous_func(site, p):
    y, x = site.pos
    a = p.a
    A, C = p.A, p.C
    return (np.array([[A(-a/2 + x, y)/a**2 + A(a/2 + x, y)/a**2, 0], [0, sin(x) + C(x, -a/2 + y)/a**2 + C(x, a/2 + y)/a**2]]))

Function generated for (-1, 0):
def _anonymous_func(site1, site2, p):
    y, x = site2.pos
    a = p.a
    A, B = p.A, p.B
    return (np.array([[-A(-a/2 + x, y)/a**2, 1.j*B(x, y)/(2*a)], [1.j*B(-a + x, y)/(2*a), 0]]))
```

```python

```
