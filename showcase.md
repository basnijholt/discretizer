# Showcase

This will be notebook showing how our stuff works

## Imports general

```python
>>> import sympy
>>> from sympy.interactive import printing
>>> printing.init_printing(use_latex='mathjax')
```

## Imports discretizer

```python
>>> from discretizer import tb_hamiltonian
>>> from discretizer import momentum_operators
```

## Defining sample expression

```python
>>> kx, ky, kz = momentum_operators
...
>>> A, B, C = sympy.symbols('A B C', commutative=False)
>>> H = sympy.Matrix([[kx*A*kx, kx*B], [B*kx, C]]); H
⎡kₓ⋅A⋅kₓ  kₓ⋅B⎤
⎢             ⎥
⎣ B⋅kₓ     C  ⎦
```

```python
>>> space_dependent = {'A', 'B'}
>>> discrete_coordinates = {'x'}
...
>>> ons, hops = tb_hamiltonian(H, space_dependent, discrete_coordinates, verbose=True,
...                            symbolic_output=True, all_hoppings=True, interpolate=False)
Discrete coordinates set to:  ['x']
```

```python
>>> ons
⎡ ⎛  a    ⎞    ⎛a    ⎞   ⎤
⎢A⎜- ─ + x⎟   A⎜─ + x⎟   ⎥
⎢ ⎝  2    ⎠    ⎝2    ⎠   ⎥
⎢────────── + ────────  0⎥
⎢     2           2      ⎥
⎢    a           a       ⎥
⎢                        ⎥
⎣          0            C⎦
```

```python
>>> hops
⎧       ⎡  ⎛  a    ⎞            ⎤        ⎡  ⎛a    ⎞         ⎤⎫
⎪       ⎢-A⎜- ─ + x⎟            ⎥        ⎢-A⎜─ + x⎟         ⎥⎪
⎪       ⎢  ⎝  2    ⎠    -ⅈ⋅B(x) ⎥        ⎢  ⎝2    ⎠   ⅈ⋅B(x)⎥⎪
⎪(-1,): ⎢────────────   ────────⎥, (1,): ⎢──────────  ──────⎥⎪
⎪       ⎢      2          2⋅a   ⎥        ⎢     2       2⋅a  ⎥⎪
⎨       ⎢     a                 ⎥        ⎢    a             ⎥⎬
⎪       ⎢                       ⎥        ⎢                  ⎥⎪
⎪       ⎢-ⅈ⋅B(-a + x)           ⎥        ⎢ⅈ⋅B(a + x)        ⎥⎪
⎪       ⎢─────────────     0    ⎥        ⎢──────────    0   ⎥⎪
⎪       ⎣     2⋅a               ⎦        ⎣   2⋅a            ⎦⎪
⎩                                                            ⎭
```

```python
>>> onsite, hoppings = tb_hamiltonian(H, space_dependent, discrete_coordinates, verbose=True)
Discrete coordinates set to:  ['x']

Function generated for (0,):
def _anonymous_func(site, p):
    x = site.pos
    C, a = p.C, p.a
    A = p.A
    return (np.array([[A(-a/2 + x)/a**2 + A(a/2 + x)/a**2, 0], [0, C]]))

Function generated for (1,):
def _anonymous_func(site1, site2, p):
    x = site2.pos
    a = p.a
    B, A = p.B, p.A
    return (np.array([[-A(a/2 + x)/a**2, 1.j*B(x)/(2*a)], [1.j*B(a + x)/(2*a), 0]]))
```

```python

```
