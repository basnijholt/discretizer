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
>>> H = sympy.Matrix([[kx*A*kx, ky*B*kx], [kx*B*ky, C]]); H
⎡kₓ⋅A⋅kₓ   k_y⋅B⋅kₓ⎤
⎢                  ⎥
⎣kₓ⋅B⋅k_y     C    ⎦
```

```python
>>> space_dependent = {'A', 'B'}
>>> discrete_coordinates = {'x', 'y'}
>>> all_hoppings=False
```

```python
>>> dh = tb_hamiltonian(H, space_dependent, discrete_coordinates, verbose=True, symbolic_output=True,
...       all_hoppings=all_hoppings, interpolate=True)
>>> dh
Discrete coordinates set to:  ['x', 'y']
⎛⎡A(x, y)   A(-a + x, y)   A(a + x, y)   ⎤  ⎧         ⎡              B(a + x, 
⎜⎢─────── + ──────────── + ───────────  0⎥, ⎪(1, -1): ⎢     0        ─────────
⎜⎢    2            2              2      ⎥  ⎪         ⎢                     2 
⎜⎢   a          2⋅a            2⋅a       ⎥  ⎪         ⎢                  4⋅a  
⎜⎢                                       ⎥  ⎨         ⎢                       
⎜⎣                 0                    C⎦  ⎪         ⎢B(x, -a + y)           
⎜                                           ⎪         ⎢────────────       0   
⎜                                           ⎪         ⎢       2               
⎝                                           ⎩         ⎣    4⋅a                

y)⎤          ⎡  A(x, y)   A(a + x, y)   ⎤          ⎡               -B(a + x, y
──⎥, (1, 0): ⎢- ─────── - ───────────  0⎥, (1, 1): ⎢      0        ───────────
  ⎥          ⎢       2           2      ⎥          ⎢                       2  
  ⎥          ⎢    2⋅a         2⋅a       ⎥          ⎢                    4⋅a   
  ⎥          ⎢                          ⎥          ⎢                          
  ⎥          ⎣           0             0⎦          ⎢-B(x, a + y)              
  ⎥                                                ⎢─────────────        0    
  ⎥                                                ⎢        2                 
  ⎦                                                ⎣     4⋅a                  

) ⎤⎫⎞
──⎥⎪⎟
  ⎥⎪⎟
  ⎥⎪⎟
  ⎥⎬⎟
  ⎥⎪⎟
  ⎥⎪⎟
  ⎥⎪⎟
  ⎦⎭⎠
```

```python
>>> onsite, hoppings = tb_hamiltonian(H, space_dependent, discrete_coordinates, verbose=True)
Discrete coordinates set to:  ['x', 'y']

Function generated for (1, 0):
def _anonymous_func(site1, site2, p):
    y, x = site2.pos
    a = p.a
    A = p.A
    return (np.array([[-A(a/2 + x, y)/a**2, 0], [0, 0]]))

Function generated for (0, 0):
def _anonymous_func(site, p):
    y, x = site.pos
    C, a = p.C, p.a
    A = p.A
    return (np.array([[A(-a/2 + x, y)/a**2 + A(a/2 + x, y)/a**2, 0], [0, C]]))

Function generated for (1, -1):
def _anonymous_func(site1, site2, p):
    y, x = site2.pos
    a = p.a
    B = p.B
    return (np.array([[0, B(a + x, y)/(4*a**2)], [B(x, -a + y)/(4*a**2), 0]]))

Function generated for (1, 1):
def _anonymous_func(site1, site2, p):
    y, x = site2.pos
    a = p.a
    B = p.B
    return (np.array([[0, -B(a + x, y)/(4*a**2)], [-B(x, a + y)/(4*a**2), 0]]))
```
