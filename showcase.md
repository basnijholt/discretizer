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
>>> from discretizer import Discretizer
>>> from discretizer import tb_hamiltonian
>>> from discretizer import momentum_operators
```

## Defining sample expression

```python
>>> kx, ky, kz = momentum_operators
...
>>> A, B, C = sympy.symbols('A B C', commutative=False)
>>> H = sympy.Matrix([[kx*A*kx +ky*A*ky, kx*B], [B*kx, C]]); H
⎡kₓ⋅A⋅kₓ + k_y⋅A⋅k_y  kₓ⋅B⎤
⎢                         ⎥
⎣       B⋅kₓ           C  ⎦
```

```python
>>> space_dependent = {'A', 'B'}
>>> discrete_coordinates = {'x', 'y'}
```

# class interface

```python
>>> tb = Discretizer(H, space_dependent, discrete_coordinates, lattice_constant=2.0, verbose=True)
Discrete coordinates set to:  ['x', 'y']

Function generated for (0, 1):
def _anonymous_func(site1, site2, p):
    y, x = site2.pos
    A = p.A
    return (np.array([[-0.25*A(x, 1.0 + y), 0], [0, 0]]))

Function generated for (1, 0):
def _anonymous_func(site1, site2, p):
    y, x = site2.pos
    B, A = p.B, p.A
    return (np.array([[-0.25*A(1.0 + x, y), 0.25*1.j*B(x, y)], [0.25*1.j*B(2.0 + x, y), 0]]))

Function generated for (0, 0):
def _anonymous_func(site, p):
    y, x = site.pos
    C = p.C
    A = p.A
    return (np.array([[0.25*A(x, -1.0 + y) + 0.25*A(x, 1.0 + y) + 0.25*A(-1.0 + x, y) + 0.25*A(1.0 + x, y), 0], [0, C]]))
```

```python
>>> tb.symbolic_hamiltonian
⎧        ⎡ ⎛     a    ⎞    ⎛   a    ⎞    ⎛  a       ⎞    ⎛a       ⎞   ⎤       
⎪        ⎢A⎜x, - ─ + y⎟   A⎜x, ─ + y⎟   A⎜- ─ + x, y⎟   A⎜─ + x, y⎟   ⎥       
⎪        ⎢ ⎝     2    ⎠    ⎝   2    ⎠    ⎝  2       ⎠    ⎝2       ⎠   ⎥       
⎪(0, 0): ⎢───────────── + ─────────── + ───────────── + ───────────  0⎥, (0, 1
⎪        ⎢       2              2              2              2       ⎥       
⎨        ⎢      a              a              a              a        ⎥       
⎪        ⎢                                                            ⎥       
⎪        ⎣                            0                              C⎦       
⎪                                                                             
⎪                                                                             
⎩                                                                             

   ⎡  ⎛   a    ⎞    ⎤          ⎡  ⎛a       ⎞            ⎤⎫
   ⎢-A⎜x, ─ + y⎟    ⎥          ⎢-A⎜─ + x, y⎟            ⎥⎪
   ⎢  ⎝   2    ⎠    ⎥          ⎢  ⎝2       ⎠   ⅈ⋅B(x, y)⎥⎪
): ⎢─────────────  0⎥, (1, 0): ⎢─────────────  ─────────⎥⎪
   ⎢       2        ⎥          ⎢       2          2⋅a   ⎥⎪
   ⎢      a         ⎥          ⎢      a                 ⎥⎬
   ⎢                ⎥          ⎢                        ⎥⎪
   ⎣      0        0⎦          ⎢ⅈ⋅B(a + x, y)           ⎥⎪
                               ⎢─────────────      0    ⎥⎪
                               ⎣     2⋅a                ⎦⎪
                                                         ⎭
```

```python
>>> tb.lat
kwant.lattice.Monatomic([[2.0, 0.0], [0.0, 2.0]], [0.0, 0.0], '')
```

```python
>>> tb.onsite, tb.hoppings
(<function _anonymous_func>,
 {HoppingKind((0, 1), kwant.lattice.Monatomic([[2.0, 0.0], [0.0, 2.0]], [0.0, 0.0], '')): <function _anonymous_func>,
  HoppingKind((1, 0), kwant.lattice.Monatomic([[2.0, 0.0], [0.0, 2.0]], [0.0, 0.0], '')): <function _anonymous_func>})
```

# function interface

```python
>>> lat, ons, hops = tb_hamiltonian(H, space_dependent, discrete_coordinates, verbose=True,
...                                 symbolic_output=True, return_conjugated_hops=True,
...                                 interpolate=True)
>>> ons, hops
Discrete coordinates set to:  ['x', 'y']
⎛⎡2⋅A(x, y)   A(x, -a + y)   A(x, a + y)   A(-a + x, y)   A(a + x, y)   ⎤  ⎧  
⎜⎢───────── + ──────────── + ─────────── + ──────────── + ───────────  0⎥, ⎪(-
⎜⎢     2             2              2             2              2      ⎥  ⎪  
⎜⎢    a           2⋅a            2⋅a           2⋅a            2⋅a       ⎥  ⎪  
⎜⎢                                                                      ⎥  ⎨  
⎜⎣                                 0                                   C⎦  ⎪  
⎜                                                                          ⎪  
⎜                                                                          ⎪  
⎝                                                                          ⎩  

       ⎡  A(x, y)   A(-a + x, y)  -ⅈ⋅B(x, y) ⎤           ⎡  A(x, y)   A(x, -a 
1, 0): ⎢- ─────── - ────────────  ───────────⎥, (0, -1): ⎢- ─────── - ────────
       ⎢       2           2          2⋅a    ⎥           ⎢       2           2
       ⎢    2⋅a         2⋅a                  ⎥           ⎢    2⋅a         2⋅a 
       ⎢                                     ⎥           ⎢                    
       ⎢    -ⅈ⋅B(-a + x, y)                  ⎥           ⎣           0        
       ⎢    ────────────────           0     ⎥                                
       ⎣          2⋅a                        ⎦                                
                                                                              

+ y)   ⎤          ⎡  A(x, y)   A(x, a + y)   ⎤          ⎡  A(x, y)   A(a + x, 
────  0⎥, (0, 1): ⎢- ─────── - ───────────  0⎥, (1, 0): ⎢- ─────── - ─────────
       ⎥          ⎢       2           2      ⎥          ⎢       2           2 
       ⎥          ⎢    2⋅a         2⋅a       ⎥          ⎢    2⋅a         2⋅a  
       ⎥          ⎢                          ⎥          ⎢                     
      0⎦          ⎣           0             0⎦          ⎢     ⅈ⋅B(a + x, y)   
                                                        ⎢     ─────────────   
                                                        ⎣          2⋅a        
                                                                              

y)  ⅈ⋅B(x, y)⎤⎫⎞
──  ─────────⎥⎪⎟
       2⋅a   ⎥⎪⎟
             ⎥⎪⎟
             ⎥⎬⎟
             ⎥⎪⎟
        0    ⎥⎪⎟
             ⎦⎪⎟
              ⎭⎠
```

```python
>>> lat, ons, hops = tb_hamiltonian(H, space_dependent, discrete_coordinates, verbose=True, lattice_constant=0.5)
Discrete coordinates set to:  ['x', 'y']

Function generated for (0, 1):
def _anonymous_func(site1, site2, p):
    y, x = site2.pos
    A = p.A
    return (np.array([[-4.0*A(x, 0.25 + y), 0], [0, 0]]))

Function generated for (1, 0):
def _anonymous_func(site1, site2, p):
    y, x = site2.pos
    B, A = p.B, p.A
    return (np.array([[-4.0*A(0.25 + x, y), 1.0*1.j*B(x, y)], [1.0*1.j*B(0.5 + x, y), 0]]))

Function generated for (0, 0):
def _anonymous_func(site, p):
    y, x = site.pos
    C = p.C
    A = p.A
    return (np.array([[4.0*A(x, -0.25 + y) + 4.0*A(x, 0.25 + y) + 4.0*A(-0.25 + x, y) + 4.0*A(0.25 + x, y), 0], [0, C]]))
```

```python

```
