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
```

## Defining sample expression

```python
>>> kx, ky, kz = discretizer.momentum_operators
>>> x, y, z = discretizer.coordinates
...
>>> A, B, C = sympy.symbols('A B C', commutative=False)
```

```python
>>> H1 = ky * A(y) * ky
>>> discretize(H1, discrete_coordinates=['y'])
⎧         ⎛  a    ⎞          ⎛  a    ⎞    ⎛a    ⎞          ⎛a    ⎞ ⎫
⎪       -A⎜- ─ + y⎟         A⎜- ─ + y⎟   A⎜─ + y⎟        -A⎜─ + y⎟ ⎪
⎪         ⎝  2    ⎠          ⎝  2    ⎠    ⎝2    ⎠          ⎝2    ⎠ ⎪
⎨(-1,): ────────────, (0,): ────────── + ────────, (1,): ──────────⎬
⎪             2                  2           2                2    ⎪
⎪            a                  a           a                a     ⎪
⎩                                                                  ⎭
```

```python
>>> H1 = kx * A(x, y) * kx
>>> discretize(H1, discrete_coordinates=['x', 'y'])
⎧           ⎛  a       ⎞            ⎛  a       ⎞    ⎛a       ⎞            ⎛a  
⎪         -A⎜- ─ + x, y⎟           A⎜- ─ + x, y⎟   A⎜─ + x, y⎟          -A⎜─ +
⎪           ⎝  2       ⎠            ⎝  2       ⎠    ⎝2       ⎠            ⎝2  
⎨(-1, 0): ───────────────, (0, 0): ───────────── + ───────────, (1, 0): ──────
⎪                 2                       2              2                    
⎪                a                       a              a                     
⎩                                                                             

     ⎞ ⎫
 x, y⎟ ⎪
     ⎠ ⎪
───────⎬
 2     ⎪
a      ⎪
       ⎭
```

```python
>>> H = sympy.Matrix([[kx*A*kx, B*kx], [kx*B, ky*C*ky+sympy.sin(x)]])
```

```python
>>> H
⎡kₓ⋅A⋅kₓ         B⋅kₓ       ⎤
⎢                           ⎥
⎣ kₓ⋅B    k_y⋅C⋅k_y + sin(x)⎦
```

```python
>>> H = substitute_functions(H, space_dependent={A: ('x', 'y'), B: ('x',), C:('y',)}); H
⎡kₓ⋅A(x, y)⋅kₓ         B(x)⋅kₓ       ⎤
⎢                                    ⎥
⎣   kₓ⋅B(x)     k_y⋅C(y)⋅k_y + sin(x)⎦
```

```python
>>> H.atoms(sympy.Symbol)
set([kₓ, k_y, x, y])
```

```python
>>> discretize(H, discrete_coordinates=['x', 'y'])
defaultdict(<function discretizer.algorithms.discretize.<locals>.<lambda>>,
            ⎧         ⎡  ⎛  a       ⎞         ⎤                                      ⎡ ⎛  
⎪         ⎢-A⎜- ─ + x, y⎟         ⎥                                      ⎢A⎜- 
⎪         ⎢  ⎝  2       ⎠   ⅈ⋅B(x)⎥                                      ⎢ ⎝  
⎪(-1, 0): ⎢───────────────  ──────⎥, (0, -1): ⎡0       0      ⎤, (0, 0): ⎢────
⎪         ⎢        2         2⋅a  ⎥           ⎢               ⎥          ⎢    
⎪         ⎢       a               ⎥           ⎢     ⎛  a    ⎞ ⎥          ⎢    
⎨         ⎢                       ⎥           ⎢   -C⎜- ─ + y⎟ ⎥          ⎢    
⎪         ⎢  ⅈ⋅B(-a + x)          ⎥           ⎢     ⎝  2    ⎠ ⎥          ⎢    
⎪         ⎢  ───────────      0   ⎥           ⎢0  ────────────⎥          ⎢    
⎪         ⎣      2⋅a              ⎦           ⎢         2     ⎥          ⎢    
⎪                                             ⎣        a      ⎦          ⎢    
⎪                                                                        ⎢    
⎩                                                                        ⎣    

a       ⎞    ⎛a       ⎞                                ⎤                      
─ + x, y⎟   A⎜─ + x, y⎟                                ⎥                      
2       ⎠    ⎝2       ⎠                                ⎥                      
───────── + ───────────                0               ⎥, (0, 1): ⎡0      0   
   2              2                                    ⎥          ⎢           
  a              a                                     ⎥          ⎢     ⎛a    
                                                       ⎥          ⎢   -C⎜─ + y
                                   ⎛  a    ⎞    ⎛a    ⎞⎥          ⎢     ⎝2    
                                  C⎜- ─ + y⎟   C⎜─ + y⎟⎥          ⎢0  ────────
                                   ⎝  2    ⎠    ⎝2    ⎠⎥          ⎢        2  
         0               sin(x) + ────────── + ────────⎥          ⎣       a   
                                       2           2   ⎥                      
                                      a           a    ⎦                      

             ⎡  ⎛a       ⎞           ⎤⎫
             ⎢-A⎜─ + x, y⎟           ⎥⎪
             ⎢  ⎝2       ⎠   -ⅈ⋅B(x) ⎥⎪
  ⎤, (1, 0): ⎢─────────────  ────────⎥⎪
  ⎥          ⎢       2         2⋅a   ⎥⎪
⎞ ⎥          ⎢      a                ⎥⎪
⎟ ⎥          ⎢                       ⎥⎬
⎠ ⎥          ⎢-ⅈ⋅B(a + x)            ⎥⎪
──⎥          ⎢────────────      0    ⎥⎪
  ⎥          ⎣    2⋅a                ⎦⎪
  ⎦                                   ⎪
                                      ⎪
                                      ⎭)
```

# Generating functions (work in progress)

```python
>>> hop = discretize(H)[1,0,0]; hop
⎡  ⎛a       ⎞           ⎤
⎢-A⎜─ + x, y⎟           ⎥
⎢  ⎝2       ⎠   -ⅈ⋅B(x) ⎥
⎢─────────────  ────────⎥
⎢       2         2⋅a   ⎥
⎢      a                ⎥
⎢                       ⎥
⎢-ⅈ⋅B(a + x)            ⎥
⎢────────────      0    ⎥
⎣    2⋅a                ⎦
```

```python
>>> from discretizer.algorithms import make_kwant_functions, make_return_string
>>> from discretizer.algorithms import value_function, assign_symbols
```

```python
>>> value = hop
>>> onsite = True
...
>>> return_string, func_symbols, const_symbols = make_return_string(value)
>>> lines = assign_symbols(func_symbols, const_symbols, onsite=onsite, discrete_coordinates=['x', 'y'])
>>> lines.append(return_string)
...
>>> f = value_function(lines, verbose=True, onsite=False)
def _anonymous_func(site1, site2, p):
    x, y = site.pos
    a = p.a
    B, A = p.B, p.A
    return (np.array([[-A(a/2 + x, y)/a**2, -I*B(x)/(2*a)], [-I*B(a + x)/(2*a), 0]]))
```

# Identifying

```python
>>> x = sympy.Symbol('x')
...
>>> A = sympy.Symbol('A')
>>> f = sympy.Function('A')
...
>>> Ax = A(x)
>>> fx = f(x)
```

```python
>>> print(isinstance(A, sympy.Symbol))
>>> print(isinstance(f, sympy.function.UndefinedFunction))
```

```python
>>> print(isinstance(Ax, sympy.Function))
>>> print(isinstance(fx, sympy.Function))
```

```python
>>> print(isinstance(Ax, sympy.function.AppliedUndef))
>>> print(isinstance(fx, sympy.function.AppliedUndef))
```

```python
>>> Ax == fx
```

# AppliedUndef may be useful to separate our functions from others

```python
>>> s = sympy.sin(x)
>>> print(isinstance(s, sympy.function.AppliedUndef))
>>> print(isinstance(s, sympy.Function))
```
