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
>>> x, y, z = discretizer.coord
...
>>> A, B, C = sympy.symbols('A B C', commutative=False)
```

```python
>>> H = sympy.Matrix([[kx*A*kx, B*kx], [kx*B, C]])
```

```python
>>> H
⎡kₓ⋅A⋅kₓ  B⋅kₓ⎤
⎢             ⎥
⎣ kₓ⋅B     C  ⎦
```

```python
>>> H = substitute_functions(H, space_dependent={A: (x, y, z), B: x}); H
⎡kₓ⋅A(x, y, z)⋅kₓ  B(x)⋅kₓ⎤
⎢                         ⎥
⎣    kₓ⋅B(x)          C   ⎦
```

```python
>>> discretize(H)
defaultdict(<function discretizer.algorithms.discretize.<locals>.<lambda>>,
            ⎧            ⎡  ⎛  a          ⎞         ⎤             ⎡ ⎛  a          ⎞    ⎛a 
⎪            ⎢-A⎜- ─ + x, y, z⎟         ⎥             ⎢A⎜- ─ + x, y, z⎟   A⎜─ 
⎪            ⎢  ⎝  2          ⎠   ⅈ⋅B(x)⎥             ⎢ ⎝  2          ⎠    ⎝2 
⎪(-1, 0, 0): ⎢──────────────────  ──────⎥, (0, 0, 0): ⎢──────────────── + ────
⎪            ⎢         2           2⋅a  ⎥             ⎢        2              
⎨            ⎢        a                 ⎥             ⎢       a               
⎪            ⎢                          ⎥             ⎢                       
⎪            ⎢   ⅈ⋅B(-a + x)            ⎥             ⎣                0      
⎪            ⎢   ───────────        0   ⎥                                     
⎪            ⎣       2⋅a                ⎦                                     
⎩                                                                             

         ⎞   ⎤             ⎡  ⎛a          ⎞           ⎤⎫
+ x, y, z⎟   ⎥             ⎢-A⎜─ + x, y, z⎟           ⎥⎪
         ⎠   ⎥             ⎢  ⎝2          ⎠   -ⅈ⋅B(x) ⎥⎪
──────────  0⎥, (1, 0, 0): ⎢────────────────  ────────⎥⎪
   2         ⎥             ⎢        2           2⋅a   ⎥⎪
  a          ⎥             ⎢       a                  ⎥⎬
             ⎥             ⎢                          ⎥⎪
            C⎦             ⎢  -ⅈ⋅B(a + x)             ⎥⎪
                           ⎢  ────────────       0    ⎥⎪
                           ⎣      2⋅a                 ⎦⎪
                                                       ⎭)
```

# Generating functions (work in progress)

```python
>>> hop = discretize(H)[1,0,0][0, 0]
```

```python
>>> hop
  ⎛a          ⎞ 
-A⎜─ + x, y, z⎟ 
  ⎝2          ⎠ 
────────────────
        2       
       a
```

```python
>>> from sympy.utilities.lambdify import lambdastr
>>> from sympy.printing.lambdarepr import LambdaPrinter
...
>>> class NumericPrinter(LambdaPrinter):
...     def _print_ImaginaryUnit(self, expr):
...         return "1.j"
```

```python
>>> def make_return_string(expr):
...
>>> #     exp
...     string = lambdastr((), hop, printer=NumericPrinter)[len('lambda : '):]
...
...     return 'return {}'.format(string)
```

```python
>>> make_return_string(H)
'return (-A(a/2 + x, y, z)/a**2)'
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
True
True
```

```python
>>> print(isinstance(Ax, sympy.Function))
>>> print(isinstance(fx, sympy.Function))
True
True
```

```python
>>> print(isinstance(Ax, sympy.function.AppliedUndef))
>>> print(isinstance(fx, sympy.function.AppliedUndef))
True
True
```

```python
>>> Ax == fx
True
```

# AppliedUndef may be useful to separate our functions from others

```python
>>> s = sympy.sin(x)
>>> print(isinstance(s, sympy.function.AppliedUndef))
>>> print(isinstance(s, sympy.Function))
False
True
```

```python
>>> s = sympy.sin
```

```python
>>> sympy.function.Function
Function
```

```python
>>> type(s)
sympy.core.function.FunctionClass
```
