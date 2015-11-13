## Assumptions
* space_dependent - these are cast to functions, not commutative
* momentum_symbols - will be specified (default to kx, ky, kz probably) and not commutative
* constant - commutative, everything else
* coordinates are always create as symbol with commutative=False
* lattice constant is always a

```python
>>> from sympy.printing.dot import dotprint
>>> from graphviz import Source
>>> graph = lambda x: Source(dotprint(x))
```

```python
>>> import sympy
>>> from sympy.interactive import printing
>>> printing.init_printing(use_latex='mathjax')
...
>>> from discretizer import coord, momentum_operators, a
>>> from discretizer import substitute_functions, derivate
```

## Defining input expression

```python
>>> A, B, C, D = sympy.symbols('A B C D', commutative=False)
>>> space_dependent = ['A', 'B']
```

```python
>>> expr = B*C*A + D; expr
B⋅C⋅A + D
```

```python
>>> expr = substitute_functions(expr, space_dependent); expr
D + B(x, y, z)⋅C⋅A(x, y, z)
```

## Calculation derivation

```python
>>> expr
D + B(x, y, z)⋅C⋅A(x, y, z)
```

```python
>>> derivate(expr, (1, 0, 0))
  - -0.5⋅ⅈ⋅B(-a + x, y, z)⋅C⋅A(-a + x, y, z)    0.5⋅ⅈ⋅B(a + x, y, z)⋅C⋅A(a + x
- ─────────────────────────────────────────── + ──────────────────────────────
                       a                                          a

, y, z)
───────
```

# Designing recursive algorithm
We must somehow handle this expression recursively
First idea may be to:
1. Expand intput expression to obtain Add-s of Mul-s
2. For every summand do:
    * split into lhs, operators, rhs
    * calculate derivative
    * pass it as input to 1.
3. gather stuff together

```python
>>> kx, ky, kz = momentum_operators
>>> Psi = sympy.Symbol('Psi', commutative=False) #sympy.Function('Psi')(*coord)
```

```python
>>> expr = expr = kx*A*kx*Psi
>>> expr = substitute_functions(expr, ['A', 'Psi']); expr
kₓ⋅A(x, y, z)⋅kₓ⋅Ψ(x, y, z)
```

```python
>>> def recursive(expr):
...     print(expr)
...     expr = sympy.expand(expr)
...
...     def do_stuff(expr):
...         lhs, operators, rhs = split_factors(expr)
...         if lhs == 1:
...             return derivate(rhs, get_powers(operators))
...         else:
...             return recursive(lhs*derivate(rhs, get_powers(operators)))
...
...     if expr.func == sympy.Mul:
...         return do_stuff(expr)
...
...     elif expr.func == sympy.Add:
...         return do_stuff(expr.args[-1]) + recursive(sympy.Add(*expr.args[:-1]))
...
...     else:
...         raise ValueError('Incorrect input', expr)
```

```python
>>> recursive(expr)
C*k_x**2*k_y*Psi + k_x*A*k_x*Psi + k_z*Psi
0
```

```python

```

```python
>>> def rec_sum(x):
...     if len(x) == 1:
...         return x[0]
...     else:
...         return x[-1] + rec_sum(x[:-1])
```

```python
>>> x = list(range(10))
```

```python
>>> sum(x)
45
```

```python
>>> rec_sum(x)
45
```

# Playing around

```python
>>> expr = kx*A*kx + C * kx**2 * ky + kz
>>> expr = expr * Psi
>>> expr = sympy.expand(expr); expr
    2                          
C⋅kₓ ⋅k_y⋅Ψ + kₓ⋅A⋅kₓ⋅Ψ + k_z⋅Ψ
```

```python
>>> graph(expr)
<graphviz.files.Source at 0x7f7a700a74a8>
```

# Generalizing

```python
>>> from discretizer import split_factors
```

```python
>>> output = []
>>> expr = sympy.expand(expr)
>>> for subexpr in expr.args:
...     output.append(split_factors(subexpr))
```

```python
>>> expr.args
⎛           2                 ⎞
⎝k_z⋅Ψ, C⋅kₓ ⋅k_y⋅Ψ, kₓ⋅A⋅kₓ⋅Ψ⎠
```

```python
>>> output
⎡             ⎛     2       ⎞               ⎤
⎣(1, k_z, Ψ), ⎝C, kₓ ⋅k_y, Ψ⎠, (kₓ⋅A, kₓ, Ψ)⎦
```

# checking powers
yeah, maybe this should be done when split is being done?

```python
>>> from discretizer import get_powers
```

```python
>>> operator = split_factors(subexpr)[1]; operator
kₓ
```

```python
>>> test_operators = [kx, kx**2, ky, kz, kx*ky, kx**2*ky, 2*kx, kx+ky]
...
>>> for operator in test_operators:
...     try:
...         print(get_powers(operator), 'from operator', operator)
...     except ValueError:
...         print('Error on', operator)
(1, 0, 0) from operator k_x
(2, 0, 0) from operator k_x**2
(0, 1, 0) from operator k_y
(0, 0, 1) from operator k_z
(1, 1, 0) from operator k_x*k_y
(2, 1, 0) from operator k_x**2*k_y
Error on 2*k_x
Error on k_x + k_y
```
