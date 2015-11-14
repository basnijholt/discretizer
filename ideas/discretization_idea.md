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
>>> from discretizer import substitute_functions
>>> from discretizer import derivate
>>> from discretizer import split_factors
>>> from discretizer import get_powers
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
>>> from discretizer import recursive
```

```python
>>> kx, ky, kz = momentum_operators
>>> Psi = sympy.Function('Psi')(*coord)
>>> A = sympy.Function('A')(*coord)
```

```python
>>> expr = expr = kx*A*kx*Psi + kx*Psi#+ kz**8*kx*ky*A*kx*Psi + ky**3*kx*ky*A*kx*Psi
>>> expr
kₓ⋅A(x, y, z)⋅kₓ⋅Ψ(x, y, z) + kₓ⋅Ψ(x, y, z)
```

```python
>>> %time recursive(expr)
CPU times: user 32 ms, sys: 4 ms, total: 36 ms
Wall time: 36.3 ms
  - -0.5⋅ⅈ⋅Ψ(-a + x, y, z)    0.5⋅ⅈ⋅Ψ(a + x, y, z)   0.25⋅A(-a + x, y, z)⋅Ψ(x,
- ───────────────────────── + ──────────────────── + ─────────────────────────
              a                        a                             2        
                                                                    a         

 y, z)   - -0.25⋅A(-a + x, y, z)⋅Ψ(-2⋅a + x, y, z)    0.25⋅A(a + x, y, z)⋅Ψ(x,
────── - ────────────────────────────────────────── + ────────────────────────
                              2                                      2        
                             a                                      a         

 y, z)   - -0.25⋅A(a + x, y, z)⋅Ψ(2⋅a + x, y, z) 
────── - ────────────────────────────────────────
                             2                   
                            a
```

# Test

```python
>>> from discretizer import tb_coefs_1d
```

```python
>>> n = 3
```

```python
>>> tb_coefs_1d(n)
⎡-0.5⋅ⅈ   1.0⋅ⅈ     -1.0⋅ⅈ   0.5⋅ⅈ⎤
⎢───────, ─────, 0, ───────, ─────⎥
⎢    3       3          3       3 ⎥
⎣   a       a          a       a  ⎦
```

```python
>>> recursive(kx**n*Psi)
0.5⋅ⅈ⋅Ψ(-2⋅a + x, y, z)   - -1.0⋅ⅈ⋅Ψ(-a + x, y, z)    1.0⋅ⅈ⋅Ψ(a + x, y, z)   -
─────────────────────── - ───────────────────────── + ──────────────────── - ─
            3                          3                        3             
           a                          a                        a              

 -0.5⋅ⅈ⋅Ψ(2⋅a + x, y, z) 
─────────────────────────
            3            
           a
```

# Playing around

```python
>>> expr = kx*A*kx + C * kx**2 * ky
>>> expr = expr * Psi
>>> expr = sympy.expand(expr); expr
    2                                             
C⋅kₓ ⋅k_y⋅Ψ(x, y, z) + kₓ⋅A(x, y, z)⋅kₓ⋅Ψ(x, y, z)
```

# Testing split_function

```python
>>> output = []
>>> expr = sympy.expand(expr)
>>> for subexpr in expr.args:
...     output.append(split_factors(subexpr))
```

```python
>>> test_operators = [kz*Psi, C*kx**2*Psi, C*kx**2*ky*Psi, ky*A*kx*B*Psi, kx, kx**2, Psi]
...
>>> tested = []
>>> output = []
>>> for operator in test_operators:
...     try:
...         output.append(split_factors(operator))
...         tested.append(operator)
...     except ValueError:
...         print('ValueError on', operator)
...     except AssertionError:
...         print('AssertionError on', operator, operator.func)
AssertionError on k_x <class 'sympy.core.symbol.Symbol'>
AssertionError on k_x**2 <class 'sympy.core.power.Pow'>
AssertionError on Psi(x, y, z) Psi
```

```python
>>> tested
⎡                    2                 2                                      
⎣k_z⋅Ψ(x, y, z), C⋅kₓ ⋅Ψ(x, y, z), C⋅kₓ ⋅k_y⋅Ψ(x, y, z), k_y⋅A(x, y, z)⋅kₓ⋅B⋅Ψ

         ⎤
(x, y, z)⎦
```

```python
>>> output
⎡                      ⎛     2            ⎞  ⎛     2                ⎞         
⎣(1, k_z, Ψ(x, y, z)), ⎝C, kₓ , Ψ(x, y, z)⎠, ⎝C, kₓ ⋅k_y, Ψ(x, y, z)⎠, (k_y⋅A(

                           ⎤
x, y, z), kₓ, B⋅Ψ(x, y, z))⎦
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
