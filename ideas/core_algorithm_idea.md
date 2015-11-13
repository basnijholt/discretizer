## Assumptions
* space_dependent - these are cast to functions, not commutative
* momentum_symbols - will be specified (default to kx, ky, kz probably) and not commutative
* constant - commutative, everything else
* coordinates are always create as symbol with commutative=False
* lattice constant is always a

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

## Substituting functions

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
>>> from sympy.printing.dot import dotprint
>>> from graphviz import Source
>>> graph = lambda x: Source(dotprint(x))
```

```python
>>> kx, ky, kz = momentum_operators
>>> Psi = sympy.Symbol('Psi', commutative=False) #sympy.Function('Psi')(*coord)
```

```python
>>> expr = kx*A*kx + C * kx**2 * ky
>>> expr = expr * Psi
>>> graph(expr)
<graphviz.files.Source at 0x7fcb55a8ca58>
```

```python
>>> expr = sympy.expand(expr); expr
    2                  
C⋅kₓ ⋅k_y⋅Ψ + kₓ⋅A⋅kₓ⋅Ψ
```

```python
>>> graph(expr)
<graphviz.files.Source at 0x7fcb55a93c18>
```

# Generalizing

```python
>>> from discretizer import split_factors
```

```python
>>> subexpr = expr.args[0]; subexpr
    2      
C⋅kₓ ⋅k_y⋅Ψ
```

```python
>>> split_factors(subexpr)
⎛     2       ⎞
⎝C, kₓ ⋅k_y, Ψ⎠
```
