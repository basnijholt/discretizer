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
<graphviz.files.Source at 0x7f2fb9d60a90>
```

```python
>>> expr = sympy.expand(expr); expr
    2                  
C⋅kₓ ⋅k_y⋅Ψ + kₓ⋅A⋅kₓ⋅Ψ
```

```python
>>> graph(expr)
<graphviz.files.Source at 0x7f2fb9d66c18>
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
>>> split_factors(subexpr)[1]
  2    
kₓ ⋅k_y
```

# checking powers
yeah, maybe this should be done when split is being done?

```python
>>> operator = split_factors(subexpr)[1]; operator
  2    
kₓ ⋅k_y
```

```python
>>> def get_powers(expr):
...     gens = [sympy.Symbol(s.name) for s in momentum_operators]
...     subs = {c: n for c,n in zip(momentum_operators, gens)}
...
...     expr = sympy.poly(operator.subs(subs), gens)
...     monomials = [(power, coef)
...         for power, coef in zip(expr.monoms(), expr.coeffs())]
...
...     if len(monomials) != 1:
...         raise ValueError('Momentum operator is not built from simple Mul-s.')
...
...     powers, coef = monomials[0]
...     if coef != 1:
...         raise ValueError('Momentum operator contains sth except momentum operators.')
...
...     return powers
```

```python
>>> test_operators = [kx, kx**2, ky, kz, kx*ky, kx**2*ky, 2*kx, kx+ky]
```

```python
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
