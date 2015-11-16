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
>>> kx, ky, kz = momentum_operators
```

```python
>>> derivate(expr, kx)
  B(-aₓ + x, y, z)⋅C⋅A(-aₓ + x, y, z)   B(aₓ + x, y, z)⋅C⋅A(aₓ + x, y, z)
- ─────────────────────────────────── + ─────────────────────────────────
                  2⋅aₓ                                 2⋅aₓ
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

### Test1

```python
>>> expr = expr = kx*A*kx*Psi; expr
kₓ⋅A(x, y, z)⋅kₓ⋅Ψ(x, y, z)
```

```python
>>> recursive(expr)
  A(-aₓ + x, y, z)⋅Ψ(x, y, z)   A(-aₓ + x, y, z)⋅Ψ(-2⋅aₓ + x, y, z)   A(aₓ + x
- ─────────────────────────── + ─────────────────────────────────── - ────────
                 2                                 2                          
             4⋅aₓ                              4⋅aₓ                           

, y, z)⋅Ψ(x, y, z)   A(aₓ + x, y, z)⋅Ψ(2⋅aₓ + x, y, z)
────────────────── + ─────────────────────────────────
      2                                2              
  4⋅aₓ                             4⋅aₓ
```

### Test2

```python
>>> expr = expr = kx*ky*Psi; expr
kₓ⋅k_y⋅Ψ(x, y, z)
```

```python
>>> recursive(expr)
Ψ(-aₓ + x, -a_y + y, z)   Ψ(-aₓ + x, a_y + y, z)   Ψ(aₓ + x, -a_y + y, z)   Ψ(
─────────────────────── - ────────────────────── - ────────────────────── + ──
        4⋅aₓ⋅a_y                 4⋅aₓ⋅a_y                 4⋅aₓ⋅a_y            

aₓ + x, a_y + y, z)
───────────────────
     4⋅aₓ⋅a_y
```

### Test3

```python
>>> expr = expr = kx**2*Psi + kx*Psi; expr
                  2           
kₓ⋅Ψ(x, y, z) + kₓ ⋅Ψ(x, y, z)
```

```python
>>> recursive(expr)
  Ψ(-aₓ + x, y, z)   Ψ(aₓ + x, y, z)   Ψ(x, y, z)   Ψ(-2⋅aₓ + x, y, z)   Ψ(2⋅a
- ──────────────── + ─────────────── - ────────── + ────────────────── + ─────
        2⋅aₓ               2⋅aₓ              2                2               
                                         2⋅aₓ             4⋅aₓ                

ₓ + x, y, z)
────────────
     2      
 4⋅aₓ
```

# Testing split_function

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
⎡                                              ⎛    2                 ⎞       
⎣(1, k_z, Ψ(x, y, z)), (C⋅kₓ, kₓ, Ψ(x, y, z)), ⎝C⋅kₓ , k_y, Ψ(x, y, z)⎠, (k_y⋅

                             ⎤
A(x, y, z), kₓ, B⋅Ψ(x, y, z))⎦
```

```python

```
