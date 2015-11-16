## Assumptions
(must be updated)

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
>>> lattice_constants = sympy.symbols('a_x a_y a_z', commutative=False)
```

```python
>>> kx, ky, kz = momentum_operators
```

```python
>>> derivate(expr, kx)
   ⎛  B(-aₓ + x, y, z)⋅C⋅A(-aₓ + x, y, z)   B(aₓ + x, y, z)⋅C⋅A(aₓ + x, y, z)⎞
-ⅈ⋅⎜- ─────────────────────────────────── + ─────────────────────────────────⎟
   ⎝                  2⋅aₓ                                 2⋅aₓ              ⎠
```

# expanding

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
>>> expr = kx*A*kx + C * kx**2 * ky + kz + ky*kx**2
>>> expr = expr * Psi
>>> graph(expr)
<graphviz.files.Source at 0x7f2f5f91ea90>
```

```python
>>> expr = sympy.expand(expr); expr
    2                           2          
C⋅kₓ ⋅k_y⋅Ψ + kₓ⋅A⋅kₓ⋅Ψ + k_y⋅kₓ ⋅Ψ + k_z⋅Ψ
```

```python
>>> graph(expr)
<graphviz.files.Source at 0x7f2f5f92a518>
```

# spliting into lhs, operators, rhs

## Check of function once it works

```python
>>> from discretizer import split_factors
```

```python
>>> test_operators = [kz*Psi, C*kx**2*Psi, C*kx**2*ky*Psi, ky*A*kx*B*Psi, kx, kx**2, Psi,
...                   Psi(*coord)]
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
```

```python
>>> tested
⎡           2        2                            2               ⎤
⎣k_z⋅Ψ, C⋅kₓ ⋅Ψ, C⋅kₓ ⋅k_y⋅Ψ, k_y⋅A⋅kₓ⋅B⋅Ψ, kₓ, kₓ , Ψ, Ψ(x, y, z)⎦
```

```python
>>> output
⎡                            ⎛    2        ⎞                                  
⎣(1, k_z, Ψ), (C⋅kₓ, kₓ, Ψ), ⎝C⋅kₓ , k_y, Ψ⎠, (k_y⋅A, kₓ, B⋅Ψ), (1, kₓ, 1), (k

                                        ⎤
ₓ, kₓ, 1), (1, 1, Ψ), (1, 1, Ψ(x, y, z))⎦
```
