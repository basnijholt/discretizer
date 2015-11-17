## Assumptions
(must be updated)

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
>>> import discretizer
>>> from discretizer import discretize_expression
>>> from discretizer import momentum_operators
>>> from discretizer import coord
```

## Defining input expression

```python
>>> A, B, C, D = sympy.symbols('A B C D', commutative=False)
>>> kx, ky, kz = momentum_operators
>>> x, y, z = coord
```

```python
>>> expr = B(x,y,z)*C*A(x,y,z) + D; expr
D + B(x, y, z)⋅C⋅A(x, y, z)
```

## Calculation derivation

```python
>>> from discretizer.algorithms import derivate
```

```python
>>> lattice_constants = sympy.symbols('a_x a_y a_z', commutative=False)
```

```python
>>> kx, ky, kz = momentum_operators
```

```python
>>> derivate(A(x, y), ky)
   ⎛  A(x, -a_y + y)   A(x, a_y + y)⎞
-ⅈ⋅⎜- ────────────── + ─────────────⎟
   ⎝      2⋅a_y            2⋅a_y    ⎠
```

# expanding

```python
>>> Psi = discretizer.algorithms.wf
```

```python
>>> expr = kx*A*kx + C * kx**2 * ky + kz + ky*kx**2
>>> expr = expr * Psi
>>> graph(expr)
<graphviz.files.Source at 0x7ff21d446780>
```

```python
>>> expr = sympy.expand(expr); expr
    2                                             2                           
C⋅kₓ ⋅k_y⋅Ψ(x, y, z) + kₓ⋅A⋅kₓ⋅Ψ(x, y, z) + k_y⋅kₓ ⋅Ψ(x, y, z) + k_z⋅Ψ(x, y, z

 
)
```

```python
>>> graph(expr)
<graphviz.files.Source at 0x7ff21d454b70>
```

# spliting into lhs, operators, rhs

## Check of function once it works

```python
>>> from discretizer.algorithms import split_factors
```

```python
>>> import numpy as np
```

```python
>>> test_operators = [kz*Psi, C*kx**2*Psi, C*kx**2*ky*Psi, ky*A*kx*B*Psi, kx, kx**2, A, Psi,
...                   3, 5.0, np.float(5), np.int(3), sympy.Integer(5), sympy.Float(5)]
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
⎡                    2                 2                                      
⎣k_z⋅Ψ(x, y, z), C⋅kₓ ⋅Ψ(x, y, z), C⋅kₓ ⋅k_y⋅Ψ(x, y, z), k_y⋅A⋅kₓ⋅B⋅Ψ(x, y, z)

        2                                       ⎤
, kₓ, kₓ , A, Ψ(x, y, z), 3, 5.0, 5.0, 3, 5, 5.0⎦
```

```python
>>> output
⎡                                              ⎛    2                 ⎞       
⎣(1, k_z, Ψ(x, y, z)), (C⋅kₓ, kₓ, Ψ(x, y, z)), ⎝C⋅kₓ , k_y, Ψ(x, y, z)⎠, (k_y⋅

                                                                              
A, kₓ, B⋅Ψ(x, y, z)), (1, kₓ, 1), (kₓ, kₓ, 1), (1, 1, A), (1, 1, Ψ(x, y, z)), 

                                                                      ⎤
(1, 1, 3), (1, 1, 5.0), (1, 1, 5.0), (1, 1, 3), (1, 1, 5), (1, 1, 5.0)⎦
```
