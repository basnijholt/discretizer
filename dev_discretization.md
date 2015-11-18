## Assumptions
(read ideas_basic_operations)

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
>>> kx, ky, kz = momentum_operators
>>> x, y, z = coord
...
>>> A, B = sympy.symbols('A B', commutative=False)
>>> A = A(x,y,z)
```

# Designing recursive algorithm
We must somehow handle this expression recursively
1. First expand
2. loop over summands to keep them separated.
3. For every summand launch recursive algorithm that discretize it (point 2 above)
short hoppings
4. gather results from summands (after the loop from 2).

```python
>>> from discretizer.algorithms import discretize_expression
>>> from discretizer.algorithms import wf
...
>>> from discretizer.algorithms import extract_hoppings
>>> from discretizer.algorithms import shortening
```

```python
>>> expr = kx**2 + kx; expr
       2
kₓ + kₓ
```

```python
>>> A/2/2
A(x, y, z)
──────────
    4
```

```python
>>> ax, ay, az = discretizer.algorithms.lattice_constants
>>> tmp = 1/ax/2
```

```python
>>> import numpy as np
```

```python
>>> tmp.subs(ax, ax/np.int64(1))
0.5
───
 aₓ
```

```python
>>> output = discretize_expression(expr)
>>> output
⎧             ⅈ      1               2                 ⅈ      1 ⎫
⎪(-1, 0, 0): ──── - ───, (0, 0, 0): ───, (1, 0, 0): - ──── - ───⎪
⎨            2⋅aₓ     2               2               2⋅aₓ     2⎬
⎪                   aₓ              aₓ                       aₓ ⎪
⎩                                                               ⎭
```

### Test1

```python
>>> expr = kx*A*kx + B + 5.0 + A; expr
5.0 + B + kₓ⋅A(x, y, z)⋅kₓ + A(x, y, z)
```

```python
>>> discretize_expression(expr)
⎧              ⎛  aₓ          ⎞                                      ⎛  aₓ    
⎪            -A⎜- ── + x, y, z⎟                                     A⎜- ── + x
⎪              ⎝  2           ⎠                                      ⎝  2     
⎨(-1, 0, 0): ───────────────────, (0, 0, 0): 5.0 + B + A(x, y, z) + ──────────
⎪                      2                                                     2
⎪                    aₓ                                                    aₓ 
⎩                                                                             

      ⎞    ⎛aₓ          ⎞               ⎛aₓ          ⎞ ⎫
, y, z⎟   A⎜── + x, y, z⎟             -A⎜── + x, y, z⎟ ⎪
      ⎠    ⎝2           ⎠               ⎝2           ⎠ ⎪
─────── + ───────────────, (1, 0, 0): ─────────────────⎬
                  2                            2       ⎪
                aₓ                           aₓ        ⎪
                                                       ⎭
```

### Test2

```python
>>> expr = kx*ky; expr
kₓ⋅k_y
```

```python
>>> discretize_expression(expr)
⎧               -1                     1                     1                
⎨(-1, -1, 0): ────────, (-1, 1, 0): ────────, (1, -1, 0): ────────, (1, 1, 0):
⎩             4⋅aₓ⋅a_y              4⋅aₓ⋅a_y              4⋅aₓ⋅a_y            

   -1    ⎫
 ────────⎬
 4⋅aₓ⋅a_y⎭
```

### Test3

```python
>>> expr = kx**2 + kx; expr
       2
kₓ + kₓ
```

```python
>>> discretize_expression(expr)
⎧             ⅈ      1               2                 ⅈ      1 ⎫
⎪(-1, 0, 0): ──── - ───, (0, 0, 0): ───, (1, 0, 0): - ──── - ───⎪
⎨            2⋅aₓ     2               2               2⋅aₓ     2⎬
⎪                   aₓ              aₓ                       aₓ ⎪
⎩                                                               ⎭
```

### Test4

```python
>>> expr = B+5; expr
5 + B
```

```python
>>> discretize_expression(expr)
{(0, 0, 0): 5 + B}
```

### Test5

```python
>>> expr = kx; expr
kₓ
```

```python
>>> discretize_expression(expr)
⎧             ⅈ               -ⅈ  ⎫
⎨(-1, 0, 0): ────, (1, 0, 0): ────⎬
⎩            2⋅aₓ             2⋅aₓ⎭
```

### Test6

```python
>>> expr = kx*ky*discretizer.algorithms.wf; expr
kₓ⋅k_y⋅Ψ(x, y, z)
```

```python
>>> try:
...     discretize_expression(expr)
>>> except AssertionError as e:
...     print("AssertionError:", e)
AssertionError: Hamiltonian should not contain Psi(x, y, z)
```

```python

```
