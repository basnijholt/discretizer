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
>>> from discretizer.algorithms import _discretize_expression as discretize_expression
>>> from discretizer import momentum_operators
>>> from discretizer import coordinates
```

## Defining input expression

```python
>>> kx, ky, kz = momentum_operators
>>> x, y, z = coordinates
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
>>> from discretizer.algorithms import extract_hoppings
>>> from discretizer.algorithms import shortening
```

```python
>>> wf = sympy.Symbol('Psi')(*coordinates)
```

```python
>>> expr = kx**2 + kx + kx*sympy.sin(x)
>>> expr = expr
>>> expr
                   2
kₓ + kₓ⋅sin(x) + kₓ
```

```python
>>> output = discretize_expression(expr, {'x', 'y', 'z'})
>>> output
⎧             ⅈ    ⅈ⋅sin(a - x)   1              2                 ⅈ    ⅈ⋅sin(
⎪(-1, 0, 0): ─── - ──────────── - ──, (0, 0, 0): ──, (1, 0, 0): - ─── - ──────
⎨            2⋅a       2⋅a         2              2               2⋅a       2⋅
⎪                                 a              a                            
⎩                                                                             

a + x)   1 ⎫
────── - ──⎪
a         2⎬
         a ⎪
           ⎭
```

### Test1

```python
>>> expr = kx*A*kx + B + 5.0 + A; expr
5.0 + B + kₓ⋅A(x, y, z)⋅kₓ + A(x, y, z)
```

```python
>>> discretize_expression(expr, {'x', 'y', 'z'})
⎧              ⎛  a          ⎞                                      ⎛  a      
⎪            -A⎜- ─ + x, y, z⎟                                     A⎜- ─ + x, 
⎪              ⎝  2          ⎠                                      ⎝  2      
⎨(-1, 0, 0): ──────────────────, (0, 0, 0): 5.0 + B + A(x, y, z) + ───────────
⎪                     2                                                    2  
⎪                    a                                                    a   
⎩                                                                             

    ⎞    ⎛a          ⎞               ⎛a          ⎞ ⎫
y, z⎟   A⎜─ + x, y, z⎟             -A⎜─ + x, y, z⎟ ⎪
    ⎠    ⎝2          ⎠               ⎝2          ⎠ ⎪
───── + ──────────────, (1, 0, 0): ────────────────⎬
               2                           2       ⎪
              a                           a        ⎪
                                                   ⎭
```

### Test2

```python
>>> expr = kx*ky; expr
kₓ⋅k_y
```

```python
>>> discretize_expression(expr, {'x', 'y', 'z'})
⎧             -1                 1                 1               -1  ⎫
⎪(-1, -1, 0): ────, (-1, 1, 0): ────, (1, -1, 0): ────, (1, 1, 0): ────⎪
⎨                2                 2                 2                2⎬
⎪             4⋅a               4⋅a               4⋅a              4⋅a ⎪
⎩                                                                      ⎭
```

### Test3

```python
>>> expr = kx**2 + kx; expr
       2
kₓ + kₓ
```

```python
>>> discretize_expression(expr, {'x', 'y', 'z'})
⎧             ⅈ    1              2                 ⅈ    1 ⎫
⎪(-1, 0, 0): ─── - ──, (0, 0, 0): ──, (1, 0, 0): - ─── - ──⎪
⎨            2⋅a    2              2               2⋅a    2⎬
⎪                  a              a                      a ⎪
⎩                                                          ⎭
```

### Test4

```python
>>> expr = B+5; expr
5 + B
```

```python
>>> discretize_expression(expr, {'x', 'y', 'z'})
{(0, 0, 0): 5 + B}
```

### Test5

```python
>>> expr = kx; expr
kₓ
```

```python
>>> discretize_expression(expr, {'x', 'y', 'z'})
⎧             ⅈ              -ⅈ ⎫
⎨(-1, 0, 0): ───, (1, 0, 0): ───⎬
⎩            2⋅a             2⋅a⎭
```

### Test6

```python
>>> wf = sympy.Symbol(discretizer.algorithms.wavefunction_name)(*coordinates)
>>> expr = kx*ky*wf; expr
kₓ⋅k_y⋅Ψ(x, y, z)
```

```python
>>> try:
...     discretize_expression(expr)
>>> except AssertionError as e:
...     print("AssertionError:", e)
```
