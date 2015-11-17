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
>>> from discretizer import discretize_expression
```

### Test1

```python
>>> expr = kx*A*kx + B + 5.0 + A; expr
5.0 + B + kₓ⋅A(x, y, z)⋅kₓ + A(x, y, z)
```

```python
>>> discretize_expression(expr)
                                                        A(-aₓ + x, y, z)⋅Ψ(x, 
B⋅Ψ(x, y, z) + A(x, y, z)⋅Ψ(x, y, z) + 5.0⋅Ψ(x, y, z) + ──────────────────────
                                                                       2      
                                                                   4⋅aₓ       

y, z)   A(-aₓ + x, y, z)⋅Ψ(-2⋅aₓ + x, y, z)   A(aₓ + x, y, z)⋅Ψ(x, y, z)   A(a
───── - ─────────────────────────────────── + ────────────────────────── - ───
                           2                                2                 
                       4⋅aₓ                             4⋅aₓ                  

ₓ + x, y, z)⋅Ψ(2⋅aₓ + x, y, z)
──────────────────────────────
               2              
           4⋅aₓ
```

### Test2

```python
>>> expr = kx*ky; expr
kₓ⋅k_y
```

```python
>>> discretize_expression(expr)
  Ψ(-aₓ + x, -a_y + y, z)   Ψ(-aₓ + x, a_y + y, z)   Ψ(aₓ + x, -a_y + y, z)   
- ─────────────────────── + ────────────────────── + ────────────────────── - 
          4⋅aₓ⋅a_y                 4⋅aₓ⋅a_y                 4⋅aₓ⋅a_y          

Ψ(aₓ + x, a_y + y, z)
─────────────────────
       4⋅aₓ⋅a_y
```

### Test3

```python
>>> expr = kx**2 + kx; expr
       2
kₓ + kₓ
```

```python
>>> discretize_expression(expr)
ⅈ⋅Ψ(-aₓ + x, y, z)   ⅈ⋅Ψ(aₓ + x, y, z)   Ψ(x, y, z)   Ψ(-2⋅aₓ + x, y, z)   Ψ(2
────────────────── - ───────────────── + ────────── - ────────────────── - ───
       2⋅aₓ                 2⋅aₓ               2                2             
                                           2⋅aₓ             4⋅aₓ              

⋅aₓ + x, y, z)
──────────────
       2      
   4⋅aₓ
```

### Test4

```python
>>> expr = B+5; expr
5 + B
```

```python
>>> discretize_expression(expr)
B⋅Ψ(x, y, z) + 5⋅Ψ(x, y, z)
```

### Test5

```python
>>> expr = kx; expr
kₓ
```

```python
>>> discretize_expression(expr)
ⅈ⋅Ψ(-aₓ + x, y, z)   ⅈ⋅Ψ(aₓ + x, y, z)
────────────────── - ─────────────────
       2⋅aₓ                 2⋅aₓ
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
