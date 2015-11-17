# Showcase

## Imports

```python
>>> import sympy
>>> from sympy.interactive import printing
>>> printing.init_printing(use_latex='mathjax')
...
>>> import discretizer
```

## Defining sample expression

```python
>>> # Helpers from discretizer to avoid misunderstanding
... kx, ky, kz = discretizer.momentum_operators
>>> x, y, z = discretizer.coord
...
>>> # Hamiltonian parameters
... A, B, C = sympy.symbols('A B C', commutative=False)
...
>>> # Hamiltonian
... H = kx * A(x, y) * kx + kx * B + C(x,y,z) ;H
kₓ⋅B + kₓ⋅A(x, y)⋅kₓ + C(x, y, z)
```

```python
>>> discretizer.discretize_expression(H)
                        ⅈ⋅B⋅Ψ(-aₓ + x, y, z)   ⅈ⋅B⋅Ψ(aₓ + x, y, z)   A(-aₓ + x
C(x, y, z)⋅Ψ(x, y, z) + ──────────────────── - ─────────────────── + ─────────
                                2⋅aₓ                   2⋅aₓ                   
                                                                              

, y)⋅Ψ(x, y, z)   A(-aₓ + x, y)⋅Ψ(-2⋅aₓ + x, y, z)   A(aₓ + x, y)⋅Ψ(x, y, z)  
─────────────── - ──────────────────────────────── + ─────────────────────── -
    2                              2                              2           
4⋅aₓ                           4⋅aₓ                           4⋅aₓ            

 A(aₓ + x, y)⋅Ψ(2⋅aₓ + x, y, z)
 ──────────────────────────────
                 2             
             4⋅aₓ
```
