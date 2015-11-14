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
>>> from discretizer import recursive
```

```python
>>> kx, ky, kz = momentum_operators
>>> Psi = sympy.Function('Psi')(*coord)
>>> A = sympy.Function('A')(*coord)
```

# Possible inputs to hopping functions as outputs from discretization

```python
>>> recursive(kx*A*kx*Psi)
0.25⋅A(-a + x, y, z)⋅Ψ(x, y, z)   - -0.25⋅A(-a + x, y, z)⋅Ψ(-2⋅a + x, y, z)   
─────────────────────────────── - ────────────────────────────────────────── +
                2                                      2                      
               a                                      a                       

 0.25⋅A(a + x, y, z)⋅Ψ(x, y, z)   - -0.25⋅A(a + x, y, z)⋅Ψ(2⋅a + x, y, z) 
 ────────────────────────────── - ────────────────────────────────────────
                2                                     2                   
               a                                     a
```

```python
>>> recursive(kx*kx*Psi)
2.0⋅Ψ(x, y, z)   - -1.0⋅Ψ(-a + x, y, z)    - -1.0⋅Ψ(a + x, y, z) 
────────────── - ─────────────────────── - ──────────────────────
       2                     2                        2          
      a                     a                        a
```

```python
>>> recursive(kx*ky*Psi)
  - -0.25⋅Ψ(-a + x, -a + y, z)    0.25⋅Ψ(-a + x, a + y, z)   0.25⋅Ψ(a + x, -a 
- ───────────────────────────── + ──────────────────────── + ─────────────────
                 2                            2                          2    
                a                            a                          a     

+ y, z)   - -0.25⋅Ψ(a + x, a + y, z) 
─────── - ───────────────────────────
                        2            
                       a
```

```python
>>> recursive( kx * A * ky * Psi ) # yeah, I know it's not hermitian but it is only one summand
  - -0.25⋅A(-a + x, y, z)⋅Ψ(-a + x, -a + y, z)    0.25⋅A(-a + x, y, z)⋅Ψ(-a + 
- ───────────────────────────────────────────── + ────────────────────────────
                         2                                            2       
                        a                                            a        

x, a + y, z)   0.25⋅A(a + x, y, z)⋅Ψ(a + x, -a + y, z)   - -0.25⋅A(a + x, y, z
──────────── + ─────────────────────────────────────── - ─────────────────────
                                   2                                          
                                  a                                          a

)⋅Ψ(a + x, a + y, z) 
─────────────────────
2
```
