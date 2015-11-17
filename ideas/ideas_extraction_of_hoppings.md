---
dictitems:
  kernelspec:
    display_name: dev2
    language: python2
    name: dev2
  language_info:
    codemirror_mode:
      name: ipython
      version: 2
    file_extension: .py
    mimetype: text/x-python
    name: python
    nbconvert_exporter: python
    pygments_lexer: ipython2
    version: 2.7.6
kernelspec:
  display_name: dev2
  language: python2
  name: dev2
language_info:
  codemirror_mode:
    name: ipython
    version: 2
  file_extension: .py
  mimetype: text/x-python
  name: python
  nbconvert_exporter: python
  pygments_lexer: ipython2
  version: 2.7.6
state:
  _allownew: true
---

```python
>>> from sympy.printing.dot import dotprint
>>> from graphviz import Source
>>> graph = lambda x: Source(dotprint(x))
```

```python
>>> import numpy
>>> import sympy
>>> from sympy.interactive import printing
>>> printing.init_printing(use_latex='mathjax')
...
>>> from discretizer import coord, momentum_operators, a
>>> from discretizer import substitute_functions
>>> from discretizer import derivate
>>> from discretizer import split_factors
>>> from discretizer import discretize_expression
```

```python
>>> lattice_constants = sympy.symbols('a_x a_y a_z')
>>> kx, ky, kz = momentum_operators
>>> Psi = sympy.Function('Psi')(*coord)
>>> A = sympy.Function('A')(*coord)
```

# Possible inputs to hopping functions as outputs from discretization

```python
>>> discretize_expression(kx*A*kx + 5)
               A(-aₓ + x, y, z)⋅Ψ(x, y, z)   A(-aₓ + x, y, z)⋅Ψ(-2⋅aₓ + x, y, 
5⋅Ψ(x, y, z) + ─────────────────────────── - ─────────────────────────────────
                              2                                 2             
                          4⋅aₓ                              4⋅aₓ              

z)   A(aₓ + x, y, z)⋅Ψ(x, y, z)   A(aₓ + x, y, z)⋅Ψ(2⋅aₓ + x, y, z)
── + ────────────────────────── - ─────────────────────────────────
                   2                                2              
               4⋅aₓ                             4⋅aₓ
```

```python
>>> discretize_expression(kx*kx)
Ψ(x, y, z)   Ψ(-2⋅aₓ + x, y, z)   Ψ(2⋅aₓ + x, y, z)
────────── - ────────────────── - ─────────────────
      2                2                    2      
  2⋅aₓ             4⋅aₓ                 4⋅aₓ
```

```python
>>> discretize_expression(kx*ky)
  Ψ(-aₓ + x, -a_y + y, z)   Ψ(-aₓ + x, a_y + y, z)   Ψ(aₓ + x, -a_y + y, z)   
- ─────────────────────── + ────────────────────── + ────────────────────── - 
          4⋅aₓ⋅a_y                 4⋅aₓ⋅a_y                 4⋅aₓ⋅a_y          

Ψ(aₓ + x, a_y + y, z)
─────────────────────
       4⋅aₓ⋅a_y
```

```python
>>> test = discretize_expression( kx * A * kx * ky + ky * kx * A * kx)
```

```python
>>> from discretizer import wf
```

```python
>>> from discretizer import extract_hoppings
```

```python
>>> hop = extract_hoppings(test)
>>> hop
⎧               ⅈ⋅A(-aₓ + x, y, z)   ⅈ⋅A(-aₓ + x, -a_y + y, z)              ⅈ⋅
⎪(-2, -1, 0): - ────────────────── - ─────────────────────────, (-2, 1, 0): ──
⎨                       2                        2                            
⎪                   8⋅aₓ ⋅a_y                8⋅aₓ ⋅a_y                        
⎩                                                                             

A(-aₓ + x, y, z)   ⅈ⋅A(-aₓ + x, a_y + y, z)              ⅈ⋅A(-aₓ + x, y, z)   
──────────────── + ────────────────────────, (0, -1, 0): ────────────────── + 
      2                       2                                  2            
  8⋅aₓ ⋅a_y               8⋅aₓ ⋅a_y                          8⋅aₓ ⋅a_y        
                                                                              

ⅈ⋅A(-aₓ + x, -a_y + y, z)   ⅈ⋅A(aₓ + x, y, z)   ⅈ⋅A(aₓ + x, -a_y + y, z)      
───────────────────────── + ───────────────── + ────────────────────────, (0, 
            2                       2                      2                  
        8⋅aₓ ⋅a_y               8⋅aₓ ⋅a_y              8⋅aₓ ⋅a_y              
                                                                              

         ⅈ⋅A(-aₓ + x, y, z)   ⅈ⋅A(-aₓ + x, a_y + y, z)   ⅈ⋅A(aₓ + x, y, z)   ⅈ
1, 0): - ────────────────── - ──────────────────────── - ───────────────── - ─
                 2                       2                       2            
             8⋅aₓ ⋅a_y               8⋅aₓ ⋅a_y               8⋅aₓ ⋅a_y        
                                                                              

⋅A(aₓ + x, a_y + y, z)                ⅈ⋅A(aₓ + x, y, z)   ⅈ⋅A(aₓ + x, -a_y + y
──────────────────────, (2, -1, 0): - ───────────────── - ────────────────────
          2                                   2                      2        
      8⋅aₓ ⋅a_y                           8⋅aₓ ⋅a_y              8⋅aₓ ⋅a_y    
                                                                              

, z)             ⅈ⋅A(aₓ + x, y, z)   ⅈ⋅A(aₓ + x, a_y + y, z)⎫
────, (2, 1, 0): ───────────────── + ───────────────────────⎪
                         2                      2           ⎬
                     8⋅aₓ ⋅a_y              8⋅aₓ ⋅a_y       ⎪
                                                            ⎭
```

```python
>>> from discretizer import shortening
```

```python
>>> shortening(hop)
⎧                  ⎛  aₓ          ⎞      ⎛  aₓ                 ⎞              
⎪               ⅈ⋅A⎜- ── + x, y, z⎟   ⅈ⋅A⎜- ── + x, -a_y + y, z⎟              
⎪                  ⎝  2           ⎠      ⎝  2                  ⎠              
⎨(-1, -1, 0): - ─────────────────── - ──────────────────────────, (-1, 1, 0): 
⎪                        2                        2                           
⎪                    2⋅aₓ ⋅a_y                2⋅aₓ ⋅a_y                       
⎩                                                                             

   ⎛  aₓ          ⎞      ⎛  aₓ                ⎞                 ⎛  aₓ         
ⅈ⋅A⎜- ── + x, y, z⎟   ⅈ⋅A⎜- ── + x, a_y + y, z⎟              ⅈ⋅A⎜- ── + x, y, 
   ⎝  2           ⎠      ⎝  2                 ⎠                 ⎝  2          
─────────────────── + ─────────────────────────, (0, -1, 0): ─────────────────
         2                        2                                   2       
     2⋅aₓ ⋅a_y                2⋅aₓ ⋅a_y                           2⋅aₓ ⋅a_y   
                                                                              

 ⎞      ⎛  aₓ                 ⎞      ⎛aₓ          ⎞      ⎛aₓ                 ⎞
z⎟   ⅈ⋅A⎜- ── + x, -a_y + y, z⎟   ⅈ⋅A⎜── + x, y, z⎟   ⅈ⋅A⎜── + x, -a_y + y, z⎟
 ⎠      ⎝  2                  ⎠      ⎝2           ⎠      ⎝2                  ⎠
── + ────────────────────────── + ───────────────── + ────────────────────────
                 2                        2                      2            
             2⋅aₓ ⋅a_y                2⋅aₓ ⋅a_y              2⋅aₓ ⋅a_y        
                                                                              

                  ⎛  aₓ          ⎞      ⎛  aₓ                ⎞      ⎛aₓ       
               ⅈ⋅A⎜- ── + x, y, z⎟   ⅈ⋅A⎜- ── + x, a_y + y, z⎟   ⅈ⋅A⎜── + x, y
                  ⎝  2           ⎠      ⎝  2                 ⎠      ⎝2        
, (0, 1, 0): - ─────────────────── - ───────────────────────── - ─────────────
                        2                        2                       2    
                    2⋅aₓ ⋅a_y                2⋅aₓ ⋅a_y               2⋅aₓ ⋅a_y
                                                                              

   ⎞      ⎛aₓ                ⎞                   ⎛aₓ          ⎞      ⎛aₓ      
, z⎟   ⅈ⋅A⎜── + x, a_y + y, z⎟                ⅈ⋅A⎜── + x, y, z⎟   ⅈ⋅A⎜── + x, 
   ⎠      ⎝2                 ⎠                   ⎝2           ⎠      ⎝2       
──── - ───────────────────────, (1, -1, 0): - ───────────────── - ────────────
                  2                                   2                      2
              2⋅aₓ ⋅a_y                           2⋅aₓ ⋅a_y              2⋅aₓ 
                                                                              

           ⎞                ⎛aₓ          ⎞      ⎛aₓ                ⎞⎫
-a_y + y, z⎟             ⅈ⋅A⎜── + x, y, z⎟   ⅈ⋅A⎜── + x, a_y + y, z⎟⎪
           ⎠                ⎝2           ⎠      ⎝2                 ⎠⎪
────────────, (1, 1, 0): ───────────────── + ───────────────────────⎬
                                 2                      2           ⎪
⋅a_y                         2⋅aₓ ⋅a_y              2⋅aₓ ⋅a_y       ⎪
                                                                    ⎭
```

```python

```
