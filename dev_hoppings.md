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
state:
  _allownew: true
---

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
```

```python
>>> kx, ky, kz = discretizer.algorithms.momentum_operators
>>> ax, ay, az = discretizer.algorithms.lattice_constants
>>> x, y, z = discretizer.algorithms.coord
...
>>> A, B = sympy.symbols('A B', commutative=False)
>>> A = A(x,y,z)
```

# Possible inputs to hopping functions as outputs from discretization

```python
>>> from discretizer import discretize_expression
```

```python
>>> from discretizer.algorithms import read_hopping_from_wf
>>> from discretizer.algorithms import extract_hoppings
>>> from discretizer.algorithms import shortening
>>> from discretizer.algorithms import wf
```

```python
>>> expr = wf.subs({x: x-2*ax, y: y-ay}); expr
Ψ(-2⋅aₓ + x, -a_y + y, z)
```

```python
>>> read_hopping_from_wf(expr)
(-2, -1, 0)
```

### Test1

```python
>>> test = discretize_expression(kx);test
ⅈ⋅Ψ(-aₓ + x, y, z)   ⅈ⋅Ψ(aₓ + x, y, z)
────────────────── - ─────────────────
       2⋅aₓ                 2⋅aₓ
```

```python
>>> hop = extract_hoppings(test); hop
⎧             ⅈ               -ⅈ  ⎫
⎨(-1, 0, 0): ────, (1, 0, 0): ────⎬
⎩            2⋅aₓ             2⋅aₓ⎭
```

```python
>>> shortening(hop)
⎧            0.5⋅ⅈ             -0.5⋅ⅈ ⎫
⎨(-1, 0, 0): ─────, (1, 0, 0): ───────⎬
⎩              aₓ                 aₓ  ⎭
```

### Test2

```python
>>> test = discretize_expression(kx*A*kx);test
A(-aₓ + x, y, z)⋅Ψ(x, y, z)   A(-aₓ + x, y, z)⋅Ψ(-2⋅aₓ + x, y, z)   A(aₓ + x, 
─────────────────────────── - ─────────────────────────────────── + ──────────
               2                                 2                            
           4⋅aₓ                              4⋅aₓ                             

y, z)⋅Ψ(x, y, z)   A(aₓ + x, y, z)⋅Ψ(2⋅aₓ + x, y, z)
──────────────── - ─────────────────────────────────
    2                                2              
4⋅aₓ                             4⋅aₓ
```

```python
>>> hop = extract_hoppings(test); hop
⎧            -A(-aₓ + x, y, z)              A(-aₓ + x, y, z)   A(aₓ + x, y, z)
⎪(-2, 0, 0): ──────────────────, (0, 0, 0): ──────────────── + ───────────────
⎨                      2                             2                  2     
⎪                  4⋅aₓ                          4⋅aₓ               4⋅aₓ      
⎩                                                                             

             -A(aₓ + x, y, z) ⎫
, (2, 0, 0): ─────────────────⎪
                       2      ⎬
                   4⋅aₓ       ⎪
                              ⎭
```

```python
>>> shortening(hop)
⎧              ⎛  aₓ          ⎞               ⎛  aₓ          ⎞    ⎛aₓ         
⎪            -A⎜- ── + x, y, z⎟              A⎜- ── + x, y, z⎟   A⎜── + x, y, 
⎪              ⎝  2           ⎠               ⎝  2           ⎠    ⎝2          
⎨(-1, 0, 0): ───────────────────, (0, 0, 0): ───────────────── + ─────────────
⎪                      2                              2                  2    
⎪                    aₓ                             aₓ                 aₓ     
⎩                                                                             

 ⎞               ⎛aₓ          ⎞ ⎫
z⎟             -A⎜── + x, y, z⎟ ⎪
 ⎠               ⎝2           ⎠ ⎪
──, (1, 0, 0): ─────────────────⎬
                        2       ⎪
                      aₓ        ⎪
                                ⎭
```

### Test3

```python
>>> test = discretize_expression(kx + A + 5);test
                                       ⅈ⋅Ψ(-aₓ + x, y, z)   ⅈ⋅Ψ(aₓ + x, y, z)
A(x, y, z)⋅Ψ(x, y, z) + 5⋅Ψ(x, y, z) + ────────────────── - ─────────────────
                                              2⋅aₓ                 2⋅aₓ
```

```python
>>> hop = extract_hoppings(test); hop
⎧             ⅈ                                          -ⅈ  ⎫
⎨(-1, 0, 0): ────, (0, 0, 0): 5 + A(x, y, z), (1, 0, 0): ────⎬
⎩            2⋅aₓ                                        2⋅aₓ⎭
```

```python
>>> shortening(hop)
⎧            0.5⋅ⅈ                                        -0.5⋅ⅈ ⎫
⎨(-1, 0, 0): ─────, (0, 0, 0): 5 + A(x, y, z), (1, 0, 0): ───────⎬
⎩              aₓ                                            aₓ  ⎭
```

```python

```
