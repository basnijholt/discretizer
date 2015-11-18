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
>>> from discretizer.algorithms import discretize_summand
>>> def inputs_for_sebastians_functions(hamiltonian):
...     """
...     """
...     assert wf not in hamiltonian.atoms(sympy.Function), \
...             "Hamiltonian should not contain {}".format(wf)
...     expression = sympy.expand(hamiltonian * wf)
...
...     if expression.func == sympy.Add:
...         summands = expression.args
...     else:
...         summands = [expression]
...
...     outputs = []
...     for summand in summands:
...         outputs.append(discretize_summand(summand))
...     return outputs
```

```python
>>> test = inputs_for_sebastians_functions(kx);test
⎡   ⎛  Ψ(-aₓ + x, y, z)   Ψ(aₓ + x, y, z)⎞⎤
⎢-ⅈ⋅⎜- ──────────────── + ───────────────⎟⎥
⎣   ⎝        2⋅aₓ               2⋅aₓ     ⎠⎦
```

```python
>>> hop = [extract_hoppings(summand) for summand in test]
>>> hop
⎡⎧             ⅈ               -ⅈ  ⎫⎤
⎢⎨(-1, 0, 0): ────, (1, 0, 0): ────⎬⎥
⎣⎩            2⋅aₓ             2⋅aₓ⎭⎦
```

```python
>>> [shortening(summand) for summand in hop]
⎡⎧             ⅈ               -ⅈ  ⎫⎤
⎢⎨(-1, 0, 0): ────, (1, 0, 0): ────⎬⎥
⎣⎩            2⋅aₓ             2⋅aₓ⎭⎦
```

### Test2

```python
>>> test = inputs_for_sebastians_functions(kx*A*kx);test
⎡A(-aₓ + x, y, z)⋅Ψ(x, y, z)   A(-aₓ + x, y, z)⋅Ψ(-2⋅aₓ + x, y, z)   A(aₓ + x,
⎢─────────────────────────── - ─────────────────────────────────── + ─────────
⎢               2                                 2                           
⎣           4⋅aₓ                              4⋅aₓ                            

 y, z)⋅Ψ(x, y, z)   A(aₓ + x, y, z)⋅Ψ(2⋅aₓ + x, y, z)⎤
───────────────── - ─────────────────────────────────⎥
     2                                2              ⎥
 4⋅aₓ                             4⋅aₓ               ⎦
```

```python
>>> hop = [extract_hoppings(summand) for summand in test]
>>> hop
⎡⎧            -A(-aₓ + x, y, z)              A(-aₓ + x, y, z)   A(aₓ + x, y, z
⎢⎪(-2, 0, 0): ──────────────────, (0, 0, 0): ──────────────── + ──────────────
⎢⎨                      2                             2                  2    
⎢⎪                  4⋅aₓ                          4⋅aₓ               4⋅aₓ     
⎣⎩                                                                            

)             -A(aₓ + x, y, z) ⎫⎤
─, (2, 0, 0): ─────────────────⎪⎥
                        2      ⎬⎥
                    4⋅aₓ       ⎪⎥
                               ⎭⎦
```

```python
>>> [shortening(summand) for summand in hop]
⎡⎧              ⎛  aₓ          ⎞               ⎛  aₓ          ⎞    ⎛aₓ        
⎢⎪            -A⎜- ── + x, y, z⎟              A⎜- ── + x, y, z⎟   A⎜── + x, y,
⎢⎪              ⎝  2           ⎠               ⎝  2           ⎠    ⎝2         
⎢⎨(-1, 0, 0): ───────────────────, (0, 0, 0): ───────────────── + ────────────
⎢⎪                      2                              2                  2   
⎢⎪                    aₓ                             aₓ                 aₓ    
⎣⎩                                                                            

  ⎞               ⎛aₓ          ⎞ ⎫⎤
 z⎟             -A⎜── + x, y, z⎟ ⎪⎥
  ⎠               ⎝2           ⎠ ⎪⎥
───, (1, 0, 0): ─────────────────⎬⎥
                         2       ⎪⎥
                       aₓ        ⎪⎥
                                 ⎭⎦
```

### Test3

```python
>>> test = inputs_for_sebastians_functions(kx + A + 5);test
⎡                 ⎛  Ψ(-aₓ + x, y, z)   Ψ(aₓ + x, y, z)⎞                      
⎢5⋅Ψ(x, y, z), -ⅈ⋅⎜- ──────────────── + ───────────────⎟, A(x, y, z)⋅Ψ(x, y, z
⎣                 ⎝        2⋅aₓ               2⋅aₓ     ⎠                      

 ⎤
)⎥
 ⎦
```

```python
>>> hop = [extract_hoppings(summand) for summand in test]
>>> hop
⎡                ⎧             ⅈ               -ⅈ  ⎫                         ⎤
⎢{(0, 0, 0): 5}, ⎨(-1, 0, 0): ────, (1, 0, 0): ────⎬, {(0, 0, 0): A(x, y, z)}⎥
⎣                ⎩            2⋅aₓ             2⋅aₓ⎭                         ⎦
```

```python
>>> [shortening(summand) for summand in hop]
⎡                ⎧             ⅈ               -ⅈ  ⎫                         ⎤
⎢{(0, 0, 0): 5}, ⎨(-1, 0, 0): ────, (1, 0, 0): ────⎬, {(0, 0, 0): A(x, y, z)}⎥
⎣                ⎩            2⋅aₓ             2⋅aₓ⎭                         ⎦
```

```python

```
