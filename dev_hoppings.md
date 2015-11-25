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

# Possible inputs to hopping functions as outputs from discretization

```python
>>> %reset -f
```

```python
>>> from discretizer import discretize
>>> import sympy
```

```python
>>> from discretizer.algorithms import wavefunction_name
>>> from discretizer.algorithms import _discretize_summand
>>> from discretizer.algorithms import coordinates
```

```python
>>> from discretizer.algorithms import read_hopping_from_wf
```

```python
>>> def inputs_for_sebastians_functions(hamiltonian, coordinates):
...     wf = sympy.Function(wavefunction_name)(*coordinates)
...     expression = sympy.expand(hamiltonian * wf)
...
...     if expression.func == sympy.Add:
...         summands = expression.args
...     else:
...         summands = [expression]
...
...     outputs = []
...     for summand in summands:
...         outputs.append(_discretize_summand(summand))
...     return outputs
```

```python
>>> x,y,z = sympy.symbols('x y z', commutative=False)
>>> ax, ay, az = sympy.symbols('a_x a_y a_z')
...
>>> wf = sympy.Function(wavefunction_name)
>>> psi = wf(x, z)
```

```python
>>> expr = psi.subs({x: x-2*ax, y: y-az}); expr
Ψ(-2⋅aₓ + x, z)
```

```python
>>> read_hopping_from_wf(expr)
(-2, 0)
```

### Test1

```python
>>> import discretizer
```

```python
>>> from discretizer.algorithms import extract_hoppings
>>> from discretizer.algorithms import shortening
```

```python
>>> kx, ky, kz = discretizer.momentum_operators
```

```python
>>> test = inputs_for_sebastians_functions(kx, coordinates=(x,));test
⎡   ⎛  Ψ(-aₓ + x)   Ψ(aₓ + x)⎞⎤
⎢-ⅈ⋅⎜- ────────── + ─────────⎟⎥
⎣   ⎝     2⋅aₓ         2⋅aₓ  ⎠⎦
```

```python
>>> hop = [extract_hoppings(summand) for summand in test]
>>> hop
⎡⎧        ⅈ          -ⅈ  ⎫⎤
⎢⎨(-1,): ────, (1,): ────⎬⎥
⎣⎩       2⋅aₓ        2⋅aₓ⎭⎦
```

```python
>>> [shortening(summand) for summand in hop]
⎡⎧        ⅈ         -ⅈ ⎫⎤
⎢⎨(-1,): ───, (1,): ───⎬⎥
⎣⎩       2⋅a        2⋅a⎭⎦
```

### Test2

```python
>>> # A = sympy.Function('A')(x, y,z)
... # test = inputs_for_sebastians_functions(kz*A*kz, coordinates=[x, y, z]);test
...
... A = sympy.Function('A')(y,z)
>>> test = inputs_for_sebastians_functions(kz*A*kz, coordinates=[y, z]);test
⎡A(y, -a_z + z)⋅Ψ(y, z)   A(y, -a_z + z)⋅Ψ(y, -2⋅a_z + z)   A(y, a_z + z)⋅Ψ(y,
⎢────────────────────── - ─────────────────────────────── + ──────────────────
⎢             2                             2                            2    
⎣        4⋅a_z                         4⋅a_z                        4⋅a_z     

 z)   A(y, a_z + z)⋅Ψ(y, 2⋅a_z + z)⎤
─── - ─────────────────────────────⎥
                       2           ⎥
                  4⋅a_z            ⎦
```

```python
>>> hop = [extract_hoppings(summand) for summand in test]
>>> hop
⎡⎧         -A(y, -a_z + z)           A(y, -a_z + z)   A(y, a_z + z)          -
⎢⎪(0, -2): ────────────────, (0, 0): ────────────── + ─────────────, (0, 2): ─
⎢⎨                   2                        2                2              
⎢⎪              4⋅a_z                    4⋅a_z            4⋅a_z               
⎣⎩                                                                            

A(y, a_z + z) ⎫⎤
──────────────⎪⎥
         2    ⎬⎥
    4⋅a_z     ⎪⎥
              ⎭⎦
```

```python
>>> [shortening(summand, discrete_coordinates=['y', 'z']) for summand in hop]
⎡⎧           ⎛     a    ⎞            ⎛     a    ⎞    ⎛   a    ⎞            ⎛  
⎢⎪         -A⎜y, - ─ + z⎟           A⎜y, - ─ + z⎟   A⎜y, ─ + z⎟          -A⎜y,
⎢⎪           ⎝     2    ⎠            ⎝     2    ⎠    ⎝   2    ⎠            ⎝  
⎢⎨(0, -1): ───────────────, (0, 0): ───────────── + ───────────, (0, 1): ─────
⎢⎪                 2                       2              2                   
⎢⎪                a                       a              a                    
⎣⎩                                                                            

 a    ⎞ ⎫⎤
 ─ + z⎟ ⎪⎥
 2    ⎠ ⎪⎥
────────⎬⎥
  2     ⎪⎥
 a      ⎪⎥
        ⎭⎦
```

### Test3

```python
>>> A = sympy.Function('A')(x, y,z)
>>> test = inputs_for_sebastians_functions(kz*A*kz, coordinates=[x, y, z]);test
⎡A(x, y, -a_z + z)⋅Ψ(x, y, z)   A(x, y, -a_z + z)⋅Ψ(x, y, -2⋅a_z + z)   A(x, y
⎢──────────────────────────── - ───────────────────────────────────── + ──────
⎢                2                                   2                        
⎣           4⋅a_z                               4⋅a_z                         

, a_z + z)⋅Ψ(x, y, z)   A(x, y, a_z + z)⋅Ψ(x, y, 2⋅a_z + z)⎤
───────────────────── - ───────────────────────────────────⎥
          2                                 2              ⎥
     4⋅a_z                             4⋅a_z               ⎦
```

```python
>>> hop = [extract_hoppings(summand) for summand in test]
>>> hop
⎡⎧            -A(x, y, -a_z + z)              A(x, y, -a_z + z)   A(x, y, a_z 
⎢⎪(0, 0, -2): ───────────────────, (0, 0, 0): ───────────────── + ────────────
⎢⎨                        2                              2                  2 
⎢⎪                   4⋅a_z                          4⋅a_z              4⋅a_z  
⎣⎩                                                                            

+ z)             -A(x, y, a_z + z) ⎫⎤
────, (0, 0, 2): ──────────────────⎪⎥
                            2      ⎬⎥
                       4⋅a_z       ⎪⎥
                                   ⎭⎦
```

```python
>>> [shortening(summand) for summand in hop]
⎡⎧              ⎛        a    ⎞               ⎛        a    ⎞    ⎛      a    ⎞
⎢⎪            -A⎜x, y, - ─ + z⎟              A⎜x, y, - ─ + z⎟   A⎜x, y, ─ + z⎟
⎢⎪              ⎝        2    ⎠               ⎝        2    ⎠    ⎝      2    ⎠
⎢⎨(0, 0, -1): ──────────────────, (0, 0, 0): ──────────────── + ──────────────
⎢⎪                     2                             2                 2      
⎢⎪                    a                             a                 a       
⎣⎩                                                                            

               ⎛      a    ⎞ ⎫⎤
             -A⎜x, y, ─ + z⎟ ⎪⎥
               ⎝      2    ⎠ ⎪⎥
, (0, 0, 1): ────────────────⎬⎥
                     2       ⎪⎥
                    a        ⎪⎥
                             ⎭⎦
```

### Test3

```python
>>> A = sympy.Function('A')(x, y)
>>> test = inputs_for_sebastians_functions(kx + A + 5, coordinates=[x, y]);test
⎡              ⎛  Ψ(-aₓ + x, y)   Ψ(aₓ + x, y)⎞                 ⎤
⎢5⋅Ψ(x, y), -ⅈ⋅⎜- ───────────── + ────────────⎟, A(x, y)⋅Ψ(x, y)⎥
⎣              ⎝       2⋅aₓ           2⋅aₓ    ⎠                 ⎦
```

```python
>>> hop = [extract_hoppings(summand) for summand in test]
>>> hop
⎡             ⎧          ⅈ            -ⅈ  ⎫                   ⎤
⎢{(0, 0): 5}, ⎨(-1, 0): ────, (1, 0): ────⎬, {(0, 0): A(x, y)}⎥
⎣             ⎩         2⋅aₓ          2⋅aₓ⎭                   ⎦
```
