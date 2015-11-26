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
>>> #test = inputs_for_sebastians_functions(kz*A*kz, coordinates=[x, y, z]);test
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
>>> #test = inputs_for_sebastians_functions(kx + A + 5, coordinates=[x, y]);test
```

```python
>>> hop = [extract_hoppings(summand) for summand in test]
>>> hop
⎡             ⎧          ⅈ            -ⅈ  ⎫                   ⎤
⎢{(0, 0): 5}, ⎨(-1, 0): ────, (1, 0): ────⎬, {(0, 0): A(x, y)}⎥
⎣             ⎩         2⋅aₓ          2⋅aₓ⎭                   ⎦
```

```python

```

```python
>>> def search_Function_in_summand(expr, path):
...     if expr.func = sympy.Mul:
...         for i in np.range(len(expr.args)):
...             if isinstance(expr.args[i], sympy.Function):
...                 path.append(i)
...                 return path
...             else:
...                 return 'None'
...     elif isinstance(expr.args[i], sympy.Function):
...         return path
...     else:
...         return 'None'
```

```python
>>> def interpolate_Function(expr):
...     path = 'None'
...     factor = 'None'
...     interupt_recursion = True
...     for i in np.arange(len(expr.args)):
...         argument = sympy.expand(expr.args[i])
...         if argument.func == sympy.Add:
...             for j in np.arange(len(argument.args)):
...                 summand = argument.args[j]
...                 if summand.func == sympy.Mul:
...                     for k in np.arange(len(summand.args)):
...                         temp = 0
...                         if summand.args[k] == sympy.Symbol('a'):
...                             temp = sympy.Mul(sympy.Mul(*summand.args[:k]),
...                                              sympy.Mul(*summand.args[k+1:]))
...                             #print(temp)
...                         if not temp == int(temp):
...                             if not factor == 'None':
...                                 interupt_recursion = False
...                             #print('found one')
...                             factor = (temp)
...                             path = np.array([i, j, k])
...     sign = np.sign(factor)
...     offsets = np.array([int(factor), sign*(int(sign * factor) + 1)])
...     weights = 1/np.abs(offsets-factor)
...     weights = weights/np.sum(weights)
...     return offsets, weights
```

```python
>>> a = sympy.Symbol('a')
>>> test = A.subs(y, y-2*a/5)#.subs(x, x+a/2)
>>> test
 ⎛     2⋅a    ⎞
A⎜x, - ─── + y⎟
 ⎝      5     ⎠
```

```python
>>> interpolate_Function(test)
(array([ 0, -1]), array([3/5, 2/5], dtype=object))
```

```python
>>> import numpy as np
```

```python
>>> test.args[1].args[1].args
(1/2, a)
```

```python
>>> test.args
(x, a + y)
```

```python
>>> int(3.5)
3
```

```python

```
