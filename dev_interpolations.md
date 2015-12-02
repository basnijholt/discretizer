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
>>> import numpy as np
```

```python
>>> def follow_path(expr, path):
...     res = expr
...     for i in np.arange(len(path)):
...         res = res.args[path[i]]
...     return res
...
>>> def interchange(expr, sub, path):
...     res = sub
...     for i in np.arange(len(path)):
...         temp = follow_path(expr, path[:-(i+1)])
...         args = list(temp.args)
...         args[path[len(path)-i-1]] = res
...         res = temp.func(*tuple(args))
...     return res
...
>>> def interpolate_Function(expr):
...     path = 'None'
...     factor = 'None'
...     change = False
...     summand_0 = 'None'
...     summand_1 = 'None'
...     res = expr
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
...                             #print('found one')
...                             factor = (temp)
...                             path = np.array([i, j, k])
...     if not factor == 'None':
...         change = True
...
...         sign = np.sign(factor)
...         offsets = np.array([int(factor), sign * (int(sign * factor) + 1)])
...         weights = 1/np.abs(offsets - factor)
...         weights = weights/np.sum(weights)
...
...         res = (  weights[0] * interchange(expr, offsets[0] * sympy.Symbol('a'), path[:-1])
...                + weights[1] * interchange(expr, offsets[1] * sympy.Symbol('a'), path[:-1]))
...
...     return sympy.expand(res), change
...
>>> def interpolate(expr):
...     change = False
...     expr = sympy.expand(expr)
...     res = expr
...
...     if isinstance(expr, sympy.Function):# and not change:
...         path = np.array([])
...         temp, change = interpolate_Function(follow_path(expr, path))
...         res = interchange(expr, temp, path)
...
...     for i in np.arange(len(expr.args)):
...         path = np.array([i])
...         if isinstance(follow_path(expr, path), sympy.Function) and not change:
...             temp, change = interpolate_Function(follow_path(expr, path))
...             res = interchange(expr, temp, path)
...
...         for j in np.arange(len(expr.args[i].args)):
...             path = np.array([i, j])
...             if isinstance(follow_path(expr, path), sympy.Function) and not change:
...                 temp, change = interpolate_Function(follow_path(expr, path))
...                 res = interchange(expr, temp, path)
...
...     if change:
...         res = interpolate(res)
...
...     return sympy.expand(res)
```

```python
>>> a, x, y, z = sympy.symbols('a x y z')
>>> A = sympy.Function('A')(x, y, z)
>>> B = sympy.Function('B')(x, y, z)
>>> test = A.subs(y, y-2*a/5)#.subs(x, x+a/2)
>>> test
 ⎛     2⋅a       ⎞
A⎜x, - ─── + y, z⎟
 ⎝      5        ⎠
```

```python
>>> interpolate_Function(test)
(3*A(x, y, z)/5 + 2*A(x, -a + y, z)/5, True)
```

```python
>>> test_2 = A*test
>>> test_2
            ⎛     2⋅a       ⎞
A(x, y, z)⋅A⎜x, - ─── + y, z⎟
            ⎝      5        ⎠
```

```python
>>> print(test)
>>> print(follow_path(test, np.array([1, 1])))
>>> interchange(test, a, np.array([1,1]))
A(x, -2*a/5 + y, z)
-2*a/5
A(x, a + y, z)
```

```python
>>> interpolate(test)
3⋅A(x, y, z)   2⋅A(x, -a + y, z)
──────────── + ─────────────────
     5                 5
```

```python
>>> test_2 = B.subs(z, z+a/2)*(A.subs(x,x+a/2) + A.subs(y, y-a/2))
>>> interpolate(test_2)
A(x, y, z)⋅B(x, y, z)   A(x, y, z)⋅B(x, y, a + z)   A(x, -a + y, z)⋅B(x, y, z)
───────────────────── + ───────────────────────── + ──────────────────────────
          2                         2                           4             

   A(x, -a + y, z)⋅B(x, y, a + z)   A(a + x, y, z)⋅B(x, y, z)   A(a + x, y, z)
 + ────────────────────────────── + ───────────────────────── + ──────────────
                 4                              4                             

⋅B(x, y, a + z)
───────────────
4
```
