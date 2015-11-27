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
...         args[path[i]] = res
...         res = temp.func(*tuple(args))
...     return res
...
>>> def interpolate_Function(expr):
...     path = 'None'
...     factor = 'None'
...     interupt = False
...     summand_0 = 'None'
...     summand_1 = 'None'
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
...         interupt = True
...
...         sign = np.sign(factor)
...         offsets = np.array([int(factor), sign * (int(sign * factor) + 1)])
...         weights = 1/np.abs(offsets - factor)
...         weights = weights/np.sum(weights)
...
...         res = (  weights[0] * interchange(expr, offsets[0] * sympy.Symbol('a'), path[:-1])
...                + weights[1] * interchange(expr, offsets[1] * sympy.Symbol('a'), path[:-1]))
...
...     return sympy.expand(res), interupt
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
(3*A(x, y)/5 + 2*A(x, -a + y)/5, True)
```

```python
>>> def interpolate(expr):
...
...     for depth in np.arange(3):
...
...             if isinstance(expr.args[i], sympy.Function):
...                 path = np.array([i])
...                 temp, change = interpolate_Function(follow_path(expr, path))
```

```python
>>> test_2 = A*test
>>> test_2
         ⎛     2⋅a    ⎞
A(x, y)⋅A⎜x, - ─── + y⎟
         ⎝      5     ⎠
```

```python
>>> interpolate(test_2)
False
True
```

```python
>>> def search()
```

```python
>>> np.linspace(1,17,17)[:-17]
array([], dtype=float64)
```

```python
>>> print(test)
>>> print(follow_path(test, np.array([1, 1])))
>>> interchange(test, a, np.array([1,1]))
A(x, -2*a/5 + y)
-2*a/5
A(x, a + y)
```

```python
>>> follow_path(test, np.array([1,1]))
-2⋅a 
─────
  5
```

```python

```
