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
...     offsets = np.array([int(factor), sign * (int(sign * factor) + 1)])
...     weights = 1/np.abs(offsets - factor)
...     weights = weights/np.sum(weights)
...     #vectoriziation seems not to work for sympy expressions^^
...     #first only make the new offsets
...     summand_0 = sympy.Mul(offsets[0]*expr.args[path[0]].args[path[1]].args[path[2]])
...     summand_1 = sympy.Mul(offsets[1]*expr.args[path[0]].args[path[1]].args[path[2]])
...     #rebuild the whole (single) argument
...     summand_0 = sympy.Add(sympy.Add(summand_0, *expr.args[path[0]].args[:path[1]]),
...                           *expr.args[path[0]].args[path[1]+1:])
...     summand_1 = sympy.Add(sympy.Add(summand_1, *expr.args[path[0]].args[:path[1]]),
...                           *expr.args[path[0]].args[path[1]+1:])
...     #rebuild the function with all arguments (also the unchanghed ones)
...     #In opposition to addition or multiplication, we have to pass
...     #all arguments to the function at once now, so we start by
...     #creating a list of the arguments that we need.
...     temp_0 = []
...     temp_1 = []
...     for i in np.arange(len(expr.args)):
...         if not i == path[0]:
...             temp_0.append(expr.args[i])
...             temp_1.append(expr.args[i])
...         else:
...             temp_0.append(summand_0)
...             temp_1.append(summand_1)
...     summand_0 = weights[0] * expr.func(*tuple(temp_0))
...     summand_1 = weights[1] * expr.func(*tuple(temp_1))
...     return sympy.expand(summand_0 + summand_1)
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
3⋅A(x, y)   2⋅A(x, -a + y)
───────── + ──────────────
    5             5
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
>>> bla=[0, 1, 2, 3]
```

```python
>>> tuple(bla)
(0, 1, 2, 3)
```

```python

```
