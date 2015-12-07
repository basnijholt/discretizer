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
>>> from discretizer.interpolation import interpolate
```

```python
>>> a, x, y, z = sympy.symbols('a x y z')
>>> A = sympy.Symbol('A')
>>> B = sympy.Symbol('B')
```

```python
>>> inp = A(x + 3*a/4)
>>> inp
 ⎛3⋅a    ⎞
A⎜─── + x⎟
 ⎝ 4     ⎠
```

```python
>>> interpolate(inp)
A(x)   3⋅A(a + x)
──── + ──────────
 4         4
```

```python

```
