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
>>> import sys, os
```

```python
>>> sys.path.append('/home/raphael/discretizer/discretizer/')
```

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
```

```python
>>> discretize_expression(kx*kx)
```

```python
>>> discretize_expression(kx*ky)
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
```

```python
>>> from discretizer import shortening
```

```python
>>> shortening(hop)
```

```python

```
