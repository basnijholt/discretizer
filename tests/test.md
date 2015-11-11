---
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
---

```python
>>> import numpy as np
>>> from collections import namedtuple
```

```python
>>> nt_test = namedtuple('nt_test', 'A P a')
```

```python
>>> par = nt_test(A = lambda x: 1.0 + (x > 1.5),
...               P = lambda x: 1.0 + (x > 1.5),
...               a = 1.0)
```

## Test I

$ H = kAk $

```python
>>> def hop_nn(xs, xt, par):
...     A = par.A
...     res = -(A(xt) + A(xs)) / (2*par.a**2)
...     return res
...
>>> def onsite(x, par):
...     A = par.A
...     res = (A(x-1)/2 + A(x) + A(x+1)/2) / (par.a**2)
```

## Test II

$H =  \left( \begin{array}{cc}
kAk & Pk \\
kP & kAk \end{array} \right) $

```python
>>> sigma_0 = np.array([[1, 0], [0, 1]])
>>> sigma_x = np.array([[0, 1], [1, 0]])
>>> sigma_y = np.array([[0, -1j], [1j, 0]])
>>> sigma_z = np.array([[1, 0], [0, -1]])
```

```python
>>> #DOUBLE CHECK SIGNS!!! THEY ARE BITCHES!!
...
... def hop_nn(xs, xt, par):
...     A = par.A
...     res = -((A(xt) + A(xs)) / (2*par.a**2) * sigma_0
...             + np.sign(xs-xt) * 1j * sigma_x * (A(xs) + A(xt)) / (4 * par.a)
...             - np.sign(xs-xt) * sigma_y * (A(xs) - A(xt)) / (4 * par.a))
...     return res
...
>>> def onsite(x, par):
...     A = par.A
...     res = (A(x-1)/2 + A(x) + A(x+1)/2) / (par.a**2) * sigma_0
...     return res
```

## Test III a

$H =  \left( \begin{array}{cc}
kAk & (k_x k_y P + P k_x k_y)/2 \\
(k_x k_y P + P k_x k_y)/2 & kAk \end{array} \right) $

```python
>>> par = nt_test(A = lambda x, y: 1.0 + (x > 1.5) + (y > 1.5),
...               P = lambda x, y: 1.0 + (x > 1.5) + (y > 1.5),
...               a = 1.0)
```

```python
>>> #SIGNS?
... def hop_x_or_y(rs, rt, par):
...     xs = rs[0]
...     xt = rt[0]
...     ys = rs[1]
...     yt = rt[1]
...     A = par.A
...     res = -(A(xs, ys) + A(xt, yt)) / (2*par.a**2) * sigma_0
...     return res
...
>>> def hop_x_and_y(rs, rt, par):
...     xs = rs[0]
...     xt = rt[0]
...     ys = rs[1]
...     yt = rt[1]
...     P = par.P
...     res =  np.sign(xs-xt) * np.sign(ys-yt) * (P(xs, ys) + P(xt, yt)) * sigma_x / (4 * par.a**2)
...     return res
```

## Test III b

$H =  \left( \begin{array}{cc}
kAk & (k_x P k_y + k_y P k_x)/2 \\
(k_x P k_y  + k_y P k_x)/2 & kAk \end{array} \right) $

```python
>>> #SIGNS?
... def hop_x_or_y(rs, rt, par):
...     xs = rs[0]
...     xt = rt[0]
...     ys = rs[1]
...     yt = rt[1]
...     A = par.A
...     res = -(A(xs, ys) + A(xt, yt)) / (2*par.a**2) * sigma_0
...     return res
...
>>> def hop_x_and_y(rs, rt, par):
...     xs = rs[0]
...     xt = rt[0]
...     ys = rs[1]
...     yt = rt[1]
...     P = par.P
```

```python

```
