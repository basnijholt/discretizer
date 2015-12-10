```python
>>> import sys
>>> import os
>>> path_to_discretizer = '~/discretizer'
>>> sys.path.append(os.path.expanduser(path_to_discretizer))
```

```python
>>> import scipy.sparse.linalg as sla
>>> import sympy
>>> from sympy.interactive import printing
>>> printing.init_printing(use_latex='mathjax')
```

```python
>>> from discretizer import Discretizer
>>> from discretizer import momentum_operators
>>> from types import SimpleNamespace
...
>>> import kwant
>>> import matplotlib.pyplot as plt
>>> %matplotlib inline
```

```python
>>> kx, ky, kz = momentum_operators
>>> V = sympy.Symbol('V')
...
>>> H = kx**2 + ky**2 - V
>>> H
       2      2
-V + kₓ  + k_y
```

```python
>>> tb = Discretizer(H, space_dependent={'V'}, verbose=True)
Discrete coordinates set to:  ['x', 'y']

Function generated for (0, 1):
def _anonymous_func(site1, site2, p):
    x, y = site2.pos
    return (-1)

Function generated for (1, 0):
def _anonymous_func(site1, site2, p):
    x, y = site2.pos
    return (-1)

Function generated for (0, 0):
def _anonymous_func(site, p):
    x, y = site.pos
    V = p.V
    return (4 - V(x, y))
```

```python
>>> tb.symbolic_hamiltonian
⎧                   4           -1           -1 ⎫
⎪(0, 0): -V(x, y) + ──, (0, 1): ───, (1, 0): ───⎪
⎨                    2            2            2⎬
⎪                   a            a            a ⎪
⎩                                               ⎭
```

```python
>>> def stadium(position):
...     (x, y) = position
...     x = max(abs(x) - 20, 0)
...     return x**2 + y**2 < 30**2
...
>>> sys = kwant.Builder()
>>> sys[tb.lat.shape(stadium, (0, 0))] = tb.onsite
...
>>> for hop, val in tb.hoppings.items():
...     sys[hop] = val
...
>>> sys = sys.finalized()
...
>>> kwant.plot(sys)
>>> plt.show()
```

```python
>>> par = SimpleNamespace(V=lambda x,y: 5e-6*(x**2 + y**2))
>>> ham = sys.hamiltonian_submatrix(args=[par], sparse=True)
>>> ev, vec = sla.eigsh(ham, k=10, which='SA')
...
>>> kwant.plotter.map(sys, abs(vec[:,0])**2)
>>> plt.show()
```
