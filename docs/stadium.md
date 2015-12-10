```python
>>> import os
>>> import sys
...
>>> cwd = os.getcwd().split('/')
>>> sys.path.append(os.path.join('/', *cwd[:cwd.index('discretizer')+1]))
```

```python
>>> from scipy.sparse.linalg import eigsh
>>> import sympy
>>> from sympy.interactive import printing
>>> printing.init_printing(use_latex='mathjax')
```

```python
>>> from discretizer import tb_hamiltonian
>>> from discretizer import momentum_operators
>>> from types import SimpleNamespace
...
>>> import kwant
>>> import matplotlib.pyplot as plt
>>> %matplotlib inline
/home/rskolasinski/environments/dev3/lib/python3.4/site-packages/kwant/solvers/default.py:18: RuntimeWarning: MUMPS is not available, SciPy built-in solver will be used as a fallback. Performance can be very poor in this case.
  "Performance can be very poor in this case.", RuntimeWarning)
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
>>> lat, ons, hops = tb_hamiltonian(H, space_dependent={'V'},
...                                 verbose=True)
Discrete coordinates set to:  ['x', 'y']

Function generated for (0, 1):
def _anonymous_func(site1, site2, p):
    y, x = site2.pos
    return (-1)

Function generated for (1, 0):
def _anonymous_func(site1, site2, p):
    y, x = site2.pos
    return (-1)

Function generated for (0, 0):
def _anonymous_func(site, p):
    y, x = site.pos
    V = p.V
    return (4 - V(x, y))
```

```python
>>> def stadium(position):
...     (x, y) = position
...     x = max(abs(x) - 20, 0)
...     return x**2 + y**2 < 30**2
...
>>> sys = kwant.Builder()
>>> sys[lat.shape(stadium, (0, 0))] = ons
...
>>> for hop, value in hops.items():
...     sys[kwant.builder.HoppingKind(hop, lat)] = value
...
>>> sys = sys.finalized()
...
>>> kwant.plot(sys)
>>> plt.show()
```

```python
>>> par = SimpleNamespace(V=lambda x,y: 5e-6*(x**2 + y**2))
>>> ham = sys.hamiltonian_submatrix(args=[par], sparse=True)
>>> ev, vec = eigsh(ham, k=10, which='SA')
...
>>> kwant.plotter.map(sys, abs(vec[:,0])**2)
>>> plt.show()
```

```python

```

```python

```
