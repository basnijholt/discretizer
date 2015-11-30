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

# Possible inputs to hopping functions as outputs from discretization

```python
>>> import sympy
>>> sympy.init_printing(use_latex='mathjax')
>>> from discretizer import discretize
```

```python
>>> from discretizer.algorithms import wavefunction_name
>>> from discretizer.algorithms import _discretize_summand
```

```python
>>> from discretizer.algorithms import read_hopping_from_wf
>>> from discretizer.algorithms import extract_hoppings
>>> from discretizer.algorithms import shortening
```

```python
>>> def get_wf(discrete_coordinates):
...     """Rturn Psi depending only on discrete_coordinates. Help function."""
...     coordinates_names = sorted(list(discrete_coordinates))
...     coordinates = [sympy.Symbol(c, commutative=False) for c in coordinates_names]
...     Psi = sympy.Function(wavefunction_name)(*coordinates)
...
...     return Psi
```

```python
>>> def test_process_on_summand(summand, discrete_coordinates):
...     out1 = _discretize_summand(summand, discrete_coordinates)
...     out2 = extract_hoppings(out1)
...     out3 = shortening(out2, discrete_coordinates)
...
...     return out1, out2, out3
```

# Support stuff

```python
>>> kx, ky, kz = sympy.symbols('k_x k_y k_z', commutative=False)
>>> x, y, z = sympy.symbols('x y z', commutative=False)
...
>>> A, B = sympy.symbols('A B', commutative=False)
```

```python
>>> discrete_coordinates = {'x', 'y'}
...
>>> Psi = get_wf(discrete_coordinates)
>>> summand = kx * Psi
```

```python
>>> summand
kₓ⋅Ψ(x, y)
```

```python
>>> out1, out2, out3 = test_process_on_summand(summand, discrete_coordinates)
```

```python
>>> out1
   ⎛  Ψ(-aₓ + x, y)   Ψ(aₓ + x, y)⎞
-ⅈ⋅⎜- ───────────── + ────────────⎟
   ⎝       2⋅aₓ           2⋅aₓ    ⎠
```

```python
>>> out2
defaultdict(int,
            ⎧          ⅈ            -ⅈ  ⎫
⎨(-1, 0): ────, (1, 0): ────⎬
⎩         2⋅aₓ          2⋅aₓ⎭)
```

```python
>>> out3
⎧          ⅈ           -ⅈ ⎫
⎨(-1, 0): ───, (1, 0): ───⎬
⎩         2⋅a          2⋅a⎭
```
