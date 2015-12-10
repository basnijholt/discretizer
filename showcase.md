# Showcase

This notebooks shows general features of interface. For examples please see:

1. [Quantum Stadium](examples/stadium.md)
2. [Edge states in HgTe](examples/qsh.md)

## Imports general

```python
>>> import sympy
>>> from sympy.interactive import printing
>>> printing.init_printing(use_latex='mathjax')
```

## Imports discretizer

```python
>>> from discretizer import Discretizer
>>> from discretizer import momentum_operators
```

## Defining sample expression

```python
>>> kx, ky, kz = momentum_operators
...
>>> A, B, C = sympy.symbols('A B C', commutative=False)
>>> H = sympy.Matrix([[kx*A*kx +ky*A*ky, kx*B], [B*kx, C]]); H
```

```python
>>> space_dependent = {'A', 'B'}
>>> discrete_coordinates = {'x', 'y'}
```

# class interface

```python
>>> tb = Discretizer(H, space_dependent, discrete_coordinates, lattice_constant=2.0, verbose=True)
```

```python
>>> tb.symbolic_hamiltonian
```

```python
>>> tb.lattice
```

```python
>>> tb.onsite, tb.hoppings
```
