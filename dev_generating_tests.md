## Assumptions
(must be updated)

* space_dependent - these are cast to functions, not commutative
* momentum_symbols - will be specified (default to kx, ky, kz probably) and not commutative
* constant - commutative, everything else
* coordinates are always create as symbol with commutative=False
* lattice constant is always a

```python
>>> from sympy.printing.dot import dotprint
>>> from graphviz import Source
>>> graph = lambda x: Source(dotprint(x))
...
>>> from sympy.interactive import printing
>>> printing.init_printing(use_latex='mathjax')
```

# Global stuff needed for tests

```python
>>> import sympy
...
>>> import discretizer
>>> from discretizer.algorithms import split_factors
```

```python
>>> kx, ky, kz = discretizer.algorithms.momentum_operators
>>> ax, ay, az = discretizer.algorithms.lattice_constants
>>> x, y, z = discretizer.algorithms.coord
...
>>> Psi = discretizer.algorithms.wf
>>> A, B = sympy.symbols('A B', commutative=False)
...
>>> ns = {'A': A, 'B': B, 'a_x': ax, 'a_y': ay, 'az': az, 'x': x, 'y': y, 'z': z}
```

# Tests for: derivate

```python
>>> from discretizer.algorithms import derivate
```

## Manual test

```python
>>> derivate(A(x), kx)
   ⎛  A(-aₓ + x)   A(aₓ + x)⎞
-ⅈ⋅⎜- ────────── + ─────────⎟
   ⎝     2⋅aₓ         2⋅aₓ  ⎠
```

## Generating automatic tests for nose

```python
>>> to_test = [
...     ('A(x)', 'kx'),
...     ('A(x)', 'ky'),
...     ('A(x)*B', 'kx'),
...     ('A(x) + B(x)', 'kx'),
...     ('A', 'kx'),
...     ('5', 'kx'),
...     ('A(x) * B(x)', 'kx'),
...     ('A(x) * B', 'kx',
...     ('Psi', 'kx')),
...     ('Psi', 'ky')
>>> ]
```

```python
>>> outputs = []
>>> for test in to_test:
...     exec("tmp = str(derivate({0[0]}, {0[1]}))".format(test))
...     out = "({}, {}): '{}',".format(test[0], test[1], tmp)
...     outputs.append(out)
...
>>> outputs = ['    ' + out for out in outputs]
>>> outputs = 'test = {\n' + "\n".join(outputs) + '\n}'
```

### This should be copied into input of corresponding test

```python
>>> print(outputs)
test = {
    (A(x), kx): '-I*(-A(-a_x + x)/(2*a_x) + A(a_x + x)/(2*a_x))',
    (A(x), ky): '0',
    (A(x)*B, kx): '-I*(-A(-a_x + x)*B/(2*a_x) + A(a_x + x)*B/(2*a_x))',
    (A(x) + B(x), kx): '-I*(-A(-a_x + x)/(2*a_x) + A(a_x + x)/(2*a_x) - B(-a_x + x)/(2*a_x) + B(a_x + x)/(2*a_x))',
    (A, kx): '0',
    (5, kx): '0',
    (A(x) * B(x), kx): '-I*(-A(-a_x + x)*B(-a_x + x)/(2*a_x) + A(a_x + x)*B(a_x + x)/(2*a_x))',
    (A(x) * B, kx): '-I*(-A(-a_x + x)*B/(2*a_x) + A(a_x + x)*B/(2*a_x))',
    (Psi, ky): '-I*(-Psi(x, -a_y + y, z)/(2*a_y) + Psi(x, a_y + y, z)/(2*a_y))',
}
```

### Checking if it works

```python
>>> exec(outputs)
```

```python
>>> for inp, out in test.items():
...     got = (derivate(*inp))
...     out = sympy.sympify(out, locals=ns)
...     assert  sympy.simplify(sympy.expand(got - out)) == 0,\
...         "Should be: derivate({0[0]}, {0[1]})=={1}. Not {2}".format(inp, out, got)
```

# Tests for: discretize_expression and discretize_summand

```python
>>> from discretizer import discretize_expression
>>> from discretizer.algorithms import discretize_summand
>>> from discretizer.algorithms import wf as Psi
```

## Manual tests

```python
>>> expr = kx*A(x); expr
kₓ⋅A(x)
```

```python
>>> discretize_summand(expr)
   ⎛  A(-aₓ + x)   A(aₓ + x)⎞
-ⅈ⋅⎜- ────────── + ─────────⎟
   ⎝     2⋅aₓ         2⋅aₓ  ⎠
```

## Generating automatic tests for nose

```python
>>> to_test = [
...     "kx * A(x)",
...     "kx * Psi",
...     "kx**2 * Psi",
...     "kx * A(x) * kx * Psi",
>>> ]
```

```python
>>> outputs = []
>>> for test in to_test:
...     exec("tmp = str(discretize_summand({}))".format(test))
...     out = "{}: '{}',".format(test, tmp)
...     outputs.append(out)
...
>>> outputs = ['    ' + out for out in outputs]
>>> outputs = 'test = {\n' + "\n".join(outputs) + '\n}'
```

### This should be copied into input of corresponding test

```python
>>> print(outputs)
test = {
    kx * A(x): '-I*(-A(-a_x + x)/(2*a_x) + A(a_x + x)/(2*a_x))',
    kx * Psi: '-I*(-Psi(-a_x + x, y, z)/(2*a_x) + Psi(a_x + x, y, z)/(2*a_x))',
    kx**2 * Psi: 'Psi(x, y, z)/(2*a_x**2) - Psi(-2*a_x + x, y, z)/(4*a_x**2) - Psi(2*a_x + x, y, z)/(4*a_x**2)',
    kx * A(x) * kx * Psi: 'A(-a_x + x)*Psi(x, y, z)/(4*a_x**2) - A(-a_x + x)*Psi(-2*a_x + x, y, z)/(4*a_x**2) + A(a_x + x)*Psi(x, y, z)/(4*a_x**2) - A(a_x + x)*Psi(2*a_x + x, y, z)/(4*a_x**2)',
}
```

### Checking if it works

```python
>>> exec(outputs)
```

```python
>>> for inp, out in test.items():
...     got = (discretize_summand(inp))
...     out = sympy.sympify(out, locals=ns)
...     assert  sympy.simplify(sympy.expand(got - out)) == 0,\
...         "Should be: discretize_summand({})=={}. Not {}".format(inp, out, got)
```

# spliting into lhs, operators, rhs
(this is already moved into tests, left just to show how it works)

```python
>>> from discretizer.algorithms import split_factors
```

```python
>>> test_operators = [kz*Psi, A*kx**2*Psi, A*kx**2*ky*Psi, ky*A*kx*B*Psi, kx, kx**2, A, Psi, 3]
...
>>> tested = []
>>> output = []
>>> for operator in test_operators:
...     try:
...         output.append(split_factors(operator))
...         tested.append(operator)
...     except ValueError:
...         print('ValueError on', operator)
...     except AssertionError:
...         print('AssertionError on', operator, operator.func)
```

```python
>>> tested
⎡                    2                 2                                      
⎣k_z⋅Ψ(x, y, z), A⋅kₓ ⋅Ψ(x, y, z), A⋅kₓ ⋅k_y⋅Ψ(x, y, z), k_y⋅A⋅kₓ⋅B⋅Ψ(x, y, z)

        2                  ⎤
, kₓ, kₓ , A, Ψ(x, y, z), 3⎦
```

```python
>>> output
⎡                                              ⎛    2                 ⎞       
⎣(1, k_z, Ψ(x, y, z)), (A⋅kₓ, kₓ, Ψ(x, y, z)), ⎝A⋅kₓ , k_y, Ψ(x, y, z)⎠, (k_y⋅

                                                                              
A, kₓ, B⋅Ψ(x, y, z)), (1, kₓ, 1), (kₓ, kₓ, 1), (1, 1, A), (1, 1, Ψ(x, y, z)), 

         ⎤
(1, 1, 3)⎦
```
