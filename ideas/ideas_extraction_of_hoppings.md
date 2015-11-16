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
>>> from discretizer import coord, momentum_operators, a
>>> from discretizer import substitute_functions
>>> from discretizer import derivate
>>> from discretizer import split_factors
>>> from discretizer import discretize_expression
```

```python
>>> kx, ky, kz = momentum_operators
>>> Psi = sympy.Function('Psi')(*coord)
>>> A = sympy.Function('A')(*coord)
```

# Possible inputs to hopping functions as outputs from discretization

```python
>>> discretize_expression(kx*A*kx*Psi)
A(-aₓ + x, y, z)⋅Ψ(x, y, z)   A(-aₓ + x, y, z)⋅Ψ(-2⋅aₓ + x, y, z)   A(aₓ + x, 
─────────────────────────── - ─────────────────────────────────── + ──────────
               2                                 2                            
           4⋅aₓ                              4⋅aₓ                             

y, z)⋅Ψ(x, y, z)   A(aₓ + x, y, z)⋅Ψ(2⋅aₓ + x, y, z)
──────────────── - ─────────────────────────────────
    2                                2              
4⋅aₓ                             4⋅aₓ
```

```python
>>> discretize_expression(kx*kx*Psi)
Ψ(x, y, z)   Ψ(-2⋅aₓ + x, y, z)   Ψ(2⋅aₓ + x, y, z)
────────── - ────────────────── - ─────────────────
      2                2                    2      
  2⋅aₓ             4⋅aₓ                 4⋅aₓ
```

```python
>>> discretize_expression(kx*ky*Psi)
  Ψ(-aₓ + x, -a_y + y, z)   Ψ(-aₓ + x, a_y + y, z)   Ψ(aₓ + x, -a_y + y, z)   
- ─────────────────────── + ────────────────────── + ────────────────────── - 
          4⋅aₓ⋅a_y                 4⋅aₓ⋅a_y                 4⋅aₓ⋅a_y          

Ψ(aₓ + x, a_y + y, z)
─────────────────────
       4⋅aₓ⋅a_y
```

```python
>>> test = discretize_expression( kx**2 * A * ky * Psi ) # yeah, I know it's not hermitian but it is only one summand
```

```python
>>> def Psi_to_hopping(Psi):
...     offset = []
...     for argument in Psi.args:
...         temp = sympy.expand(argument)
...         if temp in sympy.symbols('x y z', commutative = False):
...             offset.append(0)
...         elif temp.func == sympy.Add:
...             for arg_summands in temp.args:
...                 if arg_summands.func == sympy.Mul:
...                     if len(arg_summands.args) > 2:
...                         print('More than two factors in an argument of Psi')
...                     if not arg_summands.args[0] == sympy.symbols('a'):
...                         offset.append(arg_summands.args[0])
...                     else:
...                         offset.append(arg_summands.args[1])
...                 elif arg_summands == a:
...                     offset.append(1)
...         else:
...             print('Argument of \Psi is neither a sum nor a single space variable.')
...     return tuple(offset)
...
>>> # This function will not have the shortening included
... def extract_hoppings(expr):
...     # this line should be unneccessary in the long run. Now I want to avoid errors due to wrong formats of the input.
...     expr = sympy.expand(expr)
...
...     # The output will be stored in a dictionary, with terms for each hopping kind
...     # This format is probably not good in the long run
...     hoppings = {}
...
...     if expr.func == sympy.Add:
...         for summand in expr.args:
...             #find a way to make it readable
...             if not summand.func == sympy.Function('Psi'):
...                 for i in range(len(summand.args)):
...                     if summand.args[i].func == sympy.Function('Psi'):
...                         index = i
...                 if index < len(summand.args) - 1:
...                     print('Psi is not in the very end of the term. Output will be wrong!')
...
...                 try:
...                     hoppings[Psi_to_hopping(summand.args[-1])] += sympy.Mul(summand.args[:-1])
...                 except:
...                     hoppings[Psi_to_hopping(summand.args[-1])] = sympy.Mul(summand.args[:-1])
...             else:
...                 try:
...                     hoppings[Psi_to_hopping(summand)] += 1
...                 except:
...                     hoppings[Psi_to_hopping(summand)] = 1
...
...     else:
...         if not expr.func == sympy.Function('Psi'):
...             for i in range(len(expr.args)):
...                 if expr.args[i].func == sympy.Function('Psi'):
...                     index = i
...             if index < len(expr.args) - 1:
...                 print('Psi is not in the very end of the term. Output will be wrong!')
...
...             try:
...                 hoppings[Psi_to_hopping(expr.args[-1])] += sympy.Mul(expr.args[:-1])
...             except:
...                 hoppings[Psi_to_hopping(expr.args[-1])] = sympy.Mul(expr.args[:-1])
...         else:
...             try:
...                 hoppings[Psi_to_hopping(expr)] += 1
...             except:
...                 hoppings[Psi_to_hopping(expr)] = 1
...     return hoppings
```

```python
>>> x, y, z = coord
>>> single_term = Psi.subs(x, x+a).subs(y, y-2*a)
>>> Psi_to_hopping(single_term)
(1, -2, 0)
```

```python
>>> expr = Psi
>>> expr = A*expr.subs(x, x+2*a).subs(y, y-2*a).subs(z, z + a) + expr.subs(x, x + a)
>>> expr
A(x, y, z)⋅Ψ(2⋅a + x, -2⋅a + y, a + z) + Ψ(a + x, y, z)
```

```python
>>> hop = extract_hoppings(expr)
>>> hop
{(1, 0, 0): 1, (2, -2, 1): (A(x, y, z),)}
```

```python

```
