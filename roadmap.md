
# Roadmap of module

Maybe different name of this file?

## Requirments

discretize continues hamiltonian into tight binding one

### must do
* handle space dependent parameters mixed with differential operators
* easy and nice input

# Program input

Hamiltonian provided by user is **supposed** to be **hermitian**. Check will be done in the beginning. 
Therefore providing $H = A(x) \partial_x$ or $H = A(x) \partial_x^2$ should throw out some kind of ValueError.

Correct input should be symmetrized by user, for example $H = \frac{1}{2} \big(\partial_x A(x) + A(x) \partial_x \big)$ or $H = k_x A(x) k_x$

# Core algorithm

## Input
Easiest would bo to work internally on sympy matrices. Every element of matrix would be representent by expression build from mix of commutative and not commutative symbols, for example $\hat{H} = k_x A(x) k_x$, where $A(x)$ could stand for $-\frac{\hbar^2}{2m(x)}$.

More complete example could be:
$$
H = 
\begin{pmatrix}
k_x A(x) k_x & k_xP(x) \\
P(x)k_x & k_x B(x) k_x
\end{pmatrix}
$$

We should also allow to specify which parameters are constants

## Remarks (future improvements)
* leave space to provide ``interpolate = False`` option for midpoints
* leave space for keeping matrices as symbol up to the end

## Things to decide
* if function should work on real space coordinates or in lattice coordinates: ``f(x, y)`` or ``f(i, j)``
* default behaviour at the edges of systems (ask Michael)

# Things to do:
* write down and code test examples
* find good example of recursive algorithm
* check how sympy stores expression, how to access last one from the right
* check how Anton generate valid kwant functions

# Magnetic field through vector potential
This will go into separate module. Check issue 19.



```python

```
