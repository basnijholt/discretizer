# To do list
* shortening and reading hoppings (to be tested)
* generalazing it to work on matrices (important)
* transforming it into functions (check how Anton did it)
* writing tests, lot of tests (working on it)
* making package (almost done)
* wrapper of string input (nice, but least on priorty)

# Development guide

One should use Python 3 for development.

### Ipymd
Please use ``ipymd`` to work on markdown ``md`` instead of ``ipynb`` files.
If you use your own laptop please do:
```
pip install git+https://github.com/rossant/ipymd
```
and add to ~/.jupyter/jupyter
```
c.NotebookApp.contents_manager_class = 'ipymd.IPymdContentsManager'
```
In case it doesn't work ask Rafal. On hpc1 should work by default (dev2 tested)

### If we decide to use ``ipynb`` notebooks we may need output filters.
Let's activate them now, just in case, with:
```
sh activate_filters.sh
```

# General description

Main idea is to provide tool that will discretize Hamiltonian and prepare functions that may be used to create a kwant system.

Input is supposed to be a valid sympy expression representing hermitian Hamiltonian. In can be also a valid sympy matrix. Further development will allow providing a string input that represents Hamiltonian.

Most important, and non-trivial to implement, is support for mix of differential operators and space dependent parameters, like ``H = kx * A(x) * kx``

## Things to decide
* default behaviour at the edges of systems (ask Michael)

## Remarks (future improvements)
* leave space to provide ``interpolate = False`` option for midpoints
* leave space for keeping matrices as symbol up to the end
