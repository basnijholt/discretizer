# General description

Main idea is to provide tool that will discretize Hamiltonian and prepare functions that may be used to create a kwant system.

Input is supposed to be a valid sympy expression representing hermitian Hamiltonian. In can be also a valid sympy matrix. Further development will allow providing a string input that represents Hamiltonian.

Most important, and non-trivial to implement, is support for mix of differential operators and space dependent parameters, like ``H = kx * A(x) * kx``


#To do
* finish interface (high priority)
* write more tests
* include initial hermicity check
* include midpoints interpolation option


# To decide
* default behaviour at the edges of systems (ask Michael)
* precise handling of space coordinates


# Development guide

One should use Python 3 for development.

### ipymd
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


### ipynb
I decided to use notebooks only to provide exampels that can be directly view
with nbviewer. Therefore filters should be deactivated for now to not remove
outputs.
