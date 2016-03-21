# General description

Main idea is to provide tool that will discretize Hamiltonian and prepare functions that may be used to create a kwant system.

Input is supposed to be a valid sympy expression representing hermitian Hamiltonian. In can be also a valid sympy matrix. Further development will allow providing a string input that represents Hamiltonian.

Most important, and non-trivial to implement, is support for mix of differential operators and space dependent parameters, like ``H = kx * A(x) * kx``


# Installation (with ``pip``)
```
 pip install git+https://gitlab.kwant-project.org/r-j-skolasinski/discretizer.git
```

# Installation (with ``conda``)
```
conda install -c rskolasinski discretizer
```


# Releases info
## v0.2
    * Providing new input argument that allows to customize coordinates on
      which spatial varying parameters depend.
    * input_hamiltonian attribute now contains preprocessed Hamiltonian

## v0.1.3
    * Input Hamiltonian can be now any sympy matrix (also ImmutableMatrix)
    * Importing whole numpy into kwant function namespaces
    * Discretizer works properly even when momentum operator is defined as
      commutative

## v0.1.2
    * adding __future__ imports to provided Python 2 compatibility

## v0.1.1
    * fix discretization of non-square matrices (issue #18)
    * fix of accessing coordinates in 1D lattices (issue #19)

## v0.1.0
    * freeze of interface provided by discretizer.Discretizer


# To do
* finish interface (waiting for aproval)
* write more tests
* include initial hermicity check


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
and add to ~/.jupyter/jupyter_notebook_config.py
```
c.NotebookApp.contents_manager_class = 'ipymd.IPymdContentsManager'
```
In case it doesn't work ask Rafal. On hpc1 should work by default (dev2 tested)


### ipynb
I decided to use notebooks only to provide examples that can be directly view
with nbviewer. Therefore filters should be deactivated for now to not remove
outputs.


# License
This code was written by Rafał Skolasiński and Sebastian Rubbert.  
The license is the same as that of Kwant (2-clause BSD): http://kwant-project.org/license