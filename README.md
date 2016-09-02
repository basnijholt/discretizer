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
See [changelog](CHANGELOG).

# To do
* simplify interface
* write more tests
* include initial hermicity check


# Development guide
One should use Python 3 for development.

Notebooks with examples should be commited together with outputs so they can
be viewed easily through ``nbviewer``.

# License
This code was written by Rafał Skolasiński and Sebastian Rubbert.  
The license is the same as that of Kwant (2-clause BSD): http://kwant-project.org/license
