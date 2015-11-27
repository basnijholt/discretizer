__all__ = ['algorithms']


from .algorithms import discretize
import sympy

__version__ = '0.0.0'

momentum_operators = sympy.symbols('k_x k_y k_z', commutative=False)
coordinates = sympy.symbols('x y z', commutative=False)


def magic(hamiltonian, space_dependent=None, discrete_coordinates=None,
          verbose=False, symbolic_output=False):

    """This is only for testing of interface."""

    tmp = substitute_functions(hamiltonian, discrete_coordinates, space_dependent)
    hamiltonian, discrete_coordinates = tmp

    if verbose:
        print('Discrete coordinates set to: ', sorted(discrete_coordinates))

    discrete_hamiltonian = discretize(hamiltonian, discrete_coordinates)

    if symbolic_output:
        return discrete_hamiltonian

    tb = make_kwant_functions(discrete_hamiltonian, discrete_coordinates, verbose)
    return tb
