import numpy as np
import sympy

from .algorithms import substitute_functions
from .algorithms import discretize

from .postprocessing import offset_to_direction
from .postprocessing import make_kwant_functions
from .postprocessing import offset_to_direction

from .interpolation import interpolate_tb_hamiltonian

try:
    # normal situation
    from kwant import Builder
    from kwant import TranslationalSymmetry
    from kwant import HoppingKind
    from kwant.lattice import Monatomic
except ImportError:
    # probably run on gitlab-ci
    pass

class Discretizer(object):
    """Discretize continous Hamiltonian into its tight binding representation.

    This class provides easy and nice interface for passing models to Kwant.

    Parameters:
    -----------
    hamiltonian : sympy.Expr or sympy.Matrix instance
        Symbolic representation of a continous Hamiltonian. Momentum operators
        should be taken from ``discretizer.momentum_operators``.
    space_dependent : set of strings
        Set of parameters that will be interpreted as  as a function of
        discrete coordinates. For example ``space_dependent={'A', 'B'}``.
    discrete_coordinates : set of strings
        Set of coordinates for which momentum operators will be treated as
        differential operators. For example ``discrete_coordinates={'x', 'y'}``.
        If left as a None they will be obtained from the input hamiltonian by
        reading present coordinates and momentum operators.
    interpolate : bool
        If True all space dependent parameters in onsite and hopping will be
        interpolated to depenend only on the values at site positions.
        Default is False.
    both_hoppings_direction : bool
        If True all hoppings will be returned. For example, if set to True, both
        hoppings into (1, 0) and (-1, 0) will be returned. Default is False.
    verbose : bool
        If True additional information will be printed. Default is False.

    Attributes:
    -----------
    symbolic_hamiltonian : dictionary
        Dictionary containing symbolic result of discretization. Key is the
        direction of the hopping (zeros for onsite)
    lattice : kwant.lattice.Monatomic instance
        Lattice to create kwant system. Lattice constant is set to
        lattice_constant value.
    onsite : function
        The value of the onsite Hamiltonian.
    hoppings : dict
        A dictionary with keys being tuples of the lattice hopping, and values
        the corresponding value functions.
    discrete_coordinates : set of strings
        As in input.
    """
    def __init__(self, hamiltonian, space_dependent=None,
                 discrete_coordinates=None, lattice_constant=1,
                 interpolate=False, both_hoppings_direction=False,
                 verbose=False):
        # preprocessing
        self.input_hamiltonian = hamiltonian
        ham_func, discr_coord = substitute_functions(hamiltonian,
                                                     space_dependent,
                                                     discrete_coordinates)

        self.discrete_coordinates = discr_coord
        if verbose:
            print('Discrete coordinates set to: ', sorted(discr_coord))
            print()

        # making kwant lattice
        dim = len(discr_coord)
        self.lattice = Monatomic(lattice_constant*np.eye(dim).reshape(dim,dim))
        self.lattice_constant = lattice_constant

        # discretization
        tb_hamiltonian = discretize(ham_func, discr_coord)
        tb_hamiltonian = offset_to_direction(tb_hamiltonian, discr_coord)

        if interpolate:
            tb_hamiltonian = interpolate_tb_hamiltonian(tb_hamiltonian)

        if not both_hoppings_direction:
            keys = list(tb_hamiltonian)
            tb_hamiltonian = {k: v for k, v in tb_hamiltonian.items()
                              if k in sorted(keys)[len(keys)//2:]}

        self.symbolic_hamiltonian = tb_hamiltonian.copy()

        for key, val in tb_hamiltonian.items():
            tb_hamiltonian[key] = val.subs(sympy.Symbol('a'), lattice_constant)

        # making kwant functions
        tb = make_kwant_functions(tb_hamiltonian, discr_coord, verbose)
        self.onsite = tb.pop((0,)*len(discr_coord))
        self.hoppings = {HoppingKind(d, self.lattice): val
                         for d, val in tb.items()}

    def build(self, shape, start, symmetry=None):
        """Build Kwant's system.

        Convienient functions that simplifies building of a Kwant's system.

        Parameters:
        -----------
        shape : function
            A function of real space coordinates that returns a truth value:
            true for coordinates inside the shape, and false otherwise.
        start : 1d array-like
            The real-space origin for the flood-fill algorithm.
        symmetry : 1d array-like
            If symmetry is provided a translational invariant system will be
            built. Here symmetry should stand for a lattice vector in which
            system is translational invariant. This vector will be scalled by
            a lattice_constant passed before to ``kwant.TranslationalSymmetry``.

        Returns:
        --------
        system : kwant.Builder instance
        """
        if symmetry is None:
            sys = Builder()
        else:
            vec = self.lattice_constant * np.array(symmetry)
            sys = Builder(TranslationalSymmetry(vec))

        sys[self.lattice.shape(shape, start)] = self.onsite
        for hop, val in self.hoppings.items():
            sys[hop] = val

        return sys
