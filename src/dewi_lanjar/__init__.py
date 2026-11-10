"""dewi-lanjar: western boundary intensification of wind-driven ocean gyres.

A library and command-line tool for the barotropic vorticity equation on a
beta plane in a closed basin, from the linear Stommel (1948) and Munk (1950)
problems to inertial, eddying double gyres. Every solution comes with
diagnostics (energy and vorticity budgets, boundary-layer widths compared
with theory, transport closure, energy-concentration entropies) and every
equation has a derivation in the documentation and in
``dewi-lanjar explain``.

Quick start::

    from dewi_lanjar import Grid, Params, WindForcing, solve_steady, diagnostics

    grid = Grid(128, 128)
    p = Params(delta_s=0.04)                      # Stommel layer, 4 % of the basin width
    sol = solve_steady(grid, p, WindForcing("single-gyre"))
    d = diagnostics.compute(sol.psi, sol.zeta, grid, p, sol.forcing)
    print(diagnostics.describe(d))

Named after Dewi Lanjar, queen of the spirits of the northern Java Sea in
Javanese tradition.
"""

from . import analytic, diagnostics
from ._version import __version__
from .analytic import (
    boundary_layer_profile,
    exact_linear,
    sverdrup_interior,
    sverdrup_transport,
    theory_widths,
)
from .cases import Case, builtin_names, case_from_dict, load_case, run_case
from .floats import Floats
from .forcing import PROFILES, WindForcing
from .grid import Grid
from .params import DimensionalScales, Params, beta_at_latitude, boundary_layer_roots
from .steady import SteadySolution, solve_steady
from .transient import RunConfig, TransientResult, TransientSolver

__all__ = [
    "__version__",
    "Grid", "Params", "DimensionalScales", "WindForcing", "PROFILES",
    "beta_at_latitude", "boundary_layer_roots",
    "solve_steady", "SteadySolution",
    "TransientSolver", "RunConfig", "TransientResult", "Floats",
    "exact_linear", "sverdrup_interior", "sverdrup_transport",
    "boundary_layer_profile", "theory_widths",
    "diagnostics", "analytic",
    "Case", "load_case", "case_from_dict", "builtin_names", "run_case",
]
