# API reference

## Setup

```{eval-rst}
.. autoclass:: dewi_lanjar.Grid
   :members:

.. autoclass:: dewi_lanjar.Params
   :members:

.. autoclass:: dewi_lanjar.DimensionalScales
   :members:

.. autofunction:: dewi_lanjar.beta_at_latitude

.. autofunction:: dewi_lanjar.boundary_layer_roots

.. autoclass:: dewi_lanjar.WindForcing
   :members:
```

## Solvers

```{eval-rst}
.. autofunction:: dewi_lanjar.solve_steady

.. autoclass:: dewi_lanjar.SteadySolution

.. autoclass:: dewi_lanjar.TransientSolver
   :members:

.. autoclass:: dewi_lanjar.RunConfig

.. autoclass:: dewi_lanjar.TransientResult

.. autoclass:: dewi_lanjar.Floats
   :members:
```

## Exact and asymptotic solutions

```{eval-rst}
.. automodule:: dewi_lanjar.analytic
   :members: exact_linear, solve_mode, characteristic_roots, sverdrup_interior,
             sverdrup_transport, boundary_layer_profile, theory_widths, measure_widths
```

## Diagnostics

```{eval-rst}
.. automodule:: dewi_lanjar.diagnostics
   :members: compute, describe, kinetic_energy, enstrophy, normalized_entropy,
             gyre_latitude_index, gyre_band, tendency_terms
```

## Operators

```{eval-rst}
.. automodule:: dewi_lanjar.operators
   :members:
```

## Cases

```{eval-rst}
.. automodule:: dewi_lanjar.cases
   :members: load_case, case_from_dict, builtin_names, run_case, Case, CaseResult, apply_overrides
```
