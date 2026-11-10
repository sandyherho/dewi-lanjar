# Command-line interface

```text
dewi-lanjar [--version] <command> [options]
```

| command | purpose |
|---|---|
| `list` | list the built-in cases |
| `show CASE` | print a case as fully resolved TOML with its description |
| `init FILE [--template CASE]` | write an editable case file |
| `run CASE_OR_FILE` | run a case and write its outputs |
| `run-all [--kinds steady,transient,sweep]` | run every built-in case |
| `steady` | solve a steady linear gyre from command-line parameters |
| `scales` | scales and boundary-layer widths from dimensional ocean parameters |
| `widths` | theoretical boundary-layer roots and widths |
| `validate [--quick]` | grid convergence against the exact solution |
| `explain [TOPIC]` | print a derivation |

## `run`

```text
dewi-lanjar run CASE_OR_FILE [--out DIR] [--set SECTION.KEY=VALUE ...]
                             [--no-gif] [--no-png] [--quiet]
```

`CASE_OR_FILE` is a built-in name (`dewi-lanjar list`) or a path to a TOML
file. `--set` may be repeated and overrides any key of the case, for example
`--set physics.delta_i=0.05 --set grid.nx=256 --set floats.n=0`. Outputs go to
`DIR/<case name>/` (default `outputs/`).

## `steady`

```text
dewi-lanjar steady [--delta-s X] [--delta-m X] [--beta-factor B]
                   [--bc free-slip|no-slip] [--bc-ns free-slip|no-slip]
                   [--wind single-gyre|double-gyre|asymmetric-double-gyre|triple-gyre|cyclonic-gyre]
                   [--sign 1|-1] [--asymmetry A] [--nx N] [--ny N] [--Ly LY]
                   [--png FILE] [--why]
```

Prints the full diagnostic report and, when an exact solution exists, the
error against it. `--why` appends the derivation of the relevant balance.

## `scales`

```text
dewi-lanjar scales [--L KM] [--lat DEG | --beta B] [--r R] [--A-h A]
                   [--tau0 T] [--rho0 RHO] [--H H]
```

Shows the Sverdrup velocity, time and transport scales and the three
boundary-layer widths with every number substituted, for example

```text
dewi-lanjar scales --L 6000 --lat 35 --r 1e-6 --A-h 2000
```

## `widths`

```text
dewi-lanjar widths [--delta-s X] [--delta-m X] [--beta-factor B] [--bc ...]
```

Lists the roots of $\delta_M^3 m^3 - \delta_S m - b = 0$, marks the ones that
form the western layer, and reports the decay scale, e-folding width and
counter-current position.

## `explain`

Topics: `overview`, `nondimensionalization`, `sverdrup`, `stommel`, `munk`,
`inertial`, `fplane`, `budgets`, `numerics`, `diagnostics`, `floats`.

## Exit status

0 on success, 1 when `validate` finds an order below 1.8 or `init` would
overwrite a file, 2 on invalid input.
