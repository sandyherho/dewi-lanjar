# Case files

A case is a TOML file. Every key has a default, so a file only needs the keys
it changes. `dewi-lanjar show <case>` prints any case with all defaults filled
in, and `dewi-lanjar init my.toml --template <case>` writes one to edit.

```toml
[case]
name = "my_gyre"                 # output folder name
title = "My gyre"
kind = "transient"               # steady | transient | sweep
description = """Free text, copied into summary.txt."""

[grid]
nx = 192                         # cells in x (nodes: nx + 1)
ny = 192
Lx = 1.0                         # basin size in units of L
Ly = 1.0

[physics]
delta_s = 0.02                   # Stommel width r/(beta L)
delta_m = 0.02                   # Munk width (A_h/beta)^(1/3)/L
delta_i = 0.035                  # inertial width (U/beta)^(1/2)/L  (ignored by steady and sweep)
beta_factor = 1.0                # 1 = beta plane, 0 = f plane
bc = "no-slip"                   # western and eastern walls: free-slip | no-slip
bc_ns = "free-slip"              # southern and northern walls (default: same as bc)
from_scales = false              # true: compute the three deltas from [scales]

[wind]
profile = "single-gyre"          # single-gyre | double-gyre | asymmetric-double-gyre | triple-gyre | cyclonic-gyre
amplitude = 1.0
sign = 1                         # -1 mirrors the circulation (Southern Hemisphere)
asymmetry = 0.5                  # for asymmetric-double-gyre
seasonal_amplitude = 0.0         # periodic modulation 1 + A sin(2 pi t / P)
seasonal_period = 50.0

[run]                            # transient (n_frames also sets the steady animation length)
t_end = 200.0                    # in units of 1/(beta L)
n_frames = 160
t_record = 0.0                   # time of the first stored frame
cfl = 0.7                        # safety factor on all stability limits
dt_max = 0.5                     # optional cap on the time step
diag_per_frame = 2               # diagnostic samples per frame interval
init_noise = 0.0                 # smooth random initial vorticity (breaks symmetry)
init_seed = 0

[floats]
n = 600                          # number of Lagrangian floats (0 = none)
seed = 7
delta_i = 0.04                   # kinematic delta_I for moving floats (default: physics.delta_i, or 0.03)
laps = 2.5                       # steady animations: travel of the fastest float, in basin perimeters

[sweep]                          # kind = "sweep"
parameter = "beta_factor"        # delta_s | delta_m | beta_factor | asymmetry | amplitude
start = 0.0
stop = 1.0
n = 60
scale = "linear"                 # linear | log

[animation]
fps = 18
dpi = 90
trail = 10                       # float trail length in frames
width = 11.0                     # inches
height = 6.4

[scales]                         # optional: dimensional reporting (and from_scales)
L = 5.0e6                        # m
beta = 2.0e-11                   # 1/(m s)
r = 1.0e-6                       # 1/s
A_h = 1000.0                     # m^2/s
tau0 = 0.1                       # N/m^2
rho0 = 1025.0                    # kg/m^3
H = 1000.0                       # m
```

## Kinds

`steady`
: Sparse direct solution of the linear problem, compared with the exact
  separable solution whenever one exists (built-in wind; free-slip southern
  and northern walls when $\delta_M > 0$). The animation shows floats
  circulating through the steady gyre.

`transient`
: Time integration from rest (plus optional noise) of the full nonlinear
  equation. The animation shows the streamfunction with floats, the
  vorticity, and time series of energy, transport closure, centroid and
  entropy.

`sweep`
: A steady solution for each value of one parameter. The animation shows
  the streamfunction, the velocity section with the measured e-folding width,
  and the diagnostics as functions of the parameter.

## Choosing parameters

* Keep at least three to four grid points per boundary-layer decay scale
  (reported as `points_per_layer`). For a Munk layer the decay scale is
  $2\delta_M$; for a Stommel layer it is $\delta_S$.
* For eddying runs also watch `grid_scale_fraction`; above 0.3 refine the
  grid or raise `delta_m`.
* Inertial effects appear when $\delta_I$ approaches the frictional widths;
  free-slip walls with $\delta_I > \delta_M$ may run away into a
  basin-filling mode.
* Time is in units of $1/(\beta L)$ (about 2.8 hours for $L = 5000$ km); the
  linear spin-up takes a few Rossby transit times plus the drag time
  $1/\delta_S$.
