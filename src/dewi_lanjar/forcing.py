"""Wind-stress profiles and their curl.

Every built-in profile is a zonal wind stress ``tau_x(y)`` whose curl
``F = -d tau_x / dy`` is a finite sum of sines,
``F(y) = sum_n a_n sin(n pi y / Ly)``. The sine form satisfies the wall
conditions at ``y = 0`` and ``y = Ly`` and lets
:func:`dewi_lanjar.analytic.exact_linear` build an exact steady solution
mode by mode.

Sign convention: a negative curl in the Northern Hemisphere drives an
anticyclonic (clockwise) subtropical gyre with ``psi > 0``.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Callable, Optional

import numpy as np

#: Built-in profiles and their human-readable descriptions.
PROFILES = {
    "single-gyre": "tau = -cos(pi y/Ly): one anticyclonic subtropical gyre (Stommel 1948)",
    "double-gyre": "tau = -cos(2 pi y/Ly): subtropical gyre (south) + subpolar gyre (north)",
    "asymmetric-double-gyre": "double gyre plus a single-gyre component of weight `asymmetry`",
    "triple-gyre": "tau = -cos(3 pi y/Ly): subpolar, subtropical and tropical gyres",
    "cyclonic-gyre": "tau = +cos(pi y/Ly): one cyclonic (subpolar-type) gyre",
}


def _profile_modes(profile: str, Ly: float, asymmetry: float) -> list[tuple[int, float]]:
    """Return ``[(n, a_n), ...]`` with ``F(y) = sum a_n sin(n pi y/Ly)``.

    For ``tau = -cos(n pi y/Ly)`` the curl is
    ``F = -d tau/dy = -(n pi/Ly) sin(n pi y/Ly)``.
    """
    k = math.pi / Ly
    if profile == "single-gyre":
        return [(1, -k)]
    if profile == "cyclonic-gyre":
        return [(1, +k)]
    if profile == "double-gyre":
        return [(2, -2 * k)]
    if profile == "triple-gyre":
        return [(3, -3 * k)]
    if profile == "asymmetric-double-gyre":
        return [(2, -2 * k), (1, -asymmetry * k)]
    raise ValueError(f"unknown wind profile {profile!r}; choose from {sorted(PROFILES)}")


@dataclass(frozen=True)
class WindForcing:
    """Wind-stress curl ``F(x, y, t)`` in nondimensional units.

    ``F = sign * amplitude * F0(y) * (1 + seasonal_amplitude * sin(2 pi t / seasonal_period))``

    Parameters
    ----------
    profile:
        One of :data:`PROFILES`, or ``"custom"`` together with ``curl``.
    amplitude:
        Overall multiplier (1 means the scaling wind stress ``tau0``).
    sign:
        +1 or -1. A sign of -1 mirrors the circulation, which reproduces a
        Southern Hemisphere subtropical gyre (``beta`` stays positive).
    asymmetry:
        Weight of the single-gyre component in ``"asymmetric-double-gyre"``.
    seasonal_amplitude, seasonal_period:
        Optional periodic modulation in nondimensional time.
    curl:
        For ``profile="custom"``: a vectorized callable ``curl(x, y)``
        returning the steady curl. Exact solutions are unavailable for it.
    """

    profile: str = "single-gyre"
    amplitude: float = 1.0
    sign: float = 1.0
    asymmetry: float = 0.5
    seasonal_amplitude: float = 0.0
    seasonal_period: float = 50.0
    Ly: float = 1.0
    curl: Optional[Callable[[np.ndarray, np.ndarray], np.ndarray]] = None

    def __post_init__(self) -> None:
        if self.profile == "custom":
            if self.curl is None:
                raise ValueError("profile='custom' requires a `curl` callable")
        elif self.profile not in PROFILES:
            raise ValueError(f"unknown wind profile {self.profile!r}; choose from {sorted(PROFILES)}")
        if self.sign not in (1, -1, 1.0, -1.0):
            raise ValueError("sign must be +1 or -1")
        if self.seasonal_period <= 0:
            raise ValueError("seasonal_period must be positive")

    # --------------------------------------------------------------
    @property
    def is_sine_series(self) -> bool:
        return self.profile != "custom"

    def modes(self) -> list[tuple[int, float]]:
        """Sine-series coefficients ``[(n, a_n)]`` of the steady curl."""
        if not self.is_sine_series:
            raise ValueError("custom forcing has no sine-series representation")
        s = self.sign * self.amplitude
        return [(n, s * a) for n, a in _profile_modes(self.profile, self.Ly, self.asymmetry)]

    def modulation(self, t: float) -> float:
        if self.seasonal_amplitude == 0.0:
            return 1.0
        return 1.0 + self.seasonal_amplitude * math.sin(2 * math.pi * t / self.seasonal_period)

    def curl_steady(self, x: np.ndarray, y: np.ndarray) -> np.ndarray:
        """Steady part of the curl on arrays ``x, y`` (broadcast together)."""
        x, y = np.broadcast_arrays(np.asarray(x, float), np.asarray(y, float))
        if not self.is_sine_series:
            return self.sign * self.amplitude * np.asarray(self.curl(x, y), float)
        out = np.zeros_like(y)
        for n, a in self.modes():
            out = out + a * np.sin(n * math.pi * y / self.Ly)
        return out

    def curl_at(self, x: np.ndarray, y: np.ndarray, t: float = 0.0) -> np.ndarray:
        return self.curl_steady(x, y) * self.modulation(t)

    def tau_x(self, y: np.ndarray) -> np.ndarray:
        """Zonal wind stress whose curl is the steady ``F`` (zero-mean gauge)."""
        y = np.asarray(y, float)
        if not self.is_sine_series:
            # integrate -F along y at x = 0.5 and remove the mean
            ys = np.linspace(0, self.Ly, 2001)
            F = self.curl_steady(np.full_like(ys, 0.5), ys)
            tau = -np.concatenate([[0.0], np.cumsum(0.5 * (F[1:] + F[:-1]) * np.diff(ys))])
            tau -= tau.mean()
            return np.interp(y, ys, tau)
        out = np.zeros_like(y)
        for n, a in self.modes():
            kn = n * math.pi / self.Ly
            out = out + (a / kn) * np.cos(kn * y)
        return out

    def describe(self) -> str:
        base = PROFILES.get(self.profile, "custom curl")
        extra = []
        if self.sign < 0:
            extra.append("mirrored (Southern Hemisphere sense)")
        if self.amplitude != 1:
            extra.append(f"amplitude x{self.amplitude:g}")
        if self.seasonal_amplitude:
            extra.append(
                f"seasonal modulation {self.seasonal_amplitude:g}, period {self.seasonal_period:g}"
            )
        return base + ("" if not extra else " [" + ", ".join(extra) + "]")
