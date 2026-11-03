"""Plasma-input models with frame integration and explicit blood correction."""
from dataclasses import dataclass
import numpy as np
from scipy.integrate import cumulative_trapezoid
from scipy.optimize import least_squares
from scipy.signal import fftconvolve


def vector(x, name):
    a = np.asarray(x, dtype=float)
    if a.ndim != 1 or not len(a) or not np.all(np.isfinite(a)):
        raise ValueError(f"{name} must be a finite nonempty vector")
    return a


@dataclass(frozen=True)
class Frames:
    start: np.ndarray
    end: np.ndarray

    def __post_init__(self):
        a, b = vector(self.start, "start"), vector(self.end, "end")
        if a.shape != b.shape or np.any(a < 0) or np.any(b <= a) or np.any(a[1:] < b[:-1]-1e-10):
            raise ValueError("Frames must be ordered, nonoverlapping, positive-duration intervals")
        object.__setattr__(self, "start", a.copy())
        object.__setattr__(self, "end", b.copy())

    @property
    def mid(self): return (self.start + self.end) / 2
    @property
    def duration(self): return self.end - self.start


@dataclass(frozen=True)
class InputFunction:
    time: np.ndarray
    plasma: np.ndarray
    blood: np.ndarray | None = None

    def __post_init__(self):
        t, p = vector(self.time, "time"), vector(self.plasma, "plasma")
        if len(t) < 2 or t.shape != p.shape or t[0] != 0 or np.any(np.diff(t) <= 0) or np.any(p < 0):
            raise ValueError("Input must start at zero, increase strictly, and have nonnegative plasma")
        b = None if self.blood is None else vector(self.blood, "blood")
        if b is not None and (b.shape != t.shape or np.any(b < 0)):
            raise ValueError("Blood must match input times and be nonnegative")
        object.__setattr__(self, "time", t.copy())
        object.__setattr__(self, "plasma", p.copy())
        object.__setattr__(self, "blood", None if b is None else b.copy())


PARAMS = {"1t": ("K1", "k2"), "2t": ("K1", "k2", "k3", "k4"), "2t_irreversible": ("K1", "k2", "k3")}


def predict(frames, input_function, parameters, model="2t", vb=0., dt=.025, *, delay=0., input_tail="error"):
    """Return frame means of (1-vb)*(C1+C2)+vb*whole_blood.

    Linear input interpolation, causal trapezoidal convolution, then cumulative
    trapezoidal frame integration. Decrease dt to check numerical convergence.
    """
    if model not in PARAMS: raise ValueError("Unknown model")
    p = vector(parameters, "parameters")
    if len(p) != len(PARAMS[model]) or np.any(p < 0): raise ValueError("Invalid parameters")
    if not 0 <= vb < 1 or not np.isfinite(dt) or dt <= 0: raise ValueError("Invalid vb or dt")
    if vb and input_function.blood is None: raise ValueError("Whole blood is required for vb > 0")
    if not np.isfinite(delay) or input_tail not in ("error", "hold"):
        raise ValueError("Finite delay and input_tail='error' or 'hold' required")
    needed = max(frames.end[-1], frames.end[-1]-delay)
    if input_function.time[-1] < needed and input_tail == "error":
        raise ValueError("Input must cover scan and shifted times; explicitly select input_tail='hold' to extrapolate")
    n = int(np.ceil(frames.end[-1] / dt))
    grid = np.linspace(0, frames.end[-1], n+1); step = grid[1]
    cp = np.interp(grid-delay, input_function.time, input_function.plasma, left=0.)
    K1, k2 = p[:2]
    if model == "1t":
        h = K1 * np.exp(-k2*grid)
    else:
        k3 = p[2]; k4 = p[3] if model == "2t" else 0.
        total = k2+k3+k4
        gap = np.sqrt(max(0., total*total-4*k2*k4))
        slow, fast = (total-gap)/2, (total+gap)/2
        if gap > 1e-10:
            w = (fast-k2)/gap
            h = K1*(w*np.exp(-slow*grid)+(1-w)*np.exp(-fast*grid))
        else:
            h = K1*np.exp(-slow*grid)*(1+(slow-k2)*grid)
    ct = step * (fftconvolve(cp, h)[:n+1] - .5*(cp[0]*h + cp*h[0]))
    if vb:
        ct = (1-vb)*ct + vb*np.interp(grid-delay, input_function.time, input_function.blood, left=0.)
    integral = cumulative_trapezoid(ct, grid, initial=0.)
    return (np.interp(frames.end, grid, integral)-np.interp(frames.start, grid, integral))/frames.duration


def fit_compartment(frames, input_function, tissue, model="2t", sigma=None, vb=0., dt=.025, **options):
    """Multistart WLS with inference QC; see :mod:`petkin.fitting`."""
    from .fitting import fit_compartment as fit
    return fit(frames, input_function, tissue, model, sigma, vb, dt, **options)


def logan(frames, tissue, input_function, tstar, vb=0., **options):
    from .graphical import logan as analyze
    return analyze(frames, tissue, input_function, tstar, vb, **options)
