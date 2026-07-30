"""
llt_gui.py — Limacher Lifting-Line Solver Graphical User Interface

A Windows-friendly tkinter GUI wrapper around the lifting-line solver.
Designed for users who prefer not to edit or run Python scripts directly.

Usage
-----
    python llt_gui.py

Dependencies
------------
    numpy, matplotlib
"""

import sys, io, queue, threading, traceback
from pathlib import Path
import numpy as np

# ══════════════════════════════════════════════════════════════════════════
#  Testable utility functions (no tkinter dependency)
# ══════════════════════════════════════════════════════════════════════════


def parse_geometry_file(path):
    """
    Validate and parse a geometry file.

    Parameters
    ----------
    path : str or Path

    Returns
    -------
    y, c, th : ndarray

    Raises
    ------
    FileNotFoundError, ValueError, OSError
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f'Geometry file not found: {path}')

    try:
        data = np.loadtxt(path, skiprows=1, delimiter=',')
    except (ValueError, OSError) as e:
        raise ValueError(f'Geometry file could not be read: {e}') from e

    if data.ndim == 1:
        data = data.reshape(-1, 1)  # handle single-column edge case badly→caught below

    if data.shape[1] < 3:
        raise ValueError(
            'Expected at least three comma-separated columns: y, c, th.'
        )

    if data.shape[0] < 2:
        raise ValueError(
            'Geometry file must contain at least 2 spanwise stations.'
        )

    y = data[:, 0]
    c = data[:, 1]
    th = data[:, 2]

    if not np.all(np.isfinite(y)) or not np.all(np.isfinite(c)) or not np.all(np.isfinite(th)):
        raise ValueError(
            'Geometry file contains non-finite values in y, c, or th.'
        )

    if len(y) < 2 or np.any(np.diff(y) <= 0):
        raise ValueError(
            'Spanwise coordinates (y) must be strictly increasing.'
        )

    # Check non-zero span
    if np.abs(y[-1] - y[0]) < 1e-15:
        raise ValueError('Span (y range) must be non-zero.')

    # Negative chord
    if np.any(c < 0):
        raise ValueError('Chord values (c) must not be negative.')

    return y, c, th


def parse_airfoil_file(path):
    """
    Validate and parse an airfoil force-data file.

    Parameters
    ----------
    path : str or Path

    Returns
    -------
    alpha, cl, cd : ndarray

    Raises
    ------
    FileNotFoundError, ValueError, OSError
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f'Airfoil data file not found: {path}')

    try:
        data = np.loadtxt(path, skiprows=1, delimiter=',')
    except (ValueError, OSError) as e:
        raise ValueError(f'Airfoil data file could not be read: {e}') from e

    if data.ndim == 1:
        data = data.reshape(-1, 1)

    if data.shape[1] < 3:
        raise ValueError(
            'Expected at least three comma-separated columns: Alpha, Cl, Cd.'
        )

    if data.shape[0] < 2:
        raise ValueError(
            'Airfoil data file must contain at least 2 data rows.'
        )

    alpha = data[:, 0]
    cl = data[:, 1]
    cd = data[:, 2]

    if not np.all(np.isfinite(alpha)) or not np.all(np.isfinite(cl)) or not np.all(np.isfinite(cd)):
        raise ValueError(
            'Airfoil data contains non-finite values in Alpha, Cl, or Cd.'
        )

    if np.any(np.diff(alpha) <= 0):
        raise ValueError(
            'Airfoil data angles (Alpha) must be strictly increasing.'
        )

    return alpha, cl, cd


def build_aoa_array(min_aoa, max_aoa, inc):
    """
    Generate an inclusive angle-of-attack array.

    Parameters
    ----------
    min_aoa, max_aoa, inc : float

    Returns
    -------
    ndarray
        Inclusive angle array. Guarantees min_aoa is included.
        max_aoa is included when it falls within a tolerance.
    """
    # Number of steps, including both ends
    n_steps = int(round((max_aoa - min_aoa) / inc)) + 1
    arr = min_aoa + inc * np.arange(n_steps)

    # Clamp to max_aoa within tolerance
    tol = inc * 1e-10
    arr = arr[arr <= max_aoa + tol]

    # Ensure max_aoa is included if close
    if len(arr) > 0 and abs(arr[-1] - max_aoa) > tol and abs(arr[-1] - max_aoa) < inc * 0.5:
        arr = np.append(arr, max_aoa)

    return arr


def validate_angle_sweep(min_aoa, max_aoa, inc):
    """
    Validate angle-of-attack sweep parameters.

    Parameters
    ----------
    min_aoa, max_aoa, inc : float

    Returns
    -------
    str or None
        Error message string, or None if valid.
    """
    if not np.isfinite(min_aoa) or not np.isfinite(max_aoa) or not np.isfinite(inc):
        return 'All angle values must be finite numbers.'

    if inc <= 0:
        return 'Increment must be greater than zero.'

    if max_aoa < min_aoa:
        return 'Maximum angle must be greater than or equal to minimum angle.'

    n = int(round((max_aoa - min_aoa) / inc)) + 1
    if n > 10000:
        return (
            f'The requested sweep contains {n} angles, '
            f'exceeding the 10,000-angle limit.'
        )

    return None  # valid


# ══════════════════════════════════════════════════════════════════════════
#  Queue-based stdout redirector
# ══════════════════════════════════════════════════════════════════════════


class QueueWriter(io.StringIO):
    """
    A file-like object that writes to both an internal StringIO buffer
    and a thread-safe queue for consumption by the GUI main thread.
    """

    def __init__(self, queue):
        super().__init__()
        self._q = queue

    def write(self, s):
        super().write(s)
        self._q.put(s)

    def flush(self):
        pass


# ══════════════════════════════════════════════════════════════════════════
#  PerformanceMetrics dataclass
# ══════════════════════════════════════════════════════════════════════════


class PerformanceMetrics:
    """
    Container for aerodynamic performance metrics derived from solver output.

    Attributes are None when unavailable.
    """

    def __init__(self):
        self.zero_lift_angle = None        # deg
        self.stall_angle = None             # deg
        self.lift_slope_deg = None          # 1/deg
        self.lift_slope_rad = None          # 1/rad
        self.max_ld = None                  # dimensionless
        self.angle_max_ld = None            # deg
        self.span_efficiency = None         # dimensionless
        self.aspect_ratio = None            # dimensionless
        self.log_messages = []              # str list

    def _fmt(self, val, fmt):
        if val is None or (isinstance(val, float) and np.isnan(val)):
            return 'N/A'
        return fmt.format(val)

    @property
    def zero_lift_str(self):
        return self._fmt(self.zero_lift_angle, '{:.3f}')

    @property
    def stall_str(self):
        return self._fmt(self.stall_angle, '{:.3f}')

    @property
    def lift_slope_deg_str(self):
        return self._fmt(self.lift_slope_deg, '{:.5f}')

    @property
    def lift_slope_rad_str(self):
        return self._fmt(self.lift_slope_rad, '{:.3f}')

    @property
    def max_ld_str(self):
        return self._fmt(self.max_ld, '{:.3f}')

    @property
    def angle_max_ld_str(self):
        return self._fmt(self.angle_max_ld, '{:.3f}')

    @property
    def span_efficiency_str(self):
        return self._fmt(self.span_efficiency, '{:.4f}')

    @property
    def aspect_ratio_str(self):
        return self._fmt(self.aspect_ratio, '{:.3f}')


# ══════════════════════════════════════════════════════════════════════════
#  Aerodynamic performance analysis functions (no tkinter dependency)
# ══════════════════════════════════════════════════════════════════════════


def calculate_aspect_ratio(y_raw, c_raw):
    """
    Calculate aspect ratio from raw (unnormalised) geometry arrays.

    Parameters
    ----------
    y_raw : ndarray   Spanwise coordinate (any units, but must be consistent).
    c_raw : ndarray   Chord (same units as y_raw).

    Returns
    -------
    AR : float
    """
    span = y_raw[-1] - y_raw[0]
    y_norm = y_raw / span
    c_norm = c_raw / span
    S_norm = np.trapezoid(c_norm, y_norm)
    return 1.0 / S_norm


def estimate_zero_lift_angle(CL, AoA):
    """
    Estimate the wing zero-lift angle from CL vs root AoA.

    Parameters
    ----------
    CL : ndarray   Lift coefficient (may contain NaN).
    AoA : ndarray  Root angle of attack (deg), same length.

    Returns
    -------
    float or None   Estimated zero-lift angle (deg), or None if not bracketed.
    """
    valid = np.isfinite(CL) & np.isfinite(AoA)
    if np.sum(valid) < 2:
        return None

    CLv = CL[valid]
    AoAv = AoA[valid]

    # Check if any CL is exactly zero (within tolerance)
    tol = 1e-12
    idx = np.where(np.abs(CLv) < tol)[0]
    if len(idx) > 0:
        # Choose the one nearest zero degrees
        return AoAv[idx[np.argmin(np.abs(AoAv[idx]))]]

    # Find sign changes
    sign_changes = []
    for i in range(len(CLv) - 1):
        if CLv[i] == 0:
            continue
        if CLv[i + 1] == 0:
            continue
        if CLv[i] * CLv[i + 1] < 0:
            # Linear interpolation
            frac = -CLv[i] / (CLv[i + 1] - CLv[i])
            x0 = AoAv[i] + frac * (AoAv[i + 1] - AoAv[i])
            sign_changes.append(x0)

    if not sign_changes:
        return None

    # If multiple crossings, choose the one nearest 0 degrees
    sign_changes = np.array(sign_changes)
    return sign_changes[np.argmin(np.abs(sign_changes))]


def estimate_stall_angle(CL, AoA, alpha_zero, peak_tol=None):
    """
    Estimate the wing stall angle from CL vs root AoA.

    Parameters
    ----------
    CL : ndarray        Lift coefficient.
    AoA : ndarray       Root angle of attack (deg).
    alpha_zero : float  Zero-lift angle (deg), or None.
    peak_tol : float    Tolerance for peak plateau detection.

    Returns
    -------
    float or None       Stall angle (deg), or None if not identified.
    """
    valid = np.isfinite(CL) & np.isfinite(AoA)
    if np.sum(valid) < 3:
        return None

    CLv = CL[valid]
    AoAv = AoA[valid]

    # Determine the positive-lift branch to analyse
    if alpha_zero is not None:
        # Branch at and above alpha_zero
        mask = AoAv >= alpha_zero
        if np.sum(mask) < 2:
            # Fall back to branch containing largest positive CL
            idx_max = np.argmax(CLv)
            mask = np.zeros(len(CLv), dtype=bool)
            # Contiguous region around the max
            start = max(0, idx_max - 5)
            end = min(len(CLv), idx_max + 5)
            mask[start:end] = True
    else:
        # No zero-lift info — use branch containing largest positive CL
        if np.max(CLv) <= 0:
            return None
        idx_max = np.argmax(CLv)
        mask = np.zeros(len(CLv), dtype=bool)
        start = max(0, idx_max - 5)
        end = min(len(CLv), idx_max + 5)
        mask[start:end] = True

    branch_CL = CLv[mask]
    branch_AoA = AoAv[mask]

    if len(branch_CL) < 2:
        return None

    # Find global maximum
    idx_max = np.argmax(branch_CL)
    CL_max = branch_CL[idx_max]

    if peak_tol is None:
        peak_tol = max(0.01, 0.01 * abs(CL_max))

    # Define plateau: points within peak_tol of the maximum
    plateau_mask = np.abs(branch_CL - CL_max) <= peak_tol
    plateau_indices = np.where(plateau_mask)[0]

    # Earliest AoA in the plateau
    stall_angle = branch_AoA[plateau_indices[0]]
    stall_idx_global = np.where(AoAv == stall_angle)[0][0]

    # Check: at least one later finite point with CL below plateau
    later = np.where(AoAv > stall_angle)[0]
    if len(later) == 0:
        return None  # Max at upper boundary

    later_CL = CLv[later]
    if not np.any(later_CL < (CL_max - peak_tol)):
        return None  # No meaningful post-peak decline

    return stall_angle


def estimate_lift_slope(CL, AoA, alpha_zero, alpha_stall):
    """
    Estimate mean lift-curve slope.

    Primary: secant between zero-lift and stall.
    Fallback: OLS linear fit for non-stalling linear data.

    Parameters
    ----------
    CL : ndarray           Lift coefficient.
    AoA : ndarray          Root angle of attack (deg).
    alpha_zero : float or None   Zero-lift angle (deg).
    alpha_stall : float or None  Stall angle (deg).

    Returns
    -------
    (slope_deg, slope_rad) or None
    """
    valid = np.isfinite(CL) & np.isfinite(AoA)
    CLv = CL[valid]
    AoAv = AoA[valid]

    if len(CLv) < 2:
        return None

    # Primary: stall-based secant
    if alpha_stall is not None and alpha_zero is not None and alpha_stall > alpha_zero:
        # Find CL at stall
        stall_idx = np.argmin(np.abs(AoAv - alpha_stall))
        CL_stall = CLv[stall_idx]
        slope_deg = CL_stall / (alpha_stall - alpha_zero)
        slope_rad = slope_deg * 180.0 / np.pi
        return (slope_deg, slope_rad)

    # Fallback: OLS linear fit for non-stalling data
    if alpha_stall is None:
        n = len(CLv)
        if n < 3:
            return None
        A = np.vstack([AoAv, np.ones(n)]).T
        slope_deg, intercept = np.linalg.lstsq(A, CLv, rcond=None)[0]
        CL_fit = slope_deg * AoAv + intercept
        residuals = CLv - CL_fit
        ss_res = np.sum(residuals**2)
        ss_tot = np.sum((CLv - np.mean(CLv))**2)
        if ss_tot == 0:
            return None
        r2 = 1.0 - ss_res / ss_tot
        max_resid = np.max(np.abs(residuals))
        CL_range = np.max(CLv) - np.min(CLv)
        if CL_range > 0 and r2 >= 0.995 and max_resid <= 0.02 * CL_range:
            slope_rad = slope_deg * 180.0 / np.pi
            return (slope_deg, slope_rad)

    return None


def calculate_lift_to_drag(CL, CD):
    """
    Calculate L/D = CL/CD for each point.

    Returns NaN where CD is zero, negative, or CL/CD is non-finite.
    """
    CL = np.asarray(CL, dtype=float)
    CD = np.asarray(CD, dtype=float)
    LD = np.full_like(CL, np.nan)
    mask = np.isfinite(CL) & np.isfinite(CD) & (CD > 1e-15)
    LD[mask] = CL[mask] / CD[mask]
    return LD


def find_max_ld(CL, CD, AoA, alpha_zero):
    """
    Find maximum L/D and its root angle of attack.

    Parameters
    ----------
    CL, CD, AoA : ndarray
    alpha_zero : float or None

    Returns
    -------
    (max_ld, angle) or None
    """
    # Check if ALL section Cd values are effectively zero
    # (This is checked by the caller — we handle it here via the results)
    LD = calculate_lift_to_drag(CL, CD)

    valid = np.isfinite(LD) & np.isfinite(AoA)

    # Restrict to positive-lift region
    if alpha_zero is not None:
        pos_mask = AoA >= alpha_zero
    else:
        pos_mask = CL > 0

    valid = valid & pos_mask & np.isfinite(CL)

    if np.sum(valid) < 2:
        return None

    LDv = LD[valid]
    AoAv = AoA[valid]

    idx_max = np.argmax(LDv)
    max_ld = LDv[idx_max]
    angle = AoAv[idx_max]

    # Check if all finite LD values are at the same angle (edge case)
    return (float(max_ld), float(angle))


def calculate_induced_drag(Gamma, v, y_norm, S_norm):
    """
    Calculate induced drag coefficient from circulation and downwash.

    CDi = (2 / S_norm) * integral(Gamma * sin(v) dy_norm)

    Parameters
    ----------
    Gamma : ndarray   Normalised circulation (1-D, single angle).
    v : ndarray       Normalised downwash v/U (1-D, same length).
    y_norm : ndarray  Normalised spanwise coordinates.
    S_norm : float    Normalised planform area.

    Returns
    -------
    float or NaN
    """
    if np.any(~np.isfinite(Gamma)) or np.any(~np.isfinite(v)):
        return np.nan
    integrand = Gamma * np.sin(v)
    return (2.0 / S_norm) * np.trapezoid(integrand, y_norm)


def calculate_span_efficiency(CL, CDi, AR):
    """
    Calculate span efficiency e = CL^2 / (pi * AR * CDi).

    Returns None if inputs are not finite or meaningful.
    """
    if not np.isfinite(CL) or not np.isfinite(CDi) or not np.isfinite(AR):
        return None
    if abs(CL) < 1e-15 or CDi <= 0:
        return None
    return float(CL**2 / (np.pi * AR * CDi))


def analyze_performance(CL, CD, AoA, y_raw, c_raw, y_norm, c_norm, S_norm,
                        Gamma, v, cd_section):
    """
    Compute all aerodynamic performance metrics from solver output.

    Parameters
    ----------
    CL, CD : ndarray          Wing-level coefficients (1-D, per AoA).
    AoA : ndarray             Root angle of attack (deg).
    y_raw, c_raw : ndarray    Raw (unnormalised) geometry arrays.
    y_norm, c_norm : ndarray  Normalised geometry arrays.
    S_norm : float            Normalised planform area.
    Gamma : ndarray           Circulation (n_y x n_aoa).
    v : ndarray               Downwash (n_y x n_aoa).
    cd_section : ndarray      Section Cd from airfoil data file.

    Returns
    -------
    PerformanceMetrics
    """
    metrics = PerformanceMetrics()

    # Aspect ratio
    span = y_raw[-1] - y_raw[0]
    metrics.aspect_ratio = calculate_aspect_ratio(y_raw, c_raw)

    # Ensure arrays
    CL = np.asarray(CL, dtype=float)
    CD = np.asarray(CD, dtype=float)
    AoA = np.asarray(AoA, dtype=float)

    # Valid-data mask
    valid = np.isfinite(CL) & np.isfinite(CD) & np.isfinite(AoA)
    if np.sum(valid) < 2:
        metrics.log_messages.append('Insufficient valid data for performance analysis.')
        return metrics

    CLv = CL[valid]
    CDv = CD[valid]
    AoAv = AoA[valid]

    # ── Zero-lift angle ─────────────────────────────────────────────
    metrics.zero_lift_angle = estimate_zero_lift_angle(CL, AoA)
    if metrics.zero_lift_angle is None:
        metrics.log_messages.append(
            'Zero-lift angle not bracketed by the requested angle sweep.'
        )

    # ── Stall angle ─────────────────────────────────────────────────
    metrics.stall_angle = estimate_stall_angle(CL, AoA, metrics.zero_lift_angle)
    if metrics.stall_angle is None:
        metrics.log_messages.append(
            'Stall not observed within the requested angle sweep.'
        )

    # ── Lift slope ──────────────────────────────────────────────────
    slopes = estimate_lift_slope(CL, AoA, metrics.zero_lift_angle,
                                 metrics.stall_angle)
    if slopes is not None:
        metrics.lift_slope_deg, metrics.lift_slope_rad = slopes
        if metrics.stall_angle is None:
            metrics.log_messages.append(
                'Linear lift-curve slope; no stall detected.'
            )
    else:
        if metrics.stall_angle is not None:
            metrics.log_messages.append(
                'Lift slope unavailable: insufficient data for secant calculation.'
            )
        else:
            metrics.log_messages.append(
                'Lift slope unavailable: no stall identified and data not '
                'sufficiently linear for OLS fit.'
            )

    # ── L/D and maximum L/D ─────────────────────────────────────────
    # Check if all section Cd values are effectively zero
    cd_section = np.asarray(cd_section, dtype=float)
    all_zero_profile_drag = (
        np.all(np.isfinite(cd_section)) and np.all(np.abs(cd_section) < 1e-15)
    )

    if all_zero_profile_drag:
        metrics.max_ld = None
        metrics.angle_max_ld = None
        metrics.log_messages.append(
            'Maximum L/D is not defined meaningfully for '
            'zero-profile-drag input data.'
        )
    else:
        max_ld_result = find_max_ld(CL, CD, AoA, metrics.zero_lift_angle)
        if max_ld_result is not None:
            metrics.max_ld, metrics.angle_max_ld = max_ld_result
            # Check if max occurs at sweep boundary
            if metrics.angle_max_ld is not None:
                tol = 0.5 * np.mean(np.diff(AoA)) if len(AoA) > 1 else 1.0
                if abs(metrics.angle_max_ld - AoA[-1]) < tol:
                    metrics.log_messages.append(
                        'Maximum L/D occurs at the upper boundary of the '
                        'requested sweep; the true optimum may lie outside '
                        'the requested range.'
                    )
        else:
            metrics.log_messages.append(
                'Maximum L/D not identified.'
            )

    # ── Span efficiency ─────────────────────────────────────────────
    # Find the operating point: angle at max L/D, or fallback to
    # positive-lift pre-stall point nearest CL=0.5
    op_angle = metrics.angle_max_ld
    op_source = 'max L/D'

    if op_angle is None:
        # Fallback: find positive-lift point nearest CL=0.5
        if metrics.stall_angle is not None:
            stall_mask = AoAv <= metrics.stall_angle
        else:
            stall_mask = np.ones(len(AoAv), dtype=bool)

        pos_mask = (CLv > 0) & stall_mask
        if np.any(pos_mask):
            target = 0.5
            idx = np.argmin(np.abs(CLv[pos_mask] - target))
            op_angle = AoAv[pos_mask][idx]
            op_source = f'CL ≈ {CLv[pos_mask][idx]:.3f} (fallback)'
            metrics.log_messages.append(
                f'Span efficiency evaluated at {op_source}.'
            )

    if op_angle is not None and Gamma is not None and v is not None:
        # Find the closest AoA index
        aoa_idx = np.argmin(np.abs(AoA - op_angle))
        Gamma_op = Gamma[:, aoa_idx] if Gamma.ndim == 2 else Gamma
        v_op = v[:, aoa_idx] if v.ndim == 2 else v
        CL_op = CL[aoa_idx]
        CDi_op = calculate_induced_drag(Gamma_op, v_op, y_norm, S_norm)
        metrics.span_efficiency = calculate_span_efficiency(
            CL_op, CDi_op, metrics.aspect_ratio)

    return metrics


# ══════════════════════════════════════════════════════════════════════════
#  GUI application (tkinter)
# ══════════════════════════════════════════════════════════════════════════

# tkinter import deferred so importing llt_gui for tests doesn't require a
# display.  The 'import llt_gui' above only loads the utility functions.
# The GUI class is instantiated only when run via __main__ below.


class LLTGuiApp:
    """
    Main GUI application for the Limacher Lifting-Line Solver.
    """

    def __init__(self, root):
        self.root = root
        self._setup_complete = False

        # Application state
        self.geom_path = None
        self.airfoil_path = None
        self.last_results = None  # (CL, CD, y, Gamma, v, AoA)
        self._geom_y_raw = None
        self._geom_c_raw = None
        self._geom_y_norm = None
        self._geom_c_norm = None
        self._geom_S_norm = None
        self._cd_section = None
        self._metrics = None

        # Thread safety
        self._msg_queue = queue.Queue()
        self._worker_thread = None
        self._sim_running = False
        self._polling = True

        # Defer tkinter imports — they require the event loop
        self._init_ui()

    def _init_ui(self):
        """Set up all tkinter widgets (called after tkinter is imported)."""
        import tkinter as tk
        from tkinter import ttk, filedialog, messagebox
        from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
        import matplotlib
        matplotlib.use('TkAgg')

        self.tk = tk
        self.filedialog = filedialog
        self.messagebox = messagebox
        self.FigureCanvasTkAgg = FigureCanvasTkAgg

        self.root.title('Limacher Lifting-Line Solver')
        self.root.minsize(1200, 760)

        # ── Overall layout: 2x2 grid ─────────────────────────────────
        self.root.grid_rowconfigure(0, weight=35)    # top row 35%
        self.root.grid_rowconfigure(1, weight=65)    # bottom row 65%
        self.root.grid_columnconfigure(0, weight=40, uniform='col')   # left 40%
        self.root.grid_columnconfigure(1, weight=60, uniform='col')   # right 60%

        self._build_inputs_panel()
        self._build_geometry_panel()
        self._build_updates_panel()
        self._build_results_panel()

        # Check for default files
        self._load_default_files()

        self._setup_complete = True

    # ── Panel builders ──────────────────────────────────────────────

    def _build_inputs_panel(self):
        """Top-left: Inputs and Run controls."""
        from tkinter import ttk
        frame = ttk.LabelFrame(self.root, text='Inputs', padding=8)
        frame.grid(row=0, column=0, sticky='nsew', padx=4, pady=4)
        frame.grid_rowconfigure(4, weight=1)
        frame.grid_columnconfigure(1, weight=1)

        # Geometry file
        ttk.Label(frame, text='Wing geometry file').grid(
            row=0, column=0, sticky='w', pady=(0, 2))
        self._geom_var = self.tk.StringVar()
        geom_entry = ttk.Entry(frame, textvariable=self._geom_var, state='readonly')
        geom_entry.grid(row=0, column=1, sticky='ew', padx=(2, 4))
        ttk.Button(frame, text='Browse', command=self.browse_geometry_file).grid(
            row=0, column=2, padx=(0, 2))

        # Airfoil data file
        ttk.Label(frame, text='Airfoil data file').grid(
            row=1, column=0, sticky='w', pady=(4, 2))
        self._airfoil_var = self.tk.StringVar()
        af_entry = ttk.Entry(frame, textvariable=self._airfoil_var, state='readonly')
        af_entry.grid(row=1, column=1, sticky='ew', padx=(2, 4))
        ttk.Button(frame, text='Browse', command=self.browse_airfoil_file).grid(
            row=1, column=2, padx=(0, 2))

        # Angle sweep
        sweep_frame = ttk.Frame(frame)
        sweep_frame.grid(row=2, column=0, columnspan=3, sticky='ew', pady=(8, 2))
        ttk.Label(sweep_frame, text='Minimum angle of attack (deg)').grid(
            row=0, column=0, sticky='w', padx=(0, 4))
        self._min_aoa_var = self.tk.StringVar(value='-10')
        ttk.Entry(sweep_frame, textvariable=self._min_aoa_var, width=8).grid(
            row=0, column=1, padx=(0, 12))

        ttk.Label(sweep_frame, text='Maximum angle of attack (deg)').grid(
            row=0, column=2, sticky='w', padx=(0, 4))
        self._max_aoa_var = self.tk.StringVar(value='10')
        ttk.Entry(sweep_frame, textvariable=self._max_aoa_var, width=8).grid(
            row=0, column=3, padx=(0, 12))

        ttk.Label(sweep_frame, text='Increment (deg)').grid(
            row=0, column=4, sticky='w', padx=(0, 4))
        self._inc_var = self.tk.StringVar(value='1')
        ttk.Entry(sweep_frame, textvariable=self._inc_var, width=8).grid(
            row=0, column=5)

        # Angle count display
        self._angle_count_var = self.tk.StringVar(value='Number of requested angles: 21')
        ttk.Label(frame, textvariable=self._angle_count_var).grid(
            row=3, column=0, columnspan=3, sticky='w', pady=(4, 8))

        # Trace angle entries for live count update
        self._min_aoa_var.trace_add('write', self._update_angle_count)
        self._max_aoa_var.trace_add('write', self._update_angle_count)
        self._inc_var.trace_add('write', self._update_angle_count)

        # Run button + status
        run_frame = ttk.Frame(frame)
        run_frame.grid(row=4, column=0, columnspan=3, sticky='ew', pady=(4, 0))
        self._run_btn = ttk.Button(
            run_frame, text='Run simulation', command=self.start_simulation)
        self._run_btn.pack(side='left', padx=(0, 8))
        self._status_var = self.tk.StringVar(value='Ready')
        ttk.Label(run_frame, textvariable=self._status_var).pack(side='left')

    def _build_geometry_panel(self):
        """Bottom-left: Geometry preview."""
        from tkinter import ttk
        from matplotlib.figure import Figure
        frame = ttk.LabelFrame(self.root, text='Geometry preview', padding=4)
        frame.grid(row=1, column=0, sticky='nsew', padx=4, pady=4)
        frame.grid_rowconfigure(0, weight=1)
        frame.grid_columnconfigure(0, weight=1)

        self._geom_fig = Figure(figsize=(4, 3), dpi=100, tight_layout=True)
        self._geom_ax_chord = self._geom_fig.add_subplot(2, 1, 1)
        self._geom_ax_twist = self._geom_fig.add_subplot(2, 1, 2, sharex=self._geom_ax_chord)

        self._geom_canvas = self.FigureCanvasTkAgg(
            self._geom_fig, master=frame)
        self._geom_canvas.get_tk_widget().grid(row=0, column=0, sticky='nsew')

        self._geom_ax_chord.set_ylabel('Chord, c')
        self._geom_ax_twist.set_ylabel('Twist, th (deg)')
        self._geom_ax_twist.set_xlabel('Spanwise coordinate, y')
        self._geom_ax_chord.grid(True, alpha=0.3)
        self._geom_ax_twist.grid(True, alpha=0.3)

    def _build_updates_panel(self):
        """Top-right: Solver updates and errors."""
        from tkinter import ttk
        frame = ttk.LabelFrame(self.root, text='Solver updates', padding=4)
        frame.grid(row=0, column=1, sticky='nsew', padx=4, pady=4)
        frame.grid_rowconfigure(0, weight=1)
        frame.grid_columnconfigure(0, weight=1)

        self._log_text = self.tk.Text(frame, wrap='word', state='disabled',
                                       height=10, relief='sunken', borderwidth=2)
        self._log_text.grid(row=0, column=0, sticky='nsew')

        scrollbar = ttk.Scrollbar(frame, orient='vertical',
                                   command=self._log_text.yview)
        scrollbar.grid(row=0, column=1, sticky='ns')
        self._log_text['yscrollcommand'] = scrollbar.set

        ttk.Button(frame, text='Clear log', command=self.clear_log).grid(
            row=1, column=0, columnspan=2, pady=(4, 0))

    def _build_results_panel(self):
        """Bottom-right: Performance summary + results with notebook tabs."""
        from tkinter import ttk
        from matplotlib.figure import Figure
        frame = ttk.LabelFrame(self.root, text='Results', padding=4)
        frame.grid(row=1, column=1, sticky='nsew', padx=4, pady=4)
        frame.grid_rowconfigure(1, weight=1)
        frame.grid_columnconfigure(0, weight=1)

        # ── Scalar performance summary ──────────────────────────────
        summary_frame = ttk.Frame(frame)
        summary_frame.grid(row=0, column=0, sticky='ew', pady=(0, 4))
        summary_frame.grid_columnconfigure(1, weight=1)
        summary_frame.grid_columnconfigure(3, weight=1)

        labels_data = [
            ('Zero-lift angle, α_L=0 (deg):', 0, 0,
             'zero_lift_str', 0, 1),
            ('Stall angle, α_stall (deg):', 0, 2,
             'stall_str', 0, 3),
            ('Mean lift slope, CL_α (1/deg):', 1, 0,
             'lift_slope_deg_str', 1, 1),
            ('Mean lift slope, CL_α (1/rad):', 1, 2,
             'lift_slope_rad_str', 1, 3),
            ('Maximum L/D:', 2, 0,
             'max_ld_str', 2, 1),
            ('Angle at max L/D (deg):', 2, 2,
             'angle_max_ld_str', 2, 3),
            ('Span efficiency at max L/D, e:', 3, 0,
             'span_efficiency_str', 3, 1),
        ]

        for text, r, c, attr, vr, vc in labels_data:
            ttk.Label(summary_frame, text=text).grid(
                row=r, column=c, sticky='w', padx=(0, 2))
            var = self.tk.StringVar(value='N/A')
            ttk.Label(summary_frame, textvariable=var,
                       font=('Segoe UI', 9, 'bold')).grid(
                row=vr, column=vc, sticky='w', padx=(0, 16))
            setattr(self, f'_perf_{attr}', var)

        # ── Notebook with plot tabs ─────────────────────────────────
        notebook = ttk.Notebook(frame)
        notebook.grid(row=1, column=0, sticky='nsew')

        # Lift curve tab
        lift_frame = ttk.Frame(notebook)
        notebook.add(lift_frame, text='Lift curve')
        lift_frame.grid_rowconfigure(0, weight=1)
        lift_frame.grid_columnconfigure(0, weight=1)

        self._result_fig_lift = Figure(figsize=(5, 2.8), dpi=100, tight_layout=True)
        self._result_ax_lift = self._result_fig_lift.add_subplot(111)
        self._result_canvas_lift = self.FigureCanvasTkAgg(
            self._result_fig_lift, master=lift_frame)
        self._result_canvas_lift.get_tk_widget().grid(row=0, column=0, sticky='nsew')
        self._result_ax_lift.set_xlabel('Root angle of attack (deg)')
        self._result_ax_lift.set_ylabel('Lift coefficient, CL')
        self._result_ax_lift.grid(True, alpha=0.3)

        # Drag polar tab
        drag_frame = ttk.Frame(notebook)
        notebook.add(drag_frame, text='Drag polar')
        drag_frame.grid_rowconfigure(0, weight=1)
        drag_frame.grid_columnconfigure(0, weight=1)

        self._result_fig_drag = Figure(figsize=(5, 2.8), dpi=100, tight_layout=True)
        self._result_ax_drag = self._result_fig_drag.add_subplot(111)
        self._result_canvas_drag = self.FigureCanvasTkAgg(
            self._result_fig_drag, master=drag_frame)
        self._result_canvas_drag.get_tk_widget().grid(row=0, column=0, sticky='nsew')
        self._result_ax_drag.set_xlabel('Lift coefficient, CL')
        self._result_ax_drag.set_ylabel('Drag coefficient, CD')
        self._result_ax_drag.grid(True, alpha=0.3)

        # L/D tab
        ld_frame = ttk.Frame(notebook)
        notebook.add(ld_frame, text='L/D ratio')
        ld_frame.grid_rowconfigure(0, weight=1)
        ld_frame.grid_columnconfigure(0, weight=1)

        self._result_fig_ld = Figure(figsize=(5, 2.8), dpi=100, tight_layout=True)
        self._result_ax_ld = self._result_fig_ld.add_subplot(111)
        self._result_canvas_ld = self.FigureCanvasTkAgg(
            self._result_fig_ld, master=ld_frame)
        self._result_canvas_ld.get_tk_widget().grid(row=0, column=0, sticky='nsew')
        self._result_ax_ld.set_xlabel('Root angle of attack (deg)')
        self._result_ax_ld.set_ylabel('Lift-to-drag ratio, L/D')
        self._result_ax_ld.grid(True, alpha=0.3)

    # ── Default file loading ────────────────────────────────────────

    def _load_default_files(self):
        """Auto-select example files if found alongside llt_gui.py."""
        script_dir = Path(sys.argv[0] if getattr(sys, 'frozen', False) else __file__).parent
        geom_default = script_dir / 'elliptic_AR8.txt'
        airfoil_default = script_dir / 'thin_airfoil_data.txt'

        if geom_default.exists():
            self._select_geometry_file(str(geom_default))
        if airfoil_default.exists():
            self._select_airfoil_file(str(airfoil_default))

    # ── File browsing ───────────────────────────────────────────────

    def browse_geometry_file(self):
        path = self.filedialog.askopenfilename(
            title='Select wing geometry file',
            filetypes=[('CSV/TXT files', '*.txt *.csv'), ('All files', '*.*')])
        if path:
            self._select_geometry_file(path)

    def browse_airfoil_file(self):
        path = self.filedialog.askopenfilename(
            title='Select airfoil force-data file',
            filetypes=[('CSV/TXT files', '*.txt *.csv'), ('All files', '*.*')])
        if path:
            self._select_airfoil_file(path)

    def _select_geometry_file(self, path):
        self._geom_var.set(path)
        self.geom_path = path
        self.load_and_plot_geometry()

    def _select_airfoil_file(self, path):
        self._airfoil_var.set(path)
        self.airfoil_path = path
        self._validate_and_report_airfoil()

    # ── File validation + geometry plotting ─────────────────────────

    def load_and_plot_geometry(self):
        """Parse geometry file and update geometry preview."""
        path = self._geom_var.get()
        if not path:
            return
        try:
            y, c, th = parse_geometry_file(path)
        except (FileNotFoundError, ValueError, OSError) as e:
            self.append_log(f'Error: {e}')
            # Clear geometry plots
            self._geom_ax_chord.clear()
            self._geom_ax_twist.clear()
            self._update_geom_canvas()
            return

        self.geom_path = path
        self._geom_y_raw = y.copy()
        self._geom_c_raw = c.copy()
        # Store normalised geometry for performance analysis
        span = y[-1] - y[0]
        self._geom_y_norm = y / span
        self._geom_c_norm = c / span
        self._geom_S_norm = np.trapezoid(self._geom_c_norm, self._geom_y_norm)

        self.update_geometry_plots(y, c, th)
        self.append_log(f'Geometry file loaded: {Path(path).name}')

    def _validate_and_report_airfoil(self):
        """Parse airfoil data and log the Alpha range."""
        path = self._airfoil_var.get()
        if not path:
            return
        try:
            alpha, cl, cd = parse_airfoil_file(path)
        except (FileNotFoundError, ValueError, OSError) as e:
            self.append_log(f'Error: {e}')
            return

        self.airfoil_path = path
        self._cd_section = cd.copy()
        self.append_log(
            f'Airfoil data loaded: {Path(path).name}  '
            f'(Alpha range: {alpha[0]:.1f}° to {alpha[-1]:.1f}°, '
            f'{len(alpha)} points)'
        )

    def update_geometry_plots(self, y, c, th):
        """Redraw chord and twist preview plots."""
        self._geom_ax_chord.clear()
        self._geom_ax_twist.clear()

        self._geom_ax_chord.plot(y, c, 'o-', markersize=3, linewidth=1.2)
        self._geom_ax_chord.set_ylabel('Chord, c')
        self._geom_ax_chord.grid(True, alpha=0.3)

        self._geom_ax_twist.plot(y, th, 'o-', markersize=3, linewidth=1.2)
        self._geom_ax_twist.set_ylabel('Twist, th (deg)')
        self._geom_ax_twist.set_xlabel('Spanwise coordinate, y')
        self._geom_ax_twist.grid(True, alpha=0.3)

        self._update_geom_canvas()

    def _update_geom_canvas(self):
        self._geom_canvas.draw_idle()

    # ── Angle count update ──────────────────────────────────────────

    def _update_angle_count(self, *_):
        """Update the 'Number of requested angles' label."""
        try:
            min_aoa = float(self._min_aoa_var.get())
            max_aoa = float(self._max_aoa_var.get())
            inc = float(self._inc_var.get())
            if inc > 0 and max_aoa >= min_aoa:
                arr = build_aoa_array(min_aoa, max_aoa, inc)
                self._angle_count_var.set(f'Number of requested angles: {len(arr)}')
                return
        except (ValueError, ZeroDivisionError):
            pass
        self._angle_count_var.set('Number of requested angles: —')

    # ── Logging ─────────────────────────────────────────────────────

    def append_log(self, message):
        """Append a message to the solver updates text widget."""
        self._log_text.configure(state='normal')
        self._log_text.insert('end', message + '\n')
        self._log_text.see('end')
        self._log_text.configure(state='disabled')
        self.root.update_idletasks()

    def clear_log(self):
        self._log_text.configure(state='normal')
        self._log_text.delete('1.0', 'end')
        self._log_text.configure(state='disabled')

    # ── Simulation ──────────────────────────────────────────────────

    def start_simulation(self):
        """Validate inputs and start the solver in a worker thread."""
        if self._sim_running:
            return

        # Validate geometry file
        geom_path = self._geom_var.get()
        try:
            parse_geometry_file(geom_path)
        except (FileNotFoundError, ValueError, OSError) as e:
            self.append_log(f'Error: {e}')
            return

        # Validate airfoil file
        airfoil_path = self._airfoil_var.get()
        try:
            parse_airfoil_file(airfoil_path)
        except (FileNotFoundError, ValueError, OSError) as e:
            self.append_log(f'Error: {e}')
            return

        # Validate angle sweep
        try:
            min_aoa = float(self._min_aoa_var.get())
            max_aoa = float(self._max_aoa_var.get())
            inc = float(self._inc_var.get())
        except ValueError:
            self.append_log('Error: All angle values must be valid numbers.')
            return

        err = validate_angle_sweep(min_aoa, max_aoa, inc)
        if err is not None:
            self.append_log(f'Error: {err}')
            return

        # Everything valid — build AoA array
        AoA = build_aoa_array(min_aoa, max_aoa, inc)

        # Clear previous result plots
        self._result_ax_lift.clear()
        self._result_ax_lift.set_xlabel('Root angle of attack (deg)')
        self._result_ax_lift.set_ylabel('Lift coefficient, CL')
        self._result_ax_lift.grid(True, alpha=0.3)
        self._result_canvas_lift.draw_idle()

        self._result_ax_drag.clear()
        self._result_ax_drag.set_xlabel('Lift coefficient, CL')
        self._result_ax_drag.set_ylabel('Drag coefficient, CD')
        self._result_ax_drag.grid(True, alpha=0.3)
        self._result_canvas_drag.draw_idle()

        self._result_ax_ld.clear()
        self._result_ax_ld.set_xlabel('Root angle of attack (deg)')
        self._result_ax_ld.set_ylabel('Lift-to-drag ratio, L/D')
        self._result_ax_ld.grid(True, alpha=0.3)
        self._result_canvas_ld.draw_idle()

        # Reset performance summary
        self._metrics = None
        for attr in ['zero_lift_str', 'stall_str', 'lift_slope_deg_str',
                      'lift_slope_rad_str', 'max_ld_str', 'angle_max_ld_str',
                      'span_efficiency_str']:
            getattr(self, f'_perf_{attr}').set('N/A')

        # Log start
        self.append_log(
            f'\n--- Simulation: AoA ({min_aoa}° to {max_aoa}°, '
            f'inc {inc}°, {len(AoA)} angles) ---'
        )

        # UI state
        self._sim_running = True
        self._run_btn.configure(state='disabled')
        self._status_var.set('Running')
        self._polling = True

        # Set up message queue and stdout redirect
        self._msg_queue = queue.Queue()
        self._captured_stdout = io.StringIO()

        # Run solver in thread
        self._worker_thread = threading.Thread(
            target=self._simulation_worker,
            args=(geom_path, airfoil_path, AoA),
            daemon=True
        )
        self._worker_thread.start()

        # Start polling the message queue
        self._poll_queue()

    def _simulation_worker(self, geom_path, airfoil_path, AoA):
        """Run the solver (worker thread)."""
        from liftingline import liftingline

        # Redirect stdout to queue
        old_stdout = sys.stdout
        sys.stdout = QueueWriter(self._msg_queue)

        try:
            CL, CD, y, Gamma, v = liftingline(
                geom_path, airfoil_path, AoA=AoA, relfactor=0.01
            )
            # Put result back via queue
            self._msg_queue.put(('__RESULT__', (CL, CD, y, Gamma, v, AoA)))
        except Exception as exc:
            self._msg_queue.put(('__ERROR__', (type(exc).__name__, str(exc),
                                               traceback.format_exc())))
        finally:
            sys.stdout = old_stdout
            self._msg_queue.put('__DONE__')

    def _poll_queue(self):
        """Poll the message queue from the main thread (called via root.after)."""
        if not self._polling:
            return

        try:
            while True:
                msg = self._msg_queue.get_nowait()

                if isinstance(msg, tuple):
                    if msg[0] == '__RESULT__':
                        self._handle_simulation_success(*msg[1])
                        continue
                    elif msg[0] == '__ERROR__':
                        self._handle_simulation_failure(*msg[1])
                        continue
                elif msg == '__DONE__':
                    continue

                # Regular stdout text
                self.append_log(msg.rstrip())

        except queue.Empty:
            pass

        # Check if thread is still alive
        if self._worker_thread is not None and self._worker_thread.is_alive():
            self.root.after(100, self._poll_queue)
        else:
            self._finalize_run()

    def _finalize_run(self):
        """Clean up after simulation thread finishes."""
        self._sim_running = False
        self._worker_thread = None
        self._run_btn.configure(state='normal')
        self.root.update_idletasks()

    def _handle_simulation_success(self, CL, CD, y, Gamma, v, AoA):
        """Process successful simulation results."""
        self.last_results = (CL, CD, y, Gamma, v, AoA)

        # Check for invalid (NaN) points
        valid_mask = ~(np.isnan(CL) | np.isnan(CD))
        n_nan = np.sum(~valid_mask)

        if n_nan > 0:
            nan_angles = AoA[~valid_mask]
            self.append_log(
                f'\nWarning: {n_nan} invalid point(s) at AoA = '
                f'{", ".join(f"{a:.2f}°" for a in nan_angles)}'
            )
            self._status_var.set('Completed with invalid points')
        else:
            self._status_var.set('Completed')
            self.append_log('\nSimulation completed successfully.')

        # ── Performance analysis ────────────────────────────────────
        self._metrics = None
        if self._geom_y_raw is not None and self._cd_section is not None:
            try:
                self._metrics = analyze_performance(
                    CL, CD, AoA,
                    self._geom_y_raw, self._geom_c_raw,
                    self._geom_y_norm, self._geom_c_norm,
                    self._geom_S_norm,
                    Gamma, v, self._cd_section,
                )
            except Exception as exc:
                self.append_log(f'Performance analysis error: {exc}')

        # Update summary display
        self._update_performance_summary()

        # Log performance summary
        if self._metrics is not None:
            m = self._metrics
            self.append_log('')
            self.append_log('Performance summary')
            self.append_log('-------------------')
            zl = m.zero_lift_str
            sa = m.stall_str
            sd = m.lift_slope_deg_str
            sr = m.lift_slope_rad_str
            ml = m.max_ld_str
            am = m.angle_max_ld_str
            se = m.span_efficiency_str
            self.append_log(f'Zero-lift angle: {zl} deg')
            self.append_log(f'Stall angle: {sa} deg')
            self.append_log(f'Mean lift slope: {sd} 1/deg = {sr} 1/rad')
            self.append_log(f'Maximum L/D: {ml} at {am} deg')
            self.append_log(f'Span efficiency at maximum L/D: {se}')
            for msg in m.log_messages:
                self.append_log(f'  Note: {msg}')

        # ── Update plots ────────────────────────────────────────────
        self._update_result_plots(CL, CD, AoA)

    def _update_result_plots(self, CL, CD, AoA):
        """Update all three result plot tabs with markers."""
        LD = calculate_lift_to_drag(CL, CD)

        # Lift curve
        self._result_ax_lift.clear()
        mask = np.isfinite(CL)
        if np.any(mask):
            self._result_ax_lift.plot(
                AoA[mask], CL[mask], 'o-', markersize=4, linewidth=1.2)

        # Zero-lift marker
        if self._metrics and self._metrics.zero_lift_angle is not None:
            zl = self._metrics.zero_lift_angle
            self._result_ax_lift.axvline(
                zl, color='green', linestyle='--', alpha=0.6, linewidth=1,
                label=f'α_L=0 = {zl:.2f}°')

        # Stall marker
        if self._metrics and self._metrics.stall_angle is not None:
            sa = self._metrics.stall_angle
            # Find CL at stall
            idx = np.argmin(np.abs(AoA - sa))
            CL_stall = CL[idx]
            self._result_ax_lift.plot(
                sa, CL_stall, 'ro', markersize=8, alpha=0.8,
                label=f'Stall ≈ {sa:.1f}°')

        self._result_ax_lift.set_xlabel('Root angle of attack (deg)')
        self._result_ax_lift.set_ylabel('Lift coefficient, CL')
        self._result_ax_lift.grid(True, alpha=0.3)
        if self._result_ax_lift.get_legend_handles_labels()[0]:
            self._result_ax_lift.legend(fontsize=8)
        self._result_ax_lift.autoscale()
        self._result_canvas_lift.draw_idle()

        # Drag polar
        self._result_ax_drag.clear()
        mask = np.isfinite(CL) & np.isfinite(CD)
        if np.any(mask):
            self._result_ax_drag.plot(
                CL[mask], CD[mask], 'o-', markersize=4, linewidth=1.2)
        self._result_ax_drag.set_xlabel('Lift coefficient, CL')
        self._result_ax_drag.set_ylabel('Drag coefficient, CD')
        self._result_ax_drag.grid(True, alpha=0.3)
        self._result_ax_drag.autoscale()
        self._result_canvas_drag.draw_idle()

        # L/D ratio
        self._result_ax_ld.clear()
        mask = np.isfinite(LD)
        if np.any(mask):
            self._result_ax_ld.plot(
                AoA[mask], LD[mask], 'o-', markersize=4, linewidth=1.2)

        # Max L/D marker
        if self._metrics and self._metrics.max_ld is not None:
            ml = self._metrics.max_ld
            am = self._metrics.angle_max_ld
            self._result_ax_ld.plot(
                am, ml, 'ro', markersize=8, alpha=0.8,
                label=f'Max L/D = {ml:.1f} at {am:.1f}°')
            if self._result_ax_ld.get_legend_handles_labels()[0]:
                self._result_ax_ld.legend(fontsize=8)

        self._result_ax_ld.set_xlabel('Root angle of attack (deg)')
        self._result_ax_ld.set_ylabel('Lift-to-drag ratio, L/D')
        self._result_ax_ld.grid(True, alpha=0.3)
        self._result_ax_ld.autoscale()
        self._result_canvas_ld.draw_idle()

    def _update_performance_summary(self):
        """Update the scalar performance summary display."""
        if self._metrics is None:
            return
        m = self._metrics
        self._perf_zero_lift_str.set(m.zero_lift_str)
        self._perf_stall_str.set(m.stall_str)
        self._perf_lift_slope_deg_str.set(m.lift_slope_deg_str)
        self._perf_lift_slope_rad_str.set(m.lift_slope_rad_str)
        self._perf_max_ld_str.set(m.max_ld_str)
        self._perf_angle_max_ld_str.set(m.angle_max_ld_str)
        self._perf_span_efficiency_str.set(m.span_efficiency_str)

    def _handle_simulation_failure(self, exc_type, exc_msg, exc_tb):
        """Process simulation failure."""
        self.append_log(f'\nSimulation failed: {exc_type} — {exc_msg}')
        self.append_log(f'Traceback:\n{exc_tb}')
        self._status_var.set('Failed')

        # Clear result axes
        self._result_ax_lift.clear()
        self._result_ax_lift.set_xlabel('Root angle of attack (deg)')
        self._result_ax_lift.set_ylabel('Lift coefficient, CL')
        self._result_ax_lift.grid(True, alpha=0.3)
        self._result_canvas_lift.draw_idle()

        self._result_ax_drag.clear()
        self._result_ax_drag.set_xlabel('Lift coefficient, CL')
        self._result_ax_drag.set_ylabel('Drag coefficient, CD')
        self._result_ax_drag.grid(True, alpha=0.3)
        self._result_canvas_drag.draw_idle()

        self._result_ax_ld.clear()
        self._result_ax_ld.set_xlabel('Root angle of attack (deg)')
        self._result_ax_ld.set_ylabel('Lift-to-drag ratio, L/D')
        self._result_ax_ld.grid(True, alpha=0.3)
        self._result_canvas_ld.draw_idle()

        # Reset performance summary
        self._metrics = None
        for attr in ['zero_lift_str', 'stall_str', 'lift_slope_deg_str',
                      'lift_slope_rad_str', 'max_ld_str', 'angle_max_ld_str',
                      'span_efficiency_str']:
            getattr(self, f'_perf_{attr}').set('N/A')

    def on_close(self):
        """Clean shutdown."""
        self._polling = False
        self.root.destroy()


# ══════════════════════════════════════════════════════════════════════════
#  Entry point
# ══════════════════════════════════════════════════════════════════════════

def main():
    import tkinter as tk
    root = tk.Tk()
    app = LLTGuiApp(root)
    root.protocol('WM_DELETE_WINDOW', app.on_close)
    root.mainloop()


if __name__ == '__main__':
    main()