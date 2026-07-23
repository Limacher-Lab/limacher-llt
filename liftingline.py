"""
liftingline.py — Prandtl Lifting-Line Theory solver

Port of the MATLAB implementation (limacher-llt/liftingline.m) to Python.

Functions
---------
liftingline(geomfile, forcefile, AoA=0, relfactor=0.01)
    Main lifting-line solver. Returns CL, CD, y, Gamma, v.

LLsingle(y, c, th, alfaref, cl, cd, AoA, relfactor)
    Solve for a single angle of attack.

calcgamma(yi, ci, thi, AoA, alfa, cl, ...)
    Iterative circulation solver with relaxation.

checkmonotonic(x)
    Check if a vector is strictly monotonically increasing.
"""

import warnings

import numpy as np


# ═════════════════════════════════════════════════════════════════════════════
#  Main function
# ═════════════════════════════════════════════════════════════════════════════

def liftingline(geomfile, forcefile, AoA=0, relfactor=0.01):
    """
    [CL, CD, y, Gamma, v] = liftingline(geomfile, forcefile, AoA)

    Parameters
    ----------
    geomfile : str
        Geometry file with columns: spanwise position (normalised by span),
        chord length (normalised by span), twist angle (deg).
    forcefile : str
        Airfoil data file with columns: angle of attack (deg), lift
        coefficient, drag coefficient.  Additional columns are ignored.
    AoA : float or array_like, optional
        Wing angles of attack (deg).  Scalar or 1-D array.  Default 0.

    Returns
    -------
    CL : ndarray
        Lift coefficient for each angle in AoA.
    CD : ndarray
        Induced drag coefficient for each angle in AoA.
    y : ndarray
        Spanwise stations normalised by span.
    Gamma : ndarray
        Normalised circulation  (rows = spanwise stations,
        columns = angles of attack).
    v : ndarray
        Normalised downwash, v/U  (same shape as Gamma).
    """
    # Ensure AoA is a 1-D array
    AoA = np.atleast_1d(np.asarray(AoA, dtype=float)).ravel()

    # Read force data
    data = np.loadtxt(forcefile)
    alfaref = data[:, 0]
    cl = data[:, 1]
    cd = data[:, 2]

    # Read blade geometry data
    data = np.loadtxt(geomfile, skiprows=1)          # header row
    y = data[:, 0]
    c = data[:, 1]
    th = data[:, 2]

    # Validate geometry input
    N = len(y)
    if N < 2:
        raise ValueError(
            'Geometry file must contain at least 2 spanwise stations.'
        )
    if np.any(np.diff(y) <= 0):
        raise ValueError(
            'Spanwise coordinates must be strictly increasing.'
        )

    # Normalise chord and span coordinates by total span
    L = y[-1] - y[0]
    y = y / L
    c = c / L

    # Check that the lift curve is monotonically increasing
    _, monoInds = checkmonotonic(cl)
    if len(monoInds) < len(cl):
        print('Lift coefficient data is not monotonically increasing.')
        print('  Data will be truncated to its longest monotonic section.')

    alfaref = alfaref[monoInds]
    cl = cl[monoInds]
    cd = cd[monoInds]

    # Initialise solution arrays
    n_y = len(y)
    n_aoa = len(AoA)
    Gamma = np.zeros((n_y, n_aoa))
    v = np.zeros((n_y, n_aoa))
    CL = np.zeros(n_aoa)
    CD = np.zeros(n_aoa)

    # Solve for each angle of attack
    for ii, aoa in enumerate(AoA):
        CL[ii], CD[ii], Gamma[:, ii], v[:, ii] = (
            LLsingle(y, c, th, alfaref, cl, cd, aoa, relfactor)
        )

    return CL, CD, y, Gamma, v


# ═════════════════════════════════════════════════════════════════════════════
#  Single-angle solver
# ═════════════════════════════════════════════════════════════════════════════

def LLsingle(y, c, th, alfaref, cl, cd, AoA, relfactor):
    """
    Solve lifting-line for a single angle of attack.

    Parameters
    ----------
    y : ndarray       Normalised spanwise stations.
    c : ndarray       Normalised chord distribution.
    th : ndarray       Twist angles (deg).
    alfaref : ndarray  Reference angles of attack from airfoil data (deg).
    cl : ndarray       Lift coefficient data.
    cd : ndarray       Drag coefficient data.
    AoA : float        Wing angle of attack (deg).

    Returns
    -------
    CL : float
    CD : float
    Gamma : ndarray    Normalised circulation.
    v : ndarray        Normalised downwash (v/U).
    """
    # Iterative circulation solver
    Gamma, v, alfai, _ = calcgamma(
        y, c, th, AoA, alfaref, cl, relfactor=relfactor
    )

    # Interpolate force coefficients at converged effective AoA
    cli = np.interp(alfai, alfaref, cl)
    cdi = np.interp(alfai, alfaref, cd)

    # Force coefficients go to zero at blade ends
    cli[0] = 0.0
    cdi[0] = 0.0
    cli[-1] = 0.0
    cdi[-1] = 0.0

    # Planform area (dimensionless)
    S = np.trapz(c, y)

    # Total force coefficients (eqs. 17-18 in LLT.tex)
    CL = (1.0 / S) * np.trapz(
        c * (cli * np.cos(v) - cdi * np.sin(v)), y
    )
    CD = (1.0 / S) * np.trapz(
        c * (cli * np.sin(v) + cdi * np.cos(v)), y
    )

    return CL, CD, Gamma, v


# ═════════════════════════════════════════════════════════════════════════════
#  Iterative circulation solver
# ═════════════════════════════════════════════════════════════════════════════

def calcgamma(yi, ci, thi, AoA, alfa, cl,
              Gamma0=0, maxiter=1000, errtol=1e-6, relfactor=0.01):
    """
    Iteratively solve for the circulation distribution.

    Parameters
    ----------
    yi : ndarray       Normalised spanwise stations.
    ci : ndarray       Normalised chord.
    thi : ndarray       Twist angles (deg).
    AoA : float        Wing angle of attack (deg).
    alfa : ndarray     Reference angles of attack (deg).
    cl : ndarray       Lift coefficient data.

    Keyword Args
    ------------
    Gamma0 : float     Scaling factor for initial parabolic guess.
    maxiter : int      Maximum iterations.
    errtol : float     Convergence tolerance on Gamma.
    relfactor : float  Relaxation factor (0 < r ≤ 1).

    Returns
    -------
    Gamma : ndarray    Converged normalised circulation.
    vi : ndarray       Normalised downwash (v/U).
    alfai : ndarray    Converged effective angle of attack (deg).
    err : ndarray      Error history.
    """
    # Convert angles to radians
    alfa_rad = np.deg2rad(alfa)
    AoA_rad = np.deg2rad(AoA)
    thi_rad = np.deg2rad(thi)

    # Trailing vortex positions at midpoints between yi stations
    yj = 0.5 * (yi + np.roll(yi, -1))[:-1]

    # Initial guess: parabolic distribution, zero at ends
    Gamma = -Gamma0 * (yi - yi[0]) * (yi - yi[-1])
    Gamma[0] = 0.0
    Gamma[-1] = 0.0

    # Iterate
    err = 10.0
    history = [err]
    n_iter = 0

    while n_iter < maxiter and err > errtol:
        n_iter += 1

        # Trailing vortex strengths (Helmholtz)
        dGamma = -np.diff(Gamma)

        # Influence matrix: Y[i, j] = 1 / (4 pi (yj[j] - yi[i]))
        Yij = 1.0 / (4.0 * np.pi * (yj[np.newaxis, :] - yi[:, np.newaxis]))

        # Induced downwash
        vi = Yij @ dGamma

        # Effective angle of attack (radians)
        alfai_rad = AoA_rad + thi_rad - vi

        # Interpolate lift coefficient (extrapolate during iteration)
        cli = np.interp(np.rad2deg(alfai_rad), alfa, cl)
        cli[0] = 0.0
        cli[-1] = 0.0

        # Circulation from lift coefficient  (eq. 14 in LLT.tex)
        Gamma_new = cli * ci / 2.0
        Gamma_new[0] = 0.0
        Gamma_new[-1] = 0.0

        # Error and relaxation update
        eGamma = Gamma_new - Gamma
        err = np.max(np.abs(eGamma))
        history.append(err)

        Gamma = Gamma + relfactor * eGamma

    # Convergence report
    print(f'Total iterations: {n_iter}')
    print(f'err = {err:.6e}')
    if err > errtol:
        print(f'NOT CONVERGED WITHIN e = {errtol}')
    else:
        print(f'CONVERGED WITHIN e = {errtol}')

    # Effective angle of attack at convergence (radians → degrees)
    alfai = np.rad2deg(alfai_rad)

    # Check if converged effective AoA is within supplied airfoil data range
    interior = alfai[1:-1]
    if np.any(interior > alfa[-1]) or np.any(interior < alfa[0]):
        print('ERROR: Effective angle of attack out of range.')
        Gamma[:] = np.nan
        vi[:] = np.nan
        return Gamma, vi, alfai, history

    return Gamma, vi, alfai, history


# ═════════════════════════════════════════════════════════════════════════════
#  Utility
# ═════════════════════════════════════════════════════════════════════════════

def checkmonotonic(x):
    """
    Check whether a vector is strictly monotonically increasing.

    Parameters
    ----------
    x : array_like

    Returns
    -------
    flag : bool
        True if the entire vector is strictly increasing.
    indices : ndarray
        Indices of the longest contiguous strictly increasing section
        (1-D, matches orientation of `x`).
    """
    x = np.asarray(x)
    flag = bool(np.all(np.diff(x) > 0))

    if flag:
        indices = np.arange(len(x))
    else:
        breaks = np.where(np.diff(x) <= 0)[0]
        starts = np.concatenate([[0], breaks + 1])
        ends = np.concatenate([breaks, [len(x) - 1]])
        longest = np.argmax(ends - starts + 1)
        indices = np.arange(starts[longest], ends[longest] + 1)

    # Match orientation of input
    if isinstance(x, np.matrix) or (hasattr(x, 'shape') and x.ndim == 2
                                     and x.shape[0] > x.shape[1]):
        indices = indices[:, np.newaxis]

    return flag, indices