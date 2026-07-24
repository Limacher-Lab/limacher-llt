"""
CI validation test for the lifting-line solver.

Runs at 301 spanwise stations on an elliptic AR=8 wing with
thin-airfoil data and checks against classical theory.

Exits with code 0 if all checks pass, 1 otherwise.
"""
import sys, os, io, tempfile
import numpy as np

# Paths relative to the repo root
script_dir = os.path.dirname(os.path.abspath(__file__))
repo_dir = os.path.dirname(script_dir)
sys.path.insert(0, repo_dir)
os.chdir(repo_dir)

from liftingline import liftingline

# ── Parameters ──────────────────────────────────────────────────────────
AR = 8.0
STATIONS = 301
ALFAS = np.arange(-10, 11, 2)
FORCE_FILE = os.path.join(repo_dir, 'thin_airfoil_data.txt')

# Theoretical values
a0 = 2 * np.pi
theory_slope = a0 / (1 + a0 / (np.pi * AR))

# ── Generate geometry ──────────────────────────────────────────────────
y = np.linspace(-0.5, 0.5, STATIONS)
c_root = 4.0 / (np.pi * AR)
c = c_root * np.sqrt(1.0 - (2.0 * y)**2)
th = np.zeros_like(y)

geom_file = os.path.join(tempfile.gettempdir(), 'ci_elliptic.txt')
np.savetxt(geom_file, np.column_stack([y, c, th]),
           fmt='%.15f', delimiter=',', header='y,c,th', comments='')

# ── Run solver ─────────────────────────────────────────────────────────
old_out = sys.stdout
sys.stdout = io.StringIO()
CL, CD, _, _, _ = liftingline(geom_file, FORCE_FILE, AoA=ALFAS, relfactor=0.01)
sys.stdout = old_out

os.remove(geom_file)

# ── Checks ─────────────────────────────────────────────────────────────
exit_code = 0

# 1. Solver must converge (no NaN)
if np.any(np.isnan(CL)):
    print('FAIL: Solver returned NaN')
    exit_code = 1
else:
    print('PASS: Solver converged (no NaN)')

# 2. Zero lift at α = 0
cl0 = CL[ALFAS == 0][0]
if abs(cl0) > 1e-12:
    print(f'FAIL: CL(α=0) = {cl0:.3e}, expected < 1e-12')
    exit_code = 1
else:
    print(f'PASS: CL(α=0) = {cl0:.3e}')

# 3. Lift slope error
p = np.polyfit(np.deg2rad(ALFAS), CL, 1)
slope = p[0]
slope_err = abs(slope - theory_slope) / theory_slope * 100
if slope_err >= 0.2:
    print(f'FAIL: Lift slope error = {slope_err:.4f}%, expected < 0.2%')
    exit_code = 1
else:
    print(f'PASS: Lift slope error = {slope_err:.4f}%')

# 4. Linearity
CL_fit = np.polyval(p, np.deg2rad(ALFAS))
residuals = CL - CL_fit
ss_res = np.sum(residuals**2)
ss_tot = np.sum((CL - np.mean(CL))**2)
r2 = 1 - ss_res / ss_tot
if r2 <= 0.9999:
    print(f'FAIL: R² = {r2:.6f}, expected > 0.9999')
    exit_code = 1
else:
    print(f'PASS: R² = {r2:.6f}')

# 5. Induced drag
CD_theory = CL**2 / (np.pi * AR)
mask = np.abs(CL) > 1e-10
if np.any(mask):
    cd_errors = np.abs(CD[mask] - CD_theory[mask]) / CD_theory[mask] * 100
    max_cd_err = np.max(cd_errors)
    if max_cd_err >= 5.0:
        print(f'FAIL: Max induced drag error = {max_cd_err:.2f}%, expected < 5%')
        exit_code = 1
    else:
        print(f'PASS: Max induced drag error = {max_cd_err:.2f}%')

# ── Summary ────────────────────────────────────────────────────────────
print()
if exit_code == 0:
    print('All checks passed.')
else:
    print('Some checks failed.')

sys.exit(exit_code)