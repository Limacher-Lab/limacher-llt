"""
run_validation.py — Run full LLT validation and save results to CSV.

Does three things:
1. Convergence study: lift slope error vs number of spanwise stations
   → Output: validation_python_convergence.csv
2. Lift curve and drag polar at 301 stations
   → Output: validation_python_data.csv
3. Prints a summary to console

Usage:
    python _validation/run_validation.py
"""
import io, sys, os, tempfile
import numpy as np

script_dir = os.path.dirname(os.path.abspath(__file__))
repo_dir = os.path.dirname(script_dir)
sys.path.insert(0, repo_dir)
os.chdir(repo_dir)
from liftingline import liftingline

# ── Parameters ──────────────────────────────────────────────────────────
AR = 8.0
ALFAS = np.arange(-10, 11, 2)
RELAXATION = 0.01
STATIONS = [31, 151, 301, 601, 1001]      # for convergence study
N_DETAIL = 301                             # for detailed lift/drag data

a0 = 2 * np.pi
theory_slope = a0 / (1 + a0 / (np.pi * AR))
force_file = os.path.join(repo_dir, 'thin_airfoil_data.txt')


# ═════════════════════════════════════════════════════════════════════════
#  1. Convergence study
# ═════════════════════════════════════════════════════════════════════════

print("=" * 60)
print("  LLT Validation — Python Implementation")
print("=" * 60)

print(f"\n  AR = {AR}, AoA = {ALFAS[0]}° to {ALFAS[-1]}°, relfactor = {RELAXATION}")
print(f"\n─── Convergence Study ───────────────────────────────────")
print(f"{'Stations':>10}  {'dCL/dα (rad⁻¹)':>16}  {'Error (%)':>10}  {'Avg iters':>10}")
print("-" * 52)

conv_rows = []
for n in STATIONS:
    y = np.linspace(-0.5, 0.5, n)
    c_root = 4.0 / (np.pi * AR)
    c = c_root * np.sqrt(1.0 - (2.0 * y)**2)
    th = np.zeros_like(y)

    geom_file = os.path.join(tempfile.gettempdir(), f'elliptic_{n}.txt')
    np.savetxt(geom_file,
               np.column_stack([y, c, th]),
               fmt='%.15f', delimiter=',',
               header='y,c,th', comments='')

    # Suppress print output
    old_stdout = sys.stdout
    sys.stdout = io.StringIO()
    CL, CD, _, _, _ = liftingline(
        geom_file, force_file,
        AoA=ALFAS, relfactor=RELAXATION
    )
    captured = sys.stdout.getvalue()
    sys.stdout = old_stdout

    # Count iterations
    total_iters = 0
    n_aoa = 0
    for line in captured.split('\n'):
        if line.startswith('Total iterations:'):
            total_iters += int(line.split(':')[1].strip())
            n_aoa += 1
    avg_iters = total_iters / n_aoa if n_aoa > 0 else 0

    # Compute lift slope
    if np.any(np.isnan(CL)):
        slope = np.nan
        err = np.nan
        print(f"{n:10d}  {'—':>16}  {'—':>10}  {avg_iters:10.1f}  (did not converge)")
    else:
        p = np.polyfit(np.deg2rad(ALFAS), CL, 1)
        slope = p[0]
        err = abs(slope - theory_slope) / theory_slope * 100
        print(f"{n:10d}  {slope:16.6f}  {err:9.4f}%  {avg_iters:10.1f}")

    conv_rows.append((n, slope, err, avg_iters))
    os.remove(geom_file)

# Save convergence CSV
conv_csv = os.path.join(script_dir, 'validation_python_convergence.csv')
np.savetxt(conv_csv, np.array(conv_rows),
           fmt='%d,%.6f,%.4f,%.1f',
           delimiter=',',
           header='Stations,dCL_dalpha,Error_pct,Avg_iters',
           comments='')
print(f"\n  Saved: {conv_csv}")


# ═════════════════════════════════════════════════════════════════════════
#  2. Detailed lift-curve and drag-polar data (at N_DETAIL stations)
# ═════════════════════════════════════════════════════════════════════════

print(f"\n─── Detailed Results ({N_DETAIL} stations) ────────────────")

y_full = np.linspace(-0.5, 0.5, N_DETAIL)
c_root = 4.0 / (np.pi * AR)
c_full = c_root * np.sqrt(1.0 - (2.0 * y_full)**2)
th_full = np.zeros_like(y_full)

geom_file = os.path.join(tempfile.gettempdir(), 'elliptic_detail.txt')
np.savetxt(geom_file,
           np.column_stack([y_full, c_full, th_full]),
           fmt='%.15f', delimiter=',',
           header='y,c,th', comments='')

old_stdout = sys.stdout
sys.stdout = io.StringIO()
CL, CD, y_span, Gamma, v = liftingline(
    geom_file, force_file,
    AoA=ALFAS, relfactor=RELAXATION
)
sys.stdout = old_stdout
os.remove(geom_file)

CL_theory = theory_slope * np.deg2rad(ALFAS)
CD_theory = CL_theory**2 / (np.pi * AR)

print(f"{'AoA(°)':>8} {'CL':>10} {'CD':>10} {'CL_theory':>10} {'CD_theory':>10}")
print("-" * 50)
for i, a in enumerate(ALFAS):
    print(f"{a:8d} {CL[i]:10.6f} {CD[i]:10.6f} {CL_theory[i]:10.6f} {CD_theory[i]:10.6f}")

# Save detailed data CSV
data_csv = os.path.join(script_dir, 'validation_python_data.csv')
np.savetxt(data_csv,
           np.column_stack([ALFAS, CL, CD, CL_theory, CD_theory]),
           fmt='%d,%.10f,%.10f,%.10f,%.10f',
           delimiter=',',
           header='Alpha,CL,CD,CL_theory,CD_theory',
           comments='')
print(f"\n  Saved: {data_csv}")

print("\n" + "=" * 60)
print("  Validation complete.")
print("=" * 60)