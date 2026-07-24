"""
Convergence study: lift slope error vs number of spanwise stations.
Uses 5 station counts spanning the range 31–1001.
"""
import io, sys, os, tempfile
import numpy as np

script_dir = os.path.dirname(os.path.abspath(__file__))
repo_dir = os.path.dirname(script_dir)
sys.path.insert(0, repo_dir)
os.chdir(repo_dir)
from liftingline import liftingline

AR = 8.0
ALFAS = np.arange(-10, 11, 2)
RELAXATION = 0.01
STATIONS = [31, 151, 301, 601, 1001]

a0 = 2 * np.pi
theory = a0 / (1 + a0 / (np.pi * AR))
force_file = os.path.join(repo_dir, 'thin_airfoil_data.txt')

print(f"{'Stations':>10}  {'dCL/dα (rad⁻¹)':>16}  {'Error (%)':>10}  {'Avg iters':>10}")
print("-" * 50)

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

    # Check for NaN
    if np.any(np.isnan(CL)):
        print(f"{n:10d}  {'—':>16}  {'—':>10}  {avg_iters:10.1f}  (did not converge)")
    else:
        p = np.polyfit(np.deg2rad(ALFAS), CL, 1)
        slope = p[0]
        err = abs(slope - theory) / theory * 100
        print(f"{n:10d}  {slope:16.6f}  {err:9.4f}%  {avg_iters:10.1f}")

    os.remove(geom_file)