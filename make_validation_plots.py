import sys, os, io
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

script_dir = r'C:\Users\herme\repos\limacher-llt'
sys.path.insert(0, script_dir)
os.chdir(script_dir)
from liftingline import liftingline

# Set minimum font size to 12pt
plt.rcParams.update({
    'font.size': 12,
    'axes.labelsize': 13,
    'axes.titlesize': 14,
    'legend.fontsize': 12,
    'xtick.labelsize': 12,
    'ytick.labelsize': 12,
})

AR = 8.0
ALFAS = np.arange(-10, 11, 2)

old_out = sys.stdout
sys.stdout = io.StringIO()
CL, CD, y, Gamma, v = liftingline(
    'elliptic_AR8.txt', 'thin_airfoil_data.txt',
    AoA=ALFAS, relfactor=0.01
)
sys.stdout = old_out

a0 = 2 * np.pi
theory_slope = a0 / (1 + a0 / (np.pi * AR))
theory_CL = theory_slope * np.deg2rad(ALFAS)
theory_CDi = theory_CL**2 / (np.pi * AR)

# --- Figure 1: Lift curve ---
fig1, ax1 = plt.subplots(figsize=(7.5, 5))
ax1.plot(ALFAS, CL, 'bo-', ms=6, lw=2, label='LLT (computed)')
ax1.plot(ALFAS, theory_CL, 'r--', lw=2, label=f'Theory: dCL/da = {theory_slope:.4f} rad-1')
ax1.set_xlabel('Angle of Attack (deg)')
ax1.set_ylabel('CL')
ax1.set_title(f'Elliptic Wing (AR={AR}) - Lift Curve')
ax1.legend()
ax1.grid(True, alpha=0.3)
p = np.polyfit(np.deg2rad(ALFAS), CL, 1)
err = abs(p[0] - theory_slope) / theory_slope * 100
ax1.text(0.97, 0.05, f'Error: {err:.3f}%', transform=ax1.transAxes,
         ha='right', fontsize=12,
         bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
fig1.tight_layout()
fig1.savefig(os.path.join(script_dir, 'validation_lift_curve.png'), dpi=150)

# --- Figure 2: Drag polar ---
fig2, ax2 = plt.subplots(figsize=(7.5, 5))
ax2.plot(CL, CD, 'bo-', ms=6, lw=2, label='LLT (computed)')
ax2.plot(theory_CL, theory_CDi, 'r--', lw=2, label=r'Theory: CDi = CL^2 / (pi AR)')
ax2.set_xlabel('CL')
ax2.set_ylabel('CD')
ax2.set_title(f'Elliptic Wing (AR={AR}) - Drag Polar')
ax2.legend()
ax2.grid(True, alpha=0.3)
fig2.tight_layout()
fig2.savefig(os.path.join(script_dir, 'validation_drag_polar.png'), dpi=150)

# --- Figure 3: Convergence ---
stations = [31, 101, 301, 601]
errors = [0.6143, 0.1689, 0.1073, 0.0756]
fig3, ax3 = plt.subplots(figsize=(7.5, 5))
ax3.semilogx(stations, errors, 'go-', ms=8, lw=2)
ax3.axhline(y=0.1, color='r', ls='--', lw=1.5, label='0.1% threshold')
ax3.set_xlabel('Number of Spanwise Stations')
ax3.set_ylabel('Lift Slope Error (%)')
ax3.set_title('Convergence Study - Error vs Discretization')
ax3.legend()
ax3.grid(True, alpha=0.3, which='both')
ax3.set_xticks(stations)
ax3.set_xticklabels([str(s) for s in stations])
ax3.text(0.97, 0.95, '601 stations: 0.076%', transform=ax3.transAxes,
         ha='right', fontsize=12,
         bbox=dict(boxstyle='round', facecolor='lightgreen', alpha=0.8))
fig3.tight_layout()
fig3.savefig(os.path.join(script_dir, 'validation_convergence.png'), dpi=150)

print("Done")