"""
test_elliptic_wing.py — Verify LLT solver against classical elliptic wing theory

Tests the liftingline() solver on an untwisted elliptic planform using
thin-airfoil input data (cl = 2π·α, cd = 0).

Classical lifting-line theory predicts:
    dC_L / dα = 2π / (1 + 2/AR)                          (lift slope)
    C_Di = C_L² / (π·AR)                                   (induced drag)

Usage:
    python test_elliptic_wing.py
"""

import sys
import os

import numpy as np

# Allow running from the repo directory
script_dir = os.path.dirname(os.path.abspath(__file__))
if script_dir:
    sys.path.insert(0, script_dir)
    os.chdir(script_dir)

from liftingline import liftingline


# ═════════════════════════════════════════════════════════════════════════════
#  Parameters
# ═════════════════════════════════════════════════════════════════════════════

AR = 8.0                           # aspect ratio
N_STATIONS = 31                    # spanwise discretisation
ALFAS = np.arange(-10, 11, 2)      # angles of attack (deg)
RELAXATION = 0.01                  # relaxation factor

# Thin-airfoil 2-D lift slope
a0 = 2 * np.pi

# Theoretical 3-D lift slope for an elliptic wing (Prandtl)
THEORY_CL_ALPHA_RAD = a0 / (1 + a0 / (np.pi * AR))
THEORY_CL_ALPHA_DEG = THEORY_CL_ALPHA_RAD * np.pi / 180

# Tolerance for passing the test
TOLERANCE_PCT = 2.0                # lift slope error must be under 2%


# ═════════════════════════════════════════════════════════════════════════════
#  Run solver
# ═════════════════════════════════════════════════════════════════════════════

GEOM_FILE = "elliptic_AR8.txt"
FORCE_FILE = "thin_airfoil_data.txt"

print("=" * 60)
print("  Lifting-Line Theory — Elliptic Wing Test")
print("=" * 60)
print(f"\n  Geometry:      elliptic planform, AR = {AR}")
print(f"  Stations:      {N_STATIONS}")
print(f"  Force data:    thin airfoil (cl = 2π·α, cd = 0)")
print(f"  AoA range:     {ALFAS[0]}° to {ALFAS[-1]}°")
print(f"  Relaxation:    {RELAXATION}")
print()

CL, CD, y, Gamma, v = liftingline(
    GEOM_FILE, FORCE_FILE, AoA=ALFAS, relfactor=RELAXATION
)


# ═════════════════════════════════════════════════════════════════════════════
#  Lift slope check
# ═════════════════════════════════════════════════════════════════════════════

# Linear regression: CL vs α (radians)
p = np.polyfit(np.deg2rad(ALFAS), CL, 1)
computed_slope_rad = p[0]
computed_intercept = p[1]

error_pct = abs(computed_slope_rad - THEORY_CL_ALPHA_RAD) \
            / THEORY_CL_ALPHA_RAD * 100

print("─── Lift Slope ──────────────────────────────────────")
print(f"  Computed:      {computed_slope_rad:.4f} rad⁻¹  ({computed_slope_rad * np.pi/180:.6f} deg⁻¹)")
print(f"  Theoretical:   {THEORY_CL_ALPHA_RAD:.4f} rad⁻¹  ({THEORY_CL_ALPHA_DEG:.6f} deg⁻¹)")
print(f"  Error:         {error_pct:.4f}%")
print(f"  Intercept:     {computed_intercept:.6f}")
print(f"  Pass/Fail:     {'PASS' if error_pct < TOLERANCE_PCT else 'FAIL'}")
print()

# Detailed table
print("─── Results Table ───────────────────────────────────")
print(f"{'AoA(°)':>8} {'CL':>10} {'CD':>10} {'CL²/(π·AR)':>12}")
print("-" * 42)
for a, cl_val, cd_val in zip(ALFAS, CL, CD):
    cdp = cl_val**2 / (np.pi * AR)
    print(f"{a:8d} {cl_val:10.6f} {cd_val:10.6f} {cdp:12.6f}")


# ═════════════════════════════════════════════════════════════════════════════
#  Induced drag check
# ═════════════════════════════════════════════════════════════════════════════

print()
print("─── Induced Drag ────────────────────────────────────")
cd_errors = np.abs(CD - CL**2 / (np.pi * AR))
max_cd_error = np.max(cd_errors)
print(f"  Max |CD − CL²/(π·AR)| = {max_cd_error:.6e}")

# Summary
print()
print("─" * 60)
if error_pct < TOLERANCE_PCT:
    print(f"  ✓ PASSED: Lift slope within {TOLERANCE_PCT}% of classical theory.")
else:
    print(f"  ✗ FAILED: Lift slope error {error_pct:.2f}% exceeds {TOLERANCE_PCT}%.")
print("=" * 60)