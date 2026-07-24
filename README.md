# limacher-llt

Prandtl's Lifting-Line Theory (LLT) solver — implemented in both **MATLAB** and **Python**.

This repository provides a numerical implementation of Prandtl's classical lifting-line theory for predicting the aerodynamic forces on finite wings. It was developed for use in **ENME 570/670 Aerodynamics** at the **University of Calgary**.

## Quick Start

### MATLAB

```matlab
[CL, CD] = liftingline('elliptic_AR8.txt', 'thin_airfoil_data.txt', ...
                        -10:2:10, 'relfactor', 0.01);
plot(-10:2:10, CL)
```

### Python

```python
from liftingline import liftingline
CL, CD, y, Gamma, v = liftingline('elliptic_AR8.txt', 'thin_airfoil_data.txt',
                                   AoA=range(-10, 11, 2), relfactor=0.01)
```

## Repository Structure

| Path | Description |
|------|-------------|
| `liftingline.m` / `liftingline.py` | Main LLT solvers (MATLAB / Python) |
| `create_elliptic_wing.m/.py` | Generate elliptic wing geometry files |
| `create_straight_wing.m/.py` | Generate rectangular wing geometry files |
| `test_elliptic_wing.m/.py` | Standalone validation tests |
| `elliptic_AR8.txt` | Elliptic wing geometry (AR=8, 31 stations) |
| `thin_airfoil_data.txt` | Thin-airfoil force data ($c_l = 2\pi\alpha$, $c_d = 0$) |
| `_tests/` | CI test scripts (`test_ci.py`, `test_ci.m`) |
| `_validation/` | Validation studies, reports, and comparison data |
| `_algorithm/` | LaTeX documentation of the algorithm |
| `_common/` | Shared LaTeX support files (class, style, bibliography) |
| `.github/workflows/` | CI workflows (validate + compile LaTeX) |

## Input File Format

### Geometry file (`geomfile`)
CSV with header `y,c,th`:
- `y` — spanwise coordinate, normalised by span $b$ ($-0.5 \le y \le 0.5$)
- `c` — chord length, normalised by span $b$
- `th` — twist angle (degrees)

### Force data file (`forcefile`)
CSV with header `Alpha,Cl,Cd`:
- `Alpha` — angle of attack (degrees)
- `Cl` — section lift coefficient
- `Cd` — section drag coefficient

## Testing

Run the CI validation tests:

```bash
# Python
python _tests/test_ci.py

# MATLAB (or Octave)
run('_tests/test_ci.m')
```

Both implementations run the same 5 quality checks:
1. Solver converges (no NaN)
2. Zero lift at $\alpha = 0$ ($|C_L| < 10^{-12}$)
3. Lift slope error $< 0.2\%$ vs. classical theory
4. Linearity ($R^2 > 0.9999$)
5. Induced drag error $< 5\%$

## Validation

See `_validation/validation_report.pdf` for a full comparison of Python and MATLAB results against classical theory.

## References

Anderson, J. D. Jr. (2011). *Fundamentals of Aerodynamics*, 5th ed. McGraw-Hill.
— Ch. 5 (pp. 449–470) for lifting-line theory,
  Eq. (5.69) for elliptic-wing lift slope,
  Eq. (5.61) for induced drag.

## License

MIT — see [LICENSE](LICENSE).

## Citation

If you use this code in your work, please cite:

> Limacher, E. J. (2026). *Lifting-Line Theory Educational Codebase* (Version 1.0.0). GitHub. https://github.com/Limacher-Lab/limacher-llt

See [`CITATION.cff`](CITATION.cff) for the machine-readable citation metadata.