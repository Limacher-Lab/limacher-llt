"""
Tests for aerodynamic performance analysis functions.
"""
import sys, os
import numpy as np

# Path setup
script_dir = os.path.dirname(os.path.abspath(__file__))
repo_dir = os.path.dirname(script_dir)
sys.path.insert(0, repo_dir)

from llt_gui import (
    PerformanceMetrics,
    calculate_aspect_ratio,
    estimate_zero_lift_angle,
    estimate_stall_angle,
    estimate_lift_slope,
    calculate_lift_to_drag,
    find_max_ld,
    calculate_induced_drag,
    calculate_span_efficiency,
    analyze_performance,
)


# ── Test data ──────────────────────────────────────────────────────────

def _linear_CL(AoA, slope=0.1, zero_lift=-2.0):
    """CL = slope * (AoA - zero_lift)"""
    return slope * (np.asarray(AoA) - zero_lift)


def _make_elliptic_AR8_geom():
    """Return y_raw, c_raw for an elliptic AR=8 wing (31 stations)."""
    y = np.linspace(-0.5, 0.5, 31)
    c_root = 4.0 / (np.pi * 8)
    c = c_root * np.sqrt(1.0 - (2.0 * y)**2)
    return y, c


# ── Aspect ratio ───────────────────────────────────────────────────────

class TestCalculateAspectRatio:
    def test_elliptic_ar8(self):
        """AR should be ~8 for elliptic wing (discretisation error ~0.05)."""
        y, c = _make_elliptic_AR8_geom()
        AR = calculate_aspect_ratio(y, c)
        assert abs(AR - 8.0) < 0.1

    def test_rectangular_ar6(self):
        """Rectangular wing: AR = b^2/S = b/c (since S = b*c)."""
        y = np.array([-0.5, 0.5])
        c = np.array([1.0/6.0, 1.0/6.0])  # span=1, so b^2=1, AR = 1/S = 1/(c*1) = 6
        AR = calculate_aspect_ratio(y, c)
        assert abs(AR - 6.0) < 0.01


# ── Zero-lift angle ────────────────────────────────────────────────────

class TestEstimateZeroLiftAngle:
    def test_crossing_at_sampled_point(self):
        """CL exactly zero at a sampled AoA."""
        AoA = np.array([-5.0, -2.0, 0.0, 2.0, 5.0])
        CL = _linear_CL(AoA, slope=0.1, zero_lift=-2.0)  # CL=0 at AoA=-2
        result = estimate_zero_lift_angle(CL, AoA)
        assert result is not None
        assert abs(result - (-2.0)) < 1e-10

    def test_crossing_between_points(self):
        """Zero lift between two sampled points — linear interpolation."""
        AoA = np.array([-5.0, -3.0, 0.0, 3.0])
        CL = _linear_CL(AoA, slope=0.1, zero_lift=-2.0)
        result = estimate_zero_lift_angle(CL, AoA)
        assert result is not None
        assert abs(result - (-2.0)) < 0.1

    def test_no_crossing_all_negative(self):
        """All CL values negative — no zero crossing."""
        AoA = np.array([-10.0, -8.0, -6.0])
        CL = np.array([-1.0, -0.8, -0.6])
        result = estimate_zero_lift_angle(CL, AoA)
        assert result is None

    def test_no_crossing_all_positive(self):
        """All CL values positive — no zero crossing."""
        AoA = np.array([0.0, 2.0, 4.0])
        CL = np.array([0.1, 0.3, 0.5])
        result = estimate_zero_lift_angle(CL, AoA)
        assert result is None

    def test_picks_nearest_zero_crossing_to_zero(self):
        """Multiple zero crossings — pick nearest to 0 deg."""
        AoA = np.array([-10.0, -5.0, -2.0, 0.0, 3.0, 8.0])
        CL = np.array([-1.0, -0.5, 0.2, -0.1, 0.3, 0.8])
        result = estimate_zero_lift_angle(CL, AoA)
        assert result is not None
        # Crossing near AoA=-2 (neg to pos) and near AoA=1 (pos to neg).
        # Should pick the one nearest zero degrees.
        assert abs(result) < 3.0

    def test_nan_cl_values(self):
        """NaN values in CL should be handled gracefully."""
        AoA = np.array([-5.0, -2.0, 0.0, 2.0])
        CL = np.array([-0.5, -0.2, np.nan, 0.2])
        result = estimate_zero_lift_angle(CL, AoA)
        assert result is not None  # Should find crossing via interpolation


# ── Stall angle ────────────────────────────────────────────────────────

class TestEstimateStallAngle:
    def test_clear_stall_detected(self):
        """Interior CL maximum followed by decline -> stall."""
        AoA = np.array([0, 2, 4, 6, 8, 10, 12, 14, 16])
        CL = np.array([0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.1, 1.05, 0.9])
        result = estimate_stall_angle(CL, AoA, alpha_zero=0.0)
        assert result is not None
        assert result == 12.0  # AoA where CL peaks

    def test_no_stall_at_sweep_end(self):
        """Max CL at upper boundary — not stall."""
        AoA = np.array([0, 2, 4, 6, 8, 10])
        CL = np.array([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
        result = estimate_stall_angle(CL, AoA, alpha_zero=0.0)
        assert result is None

    def test_stall_plateau(self):
        """Plateau at peak — earliest AoA in plateau."""
        AoA = np.array([0, 2, 4, 6, 8, 10, 12, 14])
        CL = np.array([0.0, 0.2, 0.5, 0.7, 1.0, 1.0, 0.95, 0.8])
        result = estimate_stall_angle(CL, AoA, alpha_zero=0.0)
        assert result is not None
        assert result == 8.0  # Earliest in peak plateau

    def test_insufficient_points(self):
        """Too few finite points."""
        AoA = np.array([0, 2])
        CL = np.array([0.0, 0.1])
        result = estimate_stall_angle(CL, AoA, alpha_zero=0.0)
        assert result is None

    def test_no_positive_lift_branch(self):
        """All CL negative."""
        AoA = np.array([-10, -8, -6, -4])
        CL = np.array([-0.8, -0.6, -0.4, -0.2])
        result = estimate_stall_angle(CL, AoA, alpha_zero=-10.0)
        assert result is None


# ── Lift slope ─────────────────────────────────────────────────────────

class TestEstimateLiftSlope:
    def test_stall_based_secant(self):
        """Secant slope between zero-lift and stall."""
        AoA = np.array([0, 4, 8, 12, 14])
        CL = np.array([0.0, 0.4, 0.8, 1.2, 1.1])
        # CL at stall (AoA=12): 1.2
        # CL at zero_lift (AoA=0): 0.0
        # slope = 1.2 / 12 = 0.1 per deg
        slopes = estimate_lift_slope(CL, AoA, alpha_zero=0.0, alpha_stall=12.0)
        assert slopes is not None
        slope_deg, slope_rad = slopes
        assert abs(slope_deg - 0.1) < 1e-10
        assert abs(slope_rad - 0.1 * 180 / np.pi) < 1e-10

    def test_linear_fallback(self):
        """Linear data with no stall — use OLS fallback."""
        AoA = np.array([-5, -3, -1, 1, 3, 5])
        CL = _linear_CL(AoA, slope=0.08, zero_lift=-2.0)
        slopes = estimate_lift_slope(CL, AoA, alpha_zero=-2.0, alpha_stall=None)
        assert slopes is not None
        slope_deg, slope_rad = slopes
        assert abs(slope_deg - 0.08) < 0.005

    def test_poor_linear_fallback_rejected(self):
        """Non-linear data without stall — fallback rejected."""
        AoA = np.array([0, 2, 4, 6, 8, 10])
        CL = np.array([0.0, 0.3, 0.5, 0.7, 0.85, 0.9])
        slopes = estimate_lift_slope(CL, AoA, alpha_zero=0.0, alpha_stall=None)
        assert slopes is None

    def test_no_data(self):
        """No finite data."""
        slopes = estimate_lift_slope(
            np.array([np.nan, np.nan]), np.array([0, 1]),
            alpha_zero=None, alpha_stall=None)
        assert slopes is None


# ── Lift-to-drag ratio ────────────────────────────────────────────────

class TestCalculateLiftToDrag:
    def test_basic_ld(self):
        """Standard calculation."""
        CL = np.array([0.0, 0.5, 1.0])
        CD = np.array([0.01, 0.02, 0.04])
        LD = calculate_lift_to_drag(CL, CD)
        assert np.allclose(LD, [0.0, 25.0, 25.0], equal_nan=True)

    def test_zero_cd_gives_nan(self):
        """Zero CD -> NaN in L/D."""
        CL = np.array([0.5, 1.0])
        CD = np.array([0.0, 0.02])
        LD = calculate_lift_to_drag(CL, CD)
        assert np.isnan(LD[0])
        assert not np.isnan(LD[1])

    def test_negative_cd_gives_nan(self):
        """Negative CD -> NaN."""
        CL = np.array([0.5])
        CD = np.array([-0.01])
        LD = calculate_lift_to_drag(CL, CD)
        assert np.isnan(LD[0])


# ── Maximum L/D ────────────────────────────────────────────────────────

class TestFindMaxLD:
    def test_basic_max(self):
        """Standard case."""
        AoA = np.array([0, 2, 4, 6, 8, 10])
        CL = np.array([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
        CD = np.array([0.02, 0.02, 0.025, 0.035, 0.05, 0.07])
        # L/D = [0, 10, 16, 17.14, 16, 14.28]
        result = find_max_ld(CL, CD, AoA, alpha_zero=0.0)
        assert result is not None
        max_ld, angle = result
        assert abs(max_ld - 17.143) < 0.01
        assert abs(angle - 6.0) < 0.01

    def test_zero_profile_drag(self):
        """All Cd zero -> report no meaningful max."""
        AoA = np.array([0, 2, 4, 6])
        CL = np.array([0.0, 0.2, 0.4, 0.6])
        CD = np.array([0.0, 0.0, 0.0, 0.0])
        result = find_max_ld(CL, CD, AoA, alpha_zero=0.0)
        assert result is None  # No meaningful max

    def test_partial_zero_cd(self):
        """Some but not all Cd zero — may still have meaningful max."""
        AoA = np.array([0, 2, 4, 6])
        CL = np.array([0.0, 0.2, 0.4, 0.6])
        CD = np.array([0.01, 0.01, 0.0, 0.02])
        result = find_max_ld(CL, CD, AoA, alpha_zero=0.0)
        # Cd[2] is zero, so L/D[2] will be NaN
        # Should still find max among finite values
        assert result is not None
        max_ld, angle = result
        assert not np.isnan(max_ld)

    def test_max_at_sweep_end(self):
        """Max L/D at end — log boundary warning in analyze_performance."""
        AoA = np.array([0, 2, 4])
        CL = np.array([0.0, 0.2, 0.4])
        CD = np.array([0.02, 0.02, 0.02])
        result = find_max_ld(CL, CD, AoA, alpha_zero=0.0)
        assert result is not None
        max_ld, angle = result
        assert angle == 4.0  # At end


# ── Induced drag ───────────────────────────────────────────────────────

class TestCalculateInducedDrag:
    def test_elliptic_wing_symmetry(self):
        """Elliptic wing: CDi should be ~CL^2/(pi*AR)."""
        y = np.linspace(-0.5, 0.5, 31)
        span = y[-1] - y[0]
        y_norm = y / span
        AR = 8.0
        c_root = 4.0 / (np.pi * AR)
        c_norm = c_root * np.sqrt(1.0 - (2.0 * y_norm)**2)
        S_norm = np.trapezoid(c_norm, y_norm)

        CL = 0.5
        # Ideal elliptic: Gamma ~ sqrt(1-(2y)^2), v ~ constant
        Gamma = CL * c_norm / 2.0  # from solver relation
        v = CL / (np.pi * AR) * np.ones_like(y_norm)  # induced AoA in rad

        CDi = calculate_induced_drag(Gamma, v, y_norm, S_norm)
        CDi_theory = CL**2 / (np.pi * AR)
        assert abs(CDi - CDi_theory) < 0.05 * CDi_theory

    def test_finite_check(self):
        """NaN inputs should produce NaN output."""
        y_norm = np.array([-0.5, 0.5])
        S_norm = 0.1
        Gamma = np.array([np.nan, 0.0])
        v = np.array([0.1, 0.1])
        CDi = calculate_induced_drag(Gamma, v, y_norm, S_norm)
        assert np.isnan(CDi)


# ── Span efficiency ────────────────────────────────────────────────────

class TestCalculateSpanEfficiency:
    def test_elliptic_theory(self):
        """Elliptic wing should give e ~ 1."""
        e = calculate_span_efficiency(CL=0.5, CDi=0.5**2/(np.pi*8), AR=8.0)
        assert abs(e - 1.0) < 1e-10

    def test_nan_inputs(self):
        """NaN inputs -> None."""
        e = calculate_span_efficiency(CL=np.nan, CDi=0.01, AR=8.0)
        assert e is None

    def test_zero_CL(self):
        """CL near zero -> None (division not meaningful)."""
        e = calculate_span_efficiency(CL=0.0, CDi=0.01, AR=8.0)
        assert e is None

    def test_negative_CDi(self):
        """Negative CDi -> None."""
        e = calculate_span_efficiency(CL=0.5, CDi=-0.01, AR=8.0)
        assert e is None


# ── Integrated analyze_performance ─────────────────────────────────────

class TestAnalyzePerformance:
    def test_with_thin_airfoil_data(self):
        """Full pipeline with linear (inviscid) data."""
        AoA = np.arange(-10, 11, 2)
        CL = _linear_CL(AoA, slope=0.1, zero_lift=-2.0)
        CD = CL**2 / (np.pi * 8)  # induced drag only (zero profile drag)
        y_norm = np.linspace(-0.5, 0.5, 31)
        c_root = 4.0 / (np.pi * 8)
        c_norm = c_root * np.sqrt(1.0 - (2.0 * y_norm)**2)
        S_norm = np.trapezoid(c_norm, y_norm)

        Gamma_arr = np.zeros((len(y_norm), len(AoA)))
        for i in range(len(AoA)):
            Gamma_arr[:, i] = CL[i] * c_norm / 2.0
        # Elliptic wing: constant induced downwash v = CL/(pi*AR)
        AR = 8.0
        v_arr = np.zeros_like(Gamma_arr)
        for i in range(len(AoA)):
            v_arr[:, i] = CL[i] / (np.pi * AR)

        cd_section = np.array([0.0, 0.0])  # zero profile drag

        metrics = analyze_performance(
            CL, CD, AoA, y_norm, c_norm, y_norm, c_norm, S_norm,
            Gamma_arr, v_arr, cd_section
        )

        assert metrics.zero_lift_angle is not None
        assert abs(metrics.zero_lift_angle - (-2.0)) < 0.1
        assert metrics.stall_angle is None  # linear, no stall
        assert metrics.lift_slope_deg is not None  # linear fallback
        assert abs(metrics.lift_slope_deg - 0.1) < 0.005
        assert metrics.max_ld is None  # zero profile drag
        assert metrics.span_efficiency is not None
        assert abs(metrics.span_efficiency - 1.0) < 0.1

    def test_all_nan_returns_defaults(self):
        """All NaN input -> all N/A."""
        AoA = np.array([0, 2, 4])
        CL = np.array([np.nan, np.nan, np.nan])
        CD = np.array([np.nan, np.nan, np.nan])
        y_norm = np.array([-0.5, 0.5])
        c_norm = np.array([0.1, 0.1])
        S_norm = 0.1
        Gamma = np.ones((2, 3))
        v = np.ones((2, 3))
        cd_section = np.array([0.1, 0.1])

        metrics = analyze_performance(
            CL, CD, AoA, y_norm, c_norm, y_norm, c_norm, S_norm,
            Gamma, v, cd_section
        )

        assert metrics.zero_lift_angle is None
        assert metrics.stall_angle is None
        assert metrics.lift_slope_deg is None
        assert metrics.max_ld is None
        assert metrics.angle_max_ld is None
        assert metrics.span_efficiency is None


if __name__ == '__main__':
    import pytest
    import sys
    sys.exit(pytest.main([__file__, '-v']))