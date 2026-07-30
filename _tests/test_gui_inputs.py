"""
Tests for GUI input validation and angle-array generation.

These functions are extracted from the GUI module so they can be tested
without opening a tkinter window.
"""
import sys, os, tempfile, io
import numpy as np

# Path setup
script_dir = os.path.dirname(os.path.abspath(__file__))
repo_dir = os.path.dirname(script_dir)
sys.path.insert(0, repo_dir)

from llt_gui import (
    parse_geometry_file,
    parse_airfoil_file,
    build_aoa_array,
    validate_angle_sweep,
)


# ── Helpers ─────────────────────────────────────────────────────────────

def _make_geom_file(y, c, th):
    """Write a temporary geometry CSV and return the path."""
    f = tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False)
    f.write('y,c,th\n')
    for row in zip(y, c, th):
        f.write(f'{row[0]},{row[1]},{row[2]}\n')
    f.close()
    return f.name


def _make_airfoil_file(alpha, cl, cd):
    """Write a temporary airfoil-data CSV and return the path."""
    f = tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False)
    f.write('Alpha,Cl,Cd\n')
    for row in zip(alpha, cl, cd):
        f.write(f'{row[0]},{row[1]},{row[2]}\n')
    f.close()
    return f.name


def _remove(path):
    """Remove a temp file, ignoring errors."""
    try:
        os.unlink(path)
    except (OSError, PermissionError):
        pass


# ── Geometry file parsing ──────────────────────────────────────────────

class TestParseGeometryFile:
    """Tests for parse_geometry_file()."""

    def test_valid_elliptic_wing(self):
        """Parse standard elliptic geometry."""
        y = np.linspace(-0.5, 0.5, 11)
        c = 0.2 * np.sqrt(1.0 - (2.0 * y)**2)
        th = np.zeros_like(y)
        path = _make_geom_file(y, c, th)
        try:
            yp, cp, thp = parse_geometry_file(path)
            assert len(yp) == 11
            assert np.allclose(yp, y)
            assert np.allclose(cp, c)
            assert np.allclose(thp, th)
        finally:
            _remove(path)

    def test_fewer_than_three_columns(self):
        """Only two columns should raise."""
        path = tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False)
        path.write('y,c\n0.0,0.1\n0.5,0.1\n')
        path.close()
        try:
            parse_geometry_file(path.name)
            assert False, 'Expected ValueError'
        except ValueError as e:
            assert 'columns' in str(e).lower()
        finally:
            _remove(path.name)

    def test_non_increasing_y(self):
        """Non-increasing y should raise."""
        y = [0.0, 0.3, 0.2, 0.5]
        c = [0.1, 0.1, 0.1, 0.0]
        th = [0.0, 0.0, 0.0, 0.0]
        path = _make_geom_file(y, c, th)
        try:
            parse_geometry_file(path)
            assert False, 'Expected ValueError'
        except ValueError as e:
            assert 'increasing' in str(e).lower()
        finally:
            _remove(path)

    def test_fewer_than_two_rows(self):
        """Single data row should raise."""
        path = tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False)
        path.write('y,c,th\n0.0,0.1,0.0\n')
        path.close()
        try:
            parse_geometry_file(path.name)
            assert False, 'Expected ValueError'
        except ValueError:
            pass
        finally:
            _remove(path.name)

    def test_non_numeric_data(self):
        """Non-numeric data should raise."""
        path = tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False)
        path.write('y,c,th\n0.0,abc,0.0\n')
        path.close()
        try:
            parse_geometry_file(path.name)
            assert False, 'Expected ValueError'
        except (ValueError, OSError):
            pass
        finally:
            _remove(path.name)

    def test_missing_file(self):
        """Non-existent file should raise."""
        try:
            parse_geometry_file('/nonexistent/path.csv')
            assert False, 'Expected FileNotFoundError or OSError'
        except (FileNotFoundError, OSError):
            pass


# ── Airfoil file parsing ───────────────────────────────────────────────

class TestParseAirfoilFile:
    """Tests for parse_airfoil_file()."""

    def test_valid_thin_airfoil(self):
        """Parse standard thin-airfoil data."""
        alpha = np.arange(-10.0, 11.0, 2.0)
        cl = 0.109662 * alpha  # ~2π per rad
        cd = np.zeros_like(alpha)
        path = _make_airfoil_file(alpha, cl, cd)
        try:
            ap, clp, cdp = parse_airfoil_file(path)
            assert len(ap) == 11
            assert np.allclose(ap, alpha)
            assert np.allclose(clp, cl)
            assert np.allclose(cdp, cd)
        finally:
            _remove(path)

    def test_non_increasing_alpha(self):
        """Non-increasing Alpha should raise."""
        alpha = [0.0, 5.0, 3.0, 10.0]
        cl = [0.0, 0.5, 0.3, 1.0]
        cd = [0.0, 0.01, 0.01, 0.02]
        path = _make_airfoil_file(alpha, cl, cd)
        try:
            parse_airfoil_file(path)
            assert False, 'Expected ValueError'
        except ValueError as e:
            assert 'increasing' in str(e).lower()
        finally:
            _remove(path)

    def test_fewer_than_two_rows(self):
        """Single data row should raise."""
        path = tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False)
        path.write('Alpha,Cl,Cd\n0.0,0.0,0.0\n')
        path.close()
        try:
            parse_airfoil_file(path.name)
            assert False, 'Expected ValueError'
        except ValueError:
            pass
        finally:
            _remove(path.name)

    def test_fewer_than_three_columns(self):
        """Only two columns should raise."""
        path = tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False)
        path.write('Alpha,Cl\n0.0,0.0\n5.0,0.5\n')
        path.close()
        try:
            parse_airfoil_file(path.name)
            assert False, 'Expected ValueError'
        except ValueError as e:
            assert 'columns' in str(e).lower()
        finally:
            _remove(path.name)


# ── Angle-of-attack array generation ───────────────────────────────────

class TestBuildAoAArray:
    """Tests for build_aoa_array()."""

    def test_integer_sweep(self):
        """-10 to 10 in steps of 2 -> 11 points."""
        arr = build_aoa_array(-10.0, 10.0, 2.0)
        expected = np.arange(-10.0, 11.0, 2.0)
        assert len(arr) == 11
        assert np.allclose(arr, expected)

    def test_fractional_sweep(self):
        """0 to 1 in steps of 0.2 -> 6 points."""
        arr = build_aoa_array(0.0, 1.0, 0.2)
        expected = np.arange(0.0, 1.001, 0.2)
        assert len(arr) == 6
        assert np.allclose(arr, expected)

    def test_single_step_covers_range(self):
        """0 to 5 in steps of 5 -> [0, 5]."""
        arr = build_aoa_array(0.0, 5.0, 5.0)
        assert len(arr) == 2
        assert np.allclose(arr, [0.0, 5.0])

    def test_min_equals_max(self):
        """Min == max -> single element."""
        arr = build_aoa_array(3.0, 3.0, 1.0)
        assert len(arr) == 1
        assert arr[0] == 3.0


# ── Angle sweep validation ─────────────────────────────────────────────

class TestValidateAngleSweep:
    """Tests for validate_angle_sweep()."""

    def test_valid_sweep(self):
        """Standard sweep returns None."""
        err = validate_angle_sweep(-10.0, 10.0, 1.0)
        assert err is None

    def test_increment_zero(self):
        err = validate_angle_sweep(-10.0, 10.0, 0.0)
        assert err is not None
        assert 'greater than zero' in err.lower()

    def test_increment_negative(self):
        err = validate_angle_sweep(-10.0, 10.0, -1.0)
        assert err is not None
        assert 'greater than zero' in err.lower()

    def test_max_below_min(self):
        err = validate_angle_sweep(10.0, -10.0, 1.0)
        assert err is not None
        assert ('maximum' in err.lower() and 'minimum' in err.lower()) or \
               ('greater than' in err.lower())

    def test_exceeds_maximum_points(self):
        err = validate_angle_sweep(-10.0, 10.0, 0.001)
        assert err is not None
        assert 'limit' in err.lower() or '10000' in err or 'exceed' in err.lower()

    def test_nan_input(self):
        """NaN values should be rejected."""
        err = validate_angle_sweep(float('nan'), 10.0, 1.0)
        assert err is not None


# ── Run if called directly ─────────────────────────────────────────────

if __name__ == '__main__':
    import sys
    import pytest
    sys.exit(pytest.main([__file__, '-v']))