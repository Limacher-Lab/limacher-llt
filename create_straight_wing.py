"""
Create straight (rectangular) blade geometry file for a specified aspect
ratio and number of spanwise stations.

Usage:
    create_straight_wing(filename, N, AR)

The file is written as a CSV with columns: y, c, th
where y and chord are normalised by span, and twist is in degrees.
"""
import numpy as np


def create_straight_wing(filename, N, AR):
    """
    Create an untwisted rectangular planform geometry file.

    Parameters
    ----------
    filename : str   Path to the output CSV file.
    N : int          Number of spanwise stations.
    AR : float       Aspect ratio (b²/S).

    Notes
    -----
    For a rectangular wing: S = b * c, so AR = b² / (b*c) = b / c.
    With normalised span b = 1, the chord is c* = 1 / AR.
    """
    y = np.linspace(-0.5, 0.5, N)
    c = (1.0 / AR) * np.ones_like(y)
    th = np.zeros_like(y)

    data = np.column_stack([y, c, th])
    np.savetxt(filename, data, fmt='%.15f', delimiter=',',
               header='y,c,th', comments='')