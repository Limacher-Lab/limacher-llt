"""
Create elliptic blade geometry file for a specified aspect ratio and
number of spanwise stations.

Usage:
    create_elliptic_wing(filename, N, AR)

The file is written as a CSV with columns: y, c, th
where y and chord are normalised by span, and twist is in degrees.
"""
import numpy as np


def create_elliptic_wing(filename, N, AR):
    """
    Create an untwisted elliptic planform geometry file.

    Parameters
    ----------
    filename : str   Path to the output CSV file.
    N : int          Number of spanwise stations.
    AR : float       Aspect ratio (b²/S).
    """
    y = np.linspace(-0.5, 0.5, N)
    a = 4.0 / (np.pi * AR)
    c = 2.0 * a * (0.25 - y**2)**0.5
    th = np.zeros_like(y)

    data = np.column_stack([y, c, th])
    np.savetxt(filename, data, fmt='%.15f', delimiter=',',
               header='y,c,th', comments='')