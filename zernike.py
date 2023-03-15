import numpy as np
import matplotlib.pyplot as plt

"""
References:
[1] Zernike polynomials: a guide, Andre Fleck, 2011
[2] Dynamic aberration correction via spatial light modulator (SLM) for femtosecond direct laser writing: towards spherical voxels, GABRIELIUS KONTENIS, 2020
"""


###
def R_mn(m, n, r):
    """Radial function of Zernike polynomials eq.3 in [1]
    m: Angular frequency
    n: Radial order
    r: unit radius

    """

    assert r >= 0 and r <= 1
    assert type(m) == int and abs(m) <= n and m >= 0
    assert type(n) == int and n >= 0
    assert (n - m) % 2 == 0

    def summand(l, m, n, r):
        """l = 0 ... (n-m)/2"""
        numerator = pow(-1, l) * np.math.factorial(n - l) * pow(r, n - 2 * l)
        denominator = (
            np.math.factorial(l)
            * np.math.factorial(int((n + m) / 2 - l))
            * np.math.factorial(int((n - m) / 2 - l))
        )
        return numerator / denominator

    R = sum([summand(l, m, n, r) for l in range(int((n - m) / 2) + 1)])
    return R


def Z_mn(m, n, r, theta):
    """Zernike polynomial in angular coordinates, eq.2 in [1]
    m: Angular frequency
    n: Radial order
    r: unit radius
    theta: angle in rad, measured clockwise"""
    assert r >= 0 and r <= 1
    assert type(m) == int and abs(m) <= n
    assert type(n) == int and n >= 0
    assert (n - m) % 2 == 0
    assert theta <= 2 * np.pi and theta >= 0

    if m >= 0:
        return R_mn(m, n, r) * np.cos(m * theta)
    else:
        return R_mn(-m, n, r) * np.sin(m * theta)


def W(r, theta, C_mn, n_max=10):
    """Wavefront in polar coordinates, eq.1 in [1]
    r: unit radius
    theta: angle in rad, measured clockwise
    Amplitudes of Zernike polynomials, array of rank 2
    """
    assert r >= 0 and r <= 1
    assert theta <= 2 * np.pi and theta >= 0
    W_sum = 0
    for n in range(n_max + 1):
        for m in range(-n, n + 1, 2):
            if C_mn[m, n] == 0:
                continue
            else:
                W_sum += C_mn[m, n] * Z_mn(m, n, r, theta)
    return W_sum


def GetGrid(x_pix=800, y_pix=800, n_max=10):
    r_max = np.sqrt((x_pix / 2) ** 2 + (y_pix / 2) ** 2)

    grid = np.zeros((x_pix, y_pix))
    for x in range(x_pix):
        for y in range(y_pix):
            x_center = x - x_pix / 2
            y_center = y - y_pix / 2
            r = np.sqrt(x_center**2 + y_center**2) / r_max
            if y_center != 0:
                theta = np.arctan(x_center / y_center) + np.pi
                print(theta)
            else:
                print(x_center, y_center)
                theta = 1
            grid[x, y] = W(r, theta, n_max)
    return grid
