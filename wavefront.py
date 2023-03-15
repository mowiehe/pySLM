import matplotlib.pyplot as plt
import numpy as np
import zernike as zk

plt.style.use("science")

### define coefficients
n_max = 10
C_mn = np.zeros((n_max + 1, n_max + 1))

# piston, bias
C_mn[0, 0] = 0

# tip,tilt
C_mn[-1, 1] = 0
C_mn[1, 1] = 0

# astigmatism, defocus
C_mn[0, 2] = -40

# coma
C_mn[-1, 3] = 10
C_mn[1, 3] = 10

# spherical
C_mn[0, 4] = -20


### creating polar plot
r_array = np.linspace(0, 1, 1000)
theta_array = np.linspace(0, 2 * np.pi, 100)

r, theta = np.meshgrid(r_array, theta_array)
W_res = np.array(
    [
        [zk.W(r[i, j], theta[i, j], C_mn) for j in range(len(r_array))]
        for i in range(len(theta_array))
    ]
)

# -- Plot... ------------------------------------------------
fig, ax = plt.subplots(subplot_kw={"projection": "polar"})
c = ax.contourf(theta, r, W_res, cmap="Greys", levels=10)
fig.colorbar(c)

## projection
degree_set = 135
theta_index = np.argmin(abs(theta_array - degree_set / 180 * np.pi))
degree = round(theta_array[theta_index] / np.pi * 180)

n_pi = round(theta_array[theta_index] / np.pi, 2)
myfilter = theta == theta_array[theta_index]
fig2, ax2 = plt.subplots(
    subplot_kw={
        "ylabel": "W [arb.]",
        "xlabel": "r [arb.]",
        "title": f"$\\theta={degree}^\circ$",
    }
)
ax2.plot(r[myfilter], W_res[myfilter])

plt.show()
