import matplotlib.pyplot as plt
import numpy as np
import zernike as zk

plt.style.use("science")

### Define coefficients
n_max = 10
C_mn = np.zeros((n_max + 1, n_max + 1))
# piston, bias
C_mn[0, 0] = 0
# tip,tilt
C_mn[-1, 1] = 1
C_mn[1, 1] = 0
# astigmatism, defocus
C_mn[-2, 2] = 0
C_mn[0, 2] = 0
C_mn[2, 2] = 0


### creating polar plot
r_array = np.linspace(0, 1, 10)
theta_array = np.arange(0, 2 * np.pi, np.pi / 180)

r, theta = np.meshgrid(r_array, theta_array)
W_res = [
    [zk.W(r[i, j], theta[i, j], C_mn) for j in range(len(r_array))]
    for i in range(len(theta_array))
]
# W_vec = np.vectorize(W)
# W_res = W_vec(r, theta)


# -- Plot... ------------------------------------------------
fig, ax = plt.subplots(subplot_kw={"projection": "polar"})
ax.contourf(theta, r, W_res)

plt.show()