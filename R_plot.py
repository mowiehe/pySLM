import numpy as np
import matplotlib.pyplot as plt
import zernike as zk

### Fig. 3 in [1]

ms = [0, 1, 2]
ns = [[2, 4, 6, 8], [1, 3, 5, 7], [2, 4, 6, 8]]

r = np.linspace(0, 1, 100)
fig, ax = plt.subplots(3, figsize=(8, 10), subplot_kw={"xlabel": "$r$"})
for i, m in enumerate(ms):
    for n in ns[i]:

        R = np.array([zk.R_mn(m, n, ir) for ir in r])
        ax[i].plot(r, R, label=f"n = {n}")
        ax[i].set_ylabel(f"$R^{m}_n(r)$")
        ax[i].set_ylim(-1, 1)

    ax[i].legend(loc="upper left")
plt.tight_layout()
plt.show()
