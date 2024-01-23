import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('out2')
data3 = np.loadtxt('out3')

fig, axes = plt.subplots(5,1)

for ax in axes:
    ax.set_ylim(0, 4e-16)
    ax.set_xlim(0.125, 1.0)

fig.suptitle("Error (with respect to 128-bit precision Newton)")

axes[0].plot(data3[:,0], data3[:,2])
axes[0].set_title("x**(1./3.)")

axes[1].plot(data3[:,0], data3[:,3])
axes[1].set_title("Newton Solver (7 iterations)")

axes[2].plot(data3[:,0], data3[:,4])
axes[2].plot(data[:,0], data[:,4])
axes[2].set_title("Divisionless Newton Solver (MOM6 cuberoot())")

axes[3].plot(data3[:,0], data3[:,6])
axes[3].set_title("Halley Solver")

axes[4].plot(data3[:,0], data3[:,5])
axes[4].plot(data[:,0], data[:,5])
axes[4].set_title("Divisionless Halley Solver (PR #....)")

plt.show()
