import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('out3')
data3 = np.loadtxt('out3')
#data3 = np.loadtxt('out3_intel')

fig, axes = plt.subplots(6,1, figsize=(6,10))

for ax in axes:
    ax.set_ylim(0, 0.75*1.1e-16)
    #ax.set_ylim(0, 4.4e-16)
    ax.set_xlim(0.125, 1.0)

fig.suptitle("Error (with respect to 128-bit precision Newton)")

axes[0].plot(data3[:,0], data3[:,2], linewidth=0.1)
axes[0].set_title("x**(1./3.)")

axes[1].plot(data3[:,0], data3[:,3], linewidth=0.1)
axes[1].set_title("Newton x6")

axes[2].plot(data3[:,0], data3[:,4], linewidth=0.1)
axes[2].set_title("No-div Newton")

axes[3].plot(data3[:,0], data3[:,5], linewidth=0.1)
axes[3].set_title("Halley x3 + Newton")

axes[4].plot(data3[:,0], data3[:,6], linewidth=0.1)
axes[4].set_title("No-div Halley (v1)")

axes[5].plot(data3[:,0], data3[:,7], linewidth=0.1)
axes[5].set_title("No-div Halley (v2)")


for ax in axes:
    for k in [-52, -53, -54]:
        ax.axhline(2**k, linestyle="--", color='k', linewidth=0.5)

plt.subplots_adjust(hspace=0.6)

#plt.show()
plt.savefig('err.svg')
plt.savefig('err.png')
