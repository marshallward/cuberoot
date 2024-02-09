import numpy as np
import matplotlib.pyplot as plt

lw = 0.1
#lw = 1.
plot_fma = False

data_quad = np.loadtxt('err_quad.txt')
data_rel = np.loadtxt('err_root_rel.txt')


fig, axes = plt.subplots(8,2, figsize=(14,12), sharex=True)
fig.suptitle('Errors of each cubic solver.\nQuad error <--> Relative error')

for row in axes:
    #row[0].set_ylim(0, 1.5e-16)
    #row[0].set_ylim(0, 1.1e-16)
    #row[0].set_ylim(0, 4.4e-16)
    #row[0].set_ylim(0., 8.8e-16)
    #row[1].set_ylim(0, 6e-16)
    #row[0].set_ylim(0, 2e-15)

    row[0].set_ylim(0, 8e-17)
    row[1].set_ylim(0, 5.5e-16)


    for ax in row:
        ax.set_xlim(0.125, 1.0)
        #ax.set_xlim(0.125, 0.1275)

axes[0,0].set_title("x**(1./3.)")
axes[0,0].plot(data_quad[:,0], data_quad[:,2], linewidth=lw)
axes[0,1].plot(data_rel[:,0], data_rel[:,2], linewidth=lw)

axes[1,0].set_title("Newton x6")
axes[1,0].plot(data_quad[:,0], data_quad[:,3], linewidth=lw)
axes[1,1].plot(data_rel[:,0], data_rel[:,3], linewidth=lw)

axes[2,0].set_title("No-div Newton")
axes[2,0].plot(data_quad[:,0], data_quad[:,4], linewidth=lw)
axes[2,1].plot(data_rel[:,0], data_rel[:,4], linewidth=lw)

axes[3,0].set_title("Halley x3 + Newton")
axes[3,0].plot(data_quad[:,0], data_quad[:,5], linewidth=lw)
axes[3,1].plot(data_rel[:,0], data_rel[:,5], linewidth=lw)

axes[4,0].set_title("No-div Halley")
axes[4,0].plot(data_quad[:,0], data_quad[:,6], linewidth=lw)
axes[4,1].plot(data_rel[:,0], data_rel[:,6], linewidth=lw)

axes[5,0].set_title("Leroy (no bit correction)")
axes[5,0].plot(data_quad[:,0], data_quad[:,7], linewidth=lw)
axes[5,1].plot(data_rel[:,0], data_rel[:,7], linewidth=lw)

axes[6,0].set_title("cbrt_ac")
axes[6,0].plot(data_quad[:,0], data_quad[:,8], linewidth=lw)
axes[6,1].plot(data_rel[:,0], data_rel[:,8], linewidth=lw)

axes[7,0].set_title("cbrt_ac in Fortran")
axes[7,0].plot(data_quad[:,0], data_quad[:,9], linewidth=lw)
axes[7,1].plot(data_rel[:,0], data_rel[:,9], linewidth=lw)


#for row in axes:
#    for ax in row:
#        for k in [-51, -52, -53, -54]:
#            ax.axhline(2**k, linestyle="--", color='k', linewidth=0.5)

plt.subplots_adjust(hspace=0.6)

#plt.show()
plt.savefig('err.svg', bbox_inches='tight')
