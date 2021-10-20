import ruamel.yaml as yaml
import numpy as np
import matplotlib.pyplot as plt
import MatplotlibSettings

########################################################
# Loada data
data = np.loadtxt("alpha.dat")

f, (ax1, ax2) = plt.subplots(2, 1, sharex = "all", gridspec_kw = dict(width_ratios = [1], height_ratios = [3, 1]))
plt.subplots_adjust(wspace = 0, hspace = 0)

ax1.text(0.01, 0.0075, r"\textbf{$\alpha(m_e) = 1/137.036$}", fontsize = 18)
ax1.text(0.01, 0.00747, r"\textbf{$\overline{\mbox{MS}}$ NLL evolution}", fontsize = 18)
ax1.set_ylabel(r"$\displaystyle\alpha(\mu)$", fontsize = 18)
ax1.set_xscale("log")
ax1.set_xlim([0.005, 100.])
ax1.set_ylim([0.0073, 0.0076])
ax1.plot(data[:,0], data[:,1], label = r"$n_{f,\rm max} = 1$", color = "red")
ax1.plot(data[:,0], data[:,2], label = r"$n_{f,\rm max} = 3$", color = "blue")
ax1.legend(fontsize = 18)

ax2.set_xlabel(r"\textbf{$\mu$ [GeV]}")
ax2.set_ylabel(r"\textbf{Ratio}", fontsize = 14)
ax2.set_xscale("log")
ax2.set_xlim([0.005, 100.])
ax2.set_ylim([0.99, 1.035])
ax2.plot(data[:,0], data[:,1]/data[:,1], color = "red", ls = "--")
ax2.plot(data[:,0], data[:,2]/data[:,1], color = "blue", ls = "--")

plt.savefig("alpha.pdf")
plt.close()
