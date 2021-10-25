import ruamel.yaml as yaml
import numpy as np
import matplotlib.pyplot as plt
import MatplotlibSettings

########################################################
# Loada data
data = np.loadtxt("VariableNf.dat")

#plt.text(0.01, 0.0075, r"\textbf{$\alpha(m_e) = 1/137.036$}", fontsize = 18)
#plt.text(0.01, 0.00747, r"\textbf{$\overline{\mbox{MS}}$ NLL evolution}", fontsize = 18)
plt.title(r"\textbf{NLL $\overline{\mbox{MS}}$ evolution from $m_e$ to $Q=100$ GeV}")
plt.ylabel(r"\textbf{Ratio new($n_{f,\rm max}$) / ePDF}", fontsize = 18)
plt.xlabel(r"$x$")
#plt.xscale("log")
plt.xlim(0.1, 1)
plt.ylim(0.985, 1.01)
plt.xticks([])
plt.plot(data[:,0], data[:,1], label = r"$\gamma$ ($n_{f,\rm max} = 1)$", color = "red", ls = "--")
plt.plot(data[:,0], data[:,2], label = r"$e^-+e^+$ ($n_{f,\rm max} = 1)$", color = "blue", ls = "--")
plt.plot(data[:,0], data[:,3], label = r"$e^--e^+$ ($n_{f,\rm max} = 1)$", color = "green", ls = "--")

plt.plot(data[:,0], data[:,4], label = r"$\gamma$ ($n_{f,\rm max} = 3)$", color = "red", ls = "-")
plt.plot(data[:,0], data[:,5], label = r"$e^-+e^+$ ($n_{f,\rm max} = 3)$", color = "blue", ls = "-")
plt.plot(data[:,0], data[:,6], label = r"$e^--e^+$ ($n_{f,\rm max} = 3)$", color = "green", ls = "-")

plt.xticks([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1], [r"\textbf{0.1}", r"\textbf{0.2}", r"\textbf{0.3}", r"\textbf{0.4}", r"\textbf{0.5}", r"\textbf{0.6}", r"\textbf{0.7}", r"\textbf{0.8}", r"\textbf{0.9}", r"\textbf{1}"])

plt.legend(fontsize = 18, ncol = 2)

plt.savefig("VariableNf.pdf")
plt.close()
