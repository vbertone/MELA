"""
Compare MELA N-space results with numerical Mellin transforms
of the LHAPDF x-space grid.

Reads `mela_Nscan.dat` (produced by LHAPDFgrid) and computes the
corresponding Mellin moments from the LHAPDF grid for comparison.
"""
import numpy as np
import lhapdf
import mpmath as mp
import matplotlib.pyplot as plt

lhapdf.setVerbosity(0)
PDF = lhapdf.mkPDF("DeltaGluonFF")

data = []
with open("build/run/mela_Nscan.dat") as f:
    for line in f:
        if line.startswith("#"):
            continue
        parts = line.split()
        N = float(parts[0])
        Q = float(parts[1])
        gluon = float(parts[2])
        charm = float(parts[3])
        if not (np.isnan(gluon) or np.isnan(charm)):
            data.append((N, Q, gluon, charm))

data = np.array(data)
Ns = data[:, 0]
Qs = data[:, 1]
mela_gluon = data[:, 2]
mela_charm = data[:, 3]

# Compute LHAPDF Mellin transforms at each (N, Q) point
lhapdf_gluon = np.zeros(len(data))
lhapdf_charm = np.zeros(len(data))

print(f"Computing {len(data)} Mellin integrals...")
for i, (N, Q, mg, mc) in enumerate(data):
    def g_integrand(z, nn=N, qq=Q):
        return z**(nn - 1) * PDF.xfxQ(0, float(z), qq) / float(z)

    def c_integrand(z, nn=N, qq=Q):
        return z**(nn - 1) * PDF.xfxQ(4, float(z), qq) / float(z)

    lhapdf_gluon[i] = float(mp.quad(g_integrand, [PDF.xMin, PDF.xMax]))
    lhapdf_charm[i] = float(mp.quad(c_integrand, [PDF.xMin, PDF.xMax]))
print("Done.")

ratio_gluon = lhapdf_gluon / mela_gluon
ratio_charm = lhapdf_charm / mela_charm

diff_gluon_percent = (lhapdf_gluon - mela_gluon) / mela_gluon * 100
diff_charm_percent = (lhapdf_charm - mela_charm) / mela_charm * 100

print(f"\n{'Q':>10s} {'MELA_g':>14s} {'LHAPDF_g':>14s} {'diff_g %':>10s}"
      f" {'ratio_g':>10s} {'MELA_c':>14s} {'LHAPDF_c':>14s}"
      f" {'diff_c %':>10s} {'ratio_c':>10s}")
print("-" * 120)

for i in range(len(data)):
    print(f"{Qs[i]:10.2f} {mela_gluon[i]:14.7e} {lhapdf_gluon[i]:14.7e}"
          f" {diff_gluon_percent[i]:10.2f} {ratio_gluon[i]:10.4f} {mela_charm[i]:14.7e}"
          f" {lhapdf_charm[i]:14.7e} {diff_charm_percent[i]:10.2f} {ratio_charm[i]:10.4f}")

fig, axes = plt.subplots(2, 2, figsize=(12, 8))

# Top left: gluon absolute values
ax = axes[0, 0]
ax.plot(Qs, mela_gluon, "b-", label="MELA (N-space)")
ax.plot(Qs, lhapdf_gluon, "r--", label="LHAPDF (Mellin)")
ax.set_xlabel("Q [GeV]")
ax.set_ylabel("g(N=6.2, Q)")
ax.set_title("Gluon Mellin moment")
ax.legend()
ax.set_xscale("log")

# Top right: charm absolute values
ax = axes[0, 1]
ax.plot(Qs, mela_charm, "b-", label="MELA (N-space)")
ax.plot(Qs, lhapdf_charm, "r--", label="LHAPDF (Mellin)")
ax.set_xlabel("Q [GeV]")
ax.set_ylabel("c(N=6.2, Q)")
ax.set_title("Charm Mellin moment")
ax.legend()
ax.set_xscale("log")

# Bottom left: gluon ratio
ax = axes[1, 0]
ax.axhline(y=1, color="gray", linestyle=":", alpha=0.5)
ax.axhspan(0.99, 1.01, color="green", alpha=0.1, label="1% band")
ax.plot(Qs, ratio_gluon, "b.-")
ax.set_xlabel("Q [GeV]")
ax.set_ylabel("LHAPDF / MELA")
ax.set_title("Gluon ratio")
ax.set_xscale("log")
ax.legend()

# Bottom right: charm ratio
ax = axes[1, 1]
ax.axhline(y=1, color="gray", linestyle=":", alpha=0.5)
ax.axhspan(0.99, 1.01, color="green", alpha=0.1, label="1% band")
ax.plot(Qs, ratio_charm, "r.-")
ax.set_xlabel("Q [GeV]")
ax.set_ylabel("LHAPDF / MELA")
ax.set_title("Charm ratio")
ax.set_xscale("log")
ax.legend()

plt.suptitle("MELA N-space vs LHAPDF Mellin transform (N=6.2)", fontsize=14)
plt.tight_layout()
plt.savefig("compare_Nscan.pdf", dpi=350)
print("\nPlots saved to `compare_Nscan.pdf`")
