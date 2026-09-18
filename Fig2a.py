# This script reproduces the results of Fig. 2a in the manuscript 

# O. M. Rosabal, A. Azarbahram, M. Ashraf, M. Shehab, A. B. Khattak, O. L. A.
# López, and M.-S. Alouini, “Power from space: Coordinated satellite charging
# for off-grid wireless systems,” arXiv preprint arXiv:2608.25589, 2026.
# (Accepted for publication in IEEE Internet of Things Magazine)

# Version: 1.0 
# Last modified: 2026-09-18

# License: This code is licensed under MIT License. See the LICENSE file in the repository root.

# If you use this code in research resulting in a publication, please cite
# the paper above.

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Ellipse
import numpy as np
import matplotlib.pyplot as plt
from src import core

# =========================================
# simulation parameters
# =========================================
M = np.arange(2, 11)                    # num. satellites
N = np.array([9, 16])                   # num. antennas per satellite (URA)
K = 5                                   # num. devices
MCIter = 1000                           # Monte Carlo iterations

R = 6378e3                              # radius of Earth [m]
H = 500e3                               # height of the satellites
txPow = 1_000                           # per-satellite power budget
fc = 12e9                               # operating frequency
wavelength = 3e8/fc                     # wavelength
delta = wavelength/2                    # inter-antenna spacing

devSphericalCapRadius = 2e3;
devPolarAngleMax = np.arcsin(devSphericalCapRadius/R)
satPolarAngleMax = np.arcsin(10e3/(R+H))

phaseErrorStdDeg = np.array([10, 20])
phaseErrorStdRad = np.deg2rad(phaseErrorStdDeg)
nPerturbations = 100

# =========================================
# optimization loop
# =========================================
rPowNonCoherent = np.zeros((len(M), len(N)))
rPowCoherent = np.zeros((len(M), len(N)))
rPowCoherentSyncErrs = np.zeros((len(M), len(N), len(phaseErrorStdRad)))

for nIdx, n in enumerate(N):
    rng = np.random.default_rng([seed, 12345])

    for mIdx, m in enumerate(M):
        # satellites deployment
        satPos = core.satDeploy(R, H, m, satPolarAngleMax)

        for seed in np.arange(MCIter):
            # devices deployment
            devPos = core.devsDeploy(R, K, devPolarAngleMax, seed)

            h = core.channelModel(n, K, m, satPos, devPos, wavelength, delta)

            # non-coherent coordination transmission 
            precodersNonCoherent = core.nonCoherentCoordination(n, K, m, h, txPow)
            rPowNonCoherent[mIdx, nIdx] += 1/MCIter*np.mean(core.rPowFnc(n, K, m, precodersNonCoherent, h, mode="non-coherent"))

            # coherent coordination transmission
            precodersCoherent = core.coherentCoordination(n, K, m, h, txPow)
            rPowCoherent[mIdx, nIdx] += 1/MCIter*np.mean(core.rPowFnc(n, K, m, precodersCoherent, h, mode="coherent"))

            for errIdx, phaseErrorStdRadVal in enumerate(phaseErrorStdRad):

                for _ in range(nPerturbations):
                    satellitePhases = rng.normal(loc=0.0, scale=phaseErrorStdRadVal, size=m)
                
                    phaseVector = np.repeat(np.exp(1j*satellitePhases), n)
                                        
                    precodersCoherentSyncErrs = (phaseVector[:, np.newaxis]
                                                *precodersCoherent
                                                *phaseVector.conj()[np.newaxis, :])
                
                    # distributed beamforming
                    rPowCoherentSyncErrs[mIdx, nIdx, errIdx] += (1/(MCIter*nPerturbations)*
                                                                 np.mean(core.rPowFnc(n, K, m, precodersCoherentSyncErrs, h, mode="coherent")))

# =========================================
# plot results
# =========================================
axisFontSize = 15
tickFontSize = 12
legendFontSize = 14

centimeters = 1 / 2.54
fig, ax = plt.subplots(figsize=(19 * centimeters, 12 * centimeters))

ax.plot(M, rPowNonCoherent[:, 0]*1e3, "-o", color="#A2142F", markersize=8, linewidth=1.5)
ax.plot(M, rPowNonCoherent[:, 1]*1e3, "--o", color="#A2142F", markersize=8, linewidth=1.5)
ax.plot(M, rPowCoherent[:, 0]*1e3, "-s", color="#77AC30", markersize=8, linewidth=1.5)
ax.plot(M, rPowCoherent[:, 1]*1e3, "--s", color="#77AC30", markersize=8, linewidth=1.5,)
ax.plot(M, rPowCoherentSyncErrs[:, 0, 0]*1e3, "-^", color="#0072BD", markersize=8, linewidth=1.5,)
ax.plot(M, rPowCoherentSyncErrs[:, 1, 0]*1e3, "--^", color="#0072BD", markersize=8, linewidth=1.5,)
ax.plot(M, rPowCoherentSyncErrs[:, 0, 1]*1e3, "-+", color="#EDB120", markersize=8, linewidth=1.5,)
ax.plot(M, rPowCoherentSyncErrs[:, 1, 1]*1e3, "--+", color="#EDB120", markersize=8, linewidth=1.5,)

antennaLegend = ax.legend(
    handles=[
        Line2D([], [], color="black", linestyle="-", linewidth=1.5),
        Line2D([], [], color="black", linestyle="--", linewidth=1.5),
    ],
    labels=[r"$9$ antennas", r"$16$ antennas"],
    fontsize=legendFontSize,
    loc="upper left",
)
ax.add_artist(antennaLegend)

sync_legend = fig.legend(
    handles=[
        Line2D([], [], color="#77AC30", linestyle="-", linewidth=1.5),
        Line2D([], [], color="#0072BD", linestyle="-", linewidth=1.5),
        Line2D([], [], color="#EDB120", linestyle="-", linewidth=1.5),
    ],
    labels=[r"$\sigma = 0$", r"$\sigma = 10^\circ$", r"$\sigma = 20^\circ$"],
    title="synchronization error",
    fontsize=legendFontSize,
    title_fontsize=legendFontSize,
    loc="upper left",
    bbox_to_anchor=(0.0, 0.81),
    bbox_transform=ax.transAxes,
)

fig.text(0.35, 0.42,
    "coherent coordination",
    ha="center",
    va="center",
    fontsize=14,
)
ax.add_patch(
    Ellipse(
        (0.5, 0.31),
        0.025,
        0.24,
        transform=ax.transAxes,
        fill=False,
        linewidth=1,
    )
)

ax.text(
    0.70,
    0.17,
    "non-coherent coordination",
    transform=ax.transAxes,
    ha="center",
    va="center",
    fontsize=14,
)
ax.add_patch(
    Ellipse(
        (0.81, 0.12),
        0.02,
        0.08,
        transform=ax.transAxes,
        fill=False,
        linewidth=1,
    )
)

ax.grid(True)
ax.set_title("(a)", fontsize=axisFontSize)
ax.set_xlabel("number of satellites", fontsize=axisFontSize)
ax.set_ylabel("average received RF power (mW)", fontsize=axisFontSize)
ax.tick_params(axis="both", labelsize=tickFontSize)

plt.show()
