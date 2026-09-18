# This script reproduces the results of Fig. 2b in the manuscript 

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
from src import core

# =========================================
# simulation parameters
# =========================================
M = np.array([5, 10])                   # num. satellites
N = np.array([9, 16])                   # num. antennas per satellite (URA)
K = 5                                   # num. devices
Ks = 50**2                              # num. of power sensors
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

# positions of the power probes
probesPos = core.sensorsDeploy(R, Ks, devPolarAngleMax)

# =========================================
# optimization loop
# =========================================
rPowNonCoherentProbes = np.zeros((len(M), len(N), Ks))
rPowCoherentProbes = np.zeros((len(M), len(N), Ks))

for nIdx, n in enumerate(N):

    for mIdx, m in enumerate(M):
        # satellites deployment
        satPos = core.satDeploy(R, H, m, satPolarAngleMax)

        # power sensors channel coefficients
        hProbes = core.channelModel(n, Ks, m, satPos, probesPos, wavelength, delta)

        for seed in np.arange(MCIter):
            # devices deployment
            devPos = core.devsDeploy(R, K, devPolarAngleMax, seed)

            # devices channel coefficients
            h = core.channelModel(n, K, m, satPos, devPos, wavelength, delta)

            # non-coherent coordination transmission 
            precodersNonCoherent = core.nonCoherentCoordination(n, K, m, h, txPow)
            
            # coherent coordination transmission
            precodersCoherent = core.coherentCoordination(n, K, m, h, txPow)

            # spilled power at the probes
            rPowNonCoherentProbes[mIdx, nIdx, :] += (1/MCIter
                                                     *core.rPowFnc(n, Ks, m, precodersNonCoherent, hProbes, mode="non-coherent").ravel())
            rPowCoherentProbes[mIdx, nIdx, :] += (1/MCIter
                                                  *core.rPowFnc(n, Ks, m, precodersCoherent, hProbes, mode="coherent").ravel())

# empirical CCDF
nBins = 100

# memory pre-allocation
centersNonCoherent = np.zeros((len(M), len(N), nBins))
countsNonCoherent = np.zeros((len(M), len(N), nBins))

centersCoherent = np.zeros((len(M), len(N), nBins))
countsCoherent = np.zeros((len(M), len(N), nBins))

for nIdx, _ in enumerate(N):
    for mIdx, _ in enumerate(M):
        # Normalize each spatial distribution by its mean probe power.
        rPowNonCoherentNorm = (
            rPowNonCoherentProbes[mIdx, nIdx, :]
            / np.mean(rPowNonCoherentProbes[mIdx, nIdx, :])
        )
        rPowCoherentNorm = (
            rPowCoherentProbes[mIdx, nIdx, :]
            / np.mean(rPowCoherentProbes[mIdx, nIdx, :])
        )

        histogramNonCoherent, edgesNonCoherent = np.histogram(
            rPowNonCoherentNorm,
            bins=nBins,
            density=False,
        )
        countsNonCoherent[mIdx, nIdx, :] = (
            np.cumsum(histogramNonCoherent) / histogramNonCoherent.sum()
        )
        centersNonCoherent[mIdx, nIdx, :] = (
            edgesNonCoherent[:-1] + edgesNonCoherent[1:]
        ) / 2

        histogramCoherent, edgesCoherent = np.histogram(
            rPowCoherentNorm,
            bins=nBins,
            density=False,
        )
        countsCoherent[mIdx, nIdx, :] = (
            np.cumsum(histogramCoherent) / histogramCoherent.sum()
        )
        centersCoherent[mIdx, nIdx, :] = (
            edgesCoherent[:-1] + edgesCoherent[1:]
        ) / 2

# =========================================
# plot results
# =========================================
axisFontSize = 15
tickFontSize = 12
legendFontSize = 14

centimeters = 1/2.54
fig, ax = plt.subplots(figsize=(19*centimeters, 12*centimeters))

ax.plot(centersNonCoherent[0, 0, :].squeeze(),
        1 - countsNonCoherent[0, 0, :].squeeze(),
        "-",
        color="#A2142F",
        linewidth=1.5,
        )
ax.plot(centersNonCoherent[0, 1, :].squeeze(),
        1 - countsNonCoherent[0, 1, :].squeeze(),
        "--",
        color="#A2142F",
        linewidth=1.5,
        )
ax.plot(centersNonCoherent[1, 0, :].squeeze(),
    1 - countsNonCoherent[1, 0, :].squeeze(),
    "-.",
    color="#A2142F",
    linewidth=1.5,
    )
ax.plot(centersNonCoherent[1, 1, :].squeeze(),
        1 - countsNonCoherent[1, 1, :].squeeze(),
        ":",
        color="#A2142F",
        linewidth=1.5,
        )

ax.plot(centersCoherent[0, 0, :].squeeze(),
        1 - countsCoherent[0, 0, :].squeeze(),
        "-",
        color="#77AC30",
        linewidth=1.5,
        )
ax.plot(centersCoherent[0, 1, :].squeeze(),
        1 - countsCoherent[0, 1, :].squeeze(),
        "--",
        color="#77AC30",
        linewidth=1.5,
        )
ax.plot(centersCoherent[1, 0, :].squeeze(),
        1 - countsCoherent[1, 0, :].squeeze(),
        "-.",
        color="#77AC30",
        linewidth=1.5,
        )
ax.plot(centersCoherent[1, 1, :].squeeze(),
        1 - countsCoherent[1, 1, :].squeeze(),
        ":",
        color="#77AC30",
        linewidth=1.5,
        )

ax.legend(
    handles=[
        Line2D([], [], color="black", linestyle="-", linewidth=1.5),
        Line2D([], [], color="black", linestyle="--", linewidth=1.5),
        Line2D([], [], color="black", linestyle="-.", linewidth=1.5),
        Line2D([], [], color="black", linestyle=":", linewidth=1.5),
    ],
    labels=[
        r"$5$ satellites, $9$ antennas",
        r"$5$ satellites, $16$ antennas",
        r"$10$ satellites, $9$ antennas",
        r"$10$ satellites, $16$ antennas",
    ],
    fontsize=legendFontSize,
    loc="upper right",
)

fig.add_artist(
    Ellipse(
        (0.2417, 0.3510),
        0.2001,
        0.0265,
        transform=fig.transFigure,
        fill=False,
        linewidth=1,
    )
)
fig.text(
    0.2841,
    0.2877,
    "coherent coordination",
    ha="center",
    va="center",
    fontsize=14,
)

fig.add_artist(
    Ellipse(
        (0.6163, 0.5497),
        0.5287,
        0.0265,
        transform=fig.transFigure,
        fill=False,
        linewidth=1,
    )
)
fig.text(
    0.5701,
    0.4732,
    "non-coherent coordination",
    ha="center",
    va="center",
    fontsize=14,
)

ax.grid(True)
ax.set_title("(b)", fontsize=axisFontSize)
ax.set_xlabel("normalized received power", fontsize=axisFontSize)
ax.set_ylabel("CCDF", fontsize=axisFontSize)
ax.tick_params(axis="both", labelsize=tickFontSize)

plt.show()
