# This module provides plotting utilities used to reproduce the results in the
# manuscript.

# O. M. Rosabal, A. Azarbahram, M. Ashraf, M. Shehab, A. B. Khattak, O. L. A.
# López, and M.-S. Alouini, “Power from space: Coordinated satellite charging
# for off-grid wireless systems,” arXiv preprint arXiv:2608.25589, 2026.
# (Accepted for publication in IEEE Internet of Things Magazine)

# Version: 1.0
# Last modified: 2026-09-18

# License: This code is licensed under MIT License. See the LICENSE file in the repository root.

# If you use this code in research resulting in a publication, please cite
# the paper above.

import numpy as np
import matplotlib.pyplot as plt

def plotDeployment(radius,nodesPos,polarAngleMax):
    """Plots node positions and their spherical-cap deployment area.

    Args:
        radius: sphere radius [m]
        nodesPos: node Cartesian positions (nNodes, 3), [m]
        polarAngleMax: maximum polar angle [rad]

    Returns:
        None: no value returned
    """
    # Create a sphere mesh for visualization
    azimuthAngle = np.linspace(0, 2*np.pi, 60)
    polarAngle = np.linspace(0, polarAngleMax, 30)
    x = radius*np.outer(np.sin(polarAngle),np.cos(azimuthAngle))
    y = radius*np.outer(np.sin(polarAngle),np.sin(azimuthAngle))
    z = radius*np.outer(np.cos(polarAngle),np.ones_like(azimuthAngle))

    # Plot
    fig, ax = plt.subplots(subplot_kw={"projection": "3d"}, figsize=[10,10])

    # Plot sphere surface (semi-transparent)
    ax.plot_surface(x, y, z, color='lightblue', alpha=0.5, linewidth=0)

    ax.scatter(nodesPos[:,0], nodesPos[:,1], nodesPos[:,2],
            color='red', s=50, edgecolors='k')
