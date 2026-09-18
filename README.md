Power from Space: Coordinated Satellite Charging for Off-Grid Wireless Systems
==================

This code package is related to the following scientific article:

O. M. Rosabal, A. Azarbahram, M. Ashraf, M. Shehab, A. B. Khattak, O. L. A. López, and M.-S. Alouini, “Power from space: Coordinated satellite charging for off-grid wireless systems,” arXiv preprint arXiv:2608.25589, 2026. (Accepted for publication in IEEE Internet of Things Magazine)

Available at: [https://arxiv.org/pdf/2608.25589](https://arxiv.org/pdf/2608.25589)

## Abstract of Article

Satellite-enabled wireless power transfer (WPT) may be a transformative solution for charging Internet of Things (IoT) devices in off-grid scenarios where traditional technologies struggle to efficiently meet urgent energy demands. In this article, we review the advantages and limitations of microwave-based long-distance charging for satellite-enabled WPT. We then introduce our vision of coordinated space-based WPT, where multiple satellites jointly serve networks of ground devices. Potential use cases are presented highlighting application requirements. We evaluate the average received power at the target locations using two coordination schemes and perform a statistical characterization of the power spillover on undesired locations. We also shed light on the performance of inter-satellite laser WPT for different operating distances and transmit-receive apertures of the peer satellites. Moreover, we explore the integration of metasurfaces on satellite apertures and ground networks to boost energy conversion efficiency, scalability, and beam management.  Finally, we outline relevant challenges and research directions towards implementing our vision. 

## Content of Code Package

The repository contains Python and MATLAB scripts and user-defined functions required to reproduce the numerical results in the article. To run the code for the RF-WPT results, we recommend using high-performance computing resources, as the simulations may take hours or even days to complete on a personal computer. Also, we advise the use of MOSEK (via CVXPY) as the optimization solver. The code was ran using cvxpy 1.7.2 and Mosek 11.0.29.

See each file for further documentation.

## Acknowledgements

This work is partially supported by Research Council of Finland (Grants 348515 (UPRISING) and 369116 (6G Flagship)). The authors also wish to acknowledge CSC - IT Center for Science, Finland, for computational resources.
