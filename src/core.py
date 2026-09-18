# This module provides the core simulation functions used to reproduce the
# results in the manuscript.

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
import cvxpy as cp
from scipy import constants
from scipy.stats import nakagami, uniform, rayleigh, poisson
from scipy.spatial.distance import cdist

def satDeploy(R,H,M,satPolarAngleMax):
    """Deploys equally spaced satellites on a spherical ring.

    Args:
        R: earth radius [m]
        H: satellite altitude [m]
        M: number of satellites
        satPolarAngleMax: satellite-ring polar angle [rad]

    Returns:
        satPos: satellite Cartesian positions (M, 3), [m]
    """
    satAzimuthAngle = np.linspace(0, 2*np.pi, M, endpoint=False)
    satPos = np.zeros((M, 3))
    for m in range(M):
        xSats = (R+H)*np.sin(satPolarAngleMax)*np.cos(satAzimuthAngle[m])
        ySats = (R+H)*np.sin(satPolarAngleMax)*np.sin(satAzimuthAngle[m])
        zSats = (R+H)*np.cos(satPolarAngleMax)

        satPos[m, :] = np.array([xSats, ySats, zSats])

    return satPos

def devsDeploy(radius, nNodes, polarAngleMax, seed):
    """Randomly deploys nodes uniformly over a spherical cap.

    Args:
        radius: sphere radius [m]
        nNodes: number of nodes
        polarAngleMax: maximum polar angle [rad]
        seed: random-number seed

    Returns:
        nodesPos: node Cartesian positions (nNodes, 3), [m]
    """
    
    rng = np.random.default_rng(seed)

    azimuthAngle = rng.uniform(0, 2*np.pi, size=nNodes)
    u = rng.uniform(np.cos(polarAngleMax), 1, size=nNodes)
    polarAngle = np.arccos(u)

    # Spherical-to-Cartesian coordinates
    x = radius*np.sin(polarAngle)*np.cos(azimuthAngle)
    y = radius*np.sin(polarAngle)*np.sin(azimuthAngle)
    z = radius*np.cos(polarAngle)

    return np.column_stack((x, y, z))

def sensorsDeploy(radius, nSensors, polarAngleMax):
    """Deploys a deterministic grid of sensors over a spherical cap.

    Args:
        radius: sphere radius [m]
        nSensors: requested number of sensors
        polarAngleMax: maximum polar angle [rad]

    Returns:
        sensorsPos: sensor Cartesian positions (sqrt(nSensors)^2, 3), [m]
    """
 
    nSamples = int(np.sqrt(nSensors))
    polarAngle = np.arccos(np.linspace(np.cos(polarAngleMax), 1.0, nSamples))
    azimuthAngle = np.linspace(0, 2*np.pi, nSamples)

    polarAngleMx, azimuthAngleMx = np.meshgrid(polarAngle, azimuthAngle, indexing='ij')

    x = radius*np.sin(polarAngleMx)*np.cos(azimuthAngleMx)
    y = radius*np.sin(polarAngleMx)*np.sin(azimuthAngleMx)
    z = radius*np.cos(polarAngleMx)

    return np.column_stack((x.ravel(), y.ravel(), z.ravel()))

def channelModel(nAnts,nRxs,nSats,satsPos,rxsPos,wavelength,delta):
    """Computes satellite-array channel coefficients for all receivers.

    Args:
        nAnts: number of antennas per satellite
        nRxs: number of receivers
        nSats: number of satellites
        satsPos: satellite Cartesian positions (nSats, 3), [m]
        rxsPos: receiver Cartesian positions (nRxs, 3), [m]
        wavelength: carrier wavelength [m]
        delta: inter-antenna spacing [m]

    Returns:
        h: complex channel coefficients (nRxs, nSats, nAnts)
    """

    # antenna relative positions centered at the Origin
    n = np.arange(nAnts).T

    rowIdx = (n % int(np.sqrt(nAnts))) - (np.sqrt(nAnts)-1)/2
    colIdx = (n // int(np.sqrt(nAnts))) - (np.sqrt(nAnts)-1)/2

    antPos = np.column_stack((np.zeros_like(rowIdx), rowIdx, colIdx))*delta

    # channel model 
    h = np.zeros((nRxs,nSats,nAnts), dtype=np.complex128)
    for m in range(nSats):
        # compute the rotation matrix

        # nadir axis
        xHat = - satsPos[m,:]/np.linalg.norm(satsPos[m,:])

        ez = np.array([0, 0, 1])
        z = ez - (ez @ xHat)*xHat

        ey = np.array([0, 1, 0])
        y = ey - (ey @ xHat)*xHat

        if np.linalg.norm(z) > np.linalg.norm(y):
            zHat = z/np.linalg.norm(z)
            yHat = np.cross(zHat, xHat)
            yHat = yHat/np.linalg.norm(yHat)
        else:
            yHat = y/np.linalg.norm(y)
            zHat = np.cross(xHat, yHat)
            zHat = zHat/np.linalg.norm(zHat)        

        RotMx = np.column_stack((xHat, yHat, zHat))

        for k in range(nRxs):
            distSatDev = np.linalg.norm(rxsPos[k,:] - satsPos[m,:])

            # line-of-sight vector
            LoSVector = (rxsPos[k,:] - satsPos[m,:])/distSatDev
            
            # wave vector
            wVector = 2*np.pi/wavelength*(RotMx @ LoSVector)

            # channel gain
            channelGain = 1/np.sqrt(4*np.pi*distSatDev**2)

            # channel coefficients
            h[k,m,:] = np.sqrt(channelGain)*np.exp(-1j*2*np.pi/wavelength*distSatDev)*np.exp(-1j*2*np.pi/wavelength*antPos@wVector)

    return h

def rPowFnc(nAnts,nRxs,nSats,precoders,channelCoeff,mode):
    """Computes received power for coherent or non-coherent coordination.

    Args:
        nAnts: number of antennas per satellite
        nRxs: number of receivers
        nSats: number of satellites
        precoders: transmit covariance matrix or matrices (nSats*nAnts, nSats*nAnts) or (nSats, nAnts, nAnts), [W]
        channelCoeff: complex channel coefficients (nRxs, nSats, nAnts)
        mode: coordination mode

    Returns:
        rPow: received power at each receiver (nRxs, 1), [W]
    """

    rPow = np.zeros((nRxs,1), dtype=float)
    if mode == "coherent":
        # Stack per user: (K, M, N) -> (K, MN)
        h = channelCoeff.reshape(nRxs, nSats*nAnts)

        # Outer product with conjugate on the second factor: (K, MN, MN)
        H = h[..., :, None] * np.conjugate(h)[..., None, :]

        # received power constraints
        for k in range(nRxs):
            rPow[k] = np.real(np.trace(H[k,:,:] @ precoders))
            
    elif mode == "non-coherent":
        # rank-1 channel matrices
        H = channelCoeff[..., :, None] * np.conjugate(channelCoeff)[..., None, :]

        # received power computation
        for k in range(nRxs):
            for m, W in enumerate(precoders):
                rPow[k] += np.real(np.trace(H[k, m] @ W))   
            
    return rPow

def nonCoherentCoordination(nAnts,nRxs,nSats,channelCoeff,txPow):
    """Optimizes per-satellite covariance matrices for non-coherent coordination.

    Args:
        nAnts: number of antennas per satellite
        nRxs: number of receivers
        nSats: number of satellites
        channelCoeff: complex channel coefficients (nRxs, nSats, nAnts)
        txPow: per-satellite transmit-power limit [W]

    Returns:
        precoders: per-satellite covariance matrices (nSats, nAnts, nAnts), [W]
    """
    # rank-1 channel matrices
    H = channelCoeff[..., :, None] * np.conjugate(channelCoeff)[..., None, :]

    assert channelCoeff.shape == (nRxs, nSats, nAnts)
    assert H.shape == (nRxs, nSats, nAnts, nAnts)

    # variables 
    xi = cp.Variable(nonneg=True)
    W = [cp.Variable((nAnts, nAnts), hermitian=True) for _ in range(nSats)]

    constraints = []

    # transmit power constraints
    for WHatm in W:
        constraints += [WHatm >> 0]
        constraints += [cp.real(cp.trace(WHatm)) <= txPow]

    # received power constraints
    for k in range(nRxs):
        lhs = cp.Constant(0.0)
        for m, WHatm in enumerate(W):
            lhs += cp.real(cp.trace(H[k, m] @ WHatm))   
        constraints += [lhs >= xi]

    prob = cp.Problem(cp.Maximize(xi), constraints)
    prob.solve(solver=cp.MOSEK)

    # Print result.
    print("The optimal value with cooperative beamforming is", prob.value)

    WVals = [Wi.value for Wi in W]
    precoders = np.stack(WVals, axis=0)
    
    return precoders

def coherentCoordination(nAnts,nRxs,nSats,channelCoeff,txPow):
    """Optimizes a joint covariance matrix for coherent coordination.

    Args:
        nAnts: number of antennas per satellite
        nRxs: number of receivers
        nSats: number of satellites
        channelCoeff: complex channel coefficients (nRxs, nSats, nAnts)
        txPow: per-satellite transmit-power limit [W]

    Returns:
        precoder: joint covariance matrix (nSats*nAnts, nSats*nAnts), [W]
    """
    # Stack per user: (K, M, N) -> (K, MN)
    h = channelCoeff.reshape(nRxs, nSats*nAnts)

    # Outer product with conjugate on the second factor: (K, MN, MN)
    H = h[..., :, None] * np.conjugate(h)[..., None, :]

    # variables 
    xi = cp.Variable(nonneg=True)
    W = cp.Variable((nSats*nAnts, nSats*nAnts), hermitian=True)

    constraints = []

    constraints += [W >> 0]

    # transmit power constraints
    for m in range(nSats):
        startIdx = m*nAnts
        stopIdx = (m+1)*nAnts
        constraints += [cp.real(cp.trace(W[startIdx:stopIdx, startIdx:stopIdx])) <= txPow]

    # received power constraints
    for Hk in H:
        constraints += [cp.real(cp.trace(Hk @ W)) >= xi]

    prob = cp.Problem(cp.Maximize(xi), constraints)
    prob.solve(solver=cp.MOSEK)

    # Print result.
    print("The optimal value with distributed beamforming is", prob.value)

    # return the precoder matrix w/o applying solution recovery
    return W.value

def mainSimulation(nAnts,nRxs,nSats,txPow,earthRad,satsHeight,devsPolarAngleMax,satsPolarAngleMax,wavelength,delta,seed):
    """Runs both coordination optimizations for one network deployment.

    Args:
        nAnts: number of antennas per satellite
        nRxs: number of receivers
        nSats: number of satellites
        txPow: per-satellite transmit-power limit [W]
        earthRad: earth radius [m]
        satsHeight: satellite altitude [m]
        devsPolarAngleMax: maximum device polar angle [rad]
        satsPolarAngleMax: satellite-ring polar angle [rad]
        wavelength: carrier wavelength [m]
        delta: inter-antenna spacing [m]
        seed: random-number seed

    Returns:
        WOptCoordinated: per-satellite covariance matrices (nSats, nAnts, nAnts), [W]
        rPowCoordinated: mean received power [W]
        WOptDistributed: joint covariance matrix (nSats*nAnts, nSats*nAnts), [W]
        rPowDistributed: mean received power [W]
    """
    # devices deployment
    devsPos = netwDeploy(earthRad,nRxs,devsPolarAngleMax,seed)

    # satellites deployment
    satsAzimuthAngle = np.linspace(0, 2*np.pi, nSats, endpoint=False)
    satsPos = np.zeros((nSats,3))
    for m in range(nSats):
        xSats = (earthRad+satsHeight)*np.sin(satsPolarAngleMax)*np.cos(satsAzimuthAngle[m])
        ySats = (earthRad+satsHeight)*np.sin(satsPolarAngleMax)*np.sin(satsAzimuthAngle[m])
        zSats = (earthRad+satsHeight)*np.cos(satsPolarAngleMax)

        satsPos[m,:] = np.array([xSats, ySats, zSats])

    h = channelModel(nAnts,nRxs,nSats,satsPos,devsPos,wavelength,delta)

    # cooperative beamforming 
    WOptCoordinated = cooperativeTxs(nAnts,nRxs,nSats,h,txPow)
    rPowCoordinated = np.mean(rPowFnc(nAnts,nRxs,nSats,WOptCoordinated,h,mode="cooperative"))

    # distributed beamforming
    WOptDistributed = distributedTxs(nAnts,nRxs,nSats,h,txPow)
    rPowDistributed = np.mean(rPowFnc(nAnts,nRxs,nSats,WOptDistributed,h,mode="distributed"))

    return WOptCoordinated, rPowCoordinated, WOptDistributed, rPowDistributed
