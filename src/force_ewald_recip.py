# -*- coding: utf-8 -*-
"""
Created on Sat Jan  4 19:28:35 2025

@author: bkell
"""
import cupy as cp
import numba as nb
import numpy as np
import time
import os

# Set OpenBLAS to use 4 threads (replace with the desired number of threads)
os.environ["OPENBLAS_NUM_THREADS"] = "99"

def ewald_reciprocal_force(positions, charges, box_length, k_vectors, alpha):
    """
    Calculate the reciprocal space contribution to the Ewald summation forces.

    Parameters:
    positions (np.ndarray): Nx3 array of particle positions.
    charges (np.ndarray): N array of particle charges.
    box_length (float): Length of the cubic simulation box.
    k_vectors (np.ndarray): Mx3 array of reciprocal lattice vectors.
    alpha (float): Ewald parameter.

    Returns:
    np.ndarray: Nx3 array of forces on each particle.
    """
    volume = box_length ** 3
    forces = np.zeros_like(positions)

    k_norms = np.linalg.norm(k_vectors, axis=1)
    non_zero_k = k_norms != 0
    k_vectors = k_vectors[non_zero_k]
    k_norms = k_norms[non_zero_k]

    k_dot_r = np.dot(positions, k_vectors.T)
    cos_k_dot_r = cp.cos(k_dot_r)
    sin_k_dot_r = cp.sin(k_dot_r)

    structure_factor_real = cp.einsum('i,ij->j', charges, cos_k_dot_r)
    structure_factor_imag = cp.einsum('i,ij->j', charges, sin_k_dot_r)

    k_squared = k_norms ** 2
    prefactor = (4 * np.pi / volume) * np.exp(-k_squared / (4 * alpha**2)) / k_squared

    # Fix the broadcasting issue in force calculation
    force_contributions = charges[:, np.newaxis] * (
        structure_factor_real[np.newaxis, :] * sin_k_dot_r - 
        structure_factor_imag[np.newaxis, :] * cos_k_dot_r
    )
    
    # Reshape and compute the forces correctly
    forces = np.sum(
        prefactor[np.newaxis, :, np.newaxis] * 
        force_contributions[:, :, np.newaxis] * 
        k_vectors[np.newaxis, :, :],
        axis=1
    )

    return forces

# Example usage
import os
os.environ["OPENBLAS_NUM_THREADS"] = "2"
# Generate 2000 random positions and charges
np.random.seed(42)  # For reproducibility
num_particles = 1000000

# Generate random positions within the box
positions = np.random.uniform(0, 10.0, size=(num_particles, 3))

# Generate random charges that sum to zero
# First generate random charges
charges = np.random.uniform(-2, 2, size=num_particles)
# Adjust the last charge to make the sum zero
#charges[-1] = -np.sum(charges[:-1])
print(f"charge of last charge is {charges[-1]}")
box_length = 10.0
k_vectors = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
alpha = 0.5

# Start the timer
start_time = time.time()

reciprocal_forces = ewald_reciprocal_force(positions, charges, box_length, k_vectors, alpha)

# Calculate and print the elapsed time
end_time = time.time()
elapsed_time = end_time - start_time
print(f"\nExecution time: {elapsed_time:.4f} seconds")

print(f"Sum of all charges: {np.sum(charges)}")
print(f"Number of particles: {len(positions)}")
print(f"Shape of reciprocal forces: {reciprocal_forces.shape}")
print(f"Reciprocal forces from Ewald summation:\n{reciprocal_forces}")
