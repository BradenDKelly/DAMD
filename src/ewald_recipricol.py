import numpy as np
import numba as nb
import time

@nb.jit(nopython=True)
def ewald_reciprocal_contribution(positions, charges, box_length, k_vectors, alpha):
    """
    Calculate the reciprocal space contribution to the Ewald summation energy and forces.

    Parameters:
    positions (np.ndarray): Nx3 array of particle positions.
    charges (np.ndarray): N array of particle charges.
    box_length (float): Length of the cubic simulation box.
    k_vectors (np.ndarray): Mx3 array of reciprocal lattice vectors.
    alpha (float): Ewald parameter.

    Returns:
    np.ndarray: Nx3 array of forces on each particle.
    float: Reciprocal space contribution to the Ewald summation energy.
    """
    volume = box_length ** 3

    k_norms = np.linalg.norm(k_vectors, axis=1)
    non_zero_k = k_norms != 0
    k_vectors = k_vectors[non_zero_k]
    k_norms = k_norms[non_zero_k]

    k_dot_r = np.dot(positions, k_vectors.T)
    cos_k_dot_r = np.cos(k_dot_r)
    sin_k_dot_r = np.sin(k_dot_r)

    structure_factor_real = np.einsum('i,ij->j', charges, cos_k_dot_r)
    structure_factor_imag = np.einsum('i,ij->j', charges, sin_k_dot_r)

    k_squared = k_norms ** 2
    prefactor = (4 * np.pi / volume) * np.exp(-k_squared / (4 * alpha**2)) / k_squared

    structure_factor = structure_factor_real**2 + structure_factor_imag**2
    energy = np.einsum('j,j->', prefactor, structure_factor)

    force_contributions = charges[:, np.newaxis] * (structure_factor_real * sin_k_dot_r - structure_factor_imag * cos_k_dot_r)
    forces = np.einsum('ij,jk,k->ik', force_contributions, k_vectors, prefactor)

    return forces, energy

# Example usage
np.random.seed(42)
num_particles = 6000

positions = np.random.uniform(0, 10.0, size=(num_particles, 3))
charges = np.random.uniform(-2, 2, size=num_particles)
print(f"charge of last charge is {charges[-1]}")
box_length = 10.0
k_vectors = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
alpha = 0.5
n_steps = 3000000
start_t = time.time()
forces, energy = ewald_reciprocal_contribution(positions, charges, box_length, k_vectors, alpha)
end_t = time.time()
print(f"Reciprocal forces from Ewald summation: {forces}")
print(f"Reciprocal energy from Ewald summation: {energy}")
print(f"Time was {end_t - start_t}")
print(f"Time for {n_steps} steps = {(end_t - start_t)*n_steps}")
