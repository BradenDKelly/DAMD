from numba import njit, prange

@njit(parallel=True)
def ewald_real_space_contribution(positions, charges, box_length, alpha):
    """
    Calculate the real space (short-range) contribution to the Ewald summation.

    Parameters:
    positions (np.ndarray): Nx3 array of particle positions.
    charges (np.ndarray): N array of particle charges.
    box_length (float): Length of the cubic simulation box.
    alpha (float): Ewald parameter.

    Returns:
    float: Real space contribution to the Ewald summation energy.
    np.ndarray: Nx3 array of forces on each particle.
    """
    N = len(positions)
    forces = np.zeros_like(positions)
    energy = 0.0


    for i in prange(N):
        for j in range(i + 1, N):
            r_ij = positions[i] - positions[j]
            r_ij -= box_length * np.round(r_ij / box_length)  # Apply minimum image convention
            r = np.linalg.norm(r_ij)
            if r > 0:
                q_i = charges[i]
                q_j = charges[j]
                erfc_alpha_r = np.erfc(alpha * r)
                energy += q_i * q_j * erfc_alpha_r / r
                force_factor = q_i * q_j * (erfc_alpha_r + 2 * alpha * r * np.exp(-alpha**2 * r**2) / np.sqrt(np.pi)) / r**2
                forces[i] += force_factor * r_ij
                forces[j] -= force_factor * r_ij

    return energy, forces




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
    forces = np.einsum('j,jk->ik', prefactor, force_contributions[:, :, np.newaxis] * k_vectors[:, np.newaxis, :])

    return forces, energy
