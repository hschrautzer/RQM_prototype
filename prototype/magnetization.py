r"""
Module contains magnetization class of this prototype
"""
import numpy as np
import scipy.spatial as spt


class Magnetization:
    r"""
    Represents the magnetization of a lattice instance. This corresponds to CSpinLattice in Spinaker.
    """

    def __init__(self, points: np.ndarray, spins: np.ndarray) -> None:
        r"""
        Initializes the lattice. The point and spin arrays come from spinaker or from mumax. The neighbor relations
        are calculated using a kd-Tree. The data structure should be as follows:

        :param points: np.ndarray of shape: [[p1x,p1y,p1z],[p2x,p2y,p2z],[...]]
        :param spins: np.ndarray of shape: [[s1x,s1y,s1z],[s2x,s2y,s2z],[...]]
        """
        self._exchange_constant = 1  # in units per atom
        self._lattice_constant = 1.0
        self._points: np.ndarray = points
        self._spins: np.ndarray = spins
        self._kdtree = spt.KDTree(self._points)

    def energy(self) -> float:
        r"""
        This just contains nearest neighbor exchange

        :return: the energy of the configuration
        """
        energy = 0
        for i, s in enumerate(self._spins):
            _, indices = self._kdtree.query(self._points[i], k=4)
            for neigh_idx in indices:
                energy -= np.dot(s, self._spins[neigh_idx]) * self._exchange_constant
        return energy

    def gradient(self) -> np.ndarray:
        r"""
        This just contains nearest neighbor exchange

        :return: the gradient of the configuration
        """
        gradient = np.zeros_like(self._spins)
        for i, p in enumerate(self._points):
            _, indices = self._kdtree.query(p, k=4)
            for neigh_idx in indices:
                gradient[i, :] -= self._exchange_constant * self._spins[neigh_idx]
        return gradient

    def gradient_tspace_3N(self) -> np.ndarray:
        r"""

        :return:
        """
        grad = self.gradient()
        return np.asarray([grad[i] - np.dot(grad[i], s) * s for i, s in enumerate(self._spins)])
        # return grad - np.linalg.vd(grad,self._spins)*grad

    def gradient_tspace_2N(self) -> np.ndarray:
        r"""

        :return:
        """
        basis = self.basis()
        return np.asarray(
            [[np.dot(basis[:, 0, i], g), np.dot(basis[:, 1, i], g)] for i, g in enumerate(self.gradient())])

    def basis(self) -> np.ndarray:
        r"""
        Calculates the tangent space basis of the magnetic configuration
        :return: np.ndarray of shape (3,2,N)
        """
        basis = np.zeros(shape=(3, 2, len(self._spins)))
        for idx, spin in enumerate(self._spins):
            if np.abs(spin[2]) >= 0.5:
                # normalize v=(1,0,0)
                xi = np.array([1.0 - spin[0] * spin[0], -spin[0] * spin[1], -spin[0] * spin[2]])
            else:
                xi = np.array([-spin[2] * spin[0], -spin[2] * spin[1], 1 - spin[2] * spin[2]])
            xi = xi / np.linalg.norm(xi)
            eta = np.cross(xi, spin)
            eta = eta / np.linalg.norm(eta)
            basis[:, 0, idx] = xi
            basis[:, 1, idx] = eta
        return basis

    def project_to_basis(self, vec_embedding_space) -> np.ndarray:
        r"""
        Takes a 3N vector from embedding space (shape (N,3)) and projects it to tangent space (shape (N,2) by
        calculating basis[:,:,i] * vec[:,i]
        :param vec_embedding_space:
        :return:
        """
        l_basis = self.basis()

    @classmethod
    def retraction(cls, current_mag: 'Magnetization', vec_tspace: np.ndarray,
                   displacement_parameter: float = 1.0, rodriguez_threshold: float = 1.0e-5) -> 'Magnetization':
        r"""
        This is a class method and is given the current magnetization instance. It returns a new instance with the same
        properties but the spins have been rotated along a vector in tangent space of the current magnetization
        (Retraction).

        :param current_mag: Magnetization instance with current orientation of magnetic moments
        :param vec_tspace: Vector in tangent space of current configuration (3N representation)
        :param displacement_parameter: scaling of the rotation
        :param rodriguez_threshold: the threshold below which rotation will be according to taylor expansion of rod.
        :return: New magnetization instance with changed orientation of spins.
        """
        rotated_spins = np.zeros(shape=np.shape(current_mag.spins))
        for (idx, spin) in enumerate(current_mag.spins):
            disp = vec_tspace[idx, :] * displacement_parameter
            angle_i = np.linalg.norm(disp)
            if angle_i == 0.0:
                rotated_spins[idx, :] = spin
                continue
            if angle_i >= rodriguez_threshold:
                rotated_spins[idx, :] = spin * np.cos(angle_i) + disp * np.sin(angle_i) / angle_i
            else:
                rotated_spins[idx, :] = spin * np.cos(angle_i) + disp * (
                        1.0 - 1.0 / 6.0 * angle_i ** 2 + 1.0 / 120.0 * angle_i ** 4)
            rotated_spins[idx, :] = rotated_spins[idx, :] / np.linalg.norm(rotated_spins[idx, :])
        return cls(points=current_mag.points, spins=rotated_spins)

    def parallel_transport(self, vec_to_be_transported: np.ndarray, vec_tspace: np.ndarray,
                           displacement_parameter: float = 1.0) -> np.ndarray:
        pass

    @property
    def points(self) -> np.ndarray:
        r"""

        :return:
        """
        return self._points

    @property
    def spins(self) -> np.ndarray:
        r"""

        :return:
        """
        return self.spins
