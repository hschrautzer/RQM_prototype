r"""
Module contains magnetization class of this prototype
"""
import numpy as np
import scipy.spatial as spt


class magnetization:
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
        :return: the gradient of the configuration
        """
        gradient = np.zeros_like(self._spins)
        for i, p in enumerate(self._points):
            _, indices = self._kdtree.query(p, k=4)
            for neigh_idx in indices:
                gradient[i, :] -= self._exchange_constant * self._spins[neigh_idx]
        return gradient

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
