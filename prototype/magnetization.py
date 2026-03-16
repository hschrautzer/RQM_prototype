"""
Module contains magnetization class of this prototype.

Open boundaries only.
- Internally stores points/spins as flattened 3N arrays.
- Nearest-neighbor exchange only (k-nearest in real space).
- Provides tangent-space operations on (S^2)^N: projected gradient, local basis,
  retraction (exact exponential map on S^2), and parallel transport along that geodesic.
"""
from __future__ import annotations

import numpy as np
import scipy.spatial as spt


def _as_flat_3n(x: np.ndarray, *, name: str) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    if x.ndim == 1:
        if x.size % 3 != 0:
            raise ValueError(f"{name} must have length 3N, got {x.size}")
        return x.copy()
    if x.ndim == 2 and x.shape[1] == 3:
        return x.reshape(-1).copy()
    raise ValueError(f"{name} must have shape (3N,) or (N,3), got {x.shape}")


def _view_n3(x_flat: np.ndarray) -> np.ndarray:
    return x_flat.reshape((-1, 3))


def _normalize_rows(x_n3: np.ndarray, eps: float = 1e-15) -> np.ndarray:
    n = np.linalg.norm(x_n3, axis=1)
    n = np.maximum(n, eps)
    return x_n3 / n[:, None]


def project_to_tangent(m_n3: np.ndarray, v_n3: np.ndarray) -> np.ndarray:
    """Project v to the tangent space at m on S^2: v - (v·m)m (row-wise)."""
    return v_n3 - np.sum(v_n3 * m_n3, axis=1)[:, None] * m_n3


def _rodrigues_rotate(v_n3: np.ndarray, k_n3: np.ndarray, theta: np.ndarray) -> np.ndarray:
    """
    Rodrigues rotation row-wise.
    v_n3, k_n3: (N,3) with k_n3 ~ unit
    theta: (N,)
    """
    ct = np.cos(theta)[:, None]
    st = np.sin(theta)[:, None]
    kv = np.cross(k_n3, v_n3)
    k_dot_v = np.sum(k_n3 * v_n3, axis=1)[:, None]
    return v_n3 * ct + kv * st + k_n3 * k_dot_v * (1.0 - ct)


class Magnetization:
    """
    Magnetization on (S^2)^N with open boundaries.

    Stores:
      - points as flat (3N,)
      - spins  as flat (3N,)
    Provides views as (N,3) via points_n3 / spins_n3 properties.
    """

    def __init__(self, points: np.ndarray, spins: np.ndarray, *, exchange_constant: float = 1.0,
                 lattice_constant: float = 1.0, k_neighbors: int = 4, normalize_spins: bool = True,
                 rebuild_neighbors: bool = True) -> None:
        """
        :param points: (3N,) or (N,3)
        :param spins:  (3N,) or (N,3)
        :param exchange_constant: nearest-neighbor Heisenberg J
        :param lattice_constant: stored, not used directly here
        :param k_neighbors: number of neighbors per site (excluding self)
        :param normalize_spins: normalize spins onto S^2 at init
        :param rebuild_neighbors: build neighbor list now (True). If False, you must call rebuild_neighbor_list().
        """
        # In case the input arrays are already flattened this does not change these arrays.
        self._points = _as_flat_3n(points, name="points")
        self._spins = _as_flat_3n(spins, name="spins")

        if self._points.size != self._spins.size:
            raise ValueError("points and spins must have the same length 3N")

        self._exchange_constant = float(exchange_constant)
        self._lattice_constant = float(lattice_constant)
        self._k_neighbors = int(k_neighbors)
        if self._k_neighbors < 1:
            raise ValueError("k_neighbors must be >= 1")

        if normalize_spins:
            self._spins = _normalize_rows(self.spins_n3).reshape(-1)

        self._kdtree: spt.cKDTree | None = None
        self._neighbors: np.ndarray | None = None

        if rebuild_neighbors:
            self.rebuild_neighbor_list()

    # -------------------------
    # Views / basic properties
    # -------------------------
    @property
    def N(self) -> int:
        return self._spins.size // 3

    @property
    def points(self) -> np.ndarray:
        """Flat points array, shape (3N,)."""
        return self._points

    @property
    def spins(self) -> np.ndarray:
        """Flat spins array, shape (3N,)."""
        return self._spins

    @property
    def points_n3(self) -> np.ndarray:
        """Points view, shape (N,3)."""
        return _view_n3(self._points)

    @property
    def spins_n3(self) -> np.ndarray:
        """Spins view, shape (N,3)."""
        return _view_n3(self._spins)

    @property
    def exchange_constant(self) -> float:
        return self._exchange_constant

    @exchange_constant.setter
    def exchange_constant(self, value: float) -> None:
        self._exchange_constant = float(value)

    @property
    def lattice_constant(self) -> float:
        return self._lattice_constant

    @property
    def k_neighbors(self) -> int:
        return self._k_neighbors

    def set_spins(self, spins: np.ndarray, *, normalize: bool = True) -> None:
        spins_flat = _as_flat_3n(spins, name="spins")
        if spins_flat.size != self._spins.size:
            raise ValueError("new spins must have the same length 3N")
        if normalize:
            spins_flat = _normalize_rows(_view_n3(spins_flat)).reshape(-1)
        self._spins = spins_flat

    def set_points(self, points: np.ndarray, *, rebuild_neighbors: bool = True) -> None:
        points_flat = _as_flat_3n(points, name="points")
        if points_flat.size != self._points.size:
            raise ValueError("new points must have the same length 3N")
        self._points = points_flat
        if rebuild_neighbors:
            self.rebuild_neighbor_list()

    def copy(self) -> "Magnetization":
        out = Magnetization(points=self._points.copy(), spins=self._spins.copy(),
                            exchange_constant=self._exchange_constant, lattice_constant=self._lattice_constant,
                            k_neighbors=self._k_neighbors, normalize_spins=False,  # already normalized
                            rebuild_neighbors=False,  # we'll copy neighbor data if present
                            )
        if self._kdtree is not None:
            out._kdtree = self._kdtree  # safe to share (points immutable unless set_points called)
        if self._neighbors is not None:
            out._neighbors = self._neighbors.copy()
        return out

    # -------------------------
    # Neighbor list
    # -------------------------

    def rebuild_neighbor_list(self) -> None:
        """(Re)build KDTree and kNN neighbor list for open boundaries."""
        pts = self.points_n3
        self._kdtree = spt.cKDTree(pts)

        # Query k+1 and drop self
        _, idx = self._kdtree.query(pts, k=self._k_neighbors + 1)

        neigh = np.empty((self.N, self._k_neighbors), dtype=int)
        for i in range(self.N):
            row = np.atleast_1d(idx[i])
            row = row[row != i]
            if row.size < self._k_neighbors:
                raise ValueError(
                    f"Could not find {self._k_neighbors} neighbors for site {i}. "
                    "Try smaller k_neighbors or check points for duplicates / small N."
                )
            neigh[i] = row[: self._k_neighbors]
        self._neighbors = neigh

    def neighbors(self, i: int) -> np.ndarray:
        """Neighbor indices of site i, shape (k_neighbors,)."""
        if self._neighbors is None:
            raise RuntimeError("Neighbor list not built. Call rebuild_neighbor_list().")
        return self._neighbors[int(i)].copy()

    # -------------------------
    # Energy / gradients
    # -------------------------

    def energy(self) -> float:
        """
        Nearest-neighbor exchange (counted once):
            E = -J * sum_{<i,j>} s_i · s_j
        """
        if self._neighbors is None:
            raise RuntimeError("Neighbor list not built. Call rebuild_neighbor_list().")

        J = self._exchange_constant
        s = self.spins_n3

        E = 0.0
        for i in range(self.N):
            si = s[i]
            for j in self._neighbors[i]:
                if j > i:
                    E -= J * float(np.dot(si, s[j]))
        return E

    def gradient(self) -> np.ndarray:
        """
        Euclidean gradient in embedding space for exchange energy:
            ∂E/∂s_i = -J * sum_{j in N(i)} s_j
        Returns flattened (3N,).
        """
        if self._neighbors is None:
            raise RuntimeError("Neighbor list not built. Call rebuild_neighbor_list().")

        J = self._exchange_constant
        s = self.spins_n3
        g = np.zeros_like(s)

        for i in range(self.N):
            g[i] = -J * np.sum(s[self._neighbors[i]], axis=0)
        return g.reshape(-1)

    def gradient_tspace_3N(self) -> np.ndarray:
        """
        Riemannian gradient in embedding representation, flattened (3N,):
            g_R = g - (g·s) s
        """
        s = self.spins_n3
        g = _view_n3(self.gradient())
        g_tan = project_to_tangent(s, g)
        return g_tan.reshape(-1)

    # -------------------------
    # Tangent basis / coordinate transforms
    # -------------------------

    def basis(self) -> np.ndarray:
        """
        Orthonormal tangent basis for each spin.

        Returns shape (N,2,3):
            basis[i,0] and basis[i,1] are unit vectors in R^3 tangent to S^2 at spin i.
        """
        m = self.spins_n3
        B = np.zeros((self.N, 2, 3), dtype=float)

        z = np.array([0.0, 0.0, 1.0])
        x = np.array([1.0, 0.0, 0.0])

        for i, mi in enumerate(m):
            ref = x if abs(mi[2]) > 0.9 else z
            e1 = ref - np.dot(ref, mi) * mi
            n1 = np.linalg.norm(e1)
            if n1 < 1e-14:
                ref = z if ref is x else x
                e1 = ref - np.dot(ref, mi) * mi
                n1 = np.linalg.norm(e1)
            e1 /= n1
            e2 = np.cross(mi, e1)
            e2 /= np.linalg.norm(e2)
            B[i, 0] = e1
            B[i, 1] = e2
        return B

    def project_to_basis(self, vec_embedding_space: np.ndarray) -> np.ndarray:
        """
        Project vectors (flattened 3N or (N,3)) to local tangent coordinates.

        Output flattened (2N,):
          [xi_0, eta_0, xi_1, eta_1, ..., xi_{N-1}, eta_{N-1}]
        """
        v_flat = _as_flat_3n(vec_embedding_space, name="vec_embedding_space")
        if v_flat.size != self._spins.size:
            raise ValueError("vec_embedding_space must have length 3N")
        v = _view_n3(v_flat)

        B = self.basis()  # (N,2,3)
        coords = np.einsum("nai,ni->na", B, v)  # (N,2)
        return coords.reshape(-1)

    def lift_from_basis(self, vec_tangent_coords: np.ndarray) -> np.ndarray:
        """
        Convert local tangent coordinates back to embedding tangent vectors.

        Input:
          vec_tangent_coords: flattened (2N,) or (N,2)
        Output:
          flattened (3N,)
        """
        c = np.asarray(vec_tangent_coords, dtype=float)
        if c.ndim == 1:
            if c.size != 2 * self.N:
                raise ValueError("flattened tangent coords must have length 2N")
            c = c.reshape(self.N, 2)
        elif c.shape != (self.N, 2):
            raise ValueError("vec_tangent_coords must be (2N,) or (N,2)")

        B = self.basis()  # (N,2,3)
        v = np.einsum("na,nai->ni", c, B)  # (N,3)
        v = project_to_tangent(self.spins_n3, v)  # enforce tangency numerically
        return v.reshape(-1)

    def gradient_tspace_2N(self) -> np.ndarray:
        """Riemannian gradient in local coordinates, flattened (2N,)."""
        g_tan = _view_n3(self.gradient_tspace_3N())
        return self.project_to_basis(g_tan)

    def finite_difference_HX(self, X: np.ndarray, technique: str = "simple_fd", eps: float = 1.0e-6) -> np.ndarray:
        r"""
        Calculate action of the Hessian on X using finite differences

        :param X: current estimate of lowest subspace
        :param technique: the finite difference technique, currently implemented are:
            - "simple_fd": simple forward finite difference
        :param eps: displacement parameter
        :return: the action of the Hessian
        """
        shape = np.shape(X)
        dim_subspace = shape[2]
        HX = np.zeros_like(X)
        for p in range(dim_subspace):
            if technique == "simple_fd":
                HX[:, p] = self.fd_simple_HcolX(x=X[:, p], eps=eps)
            else:
                raise NotImplementedError("Finite difference procedure not yet coded.")
        return HX

    def fd_simple_HcolX(self, x: np.ndarray, eps: float) -> np.ndarray:
        r"""
        Computes finite difference action of the Hessian applied to a column of X.
        To avoid frequent basis changes we will do this computation in 3N embedding space. So the procedure is the
        following:
            - Compute the 3N representation of the gradient of the magnetization (however, this is tangent to mag.)
            - Lift the column vector x from the tangent space to the embedding space (3N)
            - Rotate the magnetization along the tangent direction that is provided by x_3N
            - Compute the gradient of the displaced magnetization
            - Transport this gradient to the tangent space of the not-displaced original magnetization
            - Compute the finite-difference of the gradients
            - Project to 2N tangent space representation

        Note: This is the naive implementation of a simple forward finite difference scheme with a fixed step. Other,
        more sophisticated schemes can be found as well in this class.

        :param x: A 2N-dimensional flattened column of X
        :param eps: The displacement parameter
        :return: The finite-difference approximation of H*x (of shape (2N))
        """
        grad_3N = self.gradient_tspace_3N()
        x_3N = self.lift_from_basis(vec_tangent_coords=x)
        mag_displaced = Magnetization.retraction(current_mag=self, vec_tspace=x_3N, displacement_parameter=eps)
        grad_3N_displaced = mag_displaced.gradient_tspace_3N()
        grad_3N_displaced_transported = Magnetization.parallel_transport(current_mag=mag_displaced,
                                                                         transport_vec=x_3N,
                                                                         vec_tspace=grad_3N_displaced,
                                                                         displacement_parameter=-1.0 * eps)
        fin_diff_3N = (grad_3N_displaced_transported - grad_3N) / eps
        Hx_2N = self.project_to_basis(vec_embedding_space=fin_diff_3N)
        return Hx_2N


    def fd_richardson_HcolX(self, x):
            #@todo
    # -------------------------
    # Retraction / parallel transport
    # -------------------------
    @classmethod
    def retraction(cls, current_mag: "Magnetization", vec_tspace: np.ndarray, displacement_parameter: float = 1.0,
                   rodrigues_threshold: float = 1.0e-8) -> "Magnetization":
        """
        Retraction using the exact exponential map on S^2, applied sitewise.

        vec_tspace: flattened (3N,) or (N,3) tangent vector at current spins.
        """
        v_flat = _as_flat_3n(vec_tspace, name="vec_tspace")
        if v_flat.size != current_mag.spins.size:
            raise ValueError("vec_tspace must have length 3N matching spins")

        m = current_mag.spins_n3
        v = _view_n3(v_flat)

        v = project_to_tangent(m, v) * float(displacement_parameter)
        theta = np.linalg.norm(v, axis=1)  # (N,)

        out = np.empty_like(m)
        small = theta < rodrigues_threshold
        large = ~small

        if np.any(large):
            th = theta[large]
            u = v[large] / th[:, None]
            out[large] = m[large] * np.cos(th)[:, None] + u * np.sin(th)[:, None]

        if np.any(small):
            th = theta[small]
            th2 = th * th
            # sin(th)/th and cos(th) series
            sinc = 1.0 - th2 / 6.0 + (th2 * th2) / 120.0
            cost = 1.0 - th2 / 2.0 + (th2 * th2) / 24.0
            out[small] = m[small] * cost[:, None] + v[small] * sinc[:, None]

        out = _normalize_rows(out)

        # Keep same geometry/neighbor settings; points unchanged
        new_obj = cls(points=current_mag.points.copy(), spins=out.reshape(-1),
                      exchange_constant=current_mag.exchange_constant, lattice_constant=current_mag.lattice_constant,
                      k_neighbors=current_mag.k_neighbors, normalize_spins=False, rebuild_neighbors=False,
                      # reuse neighbor list for same points
                      )
        # Reuse neighbor info safely (same points)
        new_obj._kdtree = current_mag._kdtree
        new_obj._neighbors = None if current_mag._neighbors is None else current_mag._neighbors.copy()
        return new_obj

    @classmethod
    def parallel_transport(cls, current_mag: "Magnetization", transport_vec: np.ndarray, vec_tspace: np.ndarray,
                           displacement_parameter: float = 1.0, rodrigues_threshold: float = 1.0e-8) -> np.ndarray:
        """
        Parallel transport a tangent vector along the geodesic induced by vec_tspace
        (sitewise on S^2). For S^2, this equals rotating the vector by the same
        Rodrigues rotation that moves the spin.

        Inputs:
          transport_vec: flattened (3N,) or (N,3), tangent at current spins
          vec_tspace:     flattened (3N,) or (N,3), tangent displacement at current spins

        Output:
          flattened (3N,), tangent at the retracted spins.
        """
        v_flat = _as_flat_3n(transport_vec, name="transport_vec")
        d_flat = _as_flat_3n(vec_tspace, name="vec_tspace")

        if v_flat.size != current_mag.spins.size or d_flat.size != current_mag.spins.size:
            raise ValueError("transport_vec and vec_tspace must have length 3N matching spins")

        m = current_mag.spins_n3
        v = _view_n3(v_flat)
        d = _view_n3(d_flat)

        v = project_to_tangent(m, v)
        d = project_to_tangent(m, d) * float(displacement_parameter)

        theta = np.linalg.norm(d, axis=1)  # (N,)
        out = np.empty_like(v)

        small = theta < rodrigues_threshold
        large = ~small

        if np.any(large):
            th = theta[large]
            u = d[large] / th[:, None]
            k = np.cross(m[large], u)
            k = _normalize_rows(k)
            out[large] = _rodrigues_rotate(v[large], k, th)

        if np.any(small):
            out_small = v[small].copy()
            md = np.cross(m[small], d[small])
            md_norm = np.linalg.norm(md, axis=1)
            ok = md_norm > 1e-15
            if np.any(ok):
                k = md[ok] / md_norm[ok][:, None]
                th = theta[small][ok]
                out_small[ok] = _rodrigues_rotate(v[small][ok], k, th)
            out[small] = out_small

        # Ensure tangency at new spins
        new_mag = cls.retraction(current_mag, d.reshape(-1), displacement_parameter=1.0,
                                 rodrigues_threshold=rodrigues_threshold)
        out = project_to_tangent(new_mag.spins_n3, out)
        return out.reshape(-1)
