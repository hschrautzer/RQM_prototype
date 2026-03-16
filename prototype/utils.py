import numpy as np

"""
Calculates L2 norm of (2N,) vector, returns (N,)
"""
def norm_n2(i_vec: np.ndarray) -> np.ndarray:
	t_vec = i_vec**2
	t_vec = np.reshape(t_vec, [2, -1])
	t_vec = np.sum(t_vec, axis=0)
	return t_vec

"""
Calculates L2 norm of (3N,) vector, returns (N,)
"""
def norm_n3(i_vec: np.ndarray) -> np.ndarray:
	t_vec = i_vec**2
	t_vec = np.reshape(t_vec, [3, -1])
	t_vec = np.sum(t_vec, axis=0)
	return t_vec

"""
Retracts a (2,N,P) array onto the Grassmaniann using QR, and discarding R. This is the projection-like retraction
"""
def GM_retraction(i_X: np.ndarray) -> np.ndarray:
	t_QR_res = np.linalg.qr(i_X)
	return t_QR_res.Q


"""
Geodesic retraction
"""
def GM_retraction_exp(X: np.ndarray, U: np.ndarray, S: np.ndarray, VT: np.ndarray, delta: float) -> np.ndarray:
	r"""
	Performs the retraction along the "geodesic" matrix exponential. Implementation of the formula:

	X(delta)	=	X*V*cos(Sigma*delta)*V^T + U*sin(Sigma*delta)*V^T
				= 	(first_term + second_term)@V^T

	:param X: the current point on the Grassmannian (2N,p)
	:param U: left-singular value matrix (2N,p)
	:param S: singular values (p)
	:param VT: right-singular values matrix (p,p)
	:param delta: displacement parameter
	:return: the updated X(delta)
	"""
	first_term = X@VT.T
	second_term = U
	for p in range(np.shape(X)[1]):
		first_term[:, p] *= np.cos(S[p]) * delta
		second_term[:, p] *= np.sin(S[p]) * delta
	return (first_term + second_term)@VT

"""
Parallel Transport

The formula for transport a vector b (2N,p) from X (2N, p) along Delta X (represented by compact SVD -> U, S, VT) with
displacement parameter delta is given by:

b' = b - [X*V*sin(S*t)+U*cos(1-S*t)]*U^T*b,

which can be written as

b' = b - TM*U^T*b

with TM beeing the transport matrix.

Since you might want to transport several vectors and not only one you only want to calculate TM once and then just
apply it.
"""


def GM_calc_transportmatrix(X: np.ndarray, U: np.ndarray, S: np.ndarray, VT: np.ndarray, delta: float) -> np.ndarray:
	r"""
	Calculate the transport matrix given by:
		TM 	= X*V*sin(S*t)+U*cos(1-S*t)
			= first_term + second_term

	:param X: the current point on the Grassmannian (2N,p)
	:param U: left-singular value matrix (2N,p)
	:param S: singular values (p)
	:param VT: right-singular values matrix (p,p)
	:param delta: displacement parameter
	:return: TM (2N, p)
	"""
	first_term = X@VT.T
	second_term = U
	for p in range(np.shape(X)[1]):
		first_term[:, p] *= np.sin(S[p]) * delta
		second_term[:, p] *= (1-np.cos(S[p])) * delta
	return first_term + second_term


def GM_parallel_transport(b: np.ndarray, U: np.ndarray, TM: np.ndarray) -> np.ndarray:
	r"""
	Apply the parallel transport matrix TM to a vector b (2N,p).
		b' = b - TM*U^T*b

	:param: Quantity that shall be transported (2N,p)
	:param: left-singular value matrix (2N,p)
	:param: Transport matrix TM (2N,p)
	:return: transport quantity b'
	"""
	return b - TM@(U.T@b)




