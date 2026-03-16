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
Retracts a (2,N,P) array onto the Grassmaniann using QR, and discarding R
"""
def GM_retraction(i_X: np.ndarray) -> np.ndarray:
	t_QR_res = np.linalg.qr(i_X)
	return t_QR_res.Q


"""
Performs a thin singular value decomposition of thin matrix X
"""