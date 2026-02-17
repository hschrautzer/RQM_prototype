r"""

"""
import numpy as np

from prototype.magnetization import magnetization

class RQM:

    def __init__(self,magnetization: magnetization, epsilon_fin_diff: float = 1.0, conv_crit: float) -> None:
        r"""

        """
        self._epsilon = epsilon_fin_diff
        self._conv_crit = conv_crit
        self.magnetization = magnetization

    def HX_findiff(self, X: np.ndarray) -> np.ndarray:
        r"""
        Calculate action of the Hessian on X using finite differences

        :param X: current estimate of lowest subspace
        :return: the action of the Hessian
        """
        shape = np.shape(X)
        dim_subspace = shape[2]
        HX = np.zeros_like(X)
        for p in range(dim_subspace):
            HX[:,:,p] = self.H_colX_findiff(x_col=X[:,:,p])

    def H_colX_findiff(self, x_col) -> np.ndarray:
        r"""
        Computes finite difference action of a column of X with hessian

        :param x_col:
        :return:
        """
        grad = self.magnetization.gradient_tspace_3N()
        mag_rotated = self.magnetization.rotate_along(vec=x_col,displacement_parameter=self._epsilon)
        grad_forward = mag_rotated.gradient_tspace_3N()
        # rotate back to tangent frame of initial configuration
        grad_forward_rotated_back = mag_rotated.parallel_transport(grad_forward,x_col,displacement_parameter=-self._epsilon)
        return (grad - grad_forward_rotated_back) / self._epsilon


    def __call__(self,v_ini: np.ndarray, N_steps: int) -> np.ndarray:
        r"""
        :param v_ini: Initial guess for evecs of lowest subspace but in 3N embedding space representation
        :param args:
        :param kwargs:
        :return:
        """
        # Dimension of lowest subspace
        dim_subspace = np.shape(v_ini)[1]


        X = np.zeros(shape=(len(self.magnetization.spins),2,dim_subspace))
        # Project initial guess to tangent space
        for p in range(dim_subspace):
            X[:,:,p] = self.magnetization.project_to_basis(vec_embedding_space=v_ini[:,:,p])

        # QR decomposition of X

        # Compute product of Hessian with each column of X (HX is of shape (N,2,p))
        HX = self.HX_findiff(X)

        # Compute X^T * HX
        rq_matrix = np.zeros(shape=(dim_subspace,dim_subspace))
        #@todo

        # Compute Rayleigh Quotient Gradient  2*HX - 2X*(X^T H X) (is again of shape (N,2,p))
        # rq_gradient =
        # @todo

        # Compute the Frobenius Norm of this gradient (used as convergence criterium)
        norm_rq = None #@todo

        for n in range(N_steps):


