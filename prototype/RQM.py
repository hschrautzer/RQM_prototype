r"""

"""
import logging
import numpy as np
from typing import Union, Tuple
from pathlib import Path
from prototype.magnetization import Magnetization

class RQM:
    r"""

    """

    @staticmethod
    def _setup_logger(name: str, verbose: int, logfilepath: Union[Path, None]) -> logging.Logger:
        r"""
        Organizes the setup of the Logger.

        :param name: The name of the logger. A good idea is to choose __name__ within the class which creates the logger
        :param verbose: The verbose-level: 10 Debug, 20: Information, 30: Warning, 40: Error
        :param logfilepath: The path defining the log-file. If None no file will be created and only the stdout will be
            used for logging.
        :return: The instance of the Logger
        """
        l_logger = logging.getLogger(name)
        if verbose in [10, 20, 30, 40]:
            formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
            console_handler = logging.StreamHandler()
            console_handler.setLevel(verbose)
            console_handler.setFormatter(formatter)
            l_logger.addHandler(console_handler)
            if logfilepath is not None:
                file_handler = logging.FileHandler(logfilepath)
                file_handler.setLevel(verbose)
                file_handler.setFormatter(formatter)
                l_logger.addHandler(file_handler)
            l_logger.setLevel(verbose)
        else:
            raise ValueError("Not a valid verbose-level")
        return l_logger

    def _test_input_vector(self, v_ini: np.ndarray) -> Tuple[int,int]:
        r"""
        Sanity check of input vector for RQM.

        :param v_ini: the input vector.
        :return: the number of atoms and the dimension of the subspace
        """
        shape_v_ini = np.shape(v_ini)
        dim_subspace = shape_v_ini[2]
        # Number of atoms
        N = shape_v_ini[0]
        if N != len(self._mag.spins):
            self._logger.error(f"Length of the input vector ({N}) does not match length of spin-array"
                               f" ({len(self._mag.spins)}")
            raise ValueError(f"Length of the input vector ({N}) does not match length of spin-array"
                             f" ({len(self._mag.spins)}")
        if shape_v_ini[1] != 3:
            self._logger.error(f"Dimensionality of input vector component is expected to be 3,"
                               f" while I got {shape_v_ini[1]}")
            raise ValueError(f"Dimensionality of input vector component is expected to be 3,"
                             f" while I got {shape_v_ini[1]}")


    def __init__(self, mag: Magnetization, simudir: Union[Path, None] = None,
                 epsilon_fin_diff: float = 1.0, conv_crit: float) -> None:
        r"""
        Initializes the Rayleigh Quotient Minimization.

        :param mag: The magnetization instance storing the points and spins.
        :param simudir: The directory where the simulation data shall be stored. If None a folder called "RQM" is
            created in the current working directory.
        """
        # The magnetization does not change throughout the RQM algorithm. We just need
        # that to estimate the Hessian using finite differences.
        self._mag: Magnetization = mag

        # Setup logger and simulation directory.
        if simudir is None:
            self._simudir = Path.cwd() / "RQM"
            self._simudir.mkdir(exist_ok=True)
        else:
            self._simudir = simudir
        self._logger = self._setup_logger(name=self.__class__.__name__, verbose=20,
                                          logfilepath=self._simudir / f"rqm.log")

        self._epsilon = epsilon_fin_diff
        self._conv_crit = conv_crit

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


    def __call__(self, v_ini: np.ndarray, N_steps: int) -> np.ndarray:
        r"""
        Calls the RQM minimization initialized with v_ini.

        :param v_ini: Initial guess for evecs of lowest subspace but in 3N embedding space representation. It's shape is
            [N,3,p], where p is the dimension of the subspace and N is the number of atoms/spins and 3 encodes the x,y,z
            component of each vector
        :return:
        """
        self._logger.info("Starting RQM...")
        N, dim_subspace = self._test_input_vector(v_ini=v_ini)

        # The 3N embedding space representation needs to the projected to 2N tangent space, since we are looking for
        # parametrisations of this 2N tangent space. If we would carry the superfluous N degrees of freedom, numerical
        # errors would accumulate.
        X = np.zeros(shape=(len(self._mag.spins),2,dim_subspace))
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





