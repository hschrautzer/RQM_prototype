import numpy as np
from dataclasses import dataclass
from prototype.magnetization import *
from prototype.utils import *
import logging
from pathlib import Path
from typing import Union

"""
Recreated from Spinaker subroutine lbfgs_step_noncurved
and its calling subroutine
Additional source: Ivanov, et al, 2021, Appendix F (https://doi.org/10.1016/j.cpc.2020.107749).
"""


class lbfgs_minimizer:
	# Hyperparameters
	N_memory: int  # Number of entries in memory
	theta_max: np.float64

	# Persistent variables
	iteration: int
	force_prev: np.ndarray
	rq_step_prev: np.ndarray
	steplength_prev: np.float64
	X_prev: np.ndarray
	rho: np.ndarray  # (d·y)**-1
	gamma: np.ndarray
	d: np.ndarray  # Difference in spin configurations in 2N
	y: np.ndarray

	# RQ variables
	rq_matrix_pxp: np.ndarray
	rq_force_2nxp: np.ndarray
	rq_gradient_norm: float

	# Volatile variables
	force: np.ndarray  # 2xN_spinsxN_modes
	rq_step: np.ndarray
	steplength: np.float64

	# Optimization parameters
	N_iter: int
	rq_grad_tol: float

	# Outputs
	result: dict

	def __init__(self,
				N_memory: int,
				theta_max: np.float64,
				N_spins: int,
				N_modes: int,
				N_iter: int,
				rq_grad_tol: float):

		self.N_memory = N_memory
		self.theta_max = theta_max
		self.N_spins = N_spins
		self.N_modes = N_modes
		self.N_iter = N_iter
		self.rq_grad_tol = rq_grad_tol
		del N_memory
		del theta_max
		del N_spins
		del N_modes
		del N_iter
		del rq_grad_tol
		# Initialization
		self.rq_step = np.zeros([2 * self.N_spins, self.N_modes])
		self.rho = np.zeros([self.N_memory])
		self.d = np.zeros([2 * self.N_spins, self.N_modes, self.N_memory])
		self.y = np.zeros([2 * self.N_spins, self.N_modes, self.N_memory])
		self.gamma = np.zeros([self.N_memory])
		self.dummy_step = np.zeros([2 * self.N_spins, self.N_modes])
		self.step_previous_dummy = np.zeros([2 * self.N_spins, self.N_modes])
		self.result = {"warnings": [], "rq_gradient_norm": []}

	def calc_step(self, force):
		# Selecting memory index for LIFO history
		ind_memory = np.mod(self.iteration, self.N_memory)
		if self.iteration == 0:
			# First iteration
			self.rq_step = force
		else:
			self.d[:, :, ind_memory] = self.steplength_prev * self.rq_step_prev
			self.y[:, :, ind_memory] = self.force_prev - force

			y_dot_d = np.dot(self.y[:, :, ind_memory].ravel(),
							 self.d[:, :, ind_memory].ravel())  # SHARP
			self.rho[ind_memory] = 1 / y_dot_d

			if self.rho[ind_memory] < 0:
				self.rq_step = force
				# Wipe memory
				return
			q = -1.0 * force

			for l in range(self.N_memory, 0, -1):  # MEGASHARP (+-)1?
				j = np.mod(l + ind_memory,
						   self.N_memory)  # (+-)1?, also, odd that this is adding memory index as an offset
				q_dot_d = np.dot(q[:, :].ravel(),
								 self.d[:, :, j].ravel())  # SHARP?
				self.gamma[j] = self.rho[j] * q_dot_d
				q[:] = q[:] - self.gamma[j] * self.y[:, :, j]

			y_dot_y = np.dot(self.y[:, :, ind_memory].ravel(),
							 self.d[:, :, ind_memory].ravel())  # SHARP

			self.dummy_step[:, :] = q[:, :] / (self.rho[ind_memory] * y_dot_y)

			for l in range(1, self.N_memory):  # (+-)1?
				if self.iteration <= self.N_memory:
					j = l
				else:
					j = np.mod(l + ind_memory, self.N_memory)  # (+-)1?
				res = np.dot(self.y[:, :, j].ravel(),
							 self.dummy_step[:, :].ravel())  # SHARP
				self.dummy_step[:, :] = self.dummy_step[:, :] \
										+ self.d[:, :, j] * (self.gamma[j] - self.rho[j] * res)
			self.rq_step = -1.0 * self.dummy_step

	def step(self, force):
		self.calc_step(force)
		theta_rms = np.dot(self.rq_step[:, :].ravel(),
						   self.rq_step[:, :].ravel())  # SHARP
		theta_rms = np.sqrt(theta_rms)
		t_res = self.theta_max / theta_rms
		if t_res < 1:
			self.steplength = t_res
		else:
			self.steplength = np.float64(1.0)

	def rq_force_calc(self, mag: Magnetization, X_2np: np.ndarray) -> tuple[np.ndarray, float]:
		t_HX_2nxp = mag.finite_difference_HX(X_2np)  # 2N_spins x N_modes

		self.rq_matrix_pxp = np.zeros([self.N_modes, self.N_modes])
		#@olafur: I adjusted it to arbitrary number of modes (before it was 2). Furthermore I think it was wrong. It
		# should be X^T not X^T^T (which I think you were doing before) (I left the line below commented out for reference)
		#rq_matrix = np.reshape(X_n2, [self.N_modes, -1]).T @ t_HX  # SHARP

		#@olafur: below I was always rq_matrix instead of self.rq_matrix. I changed that.
		self.rq_matrix_pxp = X_2np.T @ t_HX_2nxp
		rq = np.trace(self.rq_matrix_pxp)

		# Check this ↓↓↓
		rq_gradient_2nxp = t_HX_2nxp  # 2N_spins x N_modes
		rq_gradient_2nxp = rq_gradient_2nxp - X_2np @ self.rq_matrix_pxp  # 2N_spins x N_modes - 2N_spins x N_modes @ N_modes x N_modes
		rq_gradient_2nxp = 2 * rq_gradient_2nxp  # 2N_spins x N_modes

		self.rq_gradient_norm = float(np.linalg.norm(rq_gradient_2nxp, ord='fro'))
		self.result['rq_gradient_norm'].append(self.rq_gradient_norm)

		self.rq_force_2nxp = - rq_gradient_2nxp

	def minimize(self, mag: Magnetization, vec_ini: np.ndarray):
		r"""
		Apply the minimization

		:param mag: Magnetization instance, contains the magnetic configuration
		:param vec_ini: initial vector. This should be in embedding space 3N representation.
		:return:
		"""
		#@olafur: previously this was of shape (3N), I changed that to (3N, p)
		v_fin_3nxp = np.zeros([3 * self.N_spins, self.N_modes])
		rq = 0
		X_2nxp = np.zeros([2 * self.N_spins, self.N_modes])
		evals = 0
		for ind_modes in range(self.N_modes):
			X_2nxp[:, ind_modes] = mag.project_to_basis(vec_ini[:, ind_modes])
			# @olafur: we don't need that. The QR decomposition below will take care.
			#X_n2[:, ind_modes] = vec_ini[:, ind_modes] / norm_n2

		X_2nxp = GM_retraction(X_2nxp)

		self.rq_force_calc(mag, X_2nxp)

		for iter in range(self.N_iter):
			self.iteration = iter

			if self.rq_gradient_norm <= self.rq_grad_tol:
				self.result['rqm_iterations'] = iter
				self.result['status'] = "converged"
				break
			if iter >= self.N_iter:
				self.result['rqm_iterations'] = iter
				self.result['warnings'].append("RQM exceeded iterations")
				self.result['status'] = "not converged"
			self.step(self.rq_force_2nxp)
			self.steplength_prev = self.steplength
			# @olafur (I commented this out, this is supposed to be updated by parallel transport, see below)
			#self.rq_step_prev = self.rq_step
			#self.force_prev = self.force
			self.X_prev = X_2nxp

			#@olafur: the full_matrices flag was missing, we want the compact SVD
			U, S, VT = np.linalg.svd(self.rq_step, full_matrices=False)

			# Retract the configuration
			X_2nxp = GM_retraction_exp(X=X_2nxp, U=U, S=S, VT=VT, delta=float(self.steplength))
			# For efficiency the parallel transport is done in 2 steps: computing the transport matrix and applying the
			# transport matrix
			TM = GM_calc_transportmatrix(X=self.X_prev,U=U,S=S,VT=VT,delta=float(self.steplength_prev))
			self.rq_step_prev = GM_parallel_transport(b=self.rq_step,U=U,TM=TM)
			self.force_prev = GM_parallel_transport(b=self.rq_force_2nxp,U=U,TM=TM)
			for k in range(self.N_memory):
				self.d[:,:,k] = GM_parallel_transport(self.d[:,:,k],U=U,TM=TM)
				self.y[:,:,k] = GM_parallel_transport(self.y[:,:,k],U=U,TM=TM)

			# Orthogonalize (to avoid accumulating numerical errors)
			X_2nxp = GM_retraction(X_2nxp)
			# Quantities for the next iteration
			self.rq_force_calc(mag,X_2nxp)

		# Compute the Ritz-Vectors to rotate the found minimum solution of R(X) to the eigenvector representation of H
		# that we want to compute
		t_eigval, t_eigvec = np.linalg.eigh(self.rq_matrix_pxp)
		print(t_eigval)
		print(self.rq_matrix_pxp)
		# Apply Rayleigh-Ritz and represent in 3N.
		#@olafur: we have to feed this Mode-wise to the projection. I corrected that.
		X_solution = X_2nxp @ t_eigvec
		for k in range(self.N_modes):
			v_fin_3nxp[:,k] = mag.lift_from_basis(X_solution[:,k])
		self.result['eigenvalues'] = t_eigval
		return self.result, v_fin_3nxp
