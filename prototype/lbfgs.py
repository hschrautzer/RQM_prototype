import numpy as np
from dataclasses import dataclass
from magnetization import *
from utils import *
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
	N_memory: int 				# Number of entries in memory
	theta_max: np.float64

	# Persistent variables
	iteration: int
	force_prev: np.ndarray
	rq_step_prev: np.ndarray
	steplength_prev: np.float64
	X_prev: np.ndarray
	rho: np.ndarray 				# (d·y)**-1
	gamma: np.ndarray
	d: np.ndarray 					# Difference in spin configurations in 2N
	y: np.ndarray

	# Volatile variables
	force: np.ndarray # 2xN_spinsxN_modes
	rq_step: np.ndarray
	steplength: np.float64

	# Optimization parameters
	N_iter: int
	rq_grad_tol: float

	# Outputs
	result: dict

	def _init__(self,
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
			self.rq_step = np.zeros([2*self.N_spins, self.N_modes])
			self.rho = np.zeros([self.N_memory])
			self.d = np.zeros([2*self.N_spins, self.N_modes, self.N_memory])
			self.y = np.zeros([2*self.N_spins, self.N_modes, self.N_memory])
			self.gamma = np.zeros([self.N_memory])
			self.dummy_step = np.zeros([2*self.N_spins,self.N_modes])
			self.step_previous_dummy = np.zeros([2*self.N_spins, self.N_modes])

			self.result["warnings"] = []
			self.result["rq_gradient_norm"] = []

	def calc_step(self, force):
		# Selecting memory index for LIFO history
		ind_memory = np.mod(self.iteration, self.N_memory)
		if self.iteration == 0:
			# First iteration
			self.rq_step = force
		else:
			self.d[:,:, ind_memory] = self.steplength_prev * self.rq_step_prev
			self.y[:,:, ind_memory] = self.force_prev - force

			y_dot_d = np.dot(self.y[:,:,ind_memory].ravel(),
							 self.d[:,:,ind_memory].ravel()) # SHARP
			self.rho[ind_memory] = 1/y_dot_d

			if self.rho[ind_memory] < 0:
				self.rq_step = force
				# Wipe memory
				return
			q = -1.0 * force

			for l in range(self.N_memory, 0, -1): # MEGASHARP (+-)1?
				j = np.mod(l+ind_memory,self.N_memory) # (+-)1?, also, odd that this is adding memory index as an offset
				q_dot_d = np.dot(q[:,:,:].ravel(),
							 self.d[:,:,j].ravel()) #SHARP?
				self.gamma[j] = self.rho[j]*q_dot_d
				q[:] = q[:] - self.gamma[j] * self.y[:,:,j]

			y_dot_y = np.dot(self.y[:,:, ind_memory].ravel(),
							 self.d[:,:,ind_memory].ravel()) #SHARP
			
			self.dummy_step[:,:] = q[:,:,:] / (self.rho[ind_memory]* y_dot_y)
			
			for l in range(1,self.N_memory): # (+-)1?
				if self.iteration <= self.N_memory:
					j = l
				else:
					j = np.mod(l+ind_memory,self.N_memory) # (+-)1?
				res = np.dot(self.y[:,:,j].ravel(),
							 self.dummy_step[:,:].ravel()) #SHARP
				self.dummy_step[:,:] = self.dummy_step[:,:] \
								+ self.d[:,:,j] * (self.gamma[j]-self.rho[j]*res)
			self.rq_step = -1.0*self.dummy_step

	def step(self,force):
		self.calc_step(force)
		theta_rms = np.dot(self.rq_step[:,:].ravel(),
							self.rq_step[:,:].ravel()) #SHARP
		theta_rms = np.sqrt(theta_rms)
		t_res = self.theta_max/theta_rms
		if t_res < 1:
			self.steplength = t_res
		else:
			self.steplength = np.float64(1.0)

	def rq_force_calc(self, mag: Magnetization, X_n2: np.ndarray) -> tuple[np.ndarray, float]:
		t_HX = mag.finite_difference_HX(X_n2) # 2N_spins x N_modes

		rq_matrix = np.zeros([self.N_modes, self.N_modes])
		rq_matrix = np.reshape(X_n2, [2,-1]).T @ t_HX #SHARP

		rq = np.linalg.trace(rq_matrix)

		# Check this ↓↓↓
		rq_gradient = t_HX # 2N_spins x N_modes
		rq_gradient = rq_gradient - X_n2 @ rq_matrix # 2N_spins x N_modes - 2N_spins x N_modes @ N_modes x N_modes
		rq_gradient = 2 * rq_gradient # 2N_spins x N_modes

		rq_gradient_norm = float(np.linalg.norm(rq_gradient, ord='fro'))
		self.result['rq_gradient_norm'].append(rq_gradient_norm)

		rq_force = - rq_gradient
		return rq_force, rq_gradient_norm

	def minimize(self, mag: Magnetization, vec_ini: np.ndarray):
		v_fin_n3 = np.zeros([3*self.N_spins])
		rq = 0
		X_n2 = np.zeros([2*self.N_spins, self.N_modes])
		evals = 0
		for ind_modes in range(self.N_modes):
			vec_ini[:,ind_modes] = mag.project_to_basis(vec_ini[:,ind_modes])
			vec_ini[:,ind_modes] = vec_ini[:,ind_modes] / norm_n2
		
		X_n2 = GM_retraction(X_n2)

		# t_HX = mag.finite_difference_HX(X_n2) # 2N_spins x N_modes

		# rq_matrix = np.zeros([self.N_modes, self.N_modes])
		# rq_matrix = np.reshape(X_n2, [2,-1]).T @ t_HX #SHARP

		# rq = np.linalg.trace(rq_matrix)

		# # Check this ↓↓↓
		# rq_gradient = t_HX # 2N_spins x N_modes
		# rq_gradient = rq_gradient - X_n2 @ rq_matrix # 2N_spins x N_modes - 2N_spins x N_modes @ N_modes x N_modes
		# rq_gradient = 2 * rq_gradient # 2N_spins x N_modes

		rq_force, rq_gradient_norm = self.rq_force_calc(mag, X_n2)

		for iter in range(self.N_iter):
			if rq_gradient_norm <= self.rq_grad_tol:
				self.result['rqm_iterations'] = iter
				self.result['status'] = "converged"
				break
			if iter >= self.N_iter:
				self.result['rqm_iterations'] = iter
				self.result['warnings'].append("RQM exceeded iterations")
				self.result['status'] = "not converged"
			self.step(rq_force)
			self.force_prev = self.force
			self.steplength_prev = self.steplength
			self.rq_step_prev = self.rq_step

			self.X_prev = X_n2

			U, S, VT = np.linalg.svd(self.rq_step)

			X_n2 = GM_retraction_exp(X=X_n2,U=U,S=S,VT=VT,delta=self.steplength)
			# grassman retraction w. USV
			# transports w. USV






		return self.result
