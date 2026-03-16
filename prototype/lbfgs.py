import numpy as np
from dataclasses import dataclass
from magnetization import *
from utils import *

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
	force_previous: np.ndarray
	step_previous: np.ndarray
	steplength_previous: np.float64
	rho: np.ndarray 				# (d·y)**-1
	gamma: np.ndarray
	d: np.ndarray 					# Difference in spin configurations in 2N
	y: np.ndarray

	# Volatile variables
	force: np.ndarray # 2xN_spinsxN_modes
	new_step: np.ndarray
	steplength: np.float64

	def _init__(self,
				N_memory: int,
				theta_max: np.float64,
				N_spins: int,
				N_modes: int):
			
			self.N_memory = N_memory
			self.theta_max = theta_max
			self.N_spins = N_spins
			self.N_modes = N_modes
			del N_memory
			del theta_max
			del N_spins
			del N_modes
			# Initialization
			self.new_step = np.zeros([2*self.N_spins, self.N_modes])
			self.rho = np.zeros([self.N_memory])
			self.d = np.zeros([2*self.N_spins, self.N_modes, self.N_memory])
			self.y = np.zeros([2*self.N_spins, self.N_modes, self.N_memory])
			self.gamma = np.zeros([self.N_memory])
			self.dummy_step = np.zeros([2*self.N_spins,self.N_modes])
			self.step_previous_dummy = np.zeros([2*self.N_spins, self.N_modes])

	def calc_step(self, force):
		# Selecting memory index for LIFO history
		ind_memory = np.mod(self.iteration, self.N_memory)
		if self.iteration == 0:
			# First iteration
			self.new_step = force
		else:
			self.d[:,:, ind_memory] = self.steplength_previous * self.step_previous
			self.y[:,:, ind_memory] = self.force_previous - force

			y_dot_d = np.dot(self.y[:,:,ind_memory].ravel(),
							 self.d[:,:,ind_memory].ravel()) # SHARP
			self.rho[ind_memory] = 1/y_dot_d

			if self.rho[ind_memory] < 0:
				self.new_step = force
				# Wipe memory
				return
			q = -1.0 * force

			for l in range(self.N_memory, 0, -1): # MEGASHARP (+-)1?
				j = np.mod(l+ind_memory,self.N_memory) # (+-)1?, also, odd that this is adding memory index as an offset
				q_dot_d = np.dot(q[:,:,:].ravel(),
							 self.d[:,:,j].ravel()) #	SHARP?
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
			self.new_step = -1.0*self.dummy_step

	def step(self,force):
		self.calc_step(force)
		theta_rms = np.dot(self.new_step[:,:].ravel(),
							self.new_step[:,:].ravel()) #SHARP
		theta_rms = np.sqrt(theta_rms)
		t_res = self.theta_max/theta_rms
		if t_res < 1:
			self.steplength = t_res
		else:
			self.steplength = np.float64(1.0)

	def minimize(self, mag: Magnetization, vec_ini: np.ndarray, basis):
		v_fin_n3 = np.zeros([3*self.N_spins])
		rq = 0
		X_n2 = np.zeros([2*self.N_spins, self.N_modes])
		evals = 0
		for ind_modes in range(self.N_modes):
			vec_ini[:,ind_modes] = mag.project_to_basis(vec_ini[:,ind_modes])
			vec_ini[:,ind_modes] = vec_ini[:,ind_modes] / norm_n2
		
		X_n2 = GM_retraction(X_n2)
