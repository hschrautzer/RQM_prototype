import numpy as np
import pandas as pd
from pprint import pprint
import matplotlib.pyplot as plt

from pathlib import Path
import itertools as it

from prototype.magnetization import Magnetization
from prototype.lbfgs import lbfgs_minimizer

if __name__=="__main__":
	N = 2

	# Standard random number generator
	rng = np.random.default_rng()#seed=42)

	# Create NxNxN grid
	points = np.array([list(x) for x in it.product(range(N),repeat=3)])
	points = np.array(points)

	# Create random but normalized spins
	spins_n = rng.random((3*N**3,))
	for ind_spin in range(N):
		t_norm = np.linalg.norm(spins_n[3*N:3*N+3])
		spins_n[3*N:3*N+3] = spins_n[3*N:3*N+3] / t_norm

	mag = Magnetization(points=points,spins=spins_n, normalize_spins=True)

	E_grad = mag.gradient()
	
	N_epsilon = 10
	E_grad_fd_array = np.zeros((3*N**3,N_epsilon))
	epsilon_array = ((0.5)**np.arange(1,N_epsilon,1))*10**-4
	for index_epsilon in range(6):
		epsilon = epsilon_array[index_epsilon]
		recip_epsilon = 1/epsilon

		E_nm1 = np.zeros((3*N**3,))
		for ind_comp in range(E_nm1.shape[0]):
			t_spins = spins_n.copy()
			t_spins[ind_comp] = t_spins[ind_comp] - epsilon
			mag.set_spins(t_spins)
			E_nm1[ind_comp] = mag.energy()

		E_np1 = np.zeros((3*N**3,))
		for ind_comp in range(E_np1.shape[0]):
			t_spins = spins_n.copy()
			t_spins[ind_comp] = t_spins[ind_comp] + epsilon
			mag.set_spins(t_spins)
			E_np1[ind_comp] = mag.energy()

		E_grad_fd = (-0.5 * E_nm1 + 0.5*E_np1)*recip_epsilon
		E_grad_fd_array[:,index_epsilon] = E_grad_fd

	fig, ax = plt.subplots()
	ax.plot(E_grad_fd_array, ls="", marker="+")
	fig.show()

	input('Press any key to exit')

	# E_grad_err = E_grad - E_grad_fd
	# E_grad_err_abs = np.abs(E_grad_err)
	# print(f"{np.average(E_grad_err_abs)=}")
	# print(f"{np.min(E_grad_err_abs)=}")
	# print(f"{np.max(E_grad_err_abs)=}")
	
	# minimizer = lbfgs_minimizer(N_memory=3, theta_max=0.5, N_spins=mag.N, N_modes=4, N_iter=100, rq_grad_tol=1e-12)

	# # Initial vector to be random and in 3N space:
	# np.random.seed(2)
	# vec_ini_3N = np.random.rand(3*mag.N, 4)
	# results = minimizer.minimize(mag=mag, vec_ini=vec_ini_3N)
	# pprint(results)



