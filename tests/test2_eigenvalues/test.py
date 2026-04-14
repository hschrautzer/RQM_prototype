import numpy as np
import pandas as pd
from pprint import pprint
import matplotlib.pyplot as plt

from pathlib import Path
import itertools as it

from prototype.magnetization import Magnetization
from prototype.lbfgs import lbfgs_minimizer

if __name__=="__main__":
	N = 10

	# Standard random number generator
	rng = np.random.default_rng()#seed=42)

	# Create NxNxN grid
	points = np.array([list(x) for x in it.product(range(N),repeat=3)])
	points = np.array(points)

	# Create random but normalized spins
	spins = rng.random((3, N**3))
	norms = np.linalg.norm(spins, axis=0)
	for ind_row, it_row in enumerate(spins):
		spins[ind_row] = spins[ind_row] / norms

	mag = Magnetization(points=points,spins=spins)

	minimizer = lbfgs_minimizer(N_memory=3, theta_max=0.5, N_spins=mag.N, N_modes=4, N_iter=100, rq_grad_tol=1e-12)

	# Initial vector to be random and in 3N space:
	np.random.seed(2)
	vec_ini_3N = np.random.rand(3*mag.N, 4)
	results = minimizer.minimize(mag=mag, vec_ini=vec_ini_3N)
	pprint(results)

	fig, ax = plt.subplots()
	ax.plot(results["dia_eigvalues"], ls="", marker="+")
	fig.show()

	input('Press any key to exit')
