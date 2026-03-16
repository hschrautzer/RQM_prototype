import numpy as np
import pandas as pd
from pathlib import Path
from prototype.magnetization import Magnetization
from prototype.lbfgs import lbfgs_minimizer

if __name__=="__main__":
    # Load the ferromagnet configuration
    df = pd.read_csv(Path(__file__).parent / "spin_i.dat",sep=r"\s+")

    # Initialize the magnetization object.
    points = np.column_stack((df["x"].to_numpy(),df["y"].to_numpy(),df["z"].to_numpy()))
    spins = np.column_stack((df["sx"].to_numpy(),df["sy"].to_numpy(),df["sz"].to_numpy()))

    mag = Magnetization(points=points,spins=spins)

    minimizer = lbfgs_minimizer(N_memory=3, theta_max=0.5, N_spins=mag.N, N_modes=4, N_iter=1, rq_grad_tol=1e-12)

    # Initial vector to be random and in 3N space:
    vec_ini_3N = np.random.rand(3*mag.N, 4)
    vec_final = minimizer.minimize(mag=mag, vec_ini=vec_ini_3N)
    print(vec_final)

