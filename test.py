import numpy as np
import pandas as pd
from pathlib import Path
from prototype.magnetization import magnetization


df = pd.read_csv(Path.cwd() / "test_config.csv")


points = np.column_stack((df["x"].to_numpy(),df["y"].to_numpy(),df["z"].to_numpy()))
spins = np.column_stack((df["sx"].to_numpy(),df["sy"].to_numpy(),df["sz"].to_numpy()))

mag = magnetization(points=points,spins=spins)

print(mag.energy())
print(mag.gradient())
print(mag.basis())
print(mag.basis()[:,:,0])
print(mag.gradient_tspace_2N())