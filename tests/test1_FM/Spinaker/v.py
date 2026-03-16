from spinterface.core.lattice.CVectorField import CVectorField
from spinterface.visualization.lattice.CVectorFieldVisualizerPyVista import CVectorFieldVisualizerPyVista
from pathlib import Path


evec_range = range(1,5)

for k in evec_range:
    EV = CVectorField(path_vectorfield=Path.cwd() / f"evec3N_rqm_{k}.dat")
    V_EV = CVectorFieldVisualizerPyVista(field=EV,topview=True,arrow_scaler=10)
    #V_EV.show()
    V_EV.write_to_file(outpath=Path.cwd() / f"evec3N_rqm_{k}.png")