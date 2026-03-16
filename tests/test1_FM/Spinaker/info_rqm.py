from otterplot.model.cotterplotter import COtterPlotter
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
plt.rcParams['text.usetex'] = True
plt.rcParams["text.latex.preamble"] = r"\usepackage{amsmath}\usepackage{bm}"
from pathlib import Path

otter = COtterPlotter(plotmode="paperPRB",grid=(1,2))
ax = otter.axes[0]
ax_conv = otter.axes[1]
df = pd.read_csv("info_rqm.csv", sep=r"\s+")

ax.plot(df["iteration"],df["eig1"],marker="+",color="r",linestyle="None",label="lowest eval")
ax.plot(df["iteration"],df["eig2"],marker="+",color="b",linestyle="None",label="scnd. lowest eval")
ax.plot(df["iteration"],df["eig3"],marker="+",color="g",linestyle="None",label="third lowest eval")
ax.plot(df["iteration"],df["eig4"],marker="+",color="tab:pink",linestyle="None",label="fourth lowest eval")

ax_conv.semilogy(df["iteration"],df["norm_grad"],marker="+",color="k",linestyle="None")

otter.adjustlabels(ax_nr=0,xlabel="RQM iteration",ylabel="eigenvalues (eV)")
otter.adjustlabels(ax_nr=1,xlabel="RQM iteration",ylabel=r"$||\nabla R(X)||_F$")
plt.subplots_adjust(left=0.16,right=0.964,bottom=0.11,top=0.98,wspace=0.2,hspace=0.283)
plt.savefig("info_rqm.png",dpi=300)
plt.show()