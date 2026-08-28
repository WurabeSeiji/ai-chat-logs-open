#!/usr/bin/env python3
from pathlib import Path
import pandas as pd, numpy as np
import matplotlib.pyplot as plt

HERE=Path(__file__).resolve().parent
df=pd.read_csv(HERE/"N5_inflation_vs_ordering_timeseries.csv")

plt.figure(figsize=(9,5))
plt.semilogy(df["step"],np.maximum(df["H_perp"],1e-30),label="H_perp")
plt.semilogy(df["step"],np.maximum(df["four_group_error"],1e-12),label="four-group error")
plt.xlabel("step")
plt.ylabel("log scale")
plt.title("N=5: rapid decompactification vs slower geometric ordering")
plt.legend()
plt.tight_layout()
plt.savefig(HERE/"N5_inflation_vs_ordering.png",dpi=180)
