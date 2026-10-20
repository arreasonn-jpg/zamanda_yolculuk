"""Read the 4x sub4 npz checkpoint: how many steps have landed."""
import numpy as np

d = np.load("results/validated/v12_10_amr_v16f4sub4.npz")
t = d["t"]
hf = d["H_fine"]
hc = d["H_coarse"]
print("points:", len(t), " last step:", int(round(t[-1] / 2e-3)))
print("H_fine last 5:", np.round(hf[-5:], 4).tolist())
print("H_c  last 5:", np.round(hc[-5:], 4).tolist())
