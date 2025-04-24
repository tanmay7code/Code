import numpy as np
from pymoo.indicators.hv import Hypervolume

# 🔹 Load your Pareto front from a CSV file (Update the filename!)
pareto_front = np.loadtxt("f:/PYTHON/your_pareto_front.csv", delimiter=",")

# 🔹 Define a reference point (Make sure it's slightly worse than all solutions)
reference_point = np.max(pareto_front, axis=0) + 0.1  # Shift max values slightly

# 🔹 Compute hypervolume
hv = Hypervolume(ref_point=reference_point)
hypervolume_value = hv(pareto_front)
print(f" Hypervolume: {hypervolume_value}")

