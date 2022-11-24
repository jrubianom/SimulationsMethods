import numpy as np
import matplotlib.pyplot as plt

R, C2 = np.loadtxt("dragCoeff_vs_reynosldNumber.dat", unpack=True)

plt.plot(R, C2)
plt.xlabel("Reynolds number")
plt.ylabel("Drag coefficient")
plt.grid()
plt.show()
