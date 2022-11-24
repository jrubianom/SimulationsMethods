import matplotlib.pyplot as plt
import numpy as np

U, Rn, Fl = np.loadtxt("Flift_Fdrag_vs_Ufan.dat", unpack=True)

plt.plot(U, Fl)
plt.xlabel("Velocidad [celdas/click]")
plt.ylabel("Fuerza de Magnus")
plt.show()
