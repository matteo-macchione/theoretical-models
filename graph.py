import numpy as np
import matplotlib.pyplot as plt

# === Parametri da modificare ===
filename = "delta_population_avg_gsl.dat"  # Nome del file .dat

# === Lettura del file ===
# Si assume che il file sia formattato come:
# x   y1   y2   y3
data = np.loadtxt(filename)

# Separazione delle colonne
x = data[:, 0]
y1 = data[:, 1]
y2 = data[:, 2]
y3 = data[:, 3]

# === Creazione del grafico ===
plt.figure(figsize=(8, 6))
plt.plot(x, y1, label='GR', color='blue')
plt.plot(x, y2, label='f(R)', color='red')
plt.plot(x, y3, label='DGP', color='green')

# === Formattazione ===
plt.xlabel(r'$r_{\perp}/r_{500}$')
plt.ylabel(r'$\hat{\Delta}_{gz}$  [km/s]')
plt.title('')
plt.legend()
plt.grid(False)
plt.tight_layout()

# === Visualizzazione ===
plt.show()
