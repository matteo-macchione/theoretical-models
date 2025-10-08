import numpy as np
import matplotlib.pyplot as plt

# === Parametri da modificare ===
filename = "delta_population_avg_gsl.dat"  
filenameTD = "delta_TD.dat"
filenameSB = "delta_SB.dat"
filename_complete = "delta_complete_population_gsl.dat"

# === Lettura del file ===
# Si assume che il file sia formattato come:
# x   y1   y2   y3
data = np.loadtxt(filename)
dataTD = np.loadtxt(filenameTD)
dataSB = np.loadtxt(filenameSB)
data_complete = np.loadtxt(filename_complete)

# Separazione delle colonne
x = data[:, 0]
y1 = data[:, 1]
y2 = data[:, 2]
y3 = data[:, 3]
TD_GR = dataTD[:, 1]
TD_fR = dataTD[:, 2]
TD_DGP = dataTD[:, 3]
SB_GR = dataSB[:, 1]
SB_fR = dataSB[:, 2]
SB_DGP = dataSB[:, 3]
GR_complete = data_complete[:, 1]
fR_complete = data_complete[:, 2]
DGP_complete = data_complete[:, 3]


# === Creazione del grafico ===
plt.figure(figsize=(8, 6))
plt.plot(x, np.zeros_like(x), color='black', linestyle='--')  # Linea orizzontale a y=0
plt.plot(x, y1, label='GR', color='blue', linestyle='-.')
plt.plot(x, y2, label='f(R)', color='red', linestyle='-.')
plt.plot(x, y3, label='DGP', color='green', linestyle='-.')
plt.plot(x, TD_GR, label='TD GR', color='blue', linestyle='--')
plt.plot(x, TD_fR, label='TD f(R)', color='red', linestyle='--')
plt.plot(x, TD_DGP, label='TD DGP', color='green', linestyle='--')
plt.plot(x, SB_GR, label='SB GR', color='blue', linestyle=':')
plt.plot(x, SB_fR, label='SB f(R)', color='red', linestyle=':')
plt.plot(x, SB_DGP, label='SB DGP', color='green', linestyle=':')
plt.plot(x, GR_complete, label='GR combined', color='blue')
plt.plot(x, fR_complete, label='f(R) combined', color='red')
plt.plot(x, DGP_complete, label='DGP combined', color='green')

# === Formattazione ===
plt.xlabel(r'$r_{\perp}/r_{500}$')
plt.ylabel(r'$\hat{\Delta}_{gz}$  [km/s]')
plt.title('')
plt.xlim(0, 4)
plt.ylim(-30, 7)
plt.legend()
plt.grid(False)
plt.tight_layout()

# === Visualizzazione ===
plt.show()
