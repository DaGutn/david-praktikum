import numpy as np
import matplotlib.pyplot as plt
import uncertainties.unumpy as unp
from uncertainties import ufloat
from uncertainties.umath import *
from uncertainties.unumpy import (nominal_values as noms, std_devs as stds)
from scipy.optimize import curve_fit
import scipy.constants as const

def f(x, a, b):
    return x*a + b

#Plotten und fitten des Plateaus mit Fehlerbalken
U, N = np.genfromtxt("Kennlinie.dat", delimiter = "", unpack = True)
N_delta = np.sqrt(N)

U_plateau = U[4:35]
N_plateau = N[4:35]

params, cov_matrix = curve_fit(f, U_plateau, N_plateau)
U_lin = np.linspace(U_plateau[0], U_plateau[-1])

plt.plot(U, N, "bx", label = r"Messwerte mit Fehlerbalken")
plt.errorbar(U, N, xerr=0, yerr=N_delta, fmt='b')
plt.plot(U_lin, f(U_lin, *params), "r", label = r"Lineare Regression")
plt.xlabel(r"$U$ [V]")
plt.ylabel(r"$N$ [Imp]")
plt.legend(loc = "best")
plt.savefig("Plateau.png")
plt.clf()

m = ufloat(params[0], cov_matrix[0, 0])
b = ufloat(params[1], cov_matrix[1, 1])
print("Steigung und y-Achsenabschnitt der linearen Regression:")
print("{:.2u}".format(m))
print("{:.4u}".format(b))

#Plateausteigung in % pro 100V
m_plateau = m * 100**2 / f(500, *params)
#m_plateau = (1 - (f(400, *params)) / (f(500, *params))) * 100
print("Steigung in % pro 100V:")
print("{:.2u}".format(m_plateau))


#####################################################################################
#Totzeit nach der Zwei-Quellen-Methode berechnen
n_1 = ufloat(96041, np.sqrt(96041)) / 120
n_2 = ufloat(76518, np.sqrt(76518)) / 120
n_12 = ufloat(158479, np.sqrt(158479)) / 120
print("{:.3u}".format(n_1))
print("{:.3u}".format(n_2))
print("{:.3u}".format(n_12))
T = (n_1 + n_2 - n_12) / (2 * n_1 * n_2)
print("{:.2u}".format(T))


######################################################################################
#Bestimmung der Ladung pro einfallendem Teilchen
N_0 = np.array([9837, 9995, 10264, 10151, 10184, 10253, 10493, 11547])
n = unp.uarray(N_0, np.sqrt(N_0)) / 60
U, I_0 = np.genfromtxt("Zaehlrohrstrom.dat", delimiter = "", unpack = True) #I in 10^(-6) / Mikrosekunden
I = unp.uarray(I_0 * 10**(-6), 0.05 * 10**(-6))

Z = I / (const.e * n)
print(Z)

plt.errorbar(noms(U), noms(Z), xerr=0, yerr=stds(Z), fmt='bx', label=r"Freigesetzte Ladungen" )
plt.xlabel(r"U [V]")
plt.ylabel(r"Z [e]")
plt.legend(loc = "best")
plt.savefig("Ladungen.png")
plt.clf()








































