import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from uncertainties import ufloat
from uncertainties.unumpy import (nominal_values as noms,
                                  std_devs as stds)


def f(x, m, b):
    return m*x + b

#Daten einlesen
linie, spannung = np.genfromtxt("data.csv",delimiter=",",unpack=True)

D = (linie - 1) * 6
print(D)

params, covariance_matrix = curve_fit(f, spannung, D)

x_werte = np.linspace(spannung[0], spannung[-1])

plt.plot(spannung, D, "rx", label = r"Messwerte")
plt.plot(x_werte, f(x_werte, *params), "b-", label = r"Lineare Regression")
plt.xlabel(r"U / V")
plt.ylabel(r"D")
plt.legend(loc = "best")
plt.savefig("Regressiongraph.pdf")

m = ufloat(params[0], np.sqrt(covariance_matrix[0, 0]))
b = ufloat(params[1], np.sqrt(covariance_matrix[1, 1]))

print(m)
print(b)


R_innen = ufloat(10, 1)
R_außen = ufloat(15, 1)
h = ufloat(20, 1)

V = np.pi * h * (R_außen**2 - R_innen**2)
print(V)
print(stds(V) / noms(V))