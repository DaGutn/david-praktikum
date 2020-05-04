import numpy as np
import matplotlib.pyplot as plt
from uncertainties import ufloat
from scipy.optimize import curve_fit

def f(x, a, b, c):
    return (x**2)*a + b*x + c

#Wichtige Konstanten
h_eV = 4.1357*10**(-15)
d = 201.4#*10**(-12)
tau = 90*10**(-6)
lichtgeschw = 2.9979*10**(8)


#########################################################################################
#Emissionsspektrum der Kupfer-Anode gegen den Bragg-Winkel (die Wellenlänge) aufgetragen
theta_Cu, impprosekCu = np.genfromtxt("EmissionCu.txt", delimiter = "", unpack = True)

plt.plot(theta_Cu, impprosekCu, "b", label = r"Emmisionsspektrum")
plt.xlabel(r"$\theta$ [°]")
plt.ylabel(r"$n$ [Imp/s]")
plt.legend(loc = "best")
plt.savefig("EmissionCu.png")
plt.clf()


#########################################################################################
#K_alpha und K_beta -Linienenergien bestimmen
E_alpha_wahr = 8038
E_beta_wahr = 8905

lam_alpha = 2 * d * np.sin(np.radians(22.5))
E_alpha = h_eV * lichtgeschw / lam_alpha

lam_beta = 2 * d * np.sin(np.radians(20.3))
E_beta = h_eV * lichtgeschw / lam_beta

a_alpha = (E_alpha_wahr - E_alpha) / E_alpha

a_beta = (E_beta_wahr - E_beta) / E_beta

#########################################################################################
#Bestimmung der Transmission als Funktion der Wellenlänge
theta_ohne, n_ohne = np.genfromtxt("ComptonOhne.txt", delimiter = "", unpack = True)
theta_Al, n_Al = np.genfromtxt("ComptonAl.txt", delimiter = "", unpack = True)

#Zählraten korregieren wegen der Totzeit
n_ohne_kor = n_ohne/(1 - tau*n_ohne)
n_Al_kor = n_Al/(1 - tau*n_Al)

#Transmission und mit Bragg die Wellenlänge berechnen und gegeneinander auftragen
T =  n_Al_kor / n_ohne_kor
lam = 2 * d * np.sin(np.radians(theta_Al))

params, cov_matrix = curve_fit(f, lam, T)
lam_reg = np.linspace(lam[0], lam[-1])
plt.plot(lam, T, "b", label = r"Transmissionskurve")
plt.plot(lam_reg, f(lam_reg, *params), "r", label = r"Regressionsparabel")

plt.xlabel(r"$\lambda$ [pm]")
plt.ylabel(r"$T$")
plt.legend(loc = "best")
plt.savefig("Transmissionskurve.png")

#a = ufloat(params[0], np.sqrt(cov_matrix[0, 0]))
#b = ufloat(params[1], np.sqrt(cov_matrix[1, 1]))
#c = ufloat(params[2], np.sqrt(cov_matrix[2, 2]))

a = params[0]
b = params[1]
c = params[2]

print(params)


########################################################################################
#Bestimmung der Compton-Wellenlänge
I_0 = 2731
I_1 = 1180
I_2 = 1024

T_1 = I_1 / I_0
T_2 = I_2 / I_0

print(T_1)
print(T_2)

lam_1 = (-b - (b**2 + 4*a*(T_1 - c))**(1/2)) / (2 * a)
lam_2 = (-b - (b**2 + 4*a*(T_2 - c))**(1/2)) / (2 * a)
lam_C = lam_2 - lam_1

print(lam_1)
print(lam_2)
print(lam_C)

a_lam = (2.4263 - lam_C) / lam_C

print(a_lam)








