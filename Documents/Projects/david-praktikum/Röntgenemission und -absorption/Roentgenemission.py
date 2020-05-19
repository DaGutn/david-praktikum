import numpy as np
import matplotlib.pyplot as plt
from uncertainties import ufloat
from uncertainties.umath import *
from scipy.optimize import curve_fit
from scipy.signal import find_peaks

def f_alpha(x_alpha, a_alpha, b_alpha, c_alpha):
    return x_alpha**2 * a_alpha + x_alpha * b_alpha + c_alpha

def f_beta(x_beta, a_beta, b_beta, c_beta):
    return x_beta**2 * a_beta + x_beta * b_beta + c_beta

def f_sigma(x_sigma, m_sigma):
    return x_sigma * m_sigma


#Wichtige Werte
h_eV = 4.1357*10**(-15)
d = 201.4*10**(-12)
lichtgeschw = 2.9979*10**(8)
rydbergenegie = 13.6057
sommerfeldkonst = 7.297 * 10**(-3)


############################################################################################
#Überprüfung der Bragg-Bedingung
sollwinkel = 28
theta_bragg, n_bragg = np.genfromtxt("Bragg.dat", delimiter = "", unpack = True)

plt.plot(theta_bragg, n_bragg, "b.", label = r"Ereignisrate in Abhängigkeit des Winkels")
plt.xlabel(r"$\theta$ [°]")
plt.ylabel(r"$n$ [Imp/s]")
plt.legend(loc = "best")
plt.savefig("Bragg.png")
plt.clf()

peak = find_peaks(n_bragg, height = 200)
expwinkel = theta_bragg[peak[0]]
a_winkel = (sollwinkel - expwinkel) / expwinkel
print("Winkel bei dem n maximal ist (exp. Glanzwinkel 1. Ordnung):")
print(expwinkel)
print("Relative Abweichung vom Sollwinkel:")
print(a_winkel)


#############################################################################################
#Analyse der Emissionspektrums der Cu-Röntgenröhre
#Emissionskurve plotten
theta, n_Cu = np.genfromtxt("Emissionsspektrum.dat", delimiter = "", unpack = True)

plt.plot(theta, n_Cu, "b", label = r"Emissionsspektrum")
plt.xlabel(r"$\theta$ [°]")
plt.ylabel(r"$n$ [Imp/s]")
plt.legend(loc = "best")
plt.savefig("Emissionsspektrum.png")
plt.clf()

############################################################################
#K_alpha-Peak plotten und fitten
theta_alpha = theta[143:150]
n_alpha = n_Cu[143:150]
plt.plot(theta_alpha, n_alpha, "b.", label = r"$K_{\alpha}$-Linie")

theta_alpha_lin = np.linspace(theta_alpha[0] - 0.02, theta_alpha[-1] + 0.02)
params_alpha, cov_matrix_alpha = curve_fit(f_alpha, theta_alpha, n_alpha)
#a_alpha = ufloat(params_alpha[0], np.sqrt(cov_matrix_alpha[0, 0]))
#b_alpha = ufloat(params_alpha[1], np.sqrt(cov_matrix_alpha[1, 1]))
#c_alpha = ufloat(params_alpha[2], np.sqrt(cov_matrix_alpha[2, 2]))
a_alpha = params_alpha[0]
b_alpha = params_alpha[1]
c_alpha = params_alpha[2]

#peak_alpha = find_peaks(f_alpha(theta_alpha_lin, *params_alpha))
#print(peak_alpha)
#print(f_alpha(theta_alpha_lin[peak_alpha[0]], *params_alpha))

print("Winkel, Rate, Wellenlänge und Energie des Peaks der K_alpha-Linie:")
#theta_alpha_wert = -b_alpha/(2*a_alpha)
#theta_alpha_wert = ufloat(theta[145], 0.1)
theta_alpha_wert = theta[145]
n_alpha_wert = f_alpha(theta_alpha_wert, *params_alpha)
print(theta_alpha_wert)
print(n_alpha_wert)
lam_alpha = 2 * d * sin(np.pi/180 * (theta_alpha_wert))
E_alpha = h_eV * lichtgeschw / lam_alpha
print(E_alpha)
schnittpunkt_alpha_1 = (-b_alpha - (b_alpha**2 + 4*a_alpha*(n_alpha_wert/2 - c_alpha))**(1/2)) / (2 * a_alpha)
schnittpunkt_alpha_2 = (-b_alpha + (b_alpha**2 + 4*a_alpha*(n_alpha_wert/2 - c_alpha))**(1/2)) / (2 * a_alpha)
E_schnittpunkt_alpha_1 = (h_eV * lichtgeschw) / (2 * d * sin(np.pi/180 * (schnittpunkt_alpha_1)))
E_schnittpunkt_alpha_2 = (h_eV * lichtgeschw) / (2 * d * sin(np.pi/180 * (schnittpunkt_alpha_2)))
print("Halbwert des Alpha-Peaks:")
print(n_alpha_wert/2)
print("Schnittpunkte Alpha-Peak:")
print(schnittpunkt_alpha_2, E_schnittpunkt_alpha_2)
print(schnittpunkt_alpha_1, E_schnittpunkt_alpha_1)
print("Auflösungsvermögen des Alpha-Peaks:")
print(E_alpha / (E_schnittpunkt_alpha_2 - E_schnittpunkt_alpha_1))
#plt.hlines(n_alpha_wert.nominal_value/2, schnittpunkt_alpha_2.nominal_value, schnittpunkt_alpha_1.nominal_value, colors='g', linestyles='dashed', label = r"Halbwertsbreite")
plt.hlines(n_alpha_wert/2, schnittpunkt_alpha_2, schnittpunkt_alpha_1, colors='g', linestyles='dashed', label = r"Halbwertsbreite")

plt.plot(theta_alpha_lin, f_beta(theta_alpha_lin, *params_alpha), "r", label = r"Regressionsparabel")
plt.xlabel(r"$\theta$ [°]")
plt.ylabel(r"$n$ [Imp/s]")
plt.legend(loc = "best")
plt.savefig("K_alpha_Linie.png")
plt.clf()

############################################################################
#K_beta-Peak plotten und fitten
theta_beta = theta[120:127]
n_beta = n_Cu[120:127]
plt.plot(theta_beta, n_beta, "b.", label = r"$K_{\beta}$-Linie")

theta_beta_lin = np.linspace(theta_beta[0] - 0.01, theta_beta[-1] + 0.01)
params_beta, cov_matrix_beta = curve_fit(f_beta, theta_beta, n_beta)
#a_beta = ufloat(params_beta[0], np.sqrt(cov_matrix_beta[0, 0]))
#b_beta = ufloat(params_beta[1], np.sqrt(cov_matrix_beta[1, 1]))
#c_beta = ufloat(params_beta[2], np.sqrt(cov_matrix_beta[2, 2]))
a_beta = params_beta[0]
b_beta = params_beta[1]
c_beta = params_beta[2]

print("Winkel, Rate, Wellenlänge und Energie des Peaks der K_beta-Linie:")
#theta_beta_wert = -b_beta/(2*a_beta)
theta_beta_wert = theta[122]
n_beta_wert = f_beta(theta_beta_wert, *params_beta)
print(theta_beta_wert)
print(n_beta_wert)
lam_beta = 2 * d * sin(np.pi/180 * (theta_beta_wert))
E_beta = h_eV * lichtgeschw / lam_beta
print(E_beta)
schnittpunkt_beta_1 = (-b_beta - (b_beta**2 + 4*a_beta*(n_beta_wert/2 - c_beta))**(1/2)) / (2 * a_beta)
schnittpunkt_beta_2 = (-b_beta + (b_beta**2 + 4*a_beta*(n_beta_wert/2 - c_beta))**(1/2)) / (2 * a_beta)
E_schnittpunkt_beta_1 = (h_eV * lichtgeschw) / (2 * d * sin(np.pi/180 * (schnittpunkt_beta_1)))
E_schnittpunkt_beta_2 = (h_eV * lichtgeschw) / (2 * d * sin(np.pi/180 * (schnittpunkt_beta_2)))
print("Halbwert des Beta-Peaks:")
print(n_beta_wert/2)
print("Schnittpunkte und zugehörige Energien Beta-Peak:")
print(schnittpunkt_beta_2, E_schnittpunkt_beta_2)
print(schnittpunkt_beta_1, E_schnittpunkt_beta_1)
print("Auflösungsvermögen des Beta-Peaks:")
print(E_beta / (E_schnittpunkt_beta_2 - E_schnittpunkt_beta_1))
plt.hlines(n_beta_wert/2, schnittpunkt_beta_2, schnittpunkt_beta_1, colors='g', linestyles='dashed', label = r"Halbwertsbreite")

plt.plot(theta_beta_lin, f_beta(theta_beta_lin, *params_beta), "r", label = r"Regressionsparabel")
plt.xlabel(r"$\theta$ [°]")
plt.ylabel(r"$n$ [Imp/s]")
plt.legend(loc = "best")
plt.savefig("K_beta_Linie.png")
plt.clf()

#########################################################################
#Absorptionskonstanten der K,L und M-Schalen berechnen
bindungsenergie_Cu = 8980.476
ordnungszahl = 29
absorptionskonst_1 = ordnungszahl - np.sqrt((bindungsenergie_Cu) / (rydbergenegie))
absorptionskonst_2 = ordnungszahl - 2 * np.sqrt((bindungsenergie_Cu - E_alpha) / (rydbergenegie))
absorptionskonst_3 = ordnungszahl - 3 * np.sqrt((bindungsenergie_Cu - E_beta) / (rydbergenegie))
print("Die Absorptionskonstanten von Kupfer:")
print(absorptionskonst_1, absorptionskonst_2, absorptionskonst_3)


#################################################################################################################
#Absorptionskonstanten
theta_zink, n_zink = np.genfromtxt("Zink.dat", delimiter = "", unpack = True)
theta_brom, n_brom = np.genfromtxt("Brom.dat", delimiter = "", unpack = True)
theta_strontium, n_strontium = np.genfromtxt("Strontium.dat", delimiter = "", unpack = True)
theta_rubidium, n_rubidium = np.genfromtxt("Rubidium.dat", delimiter = "", unpack = True)
theta_zirkonium, n_zirkonium = np.genfromtxt("Zirkonium.dat", delimiter = "", unpack = True)
theta_gallium, n_gallium = np.genfromtxt("Gallium.dat", delimiter = "", unpack = True)

Z = np.array([30, 31, 35, 37, 38, 40])
theta_K = np.array([18.67, 17.35, 13.23, 11.79, 11.79, 9.96])
E_K_wahr = np.array([9660, 10370, 13470, 15200, 16110, 18000])
E_K = h_eV * lichtgeschw / (2 * d * np.sin(np.pi/180 * theta_K))
a_E_K = (E_K_wahr - E_K) / E_K
sigma_K = Z - np.sqrt(E_K / rydbergenegie - (sommerfeldkonst**2 * Z**4) / 4)
print(sigma_K)
print(E_K / rydbergenegie)
print(np.sqrt(E_K / rydbergenegie - (sommerfeldkonst**2 * Z**4) / 4))

plt.plot(theta_zink, n_zink, "b", label = r"Zink-Absorptionskante")
plt.xlabel(r"$\theta$ [°]")
plt.ylabel(r"$n$ [Imp/s]")
plt.legend(loc = "best")
plt.savefig("Zink.png")
plt.clf()
zink_I = (55.0 + 102)/2
print("Intensität, Winkel, Energie und ihre Abweichung vom Literaturwert der K-Kante von Zink:")
print(zink_I)
print(theta_K[0])
print(E_K[0])
print(a_E_K[0])

plt.plot(theta_gallium, n_gallium, "b", label = r"Gallium-Absorptionskante")
plt.xlabel(r"$\theta$ [°]")
plt.ylabel(r"$n$ [Imp/s]")
plt.savefig("Gallium.png")
plt.legend()
plt.clf()
gallium_I = (66 + 121)/2
print("Intensität, Winkel, Energie und ihre Abweichung vom Literaturwert der K-Kante von Gallium:")
print(gallium_I)
print(theta_K[1])
print(E_K[1])
print(a_E_K[1])

plt.plot(theta_brom, n_brom, "b", label = r"Brom-Absorptionskante")
plt.xlabel(r"$\theta$ [°]")
plt.ylabel(r"$n$ [Imp/s]")
plt.legend(loc = "best")
plt.savefig("Brom.png")
plt.clf()
brom_I = (9 + 27)/2
print("Intensität, Winkel, Energie und ihre Abweichung vom Literaturwert der K-Kante von Brom:")
print(brom_I)
print(theta_K[2])
print(E_K[2])
print(a_E_K[2])

plt.plot(theta_rubidium, n_rubidium, "b", label = r"Rubidium-Absorptionskante")
plt.xlabel(r"$\theta$ [°]")
plt.ylabel(r"$n$ [Imp/s]")
plt.legend(loc = "best")
plt.savefig("Rubidium.png")
plt.clf()
rubidium_I = (10 + 64)/2
print("Intensität, Winkel, Energie und ihre Abweichung vom Literaturwert der K-Kante von Rubidium:")
print(rubidium_I)
print(theta_K[3])
print(E_K[3])
print(a_E_K[3])

plt.plot(theta_strontium, n_strontium, "b", label = r"Strontium-Absorptionskante")
plt.xlabel(r"$\theta$ [°]")
plt.ylabel(r"$n$ [Imp/s]")
plt.legend(loc = "best")
plt.savefig("Strontium.png")
plt.clf()
strontium_I = (50 + 181)/2
print("Intensität, Winkel, Energie und ihre Abweichung vom Literaturwert der K-Kante von Strontium:")
print(strontium_I)
print(theta_K[4])
print(E_K[4])
print(a_E_K[4])

plt.plot(theta_zirkonium, n_zirkonium, "b", label = r"Zirkonium-Absorptionskante")
plt.xlabel(r"$\theta$ [°]")
plt.ylabel(r"$n$ [Imp/s]")
plt.legend(loc = "best")
plt.savefig("Zirkonium.png")
plt.clf()
zirkonium_I = (126 + 282)/2
print("Intensität, Winkel, Energie und ihre Abweichung vom Literaturwert der K-Kante von Zirkonium:")
print(zirkonium_I)
print(theta_K[5])
print(E_K[5])
print(a_E_K[5])

#Rydbergenergie bestimmen
zhochzwei = (Z - sigma_K)**2
print(zhochzwei)
plt.plot(zhochzwei, E_K, "bx", label = r"Messwerte")
sigma_K_lin = np.linspace(zhochzwei[0], zhochzwei[-1])
params_sigma, cov_matrix_sigma = curve_fit(f_sigma, zhochzwei, E_K)
plt.plot(sigma_K_lin, f_sigma(sigma_K_lin, *params_sigma), "r", label = r"Lineare Regression")
plt.xlabel(r"$Z_{eff}$")
plt.ylabel(r"$E_K$ [eV]")
plt.legend(loc = "best")
plt.savefig("Rydbergenergie.png")
plt.clf()
m = ufloat(params_sigma, cov_matrix_sigma[0, 0])
print(m)
a_rydbergenergie = (rydbergenegie - m.nominal_value) / m.nominal_value
print(a_rydbergenergie)





















































