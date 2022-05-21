import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


def model(frequency, gamma, current_max, omega_0):
    omega = 2*np.pi * frequency
    omega_0 = omega_0 * 2*np.pi
    numerator = gamma * omega * current_max
    denominator = np.sqrt(gamma**2 * omega**2 + (omega**2 - omega_0**2)**2)
    return numerator/denominator


def main():
    """"""
    data = np.loadtxt("data.csv", delimiter=",", skiprows=1)
    frequency = data[-42:-21, 0]
    current = data[-42:-21, 1]
    current_sigma = [0.1e-3]*21
    print(frequency, current)
    guessparams = [100, 0.09, 550]
    cor_parameters, pcov = curve_fit(
        model, frequency, current, p0=guessparams, sigma=current_sigma, absolute_sigma=True)
    (gamma, current_0, f_0) = cor_parameters
    ax = plt.axes()
    ax.plot(frequency, current, marker="o", linestyle="none")
    xs = np.linspace(500, 600, 1000)
    ax.plot(xs, model(xs, *cor_parameters), color="r")
    ax.set_xlabel("frequency (Hz)")
    ax.set_ylabel("current (A)")
    ax.hlines(current_0/np.sqrt(2), f_0-gamma/(4*np.pi),
              f_0+gamma/(4*np.pi), color="g", label="γ (Bandwidth)")
    ax.legend(loc="upper right")
    plt.savefig("frequency_curve.png", dpi=1200)
    cor_parameter_err = np.sqrt(np.diag(pcov))
    print(cor_parameters, cor_parameter_err)


main()
