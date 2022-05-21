import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


def model(capacitance, gamma, frequency, current_max):
    """"""
    inductance = 0.1522
    omega = frequency*2*np.pi
    numerator = gamma * omega * current_max
    denominator = np.sqrt(omega ** 2 * gamma ** 2 +
                          (omega ** 2 - 1/(inductance*capacitance))**2)
    return numerator / denominator


def main():
    """"""
    data = np.loadtxt("data.csv", delimiter=",", skiprows=1)
    capacitance = data[-21:, 2]
    current = data[-21:, 1]
    current_sigma = [0.1e-3]*21
    print(capacitance, current)
    guessparams = [160, 550, 0.09]
    cor_parameters, pcov = curve_fit(
        model, capacitance, current, p0=guessparams, sigma=current_sigma, absolute_sigma=True)

    ax = plt.axes()
    ax.plot(capacitance, current, marker="o", linestyle="none")
    xs = np.linspace(0.45e-6, 0.65e-6, 1000)
    ax.plot(xs, model(xs, *cor_parameters), color="r")
    ax.set_xlabel("capacitance (F)")
    ax.set_ylabel("current (A)")
    plt.savefig("capacitance_curve.png", dpi=1200)
    cor_parameter_err = np.sqrt(np.diag(pcov))
    print(cor_parameters, cor_parameter_err)


main()
