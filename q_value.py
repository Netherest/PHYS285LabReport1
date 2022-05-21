import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


def lorentzCurve(frequency, gamma, current_max, f_0):
    """model for current"""
    omega = frequency * 2 * np.pi
    omega_0 = f_0 * 2 * np.pi
    numerator = current_max * omega * gamma
    denominator = np.sqrt(omega**2 * gamma**2 + (omega**2 - omega_0**2)**2)
    return numerator/denominator


RESISTANCES = [0, 70, 80, 105]
GUESSPARAMS = [300, 0.09, 550]
CURRENT_SIGMA = [1e-3]*21
CAPACITANCE = 0.55e-6


def main():
    """main function"""
    data = np.loadtxt("data.csv", delimiter=",", skiprows=1)
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)
    fig.set_size_inches(10, 8)
    q_value = list()
    constants = list()
    for resistance in RESISTANCES:
        restriction = np.append(np.logical_and(data[:-21, 3] == resistance,
                                               data[:-21, 2] == CAPACITANCE), [False]*21)
        # print(restriction)
        #print(data[restriction, 0]*2*np.pi, data[restriction, 1])
        cor_parameters, pcov = curve_fit(lorentzCurve,
                                         data[restriction, 0],
                                         data[restriction, 1],
                                         p0=GUESSPARAMS, sigma=CURRENT_SIGMA,
                                         absolute_sigma=True)
        cor_parameter_err = np.sqrt(np.diag(pcov))
        print(cor_parameter_err)
        (gamma, current_max, f_0) = cor_parameters
        print(f"gamma: {gamma}\ncurrent_max:{current_max}\nf_0:{f_0}\n\n")
        q_value.append(f_0*2*np.pi/gamma)
        constants.append(cor_parameters)
    # print(constants)
    # print(constants[1])
    xs = np.linspace(500, 600, 1000)
    plt.subplots_adjust(hspace=0.3)
    ax1.set_ylim(0, 0.1)
    ax1.set_ylabel("current (A)")
    ax1.set_xlabel("frequency (Hz)")
    ax2.set_ylim(0, 0.1)
    ax2.set_ylabel("current (A)")
    ax2.set_xlabel("frequency (Hz)")
    ax3.set_ylim(0, 0.1)
    ax3.set_ylabel("current (A)")
    ax3.set_xlabel("frequency (Hz)")
    ax4.set_ylim(0, 0.1)
    ax4.set_ylabel("current (A)")
    ax4.set_xlabel("frequency (Hz)")
    ax1.set_title(r'Q={0:.2f}, $+0Ω$'.format(q_value[0]))
    ax1.plot(xs, lorentzCurve(xs, *constants[0]))
    ax2.set_title(r'Q={0:.2f}, $+70Ω$'.format(q_value[1]))
    ax2.plot(xs, lorentzCurve(xs, *constants[1]))
    ax3.set_title(r'Q={0:.2f}, $+80Ω$'.format(q_value[2]))
    ax3.plot(xs, lorentzCurve(xs, *constants[2]))
    ax4.set_title(r'Q={0:.2f}, $+105Ω$'.format(q_value[3]))
    ax4.plot(xs, lorentzCurve(xs, *constants[3]))
    plt.savefig("qvalues.png", dpi=1200)


main()
