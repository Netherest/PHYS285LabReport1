from re import A
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.optimize import curve_fit


def modelSurface(data, voltage_0, base_resistance):
    inductance = 0.1522
    resistance = base_resistance + data[2]
    return voltage_0 / (np.sqrt(resistance ** 2 + (data[0]*2 * np.pi * inductance - 1/(data[0]*2 * np.pi * data[1]))**2))


def currentSurface(data, current_max, base_resistance, constant):
    inductance = 0.1522
    gamma = (base_resistance+data[2])/inductance

    numerator = (gamma * data[0] * 2 * np.pi * current_max) * \
        np.exp(-constant*(base_resistance+data[2]))
    denominator = np.sqrt((data[0] * 2 * np.pi) ** 2 * gamma ** 2 +
                          ((data[0] * 2 * np.pi) ** 2 - 1/(inductance * data[1]))**2)
    return numerator / denominator


def main():
    data = np.loadtxt("data.csv", delimiter=",", skiprows=1)
    data[:84, 3] = data[:84, 3] + 10
    #print(data[:, 3])
    guessparams = [0.09, 25]

    current_A_sigma = np.array([0.3e-3]*len(data[:, 1]))

    cor_parameters, pcov = curve_fit(modelSurface, np.array(
        [data[:, 0], data[:, 2], data[:, 3]]), data[:, 1], p0=guessparams, sigma=current_A_sigma, absolute_sigma=True)
    cor_parameter_err = np.sqrt(np.diag(pcov))

    print(cor_parameters, cor_parameter_err)

    fig = plt.figure()
    ax = fig.add_subplot(1, 2, 1, projection='3d')
    plt.subplots_adjust(wspace=0.3)

    # Make data.
    X = np.linspace(450.0, 650.0, 500)
    Y = np.linspace(0, 300.0, 500)
    C = np.array([0.55e-6]*len(X)*len(Y))
    X, Y = np.meshgrid(X, Y)
    input_values = np.array([np.ravel(X), C, np.ravel(Y)])
    zs = modelSurface(input_values, *cor_parameters)
    Z = zs.reshape(X.shape)
    ax.plot_surface(X, Y, Z, rstride=18, cstride=18,
                    cmap='plasma', edgecolor='none', alpha=0.5)
    ax.scatter3D((data[:-21, 0]), data[:-21, 3],
                 data[:-21, 1], depthshade=False, color="r")
    ax.view_init(elev=40, azim=140)
    ax.set_xlabel("frequency (Hz)")
    ax.set_ylabel("resistance (Ohms)")
    ax.set_zlabel("current (A)")

    bx = fig.add_subplot(1, 2, 2, projection='3d')
    C = np.linspace(0.45e-6, 0.65e-6, 500)
    X = np.linspace(450.0, 650.0, 500)
    Y = np.array([0]*len(C)*len(X))
    X, C = np.meshgrid(X, C)
    input_values = np.array([np.ravel(X), np.ravel(C), Y])
    zs = modelSurface(input_values, *cor_parameters)
    Z = zs.reshape(X.shape)
    bx.plot_surface(X, C, Z, cmap='plasma', rstride=36, cstride=36, alpha=0.5)
    bx.view_init(azim=-60, elev=45)
    bx.scatter3D(data[-42:, 0], data[-42:, 2], data[-42:, 1], depthshade=False)
    bx.set_xlabel("frequency (Hz)")
    bx.set_ylabel("capacitance (μF)")
    bx.set_zlabel("current (A)")
    bx.ticklabel_format(useOffset=False, style='plain')
    bx.set_yticks([0.45e-6, 0.50e-6, 0.55e-6, 0.60e-6, 0.65e-6],
                  [0.45, 0.50, 0.55, 0.60, 0.65])
    plt.savefig("plot.png", dpi=1200)


main()
