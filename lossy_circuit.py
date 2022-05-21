import numpy as np
import matplotlib.pyplot as plt


def main():
    """"""
    data = np.loadtxt("lossy_data.csv", delimiter=",")
    frequency = data[0, :]
    current = data[1, :]*1e-3
    voltage = data[2, :]

    ax = plt.axes()
    ax.plot(frequency, current, linestyle="none", marker="o", color="r")
    #ax.set_title("Resonance and Voltage Curve for Q=5.05")
    ax.set_xlabel("Frequency (Hz)")
    ax.set_ylabel("Current (A)", color="r")
    ax.tick_params(axis='y', labelcolor="r")
    ax.set_ylim(0.0094, 0.013)
    ax2 = ax.twinx()
    ax2.plot(frequency, voltage, linestyle="none", marker="o", color="purple")
    ax2.set_ylabel("Capacitor Voltage (V)", color="purple")
    ax2.tick_params(axis='y', labelcolor="purple")
    plt.savefig("lossy.png", dpi=1200)


main()
