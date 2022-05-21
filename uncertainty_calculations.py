import numpy as np


def partial_f(f, C):
    return 2/(4 * np.pi**2 * f**3 * C)


def partial_C(f, C):
    return 1/(4 * np.pi**2 * f**2 * C**2)


def sigma_L(f, C, sigma_f, sigma_C):
    return np.sqrt(partial_f(f, C)**2 * sigma_f ** 2 + partial_C(f, C) ** 2 * sigma_C ** 2)


def main():
    """"""
    f = 550
    C = 0.55e-6
    sigma_f = 1
    sigma_C = 0.01e-6
    print(sigma_L(f, C, sigma_f, sigma_C))


main()
