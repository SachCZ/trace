import numpy as np


def density(x):
    return 0.08 * (1 - x ** 2)


def ionization(x):
    return 10 * np.ones(len(x))


def temperature(x):
    return 400 * np.ones(len(x))


def main():
    x = np.linspace(-1, 0, 50)
    np.savetxt("position.csv", x)
    np.savetxt("density.csv", density(x))
    np.savetxt("ionization.csv", ionization(x))
    np.savetxt("temperature.csv", temperature(x))


if __name__ == '__main__':
    main()
