from matplotlib import pyplot as plt
import numpy as np


def f(x, u):
    return x**2 * u - x + 2


n = 10
h = 1 / n

y = [0] * (n + 1)
x = list(np.linspace(0, 1.001, n + 1))

for i in range(1, n + 1):
    phi0 = h * f(x[i - 1], y[i - 1])
    phi1 = h * f(x[i - 1] + h / 2, y[i - 1] + phi0 / 2)
    phi2 = h * f(x[i - 1] + h / 2, y[i - 1] + phi1 / 2)
    phi3 = h * f(x[i - 1] + h, y[i - 1] + phi2)
    y[i] = y[i - 1] + 1 / 6 * (phi0 + 2 * phi1 + 2 * phi2 + phi3)


plt.plot(x, y)
plt.show()