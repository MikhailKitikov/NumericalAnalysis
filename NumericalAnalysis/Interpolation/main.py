import numpy as np
import matplotlib.pyplot as plt
from pylab import *

# значения функции (для сравнения с полиномом)
func_x = np.linspace(0.0, 1.0, 100)
func_y = [exp(-x) + x**2 for x in func_x]


# вычисление интерполяционного полинома
def lagrange(x, y, val):
    res = 0
    for j in range(len(y)):
        num = 1; denom = 1
        for i in range(len(x)):
            if i == j:
                continue
            num *= val - x[i]
            denom *= x[j] - x[i]
        res += y[j] * num / denom
    return res


# интерполяция функции по узлам Чебышева
def interpolate(a, b, nodes_cnt):
    nodes_x = np.asarray([0.5 * (a + b) + 0.5 * (b - a) * cos(pi * (2 * k - 1) / (2 * nodes_cnt))
                       for k in range(1, nodes_cnt + 1)], dtype = float)
    nodes_y = np.asarray([(np.exp(-x) + x**2) for x in nodes_x], dtype = float)
    x = np.linspace(np.min(nodes_x), np.max(nodes_x), 100)
    y = [lagrange(nodes_x, nodes_y, x_val) for x_val in x]
    return x, y


# построение графика
def draw_plot(x, y, n):
    fig = plt.figure()
    ax = plt.subplot(111)
    ax.plot(func_x, func_y, 'r-',  label = 'f(x)')
    ax.plot(x, y, 'b--', label = 'polynom')
    plt.grid(True)
    plt.title('n = {}'.format(n))
    ax.legend()
    show()


num_segments = [1, 20, 40]

for n in num_segments:
    # разделяем отрезок на n частей
    segments = np.append(np.arange(0.0, 1.0, 1.0 / n), 1.0)
    x = []; y = []

    # на каждом подотрезке интерполируем функцию
    for i in range(1, size(segments)):
        tX, tY = interpolate(segments[i], segments[i - 1], 3)
        x.extend(tX); y.extend(tY)

    # строим график
    draw_plot(x, y, n)