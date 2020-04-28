import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline


def my_func(x):
    return np.exp(-x) + x**2


def draw_cubic_spline(x, y):
    cs = CubicSpline(x, y, bc_type='natural')
    x_values = np.linspace(0.0, 1.0, 100)
    plt.plot(x_values, my_func(x_values), 'b-', label='function')
    plt.plot(x_values, cs(x_values), 'r--', label='cubic spline')
    plt.grid(True)
    plt.title('n = {}'.format(n))
    plt.legend()
    plt.show()


for n in [3, 5, 10, 20]:
    # получаем узлы
    nodes = np.linspace(0.0, 1.0, n)
    values = my_func(nodes)
    # строим кубический сплайн по n узлам
    draw_cubic_spline(nodes, values)