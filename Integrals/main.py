import numpy as np

a = 1.0
b = 2.0
eps = 1e-6


def f(x):
    return np.exp(-x) + x**2


# вычисляем интеграл на заданном числе отрезков
def integral(num_segments):
    h = (b - a) / num_segments
    res = 0.0
    res = sum(((h / 2.0) * (f(a + i * h) + f(a + i * h + h))) for i in range(num_segments))
    print(res)
    return res


num_segments = 2
prev = integral(1)
curr = integral(2)
while abs(curr - prev) > eps:
    num_segments *= 2
    prev, curr = curr, integral(num_segments)