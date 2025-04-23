#!/usr/bin/env python3

# Вариант 26

import math


a, b = 0.0, 1.0
eps = 0.001
y0 = 1.0
dy0 = 5.0
n = 100
h = (b-a) / n


def dydx(x, y1, y2):
    dy1 = y2
    dy2 = 8 * y2 - 7 * y1 + 14
    return dy1, dy2

def runge_kutta_step(x, y1, y2, h):
    k1y1, k1y2 = dydx(x, y1, y2)
    k2y1, k2y2 = dydx(x + h/2, y1 + h/2 * k1y1, y2 + h/2 * k1y2)
    k3y1, k3y2 = dydx(x + h/2, y1 + h/2 * k2y1, y2 + h/2 * k2y2)
    k4y1, k4y2 = dydx(x + h, y1 + h * k3y1, y2 + h * k3y2)

    next_y1 = y1 + h * (k1y1 + 2*k2y1 + 2*k3y1 + k4y1) / 6
    next_y2 = y2 + h * (k1y2 + 2*k2y2 + 2*k3y2 + k4y2) / 6
    return next_y1, next_y2

def exact_solution(x):
    C1 = 1
    C2 = -2
    return C1 * math.exp(7*x) + C2 * math.exp(x) + 2

x_values = list()
y_num = list()
dy_num = list()
y_exact = list()

x = a
y1 = y0
y2 = dy0

while x <= b:
    x_values.append(x)
    y_num.append(y1)
    y_exact.append(exact_solution(x))
    dy_num.append(y2)

    y1, y2 = runge_kutta_step(x, y1, y2, h)
    x += h

print('|       x         |       y\'      |       y\'\'     |       y\'точное |       delta    |')
for i in range(len(x_values)):
    print('|{:14.7f}   |{:14.7f} |{:14.7f} |{:14.7f}  |{:14.7f}  |'.format(x_values[i], y_num[i], dy_num[i], y_exact[i], abs(y_exact[i] - y_num[i])))
