#!/usr/bin/env python3

# Вариант 26

import math


a, b = 0.0, 1.0
eps = 0.001
y0 = 1.0
dy0 = 5.0
order = 2**4 - 1


def exact_solution(x):
    C1 = 1
    C2 = -2
    return C1 * math.exp(7*x) + C2 * math.exp(x) + 2


def f(y1, y2):
    return y2, 8 * y2 - 7 * y1 + 14


def runge_kutta_step(y1, y2, h):
    k1y1, k1y2 = f(y1, y2)
    k2y1, k2y2 = f(y1 + h / 2 * k1y1, y2 + h / 2 * k1y2)
    k3y1, k3y2 = f(y1 + h / 2 * k2y1, y2 + h / 2 * k2y2)
    k4y1, k4y2 = f(y1 + h * k3y1, y2 + h * k3y2)

    ny = y1 + h * (k1y1 + 2 * k2y1 + 2 * k3y1 + k4y1) / 6
    ny2 = y2 + h * (k1y2 + 2 * k2y2 + 2 * k3y2 + k4y2) / 6

    return ny, ny2

def calc():
    results = []
    h = 0.5

    for i in range(1000):
        max_error = 0
        x = a
        y1 = y0
        y2 = dy0
        current_results = []

        exact = exact_solution(x)
        current_results.append((x, y1, 0.0))

        while x < b:
            y1_h, y2_h = runge_kutta_step(y1, y2, h)

            y1_mid, y2_mid = runge_kutta_step(y1, y2, h/2)
            y1_2h, y2_2h = runge_kutta_step(y1_mid, y2_mid, h/2)

            error = abs(y1_h - y1_2h) / order
            x_new = x + h

            if x_new > b:
                break

            current_results.append((x_new, y1_2h, error))
            max_error = max(max_error, error)

            y1, y2 = y1_2h, y2_2h
            x = x_new

        results.append((h, current_results))

        if max_error < eps:
            print('|       x         |       y\'      |       error   |       y\'точное |       delta    |')
            for entry in current_results:
                x_val, approx, err = entry

                ex = exact_solution(x_val)

                print('|{:14.7f}   |{:14.7f} |{:14.7f} |{:14.7f}  |{:14.7f}  |'.format(x_val, approx, err, ex, abs(approx - ex)))
            print(f'ε = {eps}, h = {h}')
            break

        h /= 2

calc()
