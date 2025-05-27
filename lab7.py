#!/usr/bin/env python3

# Вариант 26

import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import matplotlib

"""
x_a = (x_0 + x_n) / 2 = 3.0
x_g = sqrt(x_0 * x_n) = sqrt(5.0) = 2.23606797749978969640917
x_h = 2 / (1.0 + 1/5.0) = 1.6666666666666666666

y_a = 1.77
y_g = sqrt(3.14 * 0.40) = 1.1207140580897519909542995043432
y_h = 0.709604519774011299435028248587570621468926553672316

x_0 = 1.0
x_n = 5.0

y_0 = 3.14
y_n = 0.40

z(x_a) = z(3.0) = 1.03
z(x_g) = z(2.23) = 1.40
z(x_h) = z(1.666) = 1.75

d1 = |z(x_a) - y_a| = |1.03 - 1.77| = 0.74
d2 = |1.40 - 1.120714| = 0.279286
d3 = |1.03 - 1.120714| = 0.090714
d4 = |1.40 - 1.77| = 0.37
d5 = |1.75 - 1.77| = 0.02
d6 = |1.03 - 0.70| = 0.33
d7 = |1.75 - 0.7096045| = 1.04
d8 = |1.75 - 1.1207140| = 0.62
d9 = |1.40 - 0.7096045| = 0.69
"""

def fit(x, y):
    x = 1.0 / x
    sum_u_u = np.sum(x * x)
    sum_u  = np.sum(x)
    sum_y_u = np.sum(y * x)
    sum_y  = np.sum(y)
    A = np.array([[sum_u_u, sum_u],
                  [sum_u,  len(x)]], dtype=float)
    B = np.array([sum_y_u, sum_y], dtype=float)
    a, b = np.linalg.solve(A, B)
    return a, b

def predict(x, a, b):
    return a / x + b

def rss(y_true, y_pred):
    return np.sqrt(np.sum((y_pred - y_true) ** 2)) / np.sqrt(len(y_true))

def main():
    x = np.array([1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0])
    y = np.array([3.14, 1.78, 1.62, 1.33, 1.03, 0.29, 0.36, 0.60, 0.40])

    a, b = fit(x, y)
    y_appr = predict(x, a, b)
    standard_derivation = rss(y, y_appr)

    print("Аппроксимация y = a/x + b :")
    print(f"  a = {a:.8f}")
    print(f"  b = {b:.8g}")
    print(f"  СКО = {standard_derivation:.8f}")

    xs = np.linspace(x.min(), x.max(), 300)
    ys = predict(xs, a, b)

    plt.figure(figsize=(6,4))
    plt.scatter(x, y, label="исходные данные", zorder=5)
    plt.plot(xs, ys, label=f"y = {a:.3f}/x + {b:.3f}", linewidth=2)
    plt.title("аппроксимация y = a/x + b")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    matplotlib.use('TkAgg')
    main()
