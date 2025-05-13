import numpy as np
import matplotlib
import matplotlib.pyplot as plt


A, B, C, D = 1, -1, 0, -5


def f(x):
    return A * x ** 3 + B * x ** 2 + C * x + D


def f_der(x):
    return 3 * x**2 + 2*x


def f_der_der(x):
    return 6*x - 2


intervals = [(1, 3)]
EPS = 0.001


def bisection_method(f, a, b):
    fa = f(a)
    iterations = 0
    while (b - a) / 2 > EPS:
        c = (a + b) / 2
        fc = f(c)
        if abs(fc) < EPS:
            return c, iterations
        if fa * fc < 0:
            b = c
        else:
            a = c
            fa = fc
        iterations += 1
    return (a + b) / 2, iterations


def newton_method(f, f_prime, x0):
    for iterations in range(0, 100):
        fx = f(x0)
        fpx = f_prime(x0)
        x1 = x0 - fx / fpx
        if abs(x1 - x0) < EPS:
            return x1, iterations
        x0 = x1
        iterations += 1
    return 0, -1


def main():
    print("Метод деления отрезка пополам:")
    for a, b in intervals:
        root, iters = bisection_method(f, a, b)
        print(f"[{a}, {b}]: x = {root:.5f}, итераций = {iters}")

    print()
    print("Метод Ньютона:")
    initial_guesses = [3]
    for guess in initial_guesses:
        root, iters = newton_method(f, f_der, guess)
        print(f"Начальное приближение {guess}: x = {root:.5f}, итераций = {iters}")

    x_vals = np.linspace(-2, 3, 1000)
    matplotlib.use('TkAgg')
    plt.figure(figsize=(10, 6))
    plt.plot(x_vals, f(x_vals), label='f(x)')
    plt.plot(x_vals, f_der(x_vals), '--', label="f\'(x)")
    plt.plot(x_vals, f_der_der(x_vals), ':', label="f\'\'(x)")
    plt.legend()
    plt.grid(True)
    plt.show()

if __name__ == "__main__":
    main()


