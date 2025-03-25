#!/usr/bin/env python3

# Вариант 26

import math


def f(x):
    return math.sin(x) * math.cos(x/2)


def forward_pass(d: list[float], a: list[float], b: list[float], c: list[float]) -> tuple[list[float], list[float]]:
    alpha, beta = list(), list()
    alpha.append(-c[0] / b[0])
    beta.append(d[0] / b[0])
    n = len(d)
    for i in range(1, n-1):
        tmp = b[i] + a[i-1]*alpha[i-1]
        alpha.append(-c[i] / tmp)
        beta.append((d[i] - a[i-1]*beta[i-1]) / tmp)
    beta.append(
        (d[n-1] - a[n-2]*beta[n-2]) / (b[n-1] + a[n-2]*alpha[n-2]))
    return (alpha, beta)


def backward_pass(alpha: list[float], beta: list[float]) -> list[float]:
    n = len(beta)
    x = list(range(n))
    x[n-1] = beta[n-1]
    for i in range(n-2, -1, -1):
        x[i] = alpha[i]*x[i+1] + beta[i]
    return x


def calc_table(n, a, h):
    x = list()
    y = list()
    for i in range(0, n+1):
        x.append(a+h*i)
        y.append(f(x[-1]))
    return x, y


def calc_free_coefficients(y, n, h):
    res = list()
    res.append(0)
    for i in range(1, n):
        res.append(3*(y[i+1] - 2 * y[i] + y[i-1]) / h / h)
    return res


def main():
    l = 0
    r = math.pi
    n = 32
    h = (r-l) / n

    x, y = calc_table(n, l, h)
    for i in range(len(x)):
        print('x: ', '{0:0.15f}'.format(x[i]), 'y: ', '{0:0.15f}'.format(y[i]))
    print()

    free_coeffs = calc_free_coefficients(y, n, h)

    b = [4.0 for _ in range(n)]
    c = [1.0 for _ in range(n-1)]
    a = [1.0 for _ in range(n-1)]

    alpha, beta = forward_pass(free_coeffs, a, b, c)

    c_res = backward_pass(alpha, beta)
    c_res.insert(0, 0)
    c_res.append(0)

    a_res = [y[i] for i in range(len(y))]

    b_res = [(y[i+1] - y[i]) / h - (h/3) * (c_res[i+1] + 2 * c_res[i]) for i in range(0, n)]

    d_res = [(c_res[i+1] - c_res[i]) / 3 / h for i in range(0, n-1)]
    d_res.append(-c_res[n] / 3 / h)

    print('a: ', a_res)
    print('b: ', b_res)
    print('c: ', c_res)
    print('d: ', d_res)
    print()

    y_approx = list()
    diff = list()

    for i in range(n):
        y_approx.append( a_res[i] + b_res[i] * (x[i] - x[i]) + c_res[i] * ((x[i] - x[i]) ** 2) + d_res[i] * (
                    (x[i] - x[i]) ** 3))
        diff.append(math.fabs(y_approx[i] - y[i]))

    for i in range(n):
        print('x: {0:0.7f} y: {1:0.7f} y*: {2:0.7f} d*: {3:0.7f}'.format(x[i], y[i], y_approx[i], diff[i]))
    print()

    x_star = [l+(i-1/2)*h for i in range(1, n+1)]
    y_star = [f(x) for x in x_star]

    y_approx = list()
    diff = list()

    for i in range(n):
        y_approx.append( a_res[i] + b_res[i] * (x_star[i] - x[i]) + c_res[i] * ((x_star[i] - x[i]) ** 2) + d_res[i] * (
                    (x_star[i] - x[i]) ** 3))
        diff.append(math.fabs(y_approx[i] - y_star[i]))

    for i in range(n):
        print('x*: {0:0.7f} y: {1:0.7f} y*: {2:0.7f} d*: {3:0.7f}'.format(x_star[i], y_star[i], y_approx[i], diff[i]))

main()
