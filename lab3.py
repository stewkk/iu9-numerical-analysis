#!/usr/bin/env python3

# Вариант 26

import math

EPS = 0.001

# def f(x):
#     return math.sin(x) * math.cos(x/2)

def f(x):
    return math.exp(x)

def rectangles(h, n, a, b):
    # n = int((b-a)/h)
    res = 0
    for i in range(0, n):
        x = a + (h*i)
        res += f(x + h/2)
    res *= h
    return res

def simpson(h, n, a, b):
    # n = int((b-a)/h)
    res = f(a) + f(b)
    for i in range(1, n):
        x = a + (h*i)
        if i%2 == 0:
            res += 2*f(x)
        else:
            res += 4*f(x)
    res *= h/3
    return res

def trapezoid(h, n, a, b):
    res = f(a)+f(b)
    res /= 2
    for i in range(1, n):
        res += f(a + i*h)
    return h * res

def calc(a, b, method, p):
    n = 4
    h = (b-a)/n

    iterations = 1
    ds = None
    while True:
        ds = abs(method(h, n, a, b) - method(h*2, n//2, a, b)) / ((2**p) - 1)
        if ds <= EPS:
            break
        iterations += 1
        h /= 2
        n *= 2

    return iterations, h, n, ds, method(h, n, a, b)


def main():
    a = 0
    # b = math.pi
    b = 1

    iters_rectangles, h_rect, n_rect, ds_rect, res = calc(a, b, rectangles, 2)
    print(res)
    iters_simpson, h_simpson, n_simpson, ds_simpson, res = calc(a, b, simpson, 4)
    print(res)
    iters_trapezoid, h_trapezoid, n_trapezoid, ds_trapezoid, res = calc(a, b, trapezoid, 2)
    print(res)
    I_rect = rectangles(h_rect, n_rect, a, b)
    I_simpson = simpson(h_simpson, n_simpson, a, b)
    I_trapezoid = trapezoid(h_trapezoid, n_trapezoid, a, b)

    print('|             | Прямоугольников | Симпсона  | Трапеций  |')
    print('| n           |{:10}       |{:>6}     |{:>9}  |'.format(n_rect, n_simpson, n_trapezoid))
    print('| I^*_h/2     |{:14.7f}   |{:10.7f} |{:10.7f} |'.format(I_rect, I_simpson, I_trapezoid))
    print('| R           |{:14.7f}   |{:10.7f} |{:10.7f} |'.format(ds_rect, ds_simpson, ds_trapezoid))
    print('| I^*_h/2 + R |{:14.7f}   |{:10.7f} |{:10.7f} |'.format(I_rect + ds_rect, I_simpson - ds_simpson, I_trapezoid - ds_trapezoid))


if __name__ == "__main__":
    main()
