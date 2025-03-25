#!/usr/bin/env python3


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


def main():

    d = [5.0 / 12.0, 6.0 / 12.0, 6.0 / 12.0, 5.0 / 12.0]
    b = [4.0 / 12.0, 4.0 / 12.0, 4.0 / 12.0, 4.0 / 12.0]
    c = [1.0 / 12.0, 1.0 / 12.0, 1.0 / 12.0]
    a = [1.0 / 12.0, 1.0 / 12.0, 1.0 / 12.0]

    alpha, beta = forward_pass(d, a, b, c)

    res = backward_pass(alpha, beta)

    print(res)


main()
