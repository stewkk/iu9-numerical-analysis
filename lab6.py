#!/usr/bin/env python3

# Вариант 25

EPS = 0.001
X0 = 0.0
Y0 = 0.0


def f(x, y):
    return x + 2*y + 4*((1 + x**2 + y**2)**0.5)


def dfdx(x, y):
    return 4*x / ((x**2 + y**2 + 1)**0.5) + 1


def d2fdx2(x, y):
    return 4 / ((x**2 + y**2 + 1)**0.5) - 4*(x**2) / ((x**2 + y**2 + 1)**1.5)


def dfdy(x, y):
    return 4*y / ((x**2 + y**2 + 1)**0.5) + 2


def d2fdy2(x, y):
    return 4 / ((x**2 + y**2 + 1)**0.5) - 4*(y**2) / ((x**2 + y**2 + 1)**1.5)


def d2fdxdy(x, y):
    return -(4*x*y / ((x**2 + y**2 + 1)**1.5))


def main():
    x = X0
    y = Y0
    while True:
        df_dx = dfdx(x, y)
        df_dy = dfdy(x, y)

        if max(df_dx, df_dy) <= EPS:
            break

        sqr_df_dx = df_dx**2
        sqr_df_dy = df_dy**2

        dphi_dt = -sqr_df_dx - sqr_df_dy
        d2phi_dt2 = (d2fdx2(x, y) * sqr_df_dx +
                     2 * d2fdxdy(x, y) * df_dx * df_dy +
                     d2fdy2(x, y) * sqr_df_dy)
        t = -dphi_dt / d2phi_dt2

        x -= t * df_dx
        y -= t * df_dy

    f_min = f(x, y)
    print(f"({x}, {y}) = {f_min}")
    # min{x + 2 y + 4 (1 + x^2 + y^2)^0.5}≈3.316624790355399849114932737
    # at (x, y)≈(-0.3015113445777636226468120670, -0.6030226891555272452936241339)
    analytical = 3.316624790355399849114932737
    x_analytical = -0.3015113445777636226468120670
    y_analytical = -0.6030226891555272452936241339
    print(f"analytical f({x_analytical, y_analytical}): {analytical}")
    print(f"diff f(x, y): {abs(f_min - analytical):.10f}")
    print(f"diff: ({abs(x-x_analytical):.7f}, {y-y_analytical:.7f})" )


if __name__ == "__main__":
    main()
