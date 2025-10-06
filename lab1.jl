using PyCall
pygui(:tk)
using PyPlot
using Colors

x = range(-5, stop=5, length=100)
y = range(-5, stop=5, length=100)
X = [i for i in x, j in y]
Y = [j for i in x, j in y]

Z = @. 1/(1 + X^2) + 1/(1 + Y^2)

fig = figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection="3d")

colors_blue_red = [colorant"blue", colorant"cyan", colorant"yellow", colorant"red"]
cmap_custom = PyPlot.ColorMap("blue_red", colors_blue_red)
surf = ax.plot_surface(X, Y, Z, cmap=cmap_custom, edgecolor="none")

colorbar(surf, ax=ax, shrink=0.5, aspect=5)

ax.set_xlabel("X")
ax.set_ylabel("Y")
ax.set_zlabel("Z")
ax.set_title("График функции f(x, y) = 1/(1+x²) + 1/(1+y²)")

ax.view_init(elev=30, azim=45)

show()
