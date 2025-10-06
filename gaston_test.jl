using Gaston
Gaston.config.term = "x11 font ',10' size 640,480"
Gaston.config.output = :external

t = 0:0.01:1
plot(t, sin.(2π*5*t))
