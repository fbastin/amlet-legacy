set terminal epslatex color "enhanced"

set output "gsl_e1.tex"

unset key

set label "$e^{-(x-1)^2}sin(8x)-0.89e^{-(x-1000)^2}$"
plot [-2:4]exp(-(x-1)**2)*sin(8*x)-exp(-(x-1000)**2)*0.89
