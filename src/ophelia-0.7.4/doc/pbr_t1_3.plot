set terminal epslatex color "enhanced"

set output "pbr_t1_3.tex"
set view 53,322

unset key

set label "$\frac{(x+1)^3}{3}+y$"
splot [0:5][0:5]((x+1)**3)/3+y
