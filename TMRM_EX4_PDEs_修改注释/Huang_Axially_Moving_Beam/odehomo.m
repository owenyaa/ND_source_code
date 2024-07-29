function yprime=odehomo(t,x)
global A f1 f2 w k12 k13 k22 k23 M
F=[0;0;M\[f1*cos(w*t)-k12*x(1)*x(2)^2-k13*x(1)^3;f2*cos(w*t)-k22*x(2)*x(1)^2-k23*x(2)^3]];
yprime=A*x+F;