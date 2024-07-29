function yprime=sub_ode45(t,y)
global M A beta r 
F=[0;0;-M\[r*y(1)^3;beta*y(2)^3]];
yprime=A*y+F;
