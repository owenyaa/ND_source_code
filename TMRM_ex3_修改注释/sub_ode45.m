function yprime=sub_ode45(t,y)
global M A r1 r4 f w
F=[0;0;0;0;0;-M\[r1*y(1)^3;0;0;r4*y(4)^3;0]+M\[0;0;f*sin(w*t);0;0]];
yprime=A*y+F;
