function [nu, j] = Activation(cR, cP, T,


global h R F n cR0 cP0
h = 
R = 
F = 
n = 2
cR0 = 
cP0 = 
al = .5;
g0 = -237000;    
dg = g0;  
pact = exp(-1*(dg/(R*T)));
pactr = exp(-1*(dg-g)/(R*T));
T0 = 300;
f1 = (k*T0)/h;
f2 = (k*T0)/h;
forward = cR*f1*pact;
reverse = cP*f2*pactr;
net = forward - reverse;
j1 = n*F*cR*f1*pact;
j2 = n*F*cP*f2*exp(-1(dg-g+n*F*net)/R*T);
if j1 == j2
    j0 = 1;
end
if j0 == 1
    nu = ((R*T)/(al*n*F))*log(((cR0/cR)+(cP/cP0)*j0));
    j1 = j0*exp((al*n*F*nu)/(R*T));
    j2 = j0*exp((-1*(1-al)*n*F*nu)/(R*T));
    j = j0*(j1-j2);
end
