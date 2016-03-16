%Now that we know the value of g(N) =1
%We can resolve to plot eta vs f'
eta_0 = [0,10];
q_0 = .3313;
f_0 = zeros(3,1);
f_0(1) = 0;
f_0(2) = 0;
f_0(3) = .3313;
[eta,f]=ode45(@Blasius, eta_0,f_0);
figure (1)
plot(eta,f(:,2));
grid 
title('Eta vs Ratio of Velocity to Free Stream')
xlabel('eta'); 
ylabel('u/U');
approx = f(43,2);
eta_sol = eta(43);

T = 27; %degrees C
L = 10;
n = 89;
for i = 1:n
    x(i) = (i-1)/8.8;
end
f_1 = f(:,1);
U = (10000*1.568)/5;
nu = 1.568;
Re = (U*x)/nu;
u = f(:,2).*U;
figure(2)
plot(x,u)
grid
xlabel('Length')
ylabel('Velocity Distribution')
title('Boundary Layer Velocity Distribution')

bound_thickness = (Re.^-.5).*x;
fun1 = @(x) 1-f(:,2);
fun2 = @(x) f(:,2).*(1-f(:,2));
disp_thickness = integral(fun1,0,Inf);
Mom_thickness = integral(fun2,0,Inf);
