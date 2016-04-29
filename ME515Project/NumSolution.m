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
eta_sol = eta(42);

T = 27; %degrees C
L = 10;
n = 89;
for i = 1:n
    x(i) = (i-1)/8.8;
end
f_prime = f(:,2);
nu = 1.568e-5;
U_1 = (10000*nu)/5;
Re_1 = (U_1*5)/nu;
u_1 = f(:,2).*U_1;
U_2 = (20000*nu)/10;
Re_2 = (U_2*10)/nu;
u_2 = f(:,2).*U_2;
y_1 = eta(42).*5*Re_1^-.5;
y_2 = eta(42).*10*Re_2^-.5;
y = eta.*5*Re_1^-.5;
figure(2)
plot(y,u_1)
grid
xlabel('Length')
ylabel('Velocity Distribution')
title('Boundary Layer Velocity Distribution')
figure(3)
plot(y,u_2)
grid
xlabel('Length')
ylabel('Velocity Distribution')
title('Boundary Layer Velocity Distribution')
w = (10./length(f_prime));
constant1 = 5*Re_1^-.5;
constant2 = 10*Re_2^-.5;
for i = 1:length(eta)
    disp_thickness_1(i) = constant1*(1-f_prime(i));
    disp_thickness_2(i) = constant2*(1-f_prime(i));
    Mom_thickness_1(i) =  constant1*(f_prime(i).*(1-f_prime(i)));
    Mom_thickness_2(i) =  constant2*(f_prime(i).*(1-f_prime(i)));
end
disp_thickness1 = y_1.*sum(disp_thickness_1);
disp_thickness2 = y_2.*sum(disp_thickness_2);
Mom_thickness1 = y_1.*sum(Mom_thickness_1);
Mom_thickness2 = y_2.*sum(Mom_thickness_2);
