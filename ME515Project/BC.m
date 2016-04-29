%Don't know third boundary condition for Blasius solution
N = 2000;
q = zeros(length(N),1);
z = zeros(length(N),1);
for i=1:N
    q(i)=i/4000;
    eta = [0,10]; %10 is equivalent to inf
    f = [0 0 q(i)];
    [t,x]=ode45(@Blasius,eta,f);
    z(i)=x(length(x),2); % x(length(x),2) is the value of f2(10)
end;
plot(q,z);
title('Find third BC')
grid
xlabel('q');
ylabel('L');
inf = (z(133)+z(132))/2;
q_0 = (q(133)+q(132))/2;