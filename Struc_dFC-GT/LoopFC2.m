Oxidant.Supply = linspace(.003,.038);
k = 100;
q = zeros(100,1);
T = 1023;
ASR = .25;
e2 = .6;
S2C = 2;
L = 10;
W = 10;
n = 10;
Cells = 1e5;
Pr = 15;

Currentd = zeros(n,length(q));
Current = zeros(n,length(q));
Recirc = zeros(1,length(q));
FuelOut = zeros(n,length(q));
H2out = zeros(n,length(q));
V_op = zeros(1,length(q));
P_op = zeros(1,length(q));
Avg_Currentd = zeros(length(q),1);
Avg_Current = zeros(length(q),1);

for j = 1:1:k
    Oxidant.O2 = Oxidant.Supply(j);
    [i,r,FuelFlow,H2,V,P] = FuelCell2(T,ASR,e2,S2C,Oxidant,L,W,n,Cells,Pr);
    Currentd(:,j) = i;
    Current(:,j) = mean(i)*L*W;
    Recirc(1,j) = r;
    FuelOut(:,j) = FuelFlow;
    H2out(:,j) = H2;
    V_op(:,j) = V;
    P_op(:,j) = P;
    Avg_Currentd(j,:) = sum(i)/n;
    Avg_Current(j,:) = sum(i*L*W/n);
end

figure(1)
plot(Avg_Currentd,V_op,'Linewidth',2)
ylabel('Operating Voltage')
xlabel('Average Current Density')



