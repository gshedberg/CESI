function [Tout,afr] = aflame(TXNinAir,TXNinFuel)
h = enthalpy(1200);
LHVCO = 303000; %Lower Heating Value of CO
LHVH2 = 240420; %Lower HEating Value of H2
LHVfuel = ones(100,1).*8e5;
LHVfuel = zeros(length(TXNinFuel(:,1)),1)+LHVfuel.*TXNinFuel(:,2)+LHVCO.*TXNinFuel(:,3)+LHVH2.*TXNinFuel(:,5);%LHV of anode off products

hrxn1 = 2*h(:,5)+h(:,2)-h(:,1)-1.5*h(:,7); %Ch4 + 1.5 O2 --> CO + 2 H2O
hrxn2 = h(:,3)-h(:,2)-.5*h(:,7); %CO + .5 O2 --> CO2 
hrxn3 = h(:,5)-h(:,4)-.5*h(:,7); %H2 + .5 O2 -->  H2O

Hrxn = hrxn1+hrxn2+hrxn3;


TinAir = TXNinAir(:,1);
XinAir = TXNinAir(:,2:8);
NinAir = TXNinAir(:,9);

TinFuel = TXNinFuel(:,1);
XinFuel = TXNinFuel(:,2:8);
NinFuel = TXNinFuel(:,9);

cpair = SpecHeat(TinAir,XinAir);
cpfuel = SpecHeat(TinFuel,XinFuel);

hinair = enthalpy(TinAir,XinAir,NinAir);
hinfuel = enthalpy(TinFuel,XinFuel,NinFuel);
h = enthalpy(1200);

afr = NinAir./NinFuel;

Tout =(afr./(1+afr)).*(LHVfuel./(cpfuel));

plot(afr,Tout);
