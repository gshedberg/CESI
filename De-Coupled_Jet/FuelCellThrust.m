function [Xout,Nout, W_out, n_fuel,Eff_FC] = FuelCellThrust(V_fc,Tdesign, TXNoxidant,TXfuel)

F = 96485;  
h = enthalpy(Tdesign);
hrxn1 = 2*h(:,5)-2*h(:,4)-h(:,7); %2H2 + O2 -->  2H2O
% hrxn2 = h(:,3)+h(:,4)-h(:,2)-h(:,5); %CO + H20 --> CO2 + H2

XanodeIn = TXNoxidant(:,2:8)+TXfuel(:,2:8);
NanodeIn = TXNoxidant(:,9);


I = 4*1000*F*(TXNoxidant(:,8).*TXNoxidant(:,9));  %SOFC current X # of Cells based on oxidant flow and 100% oxidant utilization
[~, H_Oxidant] = enthalpy(TXNoxidant(:,1),TXNoxidant(:,2:8),TXNoxidant(:,9));

util = V_fc*0+.70;
n_fuel = (I/(2*1000*F))./util;
[~,Hfuel] = enthalpy(TXfuel(:,1),TXfuel(:,2:8),n_fuel);

while max(abs(error)) > .001
    
    n_H2_theoretical = (I/2*1000*F))./util; %Hydrogen flow rate
  
    R1 =(I/(4000*F));

    Xout(:,1) = (XanodeIn(:,1).*NanodeIn)./Nout;
    Xout(:,2) = (XanodeIn(:,2).*NanodeIn)./Nout;
    Xout(:,3) = (XanodeIn(:,3).*NanodeIn)./Nout;                 %Mass balance
    Xout(:,4) = (XanodeIn(:,4).*NanodeIn-(2*R1))./Nout;
    Xout(:,5) = (XanodeIn(:,5).*NanodeIn+(2*R1))./Nout;
    Xout(:,6) = (XanodeIn(:,6).*NanodeIn)./Nout;
    Xout(:,7) = (XanodeIn(:,7).*NanodeIn - R1)./Nout;

    [~,Hout] = enthalpy(Tdesign, Xout,Nout);
    TanodeIn = Tdesign;
    T_error = 100;
    while abs(T_error) > .1
       [~,H_guess] = enthalpy(TanodeIn, XanodeIn, NanodeIn);
       Cp = SpecHeat(TanodeIn, XanodeIn);

       T_error = (Hfuel-H_guess)./(Cp.*NanodeIn);
       TanodeIn = TanodeIn+T_error;
    end
    Qgen = -hrxn1.*R1 - (V_fc.*I)/1000;
    Qimbalance = H_Oxidant+Hfuel+Qgen-Hout;
    
    H2error = Qimbalance./hrxn1;
    H2equiv = (I/(8000*F))./util;
    
    utilNew = util.*(1-H2error./H2equiv);
    error = util - utilNew;
    util = utilNew;
end
Nout = NanodeIn + n_fuel;
W_out = V_fc.*(I/1000);
LHV = zeros(length(util),1)+ 240420;
Eff_FC = W_out./(n_fuel.*LHV);
