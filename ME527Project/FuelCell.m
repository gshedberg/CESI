function [Xout,Nout, W_out, n_fuel,Eff_FC,util,Qimbalance,I] = FuelCell(V_fc,Tdesign, TXNoxidant,TXfuel,S2Cdesign)
F = 96485;  

% Reactions occuring in SOFC %
h = enthalpy(Tdesign);
erxn2 = .6; % Water gas shift effectiveness of reaction
erxn3 = 1; % Steam reformation of methane effectiveness of reaction
hrxn1 = 2*h(:,5)-2*h(:,4)-h(:,7); %2H2 + O2 -->  2H2O
hrxn2 = h(:,3)+h(:,4)-h(:,2)-h(:,5); %CO + H20 --> CO2 + H2 
hrxn3 = 3*h(:,4)+h(:,2)-h(:,1)-h(:,5); %CH4+H2O --> CO + 3H2

I = 4*1000*F*(TXNoxidant(:,8).*TXNoxidant(:,9));  %SOFC current X # of Cells based on oxidant flow and 100% oxidant utilization
[~, h_Oxidant] = enthalpy(TXNoxidant(:,1),TXNoxidant(:,2:8),TXNoxidant(:,9));
[H_oxidant,~] = enthalpy(TXNoxidant(:,1),TXNoxidant(:,2:8),TXNoxidant(:,9));
%Initial Assumption for reformation
util = V_fc*0+.70;
n_H2_theoretical = (I/(2*1000*F))./util;  %Hydrogen flow rate with 100% CH4 conversion
n_fuel = n_H2_theoretical./(4*TXfuel(:,2)+1*TXfuel(:,3)+1*TXfuel(:,5));
XanodeIn = zeros(length(Tdesign),7);
XanodeIn(:,1) = 1./(S2Cdesign+1);
XanodeIn(:,5) = S2Cdesign./(S2Cdesign+1);
NanodeIn = n_fuel.*(S2Cdesign+1);
RecircPerc = .3;
error = util*0+1;
while max(abs(error)) > .001 
    Under90 = (util<.9);
    util = min(util,.9);
    n_H2_theoretical = (I/(2*1000*F))./util;  %Hydrogen flow rate with 100% CH4 conversion
    n_fuel = n_H2_theoretical./(4*TXfuel(:,2)+1*TXfuel(:,3)+1*TXfuel(:,5));
    NanodeIn = n_fuel./(1-RecircPerc.*(1+2*XanodeIn(:,1)*erxn3));
    
   errorR = util*0+1;
   while max(abs(errorR))>.0001
        R3 = (XanodeIn(:,1).*NanodeIn)*erxn3;
        R2 = (XanodeIn(:,2).*NanodeIn +R3)*erxn2;
        R1 =(I/(4000*F));

        Nout = NanodeIn + 2*R3;        %Total Nout of anode 

        Xout(:,1) = (XanodeIn(:,1).*NanodeIn - R3)./Nout;
        Xout(:,2) = (XanodeIn(:,2).*NanodeIn-R2+R3)./Nout;
        Xout(:,3) = (XanodeIn(:,3).*NanodeIn+R2)./Nout;                 %Mass balance
        Xout(:,4) = (XanodeIn(:,4).*NanodeIn-(2*R1)+R2+(3*R3))./Nout;
        Xout(:,5) = (XanodeIn(:,5).*NanodeIn-R2-(R3)+(2*R1))./Nout;
        Xout(:,6) = (XanodeIn(:,6).*NanodeIn)./Nout;
        Xout(:,7) = (XanodeIn(:,7).*NanodeIn)./Nout;

        [~,hout] = enthalpy(Tdesign, Xout,Nout);
        [Hout,~] = enthalpy(Tdesign, Xout,Nout);
        H2Orecirc = Xout(:,5).*Nout.*RecircPerc;
        H2Oneeded = (n_fuel.*(TXfuel(:,2)+TXfuel(:,3)+TXfuel(:,4)).*S2Cdesign)-(n_fuel.*TXfuel(:,6));
        
        errorR = (H2Oneeded-H2Orecirc)./(Nout.*Xout(:,5));
        RecircPerc = RecircPerc + errorR;

        NanodeIn = n_fuel+Nout.*RecircPerc;
        for i = 1:1:7
            XanodeIn(:,i) = (Xout(:,i).*Nout.*RecircPerc+TXfuel(:,i+1).*n_fuel)./NanodeIn;
        end
        
        [~,h_recirc] =  enthalpy(Tdesign, Xout, Nout.*RecircPerc);
        [H_recirc,~] = enthalpy(Tdesign, Xout, Nout.*RecircPerc);
        [~,hfuel] = enthalpy(TXfuel(:,1), TXfuel(:,2:8),n_fuel);
        [Hfuel,~] = enthalpy(TXfuel(:,1), TXfuel(:,2:8),n_fuel);
        
        TanodeIn = Tdesign; 
           [~,h_guess] = enthalpy(TanodeIn, XanodeIn, NanodeIn);
           Cp = SpecHeat(TanodeIn, XanodeIn);

           T_error = (h_recirc+hfuel-h_guess)./(Cp.*NanodeIn);
           TanodeIn = TanodeIn+T_error;
   end
        
    Qreform = hrxn2.*R2+hrxn3.*R3;
    Qgen = -hrxn1.*R1- V_fc.*(I/1000);

    %energy balance: Fuel in+ oxidant in+recirculation in - reforming cooling+Qgen - flow out
    Qimbalance = h_Oxidant+hfuel+h_recirc-Qreform+Qgen-hout;
%     Qimbalance = H_oxidant-Hfuel-H_recirc+Qreform+Qgen+Hout;
    CH4error = Qimbalance./hrxn3;
    CH4eqivelant = (I/(8000*F))./util;
    
    utilNew = util.*(1-CH4error./CH4eqivelant);
    error = util-utilNew;
    for i=1:1:length(util)
        if Under90(i) ==0 && error(i)<0
            error(i) =0;
        end
    end
    util = util - error;  %Updating util     
end
Nout = Nout.*(1-RecircPerc);

W_out = V_fc.*(I/1000);          %Total power out based on Voltage and current
LHV = zeros(length(util),1)+8e5;
Eff_FC = W_out./(n_fuel.*LHV);









