function [V_fc,E0] = nernst(Pr,T,V_loss)
Ru = 8.314;
F = 96485;
%Electrochemical reaction: %2H2 + O2 -->  2H2O

%Ideal Nernst Potential of SOFC%
[h,~] = enthalpy(T);
[s,~] = entropy(T);

E0 = -((h(:,5)-s(:,5).*T)-(h(:,4)-s(:,4).*T)-.5*(h(:,7)-s(:,7).*T))/(2*F);

%Actual Voltage calculation dependent on the Pessure Ratio, TempDesign and
%Composition treated as constant
V_nernst = E0-(((Ru.*T)./(4*F)).*log((1./Pr)*(.15^.5/(.30^.5*1))));

V_fc = V_nernst - V_loss;

