function SOFCNernst(nCurrent,T,P,AnodeH2_H2OAvg,O2,block)
 F Ru Tags
nodes = length(nCurrent);
H2 = AnodeH2_H2OAvg(1:nodes,1);
H2O = AnodeH2_H2OAvg(1:nodes,2);

[h,~] = enthalpy(T);
[s,~] = entropy(T);
E0 = -((h(:,5)-s(:,5).*T)-(h(:,4)-s(:,4).*T)-.5*(h(:,7)-s(:,7).*T))/(2*F);
% frac = max(.0001,sqrt(O2).*H2./H2O); %prevent log(0) by always assuming at least .02% H2
frac = sqrt(O2).*H2./H2O; 
V_Nernst = E0 + Ru*T/(2*F).*log(abs(sqrt(P/101.325)*frac));
currentDen = nCurrent/block.A_Node;
frac2 = O2-Ru*T/(4*F).*currentDen*block.t_Cath/block.DeffO2/(P*1000);
V_activation = Ru*T/(4*F*block.alpha).*log(currentDen/block.Io/(P/101.325)./frac2);
V_concentration =0*nCurrent;
V_ohmic = currentDen*block.t_Membrane.*T/block.ElecConst./(exp(-block.deltaG./(Ru*T)));

Tags.(block.name).LocalCurrentDensity = nCurrent/block.A_Node/1e4;%current in A/cm2
Tags.(block.name).nCurrent = nCurrent;
Tags.(block.name).nVoltage = V_Nernst-V_activation-V_concentration-V_ohmic;
Tags.(block.name).nPower = nCurrent.*Tags.(block.name).nVoltage;% power per node in W
Tags.(block.name).Power = sum(Tags.(block.name).nPower)*block.Cells/1000;%power for stack in kW
Tags.(block.name).LocalNernst = V_Nernst;
Tags.(block.name).LocalActivation = V_activation;
Tags.(block.name).LocalConcentration = V_concentration;
Tags.(block.name).LocalOhmic = V_ohmic;
