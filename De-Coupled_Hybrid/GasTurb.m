function [Wnet, Wt, Wc, Eta] = GasTurb(Pr, Nt, Nc)
global Cv Cp
oneVec = ones(1,length(Pr));
Cv  = .7180;
Cp = 1.0050;
gam = Cp/Cv;
P1 = 1;
T1 = 300;
PrMat = Pr'*oneVec;
persistent P_feed
    P_feed = PrMat'*P1;
    
persistent T2s
     T2s = PrMat.^(1-1/gam)*T1;
persistent T2
     T2 = (T2s-T1)./Nc+T1;
Wc = Cp*(T2-T1);


T3Mat = T3*oneVec;
T4s = PrMat.^(-1+1/gam).*T3Mat;
T4 = (T4s-T3Mat).*Nt+T3Mat;
Wt = Cp*(T3Mat-T4);
Wnet = Wt-Wc;
Eta = Wnet./((T3Mat-T2)*Cp);
end


