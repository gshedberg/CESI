function [psi_s_destroyed,Efficiency]= exergy(varargin)
TXNin = varargin{1};
Win = varargin{2};
Qin = varargin{3};
TXNout = varargin{4};
Wout = varargin{5};
Qout = varargin{6};
if length(varargin) == 9
   psi_s_in_1 = zeros(length(TXNin),7);
   psi_s_in_2 = zeros(length(TXNin),7);
   psi_s_out = zeros(length(TXNin),7);
   psi_s_destroyed = zeros(length(TXNin),7);
    TXN_2 = varargin{7};
    W_2 = varargin{8};
    Q_2 = varargin{9};
    T_0 = ones(length(TXNin(:,1)),1).*300;
    [~,H_0_1] = enthalpy(T_0,TXNin(:,2:8),TXNin(:,9));
    [h_s_0,~] = enthalpy(T_0);
    [~,S_0_1] = entropy(T_0,TXNin(:,2:8),TXNin(:,9));
    [s_s_0,~] = entropy(T_0);
    [~,H_0_2] = enthalpy(T_0,TXN_2(:,2:8),TXN_2(:,9));
    [~,S_0_2] = entropy(T_0,TXN_2(:,2:8),TXN_2(:,9));
    [h_s_0_2,~] = enthalpy(T_0);
    [s_s_0_2,~] = entropy(T_0);

    Tin_1 = TXNin(:,1);
    Xin_1 = TXNin(:,2:8);
    Nin_1 = TXNin(:,9);
    Tin_2 = TXN_2(:,1);
    Xin_2 = TXN_2(:,2:8);
    Nin_2 = TXN_2(:,9);
    
    [~,H_1] = enthalpy(Tin_1,Xin_1,Nin_1);
    [~,S_1] = entropy(Tin_1,Xin_1,Nin_1);
    [~,H_2] = enthalpy(Tin_2,Xin_2,Nin_2);
    [~,S_2] = entropy(Tin_2,Xin_2,Nin_2);
    [h_s_in_1,~] = enthalpy(Tin_1);
    [s_s_in_1,~] = entropy(Tin_1);
    [h_s_in_2,~] = enthalpy(Tin_2);
    [s_s_2_in,~] = entropy(Tin_2);
    
    Tout = TXNout(:,1);
    Xout = TXNout(:,2:8);
    Nout = TXNout(:,9);
    
    [~,H_out] = enthalpy(Tout,Xout,Nout);
    [~,S_out] = entropy(Tout,Xout,Nout);
    [h_s_out,~] = enthalpy(Tout);
    [s_s_out,~] = entropy(Tout);

    psi_in = ((H_1-H_0_1)-T_0.*(S_1-S_0_1))+((H_2-H_0_2)-T_0.*(S_2-S_0_2));
    for i = 1:7
        psi_s_in_1(:,i) = (Nin_1.*Xin_1(:,i)).*((h_s_in_1(i)-h_s_0(i)) - T_0.*(s_s_in_1(i) - s_s_0(i)));
        psi_s_in_2(:,i) = (Nin_2.*Xin_2(:,i)).*((h_s_in_2(i)-h_s_0_2(i)) - T_0.*(s_s_2_in(i) - s_s_0_2(i)));
        psi_s_out(:,i) = (Nout.*Xout(:,i)).*((h_s_out(i)-h_s_0(i)) - T_0.*(s_s_out(i) - s_s_0(i))+((h_s_out(i)-h_s_0(i)) - T_0.*(s_s_out(i) - s_s_0(i))));
        psi_s_destroyed(:,i) = ((psi_s_in_1(:,i)+Win+Qin)+psi_s_in_2(:,i)+W_2+Q_2) - psi_s_out(:,i)-Wout-Qout;
        Efficiency(:,i) = 1- psi_s_destroyed(:,i)./((psi_s_in_1(:,i)+Win+Qin)+psi_s_in_2(:,i)+W_2+Q_2);
    end
else
    psi_s_in = zeros(length(TXNin),7);
    psi_s_out = zeros(length(TXNin),7);
    psi_s_destroyed = zeros(length(TXNin),7);
    T_0 = ones(length(TXNin(:,1)),1).*300;
    [~,H_0] = enthalpy(T_0,TXNin(:,2:8),TXNin(:,9));
    [hs_0,~] = enthalpy(T_0);
    [~,S_0] = entropy(T_0,TXNin(:,2:8),TXNin(:,9));
    [ss_0,~] = entropy(T_0);

    Tin = TXNin(:,1);
    Xin = TXNin(:,2:8);
    Nin = TXNin(:,9);

    [~,H_1] = enthalpy(Tin,Xin,Nin);
    [hs_1,~] = enthalpy(Tin);
    [~,S_1] = entropy(Tin,Xin,Nin);
    [ss_1,~] = entropy(Tin);
    Tout = TXNout(:,1);
    Xout = TXNout(:,2:8);
    Nout = TXNout(:,9);
    [~,H_2] = enthalpy(Tout,Xout,Nout);
    [hs_2,~] = enthalpy(Tout);
    [~,S_2] = entropy(Tout,Xout,Nout);
    [ss_2,~] = entropy(Tout);
    for i = 1:7
        psi_s_in(:,i) = (Nin.*Xin(:,i)).*((hs_1(i)-hs_0(i)) - T_0.*(ss_1(i) - ss_0(i)));
        psi_s_out(:,i) = (Nout.*Xout(:,i)).*((hs_2(i)-hs_0(i)) - T_0.*(ss_2(i) - ss_0(i)));
        psi_s_destroyed(:,i) = psi_s_in(:,i)+Win+Qin - psi_s_out(:,i)-Wout-Qout;
        Efficiency(:,i) = 1- (psi_s_destroyed(:,i)./(psi_s_in(:,i)+Win+Qin));
    end
end