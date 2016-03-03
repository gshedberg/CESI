n = 100;
for i=1:n 
    Efficiency = zeros(i,n);
    Eff_FC = zeros(i,n);
    Eff_GT = zeros(i,n);
    W_net = zeros(i,n);
    W_fc = zeros(i,n);
    W_gt = zeros(i,n);
    T_out = zeros(i,n);
    N_out = zeros(i,n);
    V_fc = zeros(i,n);
    R_actual = zeros(i,n);
    FC_util = zeros(i,n);
    E0 = zeros(i,n);
    Qimbalance = zeros(i,n);
    recovery = linspace(.05,.7,(1/n));
    VaryHybrid_Final
    
end