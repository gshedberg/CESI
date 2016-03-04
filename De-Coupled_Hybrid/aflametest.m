j=1;
for phi=0.5:0.1:1
    phi_store(j) = phi ;
    
    T(1) = 1000 ; %[K] Initial flame temperature guess
    del_T = 10 ;
    T_ad(1) = 1000 ; %[K] Initial Adiabatic flame temperature setting
    i = 1 ; %[] counter
    T_inf = 298 ;
    
    while T(i) < T_ad(i)
        
        del_Hc = 50021 ; %kJ/kg Heat of combustion for methane: Combustion Science and Engineering pg 990
        n_CH4       = 0         ;
        % number of mols of reactants
        
        %phi = 1
        n_CO2       = 1*44         ; %[kg]
        n_H2O       = 2*18         ; %[kg]
        n_N2        = 3.76*2*28    ; %[kg]
        n_O2        = 0         ;
        
        n_PMMA      = 0         ;
        
        %Specific heats
        
        c_p_CO2     = -6.40953e-22*T(i)^6 +1.35711e-17*T(i)^5 -1.14921e-13*T(i)^4 +5.00001e-10*T(i)^3 -1.19628e-6*T(i)^2 +1.56371e-3*T(i) +4.71107e-1 ; %[kJ/(kg K)] 175K-6000K Engineering Toolbox
        
        c_p_H2O     = 1.24400e-21*T(i)^6   -2.52232e-17*T(i)^5 +1.98003e-13*T(i)^4 -7.41265e-10*T(i)^3 +1.23513e-6*T(i)^2 -2.26210e-4*T(i) +1.84352e0 ; %[KJ/(kg K)] 175K-6000K Engineering Toolbox
        
        c_p_N2      = 6.68020e-22*T(i)^6  -1.28835e-17*T(i)^5 +9.47330e-14*T(i)^4 -3.25777e-10*T(i)^3  +4.92137e-7*T(i)^2 -1.29658e-4*T(i) +1.04174e0 ; %[KJ/(kg K)] 175K-6000K Engineering Toolbox
        
        c_p_O2      = 2e-22*T(i)^6  -3e-18*T(i)^5 +2e-14*T(i)^4 -2e-11*T(i)^3 -9e-8*T(i)^2 +0.0003*T(i) +0.83   ; %[KJ/(kg K)] 175K-6000K Engineering Toolbox
        c_p_CH4     = 1e-20*T(i)^6  -2e-16*T(i)^5 +1e-12*T(i)^4 -4e-09*T(i)^3  +4e-6*T(i)^2 +0.0009*T(i) +1.7695 ; %[KJ/(kg K)] Combustion Science and Engineering pg 1012-1013
        c_p_PMMA = 0 ;
        
        % T_ad(i+1) = T_inf + del_Hc/(n_CO2*c_p_CO2 +n_H2O*c_p_H2O +n_N2*c_p_N2 +n_O2*c_p_O2 +n_CH4*c_p_CH4 + n_PMMA*c_p_PMMA) ;
        T_ad(i+1) = T_inf + (del_Hc*16)/(n_CO2*c_p_CO2 +n_H2O*c_p_H2O +n_N2*c_p_N2 ) ;
        
        T(i+1) = T(i) + del_T ;
        i=i+1 ;
        
        if i > 10000
            break
        end
        
    end
    
    format short
    T_ad_soln(j) = T_ad(end-1) ;
    
    T_ad_stoic_turns = 2226 ; %[K] adiabatic flame temperature of stoiciometric methane from Turns
    
    j=j+1 ;
end