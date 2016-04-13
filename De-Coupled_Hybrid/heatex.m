function [T_hot,T_cold] = heatex(T_hot_in,X_hot_in,N_hot,T_cold_in,X_cold_in,N_cold,P_h, P_c,Ofrac_VT,Ofrac_recirc,HX_dTmin)
%This function simulates a heat exchanger 
%Marshall's code%
?
[~,H_ColdIn] = enthlapy(T_cold_in,X_cold_in, N_cold);                  %Enthalpy (kJ/kg) of the cold inlet
[~,H_HotIn] = enthalpy(T_hot_in,X_hot_in,N_hot);                       %Enthalpy (kJ/kg) of the hot outlet
Delta_c = N_cold.*(enthalpy((T_hot_in-HX_dTmin),X_cold_in,N_cold)-H_ColdIn); %Difference in enthalpys (kJ/kg) of the cold side
Delta_h = N_hot.*(H_HotIn-enthalpy((T_cold_in+HX_dTmin),X_hot_in,N_cold));  %Difference in enthalpys (kJ/kg) of the hot side
?
T_error = ones(length(Delta_c),1);                                  %Presetting T_error to be a column vector
T_cold = ones(length(T_hot_in),1);                                           %Presetting T_ho to be a column vector
T_hot = ones(length(T_hot_in),1);                                           %Presetting T_co to be a column vector
H_guess = ones(length(T_hot_in),1);
Tnew =  ones(length(T_hot_in),1);
deltaT =  ones(length(T_hot_in),1);
A = linspace(1, length(T_cold),length(T_cold))';
    tol = 1e-3;
if Delta_c > Delta_h                                               %If the difference in enthalpy on the cold side is greater than the difference on the hot side
    Q = Delta_h;                                                    %Setting heat flux (kJ/kg) to the enthalpy difference on hot side
    T_hot = T_cold_in+HX_dTmin;                                            %?Since Delta_h is smaller, T_ci must be hotter than T_hi; so T_ho is now T_ci plus HX_dTmin
    T_cold = T_hot_in-HX_dTmin;                                            %?Since Delta_h is smaller, T_hi must be colder than T_ci; so T_co is now T_ci minus HX_dTmin
    H_cold = H_ColdIn+Q./N_cold;                                     %Enthalpy of cold inlet plus heat flux multiplied by the mass flow rate (kg/s) of the cold side
    Cp = SpecHeat(T_cold,X_cold_in);                                %Specific heat (kJ/kg) based on temperature and pressure (kPa) of the cold outlet
?
    while ~isempty(A)                                     %Iterate until the error is smaller than 1e-3
        H_guess(A) = enthalpy(T_cold(A),X_cold_in(A),N_cold(A));                       %Guess of enthalpy based on new cold temperature and cold pressure
        T_error(A) = (H_guess(A) - H_cold(A))/Cp(A).*N_cold(A);                             %Adjusting the error in temperature based on known enthalpy and specific heat of the cold side
        T_cold(A) = T_cold(A)-0.3.*T_error(A);                                   %Subtraction of a portion of the T_error from cold outlet temp to get closer to the actual temp
        A = find(abs(T_error)>tol);
    end
else                                                                 %If the difference in enthalpy on the hot side is greater than the difference on the cold side
    Q = Delta_c;                                                    %Setting heat flux (kJ/kg) to the enthalpy difference on hot side
    
    T_cold = T_hot_in-HX_dTmin;                                    %?Since Delta_h is larger, T_hi must be hotter than T_ci; so T_co is now T_ci minus HX_dTmin
    dTcold = T_cold-T_cold_in;
    dThot = dTcold.*Delta_c./Delta_h;
    T_hot = T_hot_in-dThot;                                    %?Since Delta_h is larger, T_ci must be colder than T_hi; so T_ho is now T_ci plus HX_dTmin
    H_hot = H_HotIn-Q./N_hot;                                       %Enthalpy of hot inlet minus heat flux multiplied by the mass flow rate (kg/s) of the cold side
    Cp = (.25*SpecHeat(T_hot_in,X_hot_in)+.75*SpecHeat(T_hot,X_hot_in));                                %Specific heat (kJ/kg) based on temperature and pressure (kPa) of the hot outlet 
    deltaT = Q./(N_hot.*Cp);
    Tnew = max(T_hot_in-deltaT,T_cold_in+2*HX_dTmin);
    T_hot = .5*T_hot+0.5*Tnew;  
%     T_ho = TemperatureHydrogen(H_hot,P_h,T_ho);
    while ~isempty(A)                                       %Iterate until the error is smaller than 1e-3
                                         %Subtraction of a portion of the T_error from hot outlet temp to get closer to the actual temp
        H_guess(A) = EnthalpyHydrogen(T_hot(A), P_h(A));                       %Guess of enthalpy based on new cold temperature and cold pressure
        T_error(A) = (H_guess(A) - H_hot(A))./Cp(A);                             %Adjusting the error in temperature based on known enthalpy and specific heat of the cold side
        T_hot(A) = T_hot(A)-0.05.*T_error(A);                                   %Subtraction of a portion of the T_error from cold outlet temp to get closer to the actual temp
        A = find(abs(T_error)>tol);
    end
end
?