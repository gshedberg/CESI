function [T_out] = heatex2(n_air, T_in)
[T2] = heatex(n_air, T_in);
x_in = [0 0 0 0 0 0 1];
H_out = n_air*(enthalpy(T2) - enthalpy(T_in));
T_guess = 400;
T_error = 100;
while T_error > .1
    H_guess = enthalpy(T_guess, x_in, n_air);
    Cp = SpecHeat(T_guess, x_in);
    
    T_error = (H_guess - H_out)/(Cp*n_air);
    T_out = T_guess+T_error;
end
T_out = sum(T_out(1:7)/7);


