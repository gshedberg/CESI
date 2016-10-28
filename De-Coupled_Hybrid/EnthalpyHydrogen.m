function H = EnthalpyHydrogen(varargin)
% This function calculates the Enthalpy of cryogenic hydrogen using a
% Using a MatLab function named refpropm.m that calls a program called RefProp
% Input a temperature (K) and pressure (kPa) it calculates enthalpy (J/kg) of equilibrium hydrogen
% Input a temperature (K), Pressure (kPa), and ortho fraction and calculates properties of the mixture

T = varargin{1};                                                    %Temperature input (K)
P = varargin{2};                                                    %Pressure input (kPa)
    
if length(varargin) == 2      
    Ofrac = OP_Equilib(T);                                          %Equilibrium ortho-parahydrogen at a given temperature (K)
elseif length(varargin) == 3 
    Ofrac = varargin{3};                                            %Specific ortho fraction
end
    
if length(varargin) == 2 || length(varargin) == 3 
    Hp = refproparray('H','T',T,'P',P,'PARAHYDROGEN')/1000';             %Enthalpy (kJ/kg) of para hydrogen (Hp)
    Ho = zeros(length(Hp),1);
    A = logical((T < 34.59624).*(P >= 10000)); %handles cases when some entries exceed values, and others don't
    Ho(A) = Hp(A) - (T(A)-20)*.2403861 + 702.4976; %At 34.59624 degrees the difference between ortho and para is 3.5087kJ/kg + the conversion energy
    if nnz(~A)>0
        Ho(~A) = refproparray('H','T',T(~A),'P',P(~A),'ORTHOHYDROGEN')/1000' + 702.4976; %Enthalpy (kJ/kg) of ortho hydrogen (Ho) plus added correction constant
    end
elseif length(varargin) == 4 
    Ofrac = varargin{3};                                                   %Specific ortho fraction
    Quality = varargin{4}; 
    Ho_g = refproparray('H','P',P,'Q',1,'ORTHOHYDROGEN')/1000' + 702.4976; %Enthalpy (kJ/kg) of ortho hydrogen (Ho) plus added correction constant
    Ho_f = refproparray('H','P',P,'Q',0,'ORTHOHYDROGEN')/1000' + 702.4976; %Enthalpy (kJ/kg) of ortho hydrogen (Ho) plus added correction constant
    Ho = Ho_f + Quality.*(Ho_g - Ho_f);
    Hp_g = refproparray('H','P',P,'Q',1,'PARAHYDROGEN')/1000';             %Enthalpy (kJ/kg) of para hydrogen (Hp)
    Hp_f = refproparray('H','P',P,'Q',0,'PARAHYDROGEN')/1000';             %Enthalpy (kJ/kg) of para hydrogen (Hp)
    Hp = Hp_f + Quality.*(Hp_g - Hp_f);
end
H = Ofrac.*Ho + (1-Ofrac).*Hp;                                      %Enthalpy (kJ/kg) based on ortho fraction (Ofrac)

% %% where does 702.5976 come from (plot the following
% O_P_Conversion = [0.7033	0.703	0.7027	0.7034	0.7032	0.7029	0.7036	0.7033	0.7031	0.7028	0.7035	0.7032	0.703	0.7027	0.7034	0.7031	0.7028	0.7036	0.7033	0.703	0.7027	0.7034	0.7031	0.7028	0.7035	0.7031	0.7028	0.7024	0.703	0.7026	0.7022	0.7017	0.7022	0.7017	0.7011	0.7014	0.7007	0.7	0.6991	0.6993	0.6983	0.6972	0.6971	0.6959	0.695	0.693	0.693	0.691	0.689	0.688	0.686	0.684	0.683	0.681	0.678	0.677	0.674	0.671	0.669	0.666	0.663	0.66	0.657	0.654	0.65	0.646	0.643	0.638	0.635	0.63	0.627	0.622	0.618	0.614	0.608	0.604	0.599	0.594	0.589	0.584	0.578	0.573	0.568	0.562	0.557	0.551	0.546	0.539	0.533	0.528	0.522	0.516	0.51	0.504	0.497	0.491	0.486	0.479	0.473	0.467	0.46	0.455	0.448	0.442	0.435	0.429	0.423	0.418	0.411	0.405	0.399	0.393	0.387	0.381	0.375	0.369	0.362	0.356	0.35	0.345	0.339	0.333	0.328	0.322	0.317	0.311	0.306	0.3	0.295	0.29	0.284	0.279	0.274	0.269	0.264	0.259	0.254	0.249	0.245	0.24	0.235	0.23	0.227	0.222	0.217	0.213	0.209	0.204	0.201	0.196	0.192	0.189	0.184	0.181	0.178	0.173	0.17	0.166	0.163	0.16	0.155	0.152	0.149	0.146	0.143	0.139	0.137	0.134	0.131	0.128	0.125	0.122	0.12	0.117	0.114	0.112	0.109	0.107	0.105	0.102	0.1	0.097	0.096	0.093	0.091	0.088	0.087	0.085	0.083	0.081	0.079	0.077	0.075	0.074	0.072	0.071	0.068	0.067	0.065	0.064	0.063	0.061	0.06	0.058	0.057	0.056	0.054	0.053	0.052	0.051	0.049	0.048	0.047	0.046	0.044	0.043	0.042	0.041	0.04	0.04	0.038	0.037	0.036	0.035	0.035	0.034	0.033	0.032	0.032	0.031	0.03	0.029	0.029	0.028	0.027	0.027	0.026	0.025	0.025	0.024	0.023	0.023	0.022	0.022	0.021	0.02	0.02	0.019	0.019	0.018	0.017	0.018	0.017	0.017	0.016	0.016	0.015	0.015	0.014	0.015	0.014	0.013	0.013	0.012	0.012	0.012	0.012	0.012	0.011	0.011	0.011	0.011	0.01	0.01	0.009	0.01	0.009	0.009	0.008	0.008	0.009	0.008	0.008	0.007	0.008	0.007	0.007	0.007	0.007	0.007	0.006	0.006	0.006	0.006	0.006	0.005	0.006	0.005	0.005	0.005];
% T = linspace(15,400,300);
% Hp = refproparray('H','T',T,'P',101.325,'PARAHYDROGEN')/1e6;
% Ho = refproparray('H','T',T,'P',101.325,'ORTHOHYDROGEN')/1e6;
% figure(1) % shows refprop enthalpy values for Ortho and Para
% plot(T,Hp)
% hold on
% plot(T,Ho,'r')
% figure(2) %shows the difference between refprop ortho and para with this constant added to ortho, compared to a vector of the conversion enthalpy
% plot(T,Ho+.7024976-Hp)
% hold on
% plot(T,O_P_Conversion,'r')
