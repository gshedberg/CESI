function [h,s] = enthalpy_refprop(varargin) % enthalpy (h) and entropy (s) calculated based on Refprop parameters T, P for a given flow composition X
T = varargin{1}; %Kelvin
P = varargin{2}; %kPa
X = varargin{3};
N = varargin{4}; %kmol/sec
% if length(varargin) > 4
%     S = varargin{5};
%     %Enthalpies kJ/kg
%     h_CH4 = refproparray('H','P',P,'S',S,'Methane')/1000;
%     h_CO = refproparray('H','P',P,'S',S,'CO')/1000;
%     h_CO2 = refproparray('H','P',P,'S',S,'CO2')/1000;
%     h_H2 = refproparray('H','P',P,'S',S,'Hydrogen')/1000;
%     h_H2O = refproparray('H','P',P,'S',S,'Water')/1000;
%     h_N2 = refproparray('H','P',P,'S',S,'Nitrogen')/1000;
%     h_O2 = refproparray('H','P',P,'S',S,'Oxygen')/1000;
% else
%   Enthalpies kJ/kg
h_CH4 = refproparray('H','T',T,'P',P,'Methane')/1000;
h_CO = refproparray('H','T',T,'P',P,'CO')/1000;
h_CO2 = refproparray('H','T',T,'P',P,'CO2')/1000;
h_H2 = refproparray('H','T',T,'P',P,'Hydrogen')/1000;
h_H2O = refproparray('H','T',T,'P',P,'Water')/1000;
h_N2 = refproparray('H','T',T,'P',P,'Nitrogen')/1000;
h_O2 = refproparray('H','T',T,'P',P,'Oxygen')/1000;

%MolarMasses (kg/kmol)
m_CH4 = refpropm('M','T',T,'P',P,'Methane'); 
m_CO = refpropm('M','T',T,'P',P,'CO');
m_CO2 = refpropm('M','T',T,'P',P,'CO2');
m_H2 = refpropm('M','T',T,'P',P,'Hydrogen');
m_H2O = refpropm('M','T',T,'P',P,'Water');
m_N2 = refpropm('M','T',T,'P',P,'Nitrogen');
m_O2 = refpropm('M','T',T,'P',P,'Oxygen');
%-------------------------------------------------------------------------%
n = length(T);
h = zeros(n,7);
s = zeros(n,7);
%Mass Flow
m = zeros(n,7);
%Entropies kJ/kg-K
s_CH4 = refproparray('S','T',T,'P',P,'Methane')/1000;
s_CO = refproparray('S','T',T,'P',P,'CO')/1000;
s_CO2 = refproparray('S','T',T,'P',P,'CO2')/1000;
s_H2 = refproparray('S','T',T,'P',P,'Hydrogen')/1000;
s_H2O = refproparray('S','T',T,'P',P,'Water')/1000;
s_N2 = refproparray('S','T',T,'P',P,'Nitrogen')/1000;
s_O2 = refproparray('S','T',T,'P',P,'Oxygen')/1000;
for i=1:1:n
    h(i,1) = h_CH4(i).*X(i,1);
    h(i,2) = h_CO(i).*X(i,2);
    h(i,3) = h_CO2(i).*X(i,3);
    h(i,4) = h_H2(i).*X(i,4);
    h(i,5) = h_H2O(i).*X(i,5);
    h(i,6) = h_N2(i).*X(i,6);
    h(i,7) = h_O2(i).*X(i,7); 

    s(i,1) = s_CH4(i).*X(i,1); 
    s(i,2) = s_CO(i).*X(i,2); 
    s(i,3) = s_CO2(i).*X(i,3); 
    s(i,4) = s_H2(i).*X(i,4); 
    s(i,5) = s_H2O(i).*X(i,5); 
    s(i,6) = s_N2(i).*X(i,6); 
    s(i,7) = s_O2(i).*X(i,7);
    
    m(i,1) = m_CH4; 
    m(i,2) = m_CO; 
    m(i,3) = m_CO2; 
    m(i,4) = m_H2; 
    m(i,5) = m_H2O; 
    m(i,6) = m_N2; 
    m(i,7) = m_O2;
end
m_dot = sum(m.*X,2).*N;  %kg/sec
h = sum(h,2).*m_dot; %kJ/sec
s = sum(s,2); %kJ/kg-K