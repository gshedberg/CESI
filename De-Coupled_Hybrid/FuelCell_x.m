function [i_x, E_x] = CurrentDist(V_design,Current,S2C,nodes)
F = 96485;
W = .0045 %Width of cells
i = Current/nodes;
for k = 1:1:nodes
    X_h2(k) = 1+e_wgs/3 - (2.*n_O2*r)/(3*n_CH4) - (W.*(1-r))./

