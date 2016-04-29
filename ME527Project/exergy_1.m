function [X, X_s] = exergy(varargin)
T = varargin{1};
T0 = 300;
P0 = 100;
X =[];
if length(varargin)>1
    X = varargin{2};
    N = varargin{3};
end
[h0, h_s0] = enthalpy(T0)