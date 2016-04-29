% Blasius Equation = f*f''+2*f''' = 0
%B.C's eta = 0; f=f'=0
%      eta =inf; f'=1

%Reduce 3rd order ODE to 1st order ODE
%f'=g;g'=h;h'=-f*h/2
%f(1) = 0; g(1) = 0; g(inf) = 1; f'(inf)-1=0

function df = Blasius(eta,f)
df = zeros(3,1);
df(1) = f(2);           
df(2) = f(3); 
df(3) = -f(1)*f(3)/2; 
end
