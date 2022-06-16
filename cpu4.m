function dAdt = cpu4(t,A,n,S,f,param)
% This function defines the linear system
% $A' = G(A)$
% Here we are donsidering memory decay as well.
theta = linspace(-pi,pi,n+1)';
theta = theta(1:end-1); 

% Unpack parameters
alpha = param(1);
beta = param(2);
md = param(3);
M = param(4);

% Unpack A
ith = 1:n;
izL = n+1:2*n;
izR = 2*n+1:3*n;

th = A(ith);
zL = A(izL);
zR = A(izR);

u = (cos(theta-th)+1).^5;
% Define the vectorial field
dAdt(ith,1) = S(t);
dAdt(izL,1) = alpha*subplus(-S(t))-beta*u-md*zL+M*eye(n)*f(zL);
dAdt(izR,1) = alpha*subplus(S(t))-beta*u-md*zR+M*eye(n)*f(zR);
end

