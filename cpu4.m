function dAdt = cpu4(t,A,n,S,f,M,param)
% This function defines the linear system
% $A' = G(A)$
% Here we are donsidering memory decay as well.
theta = linspace(-pi,pi,n+1)';
theta = theta(1:end-1); 

% Unpack parameters
alpha = param(1);
beta = param(2);
md = param(3);
m = param(4);

% Unpack A
ith = 1:n;
iu  = n+1:2*n;
izL = 2*n+1:3*n;
izR = 3*n+1:4*n;

th = A(ith);
u  = A(iu);
zL = A(izL);
zR = A(izR);

% Define the vectorial field
dAdt(ith,1) = S(t);
dAdt(iu,1)  = -u+M*f(u)+(cos(theta-th)+1).^5;
dAdt(izL,1) = alpha*subplus(-S(t))-beta*u-md*zL+m*eye(n)*f(zL);
dAdt(izR,1) = alpha*subplus(S(t))-beta*u-md*zR+m*eye(n)*f(zR);
end