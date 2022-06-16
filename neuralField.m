function dudt = neuralField(t,u,n,f,M,IApp)
% This function defines the linear system
% $U' = -U + Mf(U)
% Both matrix M and function f are given as parameters.

% Define the discretization of x
x = linspace(-pi,pi,n+1)';
x = x(1:end-1); 
deltax = 2*pi/n;

%IApp = @(z,t) 8*(z >= (x(2)) & z<=(x(5)+deltax/2))*(t>=25); 
%IApp = @(z,t) 8*(z>=-pi+t/20&z<-pi+t/20+7*deltax);
%IApp = @(x,t) 2*(cos(x-0.08*t)+1).^2;

% Define the vectorial field
dudt = -u+M*f(u)+(IApp(x,t));

end

