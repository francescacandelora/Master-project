function [fig,usol,T] = driverjump(p,n,Iapp,u2)

mu = p(1);      
h = p(2);       
a1 = p(3);     
b1 = p(4);      
a2 = p(5);      
b2 = p(6);

% discretization of
% $x \in [-\pi,\pi]$
x = linspace(-pi,pi,n+1)';
x = x(1:end-1); 
% now $x$ has 16 components, distanced by
deltax = 2*pi/n;

% interval of integration
trange = [0:0.02:20]; %[0:2:1000]

% Initial condition. It's a function, as U is a vector with 
% $U_i = U(x_i,t)$
if nargin<4
    u0 = (Iapp(x,0)); 
else 
    u0 = u2(x);
end

% define the function f
f = @(u) 1./(1+exp(-mu*(u-h)));

% define w_ex, needed for M, where w_ex is the periodic extension of w on  
% the interval [-3pi,3pi]
w = @(z) a1*exp(-b1*abs(z))-a2*exp(-b2*abs(z));
w_ex = @(z) w(z+2*pi).*(heaviside(-z-pi))+w(z).*(heaviside(z+pi)-heaviside(z-pi))+w(z-2*pi).*heaviside(z-pi);

% Define M as a matrix of differences of w computed on the x grid
M = zeros(n,n);
for j = 1:n
    for k = 1:n
        M(j,k) = w_ex(x(j)-x(k))*deltax;
    end
end


% Now call the function for the vectorial field of the equation we wish to
% solve
nf = @(t,u) neuralField(t,u,n,f,M,Iapp);

% Use ode45 for solutions
[tsol,usol] = ode45(nf,trange,u0);
x = [x;pi];
usol = [usol,usol(:,1)];

[X,T] = meshgrid(x,tsol); %I need this line because T is an output



fig = figure;
%size(usol)
% surf(T,X,usol);

pcolor(T,X,usol)

xlabel('t');ylabel('x');zlabel('u');
shading flat;
colormap(flipud(bone))
%colormap(brewermap(9,'Blues'))
colorbar; % caxis o clim ([a b]) forces black to be at b and white at a, so I have the same values for same colors in different images.
%caxis([-15, 25]);
hold on;


pva = angle(usol(:,1:16)*exp(i*x(1:16,:)));

plot(linspace(T(1),T(end),length(pva)),pva,'.','Color','#D95319','Linewidth',2); 


% hold off;
end


