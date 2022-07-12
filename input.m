%% ora provo ad aggiungere IApp

clear
close all
clc

n = 16; %16

mu = 10;   %10   
h = 3;       %3
a1 = 30;     %30    % 20
b1 = 5;      %5   $ 4.5
a2 = 5.5;     %5.5   % 21
b2 = 1.5;    %0.1     % 5

p = [mu,h,a1,b1,a2,b2];

x = linspace(-pi,pi,n+1)';
x = x(1:end-1); 
% % now $x$ has 16 components, distanced by
% deltax = 2*pi/n;

% Strength of input and/or init cond
s = 8;
% exponent of input and/or init cond
pa = 2;

%Iapp = @(z,t) s*(cos(z-0.004*t)+1).^pa;

% The following velocity varies very quickly
c = @(t) -1+0.2*sin(4*t)+0.2*sin(2*t)+1.5*sin(0.5*t)+20*sin(0.005*t);

% % I can use a velocity that slowly oscillates:
% c = @(t) 0.8*sin(0.1*t)+0.1*sin(t)+0.0005*sin(10*t); % 0.0001 0.01 0.0005
Iapp = @(z,t) s*(cos(z-c(t))+1).^pa;

% % Here I make some jumps with the above chosen velocity. 
% % For the local model I multiply by 3. For the global by 8. 
% Iapp = @(z,t) s*(cos(z-c(t))+1).^pa.*(t<=6) +(s*(cos(2+z-c(t-5))+1).^pa).*(t>=7&t<=12)+s*(cos(z-c(t))+1).^pa.*(t>=16);
% Iapp = @(z,t) s*(cos(z-c(t))+1).^pa.*(t<=10);

% 
% c = @(t) 10*exp(-t/20);
% Iapp = @(z,t) s*(cos(z/2+(-c(t))+1)).^pa; %come sopra


% Now this function evaluated at x and 0 becomes the initial condition,
% and I make the red line depend on the maximum values of x. So I only need
% to change it here.
% If I want the system with no input I just write @(z,t) 0 in the next line. but
% then I needto give an initial condition only depending on x.


u0 = @(z) 0.5*(cos(z)+1).^3; %usual initial condition.
% u0 = @(z) 8*(cos(z-1.5)+1).^2;

%% Global model
[fig,usol,T] = driverinput(p,n,Iapp, u0);
% [fig,usol,T] = driverfuncIAppcos(p,n,@(z,t) 0, u0);

%% Local model
% [fig,usol,T] = driverlocalfig(p,n,Iapp,u0);

%what I want to do now is to reproduce a PVA andthe maximum value of the 
% input function. 

z = linspace(-pi,pi,10000)';
T2 =  linspace(T(1),T(end),10000);
[m,zmax] = max(Iapp(z,T2));
pva = angle(usol(:,1:16)*exp(i*x(1:16,:)));

figure()
plot(T2,wrapToPi(z(zmax)),'.');
axis([T(1) T(end) -pi pi])
%axis tight;
hold on;
plot(linspace(T(1),T(end),length(pva)),pva,'.'); 
%axis tight;
hold on;
legend('input','PVA')
pbaspect([2 1 1])

xlabel('t');ylabel('\theta');zlabel('u');
