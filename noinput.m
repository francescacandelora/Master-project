%% No input

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

c = @(t) 10*exp(-t/20);
Iapp = @(z,t) 15*(cos(z/2+(-c(t))+1)).^pa; %come sopra
Iapp = @(z,t) 0*Iapp(z,t);


u0 = @(z) 0.5*(cos(z)+1).^3; %usual initial condition.

%% Global model
[fig,usol,T] = drivernoinput(p,n,Iapp, u0);


