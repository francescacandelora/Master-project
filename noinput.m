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

theta = linspace(-pi,pi,n+1)';
theta = theta(1:end-1); 
% % now $x$ has 16 components, distanced by
% deltax = 2*pi/n;

% Strength of input and/or init cond
s = 8;
% exponent of input and/or init cond
pa = 2;

c = @(t) 10*exp(-t/20);
Iapp = @(z,t) 15*(cos(z/2+(-c(t))+1)).^pa; %come sopra
Iapp = @(z,t) 0*Iapp(z,t);


% u0 = @(z) 0.5*(cos(z-1.5)+1).^3; %usual initial condition.
u0 = @(z) 0*z + 1.5 + 0.5*rand(n,1);


%% Global model
[fig,usol,T] = drivernoinput(p,n,Iapp, u0);

%%
statePlot = figure(); plot(theta,u0(theta)); xlabel('\theta'); ylabel('u0(\theta)');
pause;

for it = 1:length(T)/4
  plot(theta,usol(it,1:end-1)'); xlabel('x'); ylabel('u0(x)');
  %title([' t = ' num2str(t(it))]);
  drawnow;
end
%%
% theta =  [theta;pi];
% figure();
% plot(theta,u0(theta),'LineWidth',2)
% hold on;
% plot(theta,usol(end,:),'LineWidth',2,'Color','#D95319')
% title('Plot of U_0 and U(20)');
% legend({'U_0','U(20)'})