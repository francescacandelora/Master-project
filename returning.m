%% Returning file with a unique plot for population (with S1 and S2)

close all, clear all, clc;

angles = linspace(0.1,2*pi,10);
velocity = linspace(pi/15,6*pi/15,30);

n = 8;

% Define the discretization of space
theta = linspace(-pi,pi,n+1)';
theta = theta(1:end-1); 
% now $\theta$ has 8 components, distanced by
deltatheta = 2*pi/n;

% Time interval
ang = angles(4);
% vel = velocity(30);
vel = 1.26;
% trange = [0:0.02:1.66];
trange = [0:0.02:7];

% Input
% S1 = @(t) vel*(t<=ang/vel);
S1 = @(t) sin(t)+0.5*sin(0.2*(t));

% TB1 and CPU4 population
alpha = 20; %10
beta = 0.01; %1/30
md = 0.1;
M = 10;
param = [alpha,beta,md,M];

mu = 10;
h =1;

f = @(z) 1./(1+exp(-mu.*(z-h)));

dAdt = @(t,A) cpu4(t,A,n,S1,f,param);

th0 = [0*theta];
zL0 = [0*theta];
zR0 = [0*theta];

A0 = [th0;zL0;zR0];

[tsol1,Asol] = ode45(dAdt,trange,A0);

ith = 1:n;
izL = n+1:2*n;
izR = 2*n+1:3*n;

thsol1 = Asol(:,ith);
zLsol1 = Asol(:,izL);
zRsol1 = Asol(:,izR);

theta = [theta;pi];
thsol1 = [thsol1,thsol1(:,1)];

[TH,T] = meshgrid(theta,tsol1); 

% CPU1 population
usurf1 = ((cos(theta-thsol1')+1).^5)';
% usurf = ((cos(theta-S(0.1*T').*T')+1).^5)';
vL1 = zRsol1-circshift(usurf1(:,1:end-1),1,2);
vR1 = zLsol1-circshift(usurf1(:,1:end-1),-1,2);

%% Integral
% We now want to compute S(t) as a result of the steering.
% S(t) = \int_{-\pi}^{\pi} v_R(\theta,t)d\theta-\int_{-\pi}^{\pi} v_L(\theta,t)d\theta
% Since we have discrete constants I can also sum the values od v_R and
% v_L, does not really matter
Steering_left = sum(vL1(end,:))
Steering_right = sum(vR1(end,:))
Steering1 = sum(vR1(end,:))-sum(vL1(end,:))

%% Return
% Here we have no longer S(t) as an input, so what we do is we let the fly
% free to come back to initial position
    
trange2 = [trange(end):0.02:20];
% trange2 = [trange(end):0.02:10];

th0 = thsol1(end,1:end-1);
zL0 = zLsol1(end,:);
zR0 = zRsol1(end,:);

A0 = [th0';zL0';zR0'];

dAdt = @(t,A) cpu4return(t,A,n,vel,f,param);
    
[tsol2,Asol] = ode45(dAdt,trange2,A0);

    ith = 1:n;
    izL = n+1:2*n;
    izR = 2*n+1:3*n;

    thsol2 = Asol(:,ith);
    zLsol2 = Asol(:,izL);
    zRsol2 = Asol(:,izR);
    
    thsol2 = [thsol2,thsol2(:,1)];
    %%

    [TH,T] = meshgrid(theta,tsol2);

    % CPU1 population
    usurf2 = ((cos(theta-thsol2')+1).^5)';
    % usurf = ((cos(theta-S(0.1*T').*T')+1).^5)';
    vL2 = zRsol2-circshift(usurf2(:,1:end-1),1,2);
    vR2 = zLsol2-circshift(usurf2(:,1:end-1),-1,2);
%% Make final plots
tsol = [tsol1;tsol2];
usurf = [usurf1;usurf2];
zRsol = [zRsol1;zRsol2];
zLsol = [zLsol1;zLsol2];
zLsol = [zLsol,zLsol(:,1)];
zRsol = [zRsol,zRsol(:,1)];
vR = [vR1;vR2];
vL = [vL1;vL2];
%%
[TH,T] = meshgrid(theta,tsol); 
% First plot: u. I want to see a bump moving over time.
figure();

% surf(T,TH,usurf);
pcolor(T,TH,usurf);
shading flat;
colormap(flipud(bone))
axis tight;
xlabel('t');ylabel('\theta');zlabel('u');
title('Plot of u(\theta,t)')
pbaspect([2 1 1])
pva = angle(usurf(:,1:9)*exp(i*theta(1:9)));
hold on;
plot(linspace(T(1),T(end),length(pva)),pva,'.','Linewidth',2,'Color','#D95319'); 
%%
% Second plot: zL. I want to see a surface that changes

figure()

surf(T,TH,zLsol);
shading flat;
colormap(flipud(bone))
axis tight;
xlabel('t');ylabel('\theta');zlabel('z_L');
title('Plot of z_L(\theta,t)')

% Third plot: zR. I want to see a surface that changes
figure()
    
surf(T,TH,zRsol);
shading flat;
colormap(flipud(bone))
axis tight;
xlabel('t');ylabel('\theta');zlabel('z_R');
title('Plot of z_R(\theta,t)')

% Fourth plot: vL. I'll first plot a surface but then I am only interested
% into the last moment in time for this. 
[X,T] = meshgrid([8:16],tsol); 
figure()
surf(T,X,[vL,vL(:,1)]);
shading flat;
colormap(flipud(bone))
axis tight;
xlabel('t');ylabel('\theta');zlabel('v_L');
title('Plot of v_L(\theta,t)')

% Fifth plot: vR. I'll first plot a surface but then I am only interested
% into the last moment in time for this. 
[X,T] = meshgrid([2:10],tsol); 
figure()
surf(T,X,[vR,vR(:,1)]);
shading flat;
colormap(flipud(bone))
axis tight;
xlabel('t');ylabel('\theta');zlabel('v_R');
title('Plot of v_R(\theta,t)')

Steering = sum(vR(end,:))-sum(vL(end,:))
stvec = sum(vR,2)-sum(vL,2);
figure();plot(tsol,stvec,'LineWidth',1.5)
figure();
plot(tsol,tanh(stvec/100),'LineWidth',1.5)

%% Bar charts t=0 

% zR
figure()
subplot(2,2,1)
bar(1:8,zRsol1(1,:),'FaceColor','#D95319')
ylim([-40 120])
title('zR(\theta,t=0)')

% ZL

subplot(2,2,2)
bar(9:16,zLsol1(1,:),'FaceColor','#EDB120')
ylim([-40 120])
title('zL(\theta,t=0)')

% vL
subplot(2,2,3)
bar(8:15,vL1(1,:))
ylim([-40 120])
title('vL(\theta,t=0)')

% vR
subplot(2,2,4)
bar(2:9,vR1(1,:),'FaceColor',	'#4DBEEE')
ylim([-40 120])
title('vR(\theta,t=0)')

%% Bar charts when sitmulus removed

% zR
figure()
subplot(2,2,1)
bar(1:8,zRsol1(end,:),'FaceColor','#D95319')
ylim([-40 120])
title('zR(\theta,t=7)')

% ZL

subplot(2,2,2)
bar(9:16,zLsol1(end,:),'FaceColor','#EDB120')
ylim([-40 120])
title('zL(\theta,t=7)')

% vL
subplot(2,2,3)
bar(8:15,vL1(end,:))
ylim([-40 120])
title('vL(\theta,t=7)')

% vR
subplot(2,2,4)
bar(2:9,vR1(end,:),'FaceColor',	'#4DBEEE')
ylim([-40 120])
title('vR(\theta,t=7)')

%% Bar charts whent=t(end)
t=20;
% zR
figure()
subplot(2,2,1)
bar(1:8,zRsol2(end,:),'FaceColor','#D95319')
ylim([-40 120])
title('zR(\theta,t=20)')

% ZL

subplot(2,2,2)
bar(9:16,zLsol2(end,:),'FaceColor','#EDB120')
ylim([-40 120])
title('zL(\theta,t=20)')

% vL
subplot(2,2,3)
bar(8:15,vL2(end,:))
ylim([-40 120])
title('vL(\theta,t=20)')

% vR
subplot(2,2,4)
bar(2:9,vR2(end,:),'FaceColor',	'#4DBEEE')
ylim([-40 120])
title('vR(\theta,t=20)')