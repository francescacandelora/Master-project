%% Returning but with u neural field instead of travelling wave.
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
m = 10;
param = [alpha,beta,md,m];

a = 2;
c = 0.5;

mu = 10;
h =1;

f = @(z) 1./(1+exp(-mu.*(z-h)));

w = @(z) a*exp(-abs(z))-c;
w_ex = @(z) w(z+2*pi).*(heaviside(-z-pi))+w(z).*(heaviside(z+pi)-heaviside(z-pi))+w(z-2*pi).*heaviside(z-pi);
M = zeros(n,n);
for j = 1:n
    for k = 1:n
        M(j,k) = w_ex(theta(j)-theta(k))*deltatheta;
    end
end

dAdt = @(t,A) cpu4u(t,A,n,S1,f,M,param);

th0 = [0*theta];
u0  = 0.5*(cos(theta)+1).^3;
zL0 = [0*theta];
zR0 = [0*theta];

A0 = [th0;u0;zL0;zR0];

[tsol1,Asol] = ode45(dAdt,trange,A0);

ith = 1:n;
iu  = n+1:2*n;
izL = 2*n+1:3*n;
izR = 3*n+1:4*n;

thsol1 = Asol(:,ith);
usol1  = Asol(:,iu);
zLsol1 = Asol(:,izL);
zRsol1 = Asol(:,izR);

theta = [theta;pi];
thsol1 = [thsol1,thsol1(:,1)];

[TH,T] = meshgrid(theta,tsol1); 

usol1 = [usol1,usol1(:,1)];
% CPU1 population
vL1 = zRsol1-circshift(usol1(:,1:end-1),1,2);
vR1 = zLsol1-circshift(usol1(:,1:end-1),-1,2);

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
u0  = usol1(end,1:end-1);
zL0 = zLsol1(end,:);
zR0 = zRsol1(end,:);

A0 = [th0';u0';zL0';zR0'];

dAdt = @(t,A) cpu4returnu(t,A,n,vel,f,M,param);
    
[tsol2,Asol] = ode45(dAdt,trange2,A0);

ith = 1:n;
iu  = n+1:2*n;
izL = 2*n+1:3*n;
izR = 3*n+1:4*n;

thsol2 = Asol(:,ith);
usol2  = Asol(:,iu);
zLsol2 = Asol(:,izL);
zRsol2 = Asol(:,izR);
    
thsol2 = [thsol2,thsol2(:,1)];
    %%

[TH,T] = meshgrid(theta,tsol2);

usol2 = [usol2,usol2(:,1)];
% CPU1 population
vL2 = zRsol2-circshift(usol2(:,1:end-1),1,2);
vR2 = zLsol2-circshift(usol2(:,1:end-1),-1,2);

%% Make final plots
tsol = [tsol1;tsol2];
usol = [usol1;usol2];
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
pcolor(T,TH,usol);
shading flat;
colormap(flipud(bone))
axis tight;
xlabel('t');ylabel('\theta');zlabel('u');
title('Plot of u(\theta,t)')
pbaspect([2 1 1])
pva = angle(usol(:,1:9)*exp(i*theta(1:9)));
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