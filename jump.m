clear
close all
clc

n = 16; %16

mu = 10;   %10   
h = 3;       %3
a1 = 30;     %30
b1 = 5;      %5
a2 = 5.5;     %5.5
b2 = 1.5;    %0.1


p = [mu,h,a1,b1,a2,b2];

x = linspace(-pi,pi,n+1)';
x = x(1:end-1); 
% % now $x$ has 16 components, distanced by
% deltax = 2*pi/n;

%Iapp = @(z,t) 5*(cos(z)+1).^2;
jumps = [-5/6*pi,-4/6*pi,-3/6*pi,-2/6*pi,2/6*pi,3/6*pi,4/6*pi,5/6*pi]
c = @(t) -1+0.2*sin(4*t)+0.2*sin(2*t)+1.5*sin(0.5*t)+20*sin(0.005*t);
%Iapp = @(z,t) 8*(cos(z-c(t))+1).^2.*(t<=10&t>1)+(5*(cos(z-c(t)-ju)+1).^2).*(t>10);
% Iapp = @(z,t) 8*(cos(z-c(t))+1).^2.*(t>=5&t<10);

% This is to see all possible jumps
for j=8:8
    % HereI added the first 4 randomply permuted jumps
    % Iapp = @(z,t) Iapp(z,t)+8*(cos(z-c(t)-jumps(i))+1).^2.*(t>=5*(j+1)&t<5*(j+2)); 
    Iapp = @(z,t) 8*(cos(z-c(t))+1).^2.*(t<=10)+(8*(cos(z-c(t)-jumps(j))+1).^2).*(t>10);
    u0 = @(z) 0.5*(cos(z)+1).^3;

    [fig,usol,T] = driverjump(p,n,Iapp,u0);
    
%     [m,umax] = max(usol');
    z = linspace(-pi,pi,10000)';
    T2 =  linspace(T(1),T(end),10000);
    [m,zmax] = max(Iapp(z,T2));

    figure()
    pva = angle(usol(:,1:16)*exp(i*x(1:16,:)));
    
   
    plot(T2,wrapToPi(z(zmax)),'.');
    hold on;
    plot(linspace(T(1),T(end),length(pva)),pva,'.','Linewidth',2); 
%     axis([T(1) T(end) -pi pi])
    axis tight;
    
    
    [m,umax] = max(usol');
    plot(linspace(T(1),T(end),length(umax)),wrapToPi(x(umax)),'.','Color','#77AC30','Linewidth',2); 
    legend('input','PVA','max')
    hold off;
    pbaspect([2 1 1])
    xlabel('t');ylabel('\theta');zlabel('u');
end


