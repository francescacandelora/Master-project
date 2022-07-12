function [result,umax] = bifurcation(s,p,w)
    % j is going to be fixed. 
    % This function creates an input with strength s, width w (given by p)
    % and makes a jump of j

    
    n = 16; %16

    mu = 10;   %10   
    h = 3;       %3
    a1 = 30;     %30
    b1 = 5;      %5
    a2 = 5.5;     %5.5
    b2 = 1.5;    

    par = [mu,h,a1,b1,a2,b2];

    x = linspace(-pi,pi,n+1)';
    x = x(1:end-1); 
    deltax = 2*pi/n;

    c = @(t) 0.001*sin(0.1*t)+0.01*sin(t)+0.0005*sin(10*t);
    %% HERE 1/2
    % 5/6*pi = 150°
    jump = 5*pi/6;
    Iapp = @(z,t) (cos(z-c(t))+1).^p.*(t<=8)+ (cos(z-c(t)-jump)+1).^p.*(t>8); 
    % Now I normalize Iapp, so that its maximum is 1, and I multiply it by
    % s.
    Iapp = @(z,t) s*Iapp(z,t)/Iapp(0,0); 

    u0 = @(z) 0.5*(cos(z)+1).^3; %usual initial condition.
    [usol,T] = driverjumpbif(par,n,Iapp,u0);
%     [fileID,usol,T] = driverjumpfig(par,n,Iapp,u0);
    

    % Draw the graph to see if the solution jumps. 
    
    [m,umax] = max(usol');
%     pva = angle(usol(:,1:16)*exp(i*x(1:16,:)));
%     result = max(abs(diff(pva)));
    result = max(abs(diff(x(umax))));
    result = min(result/deltax,2);
    
%     figure()
%     %subplot(4,1,1)
%     T2 = linspace(T(1),T(end),length(umax))
%     plot(linspace(T(1),T(end),length(umax)),wrapToPi(x(umax)),'.'); 
%     axis([T(1) T(end) -pi pi])
%     %axis tight;
%     hold on;
%     plot(T2,wrapToPi(z(zmax)),'.');
%     %axis tight;
%     hold on;
%     pbaspect([4 1 1])
%     title("s="+s+"width="+w)
    
    %% HERE 2/2
    fileID = fopen('recorder150.txt','a'); 
    % Using 'w' rewrites, using 'a' appends.
    formatSpec = '%3.2f %3.2f %f\n';
    
    fprintf(fileID,formatSpec,s,w,result);
    fclose(fileID);
end

