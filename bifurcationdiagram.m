 %% Try to see if bifurcation60.m works
close all

N = 20;
%% This is a global model. 
% The strength goes from 1 to 13 and the width from 0 to pi. 
% So I am going to copmute 100 x 100 values.
% 
strength = linspace(2,18,N);
widths = linspace(0.2,3,N);
ps = [199.5504,computeexponents(widths(2:end))]; %199.5504   30.7910   12.0015    6.3977    4.0054    2.7690    2.0490    1.5941    1.2894    1.0761

%%
for j=1:N
    for k=1:N
        close all
        bifurcation(strength(j),ps(k),widths(k));
    end
end

%%
fileID = fopen('recorder150.txt','r');
readglobal(fileID,N,strength,widths)


