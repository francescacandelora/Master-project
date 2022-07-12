function readglobal(fileID,N,strength,widths)
    %% Read global file
    formatSpec = '%f %f %f';
    A = fscanf(fileID,formatSpec,[3 N^2])';

    % A(i,1) and A(i,2) are exchanged because we want w on the x axis and s on
    % the y axis. 
    % x = [0.2 pi/2];
    % y = [1 13/2];
    C = reshape(A(:,3),N,N);
    % imagesc(x,y,C')
    y = linspace(strength(1),strength(end),N);
    x = linspace(widths(1),widths(end),N);
    pcolor(x,y,C')
    shading flat
    xlabel('width of input');
    ylabel('strength of input');
    caxis([0 2]);
    colormap spring
    text(2.7,16,'Flow');
    text(1.5,2.5,'No effect');
    text(0.5,10,'Jump');
    title('150° jumps')



end