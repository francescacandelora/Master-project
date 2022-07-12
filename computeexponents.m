function p = computeexponents(widths)
    p = [];
    x0 = 1.1;
    for i=length(widths):-1:1
        fun = @(p) width(p,widths(i));
        p = [fsolve(fun,x0),p];
        x0 = p(1);
    end
end