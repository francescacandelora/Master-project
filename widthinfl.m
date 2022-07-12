function w = widthinfl(p)
    % Ths function outputs the FWHM of a function defined as Iapp, given
    % the exponent p, defined on [-pi,pi]
    w=zeros(length(p),1);
    for i=1:length(p)
        Iapp = @(z,t) (cos(z)+1).^(p(i)); 
        syms x
        Iapp1= Iapp(x,0);
        Iapp2=diff(Iapp1,2);
        inflec_pt = sort(double(solve(Iapp2,'MaxDegree',2)));
        w(i)=double(inflec_pt(2)-inflec_pt(1));
    end
end

