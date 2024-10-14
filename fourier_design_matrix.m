function [X] = fourier_design_matrix(n,k,opt)

if (opt == 0)
    
    X=[ones(n,1) linspace(-1,1,n)'];
    
    for j=1:k
        rt=cos((2*pi*j*n^-1).*[1:n])';
        it=sin((2*pi*j*n^-1).*[1:n])';
        X=[X rt it];
    end
    
else
    
    t=linspace(0,1,n);
    X=[];
    for j=1:k
        rt=0.5-0.5.*cos(2*pi*j.*t)';
        it=sin(2*pi*j.*t)';
        X=[X rt it];
    end
    
end