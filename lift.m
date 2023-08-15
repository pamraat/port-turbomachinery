function [cl, eps] = lift(beta2, gamma)
    i0 = 0;
    is = 18*pi/180;
    clmax = 1.4;
    eps = 0.04;
    
    if(abs(beta2 - gamma) > 1.25*is)
        cl = 0;
    else
        cl = clmax*sin(pi/2*(gamma - beta2 + i0)/(is - i0));
    end
end