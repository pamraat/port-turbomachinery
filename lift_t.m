function [cl, eps] = lift_t(beta4, gamma)
    i0 = 0;
    is = 18*pi/180;
    clmax = 1.2;
    eps = 0.04;
    if(abs(beta4 - gamma) > 1.25*is)
        cl = 0;
    else
        cl = clmax*sin(pi/2*(beta4 - gamma + i0)/(is - i0));
    end
    if (cl<0)
        cl = 0;
    end
end