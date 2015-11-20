function a_n = getMaxREweights(order)
%GETMAXREWEIGHTS Summary of this function goes here
%   Detailed explanation goes here

    a_n = [];
    for n=0:order
        temp = legendre(n, cos(137.9*(pi/180)/(order+1.51)));
        a_n = [a_n; ones(2*n+1,1)*temp(1)];
    end

end

