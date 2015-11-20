function a_n = getMaxREweights(order)
%GETMAXREWEIGHTS Returns the order weights for encoding or decoding that maximise energy vectors.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Archontis Politis, 15/11/2015
%   archontis.politis@aalto.fi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    a_n = [];
    for n=0:order
        temp = legendre(n, cos(137.9*(pi/180)/(order+1.51)));
        a_n = [a_n; ones(2*n+1,1)*temp(1)];
    end

end

