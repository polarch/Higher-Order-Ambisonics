function a_n = getMaxREweights(order)
%GETMAXREWEIGHTS Returns the order weights for encoding or decoding that maximise energy vectors.
% GETMAXREWEIGHTS computes the per-order weights that need to be applied in
% order for a HOA decoder (or encoder) to exhibit a maximum-norm energy
% vector for a certain direction-of-arrival of the sound. The relation is
% based on the one found in
%
%   Zotter, F., Frank, M. (2012). All-Round Ambisonic Panning and Decoding. 
%   Journal of the Audio Engineering Society, 60(10), 807:820.
%
% Inputs:   
%   order:   order of the HOA decoding.
%
% Outputs:
%   a_n:    [(order+1)^2] vector of weights
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

