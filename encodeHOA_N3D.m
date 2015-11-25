function hoasig = encodeHOA_N3D(order, signals, src_directions)
%ENCODEHOA_N3D Encode a number of sound sources in HOA signals.
% ENCODEHOA_N3D encodes a number of signals coming from certain directions,
% to ideal HOA signals. Essentially, this corresponds to the signals
% multiplied with the spherical harmonic values for the source directions,
% up to a specified order. Orthonormalized (N3D) spherical harmonics are
% assumed.
%
% Inputs:   
%   order:      maximum order for the HOA encoding.
%   signals:    matrix of [L x K] signal values, where L is the length of
%               the signals and K is the number of them.
%   src_dirs:   source directions in [azi1 elev1; azi2 elev2;... ; aziK elevK]
%               convention, in degrees
%
% Outputs:
%   hoasig:     [L x (order+1)^2] HOA signals
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Archontis Politis, 15/11/2015
%   archontis.politis@aalto.fi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% create the encoding matrix
E = getRSH(order, src_directions);

% encode to HOA signals
hoasig = signals * E.';
    
end
