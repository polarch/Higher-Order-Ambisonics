function WXYZ = encodeBformat(signals, src_directions)
%ENCODE_BFORMAT Encode a number of sound sources in ideal B-format.
% ENCODE_BFORMAT encodes a number of signals coming from certain directions,
% to ideal B-format signals. The traditional sqrt(2) factor is assumed to
% be on the dipoles.
%
% Inputs:   
%   signals:    matrix of [L x K] signal values, where L is the length of
%               the signals and K is the number of them.
%   src_dirs:   source directions in [azi1 elev1; azi2 elev2;... ; aziK elevK]
%               convention, in degrees
%
% Outputs:
%   WXYZ:       [L x 4] B-format signals
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Archontis Politis, 15/11/2015
%   archontis.politis@aalto.fi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% number of signals
N_src = size(signals, 2);

% azimuth and elevation converted to radians
src_az = src_directions(:, 1)*pi/180;
src_el = src_directions(:, 2)*pi/180;

% create the b-format encoding matrix
gw = ones(N_src,1);
gx = cos(src_az).*cos(src_el);
gy = sin(src_az).*cos(src_el);
gz = sin(src_el);

Gwxyz = [gw/sqrt(2) gx gy gz];

% encode to B-format
WXYZ = signals * Gwxyz;
    
end
