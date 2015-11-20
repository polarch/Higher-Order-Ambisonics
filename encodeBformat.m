function WXYZ = encodeBformat(signals, src_directions)
%ENCODE_BFORMAT Encode a number of sound sources in ideal B-format.
%   ENCODE_BFORMAT encodes a number of N signals, given as a matrix of
%   N signal columns, to 'classic' B-format signals W, X, Y, Z. The
%   src_directions are the azimuth and elevations of each signal and it
%   should be either an Nx2 matrix, specifying the direction of each
%   signal, or a vector [azi elev], in case all the signals are encoded
%   with a single direction.

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

Gwxyz = [gw sqrt(2)*[gx gy gz]];

% encode to B-format
WXYZ = signals * Gwxyz;
    
end
