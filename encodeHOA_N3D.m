function hoasig = encodeHOA_N3D(order, signals, src_directions)
%ENCODE_BFORMAT Encode a number of sound sources in ideal B-format.
%   ENCODE_BFORMAT encodes a number of N signals, given as a matrix of
%   N signal columns, to 'classic' B-format signals W, X, Y, Z. The
%   src_directions are the azimuth and elevations of each signal and it
%   should be either an Nx2 matrix, specifying the direction of each
%   signal, or a vector [azi elev], in case all the signals are encoded
%   with a single direction.

% create the encoding matrix
E = getRSH(order, src_directions);

% encode to HOA signals
hoasig = signals * E.';
    
end
