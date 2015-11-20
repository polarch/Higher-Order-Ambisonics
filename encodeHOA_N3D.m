function hoasig = encodeHOA_N3D(order, signals, src_directions)
%ENCODEHOA_N3D Encode a number of sound sources in HOA signals.
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
