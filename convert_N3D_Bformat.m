function outsig = convert_N3D_Bformat(sig, type)
%CONVERT_N3D_BFORMAT Converts between 1st order N3D/ACN signals to B-format
% CONVERT_N3D_BFORMAT Converts between 1st-order orthonormalized (N3D) 
% signals, with ACN channel ordering, and traditional B-format.
%
% Inputs:
%   insig:  The N3D/ACN or B-format signals
%   type:   'b2n' to convert from B-format to N3D/ACN, or 'n2b' to do the
%           opposite
%
% Outputs:
%   outsig:  the converted signals
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Archontis Politis, 15/11/2015
%   archontis.politis@aalto.fi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch lower(type)
    case 'b2n'
        BFsig = sig;
        % transformation matrix from B-format to N3D (without the sqrt(2)
        % convention, just pure monopole and dipoles)
        T_BFtoN3D =  [sqrt(1/4/pi)  0               0               0;
            0             0               sqrt(3/4/pi)    0;
            0             0               0               sqrt(3/4/pi);
            0             sqrt(3/4/pi)    0               0];
        
        % remove sqrt(2) scaling convention (assuming it applied to dipoles)
        BFsig(:,2:4) = BFsig(:,2:4)/sqrt(2);
        
        % apply to signals (assuming column signals)
        N3Dsig = BFsig * T_BFtoN3D.';
        outsig = N3Dsig;
        
    case 'n2b'
        N3Dsig = sig;
        % transformation matrix from N3D to B-format (without the sqrt(2)
        % convention, just pure monopole and dipoles)
        T_N3DtoBF =  [sqrt(4*pi)    0               0               0;
            0             0               0               sqrt(4*pi/3);
            0             sqrt(4*pi/3)    0               0;
            0             0               sqrt(4*pi/3)    0];
        
        % apply to signals (assuming column signals)
        BFsig = N3Dsig(:,1:4) * T_N3DtoBF.';
        
        % apply sqrt(2) scaling convention (at dipoles)
        BFsig(:,2:4) = sqrt(2)*BFsig(:,2:4);
        outsig = BFsig;        
end
