function outsig = convert_N3D_FuMa(sig, type)
%CONVERT_N3D_FuMa Converts between 1st-3rd order N3D/ACN signals to FuMa
% CONVERT_N3D_BFORMAT Converts between up-to 3rd-order orthonormalized (N3D) 
% signals, with ACN channel ordering, and FuMa normalization and ordering.
%
% Inputs:
%   insig:  The N3D/ACN or FuMa signals
%   type:   'fuma2n' to convert from FuMa to N3D/ACN, or 'n2fuma' to do the
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

% map from FuMa to ACN
map_fuma2acn = [1  3  4  2  9  7  5  6  8  16  14  12  10  11  13  15];
%               W  Y  Z  X  V  T  R  S  U   Q   O   M   K   L   N   P

% FuMa to N3D normalizations (after channel mapping to ACN)
norm_fuma2n3d = [   sqrt(2)      % W
                    sqrt(3)      % Y
                    sqrt(3)      % Z
                    sqrt(3)      % X
                    sqrt(15)/2   % V
                    sqrt(15)/2   % T
                    sqrt(5)      % R
                    sqrt(15)/2   % S
                    sqrt(15)/2   % U
                    sqrt(35/8)   % Q
                    sqrt(35)/3   % O
                    sqrt(224/45) % M
                    sqrt(7)      % K
                    sqrt(224/45) % L
                    sqrt(35)/3   % N
                    sqrt(35/8)]; % P

% invert mapping from ACN to FuMa
[~,map_acn2fuma] = sort(map_fuma2acn);

% N3D to FuMa normalizations (after channel mapping to FuMa)
norm_n3d2fuma = 1./norm_fuma2n3d;
norm_n3d2fuma = norm_n3d2fuma(map_acn2fuma);

% limit operations up to 3rd-order if input signals or of higher order
order = sqrt(size(sig,2))-1;
if order>3
    sig = sig(1:16,:);
    order = 3;
    warning('the input signals are of HOA order higher than 3, limiting to order 3 for FuMa convention')
end
nSH = (order+1)^2;
                
switch lower(type)
    case 'fuma2n'
        FUMAsig = sig;
        % re-order channels to ACN
        N3Dsig = FUMAsig(:,map_fuma2acn(1:nSH));
        % re-normalize to N3D
        N3Dsig = N3Dsig*diag(norm_fuma2n3d(1:nSH));
        outsig = N3Dsig;        
        
    case 'n2fuma'
        N3Dsig = sig;
        % re-order channels to FuMa
        FUMAsig = N3Dsig(:,map_acn2fuma(1:nSH));
        % re-normalize to FuMa
        FUMAsig = FUMAsig*diag(norm_n3d2fuma(1:nSH));
        outsig = FUMAsig;      
end
