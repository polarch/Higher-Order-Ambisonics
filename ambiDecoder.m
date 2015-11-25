function [D, order] = ambiDecoder(ls_dirs, method, rE_WEIGHT, order)
%AMBIDECODER Returns a HOA decoding matrix for a loudspeaker setup.
% AMBIDECODER computed a HOA decoding matrix for a loudspeaker setup, based
% on one of five implemented decoder designs. The methods are:
%   - Sampling decoder (SAD)
%   - Mode-matching decoder (MMD)
%   - Energy-preserving decoder (EPAD)
%   - All-round ambisonic panning (ALLRAD)
%   - Constant Angular Spread Decoder (CSAD)
%
%   The two first ones are traditional design approaches. The last three
%   ones are recently proposed design methods that are more flexible and
%   more psychoacoustically motivated. For references, check the
%   TEST_AMBI_SCRIPT.m. Additionally, the so-called max-rE weighting can be
%   enabled for all of the above decoders apart from CSAD. The max-rE
%   weighting maximises the norm of the energy vector for all decoding
%   directions, which in ambisonic literature is considered to reduce 
%   localization blur.
%
% Inputs:   
%   ls_dirs: speaker directions in [azi1 elev1; azi2 elev2;... ; aziL elevL]
%            convention, in degrees, for L loudspeakers in the layout
%   method:  {'SAD','MMD','EPAD','ALLRAD','CSAD'} for any one of the 
%            respective decoding methods
%   rE_WEIGHT:  Enable max-rE weighting. If CSAD is chosen, then this does
%               nothing.
%   order:      Ambisonic order of the array. If it is left blank, then an
%               equivalent ambisonic order is computed by
%               getLayoutAmbisonicOrder().
%
% Outputs:
%   D:      The [L x (order+1)^2] ambisonic decoding matrix
%   order:  The order of the layout, useful to see what it is if it is not 
%           defined and it is computed internally.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Archontis Politis, 15/11/2015
%   archontis.politis@aalto.fi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% find speaker triplets
ls_num = size(ls_dirs,1);

if ~exist('order','var')
    % compute ambisonic equivalent order (from Zotter & Frank)
    order = getLayoutAmbisonicOrder(ls_dirs);
end

if ~exist('rE_WEIGHT','var')
    rE_WEIGHT = 1;
end

if rE_WEIGHT
    a_n = getMaxREweights(order);
end

% compute real SH matrix for layout
Y_ls = getRSH(order, ls_dirs);

switch lower(method)
    case 'sad'
        D = (4*pi)/ls_num * Y_ls.';
    case 'mmd'
        D = pinv(Y_ls);
    case 'epad'
        [U,S,V] = svd(Y_ls);
        S_trunc = S(1:(order+1)^2,1:(order+1)^2);
        V_trunc = V(:,1:(order+1)^2);
        D = (4*pi)/ls_num * V_trunc*U.';
    case 'allrad'
        D = allrad(ls_dirs, order);        
    case 'csad'
        D = csad(ls_dirs, order);
end

% apply rE weights if defined
if rE_WEIGHT && ~isequal(method,'csad')
    D = D*diag(a_n);
end

end
