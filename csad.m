function D_csad = csad(ls_dirs, order)
%CSAD Implements the Constant Angular Spread HOA decoding of Epain, Jin & Zotter.
% CSAD implements the Constant Angular Spread HOA decoding published by
% Epain, Jin & Zotter in
%
%   Epain, N., Jin, C. T., Zotter, F. (2014). Ambisonic Decoding With Constant 
%   Angular Spread. Acta Acustica United with Acustica, 100, 928:936.
%
% This ambisonic decoder is energy preserving, has very low directional error
% in terms of energy vectors, and constant angular spread (constant magnitude 
% of energy vectors). The code requires the VBAP/VBIP library and the
% Spherical Harmonic transform library (for VBAP/VBIP and t-Designs accordingly)
% found in:
%
% <https://github.com/polarch/Vector-Base-Amplitude-Panning>
% <https://github.com/polarch/Spherical-Harmonic-Transform>
%
% Note that the implementation here is a slightly simplified version of the
% one presented in the reference. The difference is on the way the
% spread sources are generated: in the reference a dense grid of spread sources
% are generated around each evaluation point, they are then weighted 
% radially with an angular window, that is adjusted to meet the spread
% constraint. While in this code, a fixed number of symmetric sources are
% generated around the evaluation direction and their angular distance is
% adjusted till the spread constraint is met. Even though the first
% approach is more elaborate, the simpler approach here seems to deliver
% results almost cloe to the ones in the reference. For more information,
% see the publication above.
%
% Inputs:   
%   ls_dirs: speaker directions in [azi1 elev1; azi2 elev2;... ; aziL elevL]
%            convention, in degrees
%   order:   order of the HOA decoding matrix. For an irregular speaker
%            layout, an well'beahved decoding order can be found by the
%            getLayoutAmbisonicOrder() function.
%
% Outputs:
%   D_csad:  [L x (order+1)^2] decoding matrix
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Archontis Politis, 15/11/2015
%   archontis.politis@aalto.fi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ls_dirs_rad = ls_dirs*pi/180;
[U_ls(:,1), U_ls(:,2), U_ls(:,3)] = sph2cart(ls_dirs_rad(:,1), ls_dirs_rad(:,2), 1);
% get a dense grid of evaluation directions
[~, tdirs_rad] = getTdesign(4*order);
tdirs = tdirs_rad*180/pi;
% Compute a dense VBIP gain table, with 5deg resolution, compute energy 
% vector norms, and from that find the maximum spread value
Gvbip = getGainTable(ls_dirs, [5 5], 0, 'vbip');
E = sum(Gvbip.^2,2);
rE = (Gvbip.^2 * U_ls) ./ (E*ones(1,3));
rE_mag = sqrt(sum(rE.^2,2));
spread = 2*acos(2*rE_mag-1);
max_spread = max(spread);
% get triangulation and loudspeaker matrix inversions to reuse on the
% optimization
ls_groups = findLsTriplets(ls_dirs);
ls_invMtx = invertLsMtx(ls_dirs,ls_groups);

ls_num = size(ls_dirs,1);
vs_num = size(tdirs,1);
% initialize MDIP gain table
Gmdip = zeros(ls_num, vs_num);
% define the optimization function, and run for all evaluation directions
for nd = 1:size(tdirs,1);
    disp(['Solving ' num2str(nd) '/' num2str(size(tdirs,1)) ' panning directions'])
    dir = tdirs(nd,:);
    function2minimise = @(x) findVbipGains(dir, ls_groups, ls_invMtx, x, U_ls) - max_spread*180/pi;
    
    opti_spread(nd) = fzero(function2minimise,max_spread*180/pi);
    Gmdip(:,nd) = vbip(dir, ls_groups, ls_invMtx, opti_spread(nd), 32, 4);
end

% approximate the MDIP table with a pseudo-inverse - that's the CSAD
% decoding matrix
Y_td = getRSH(order, tdirs);
D_csad = Gmdip * Y_td' * inv(Y_td*Y_td');

end

% function to return VBIP gains for a certain spread value
function spread_deg = findVbipGains(dir, ls_groups, ls_invMtx, target_spread_deg, U_ls)

G = vbip(dir, ls_groups, ls_invMtx, target_spread_deg, 32, 4);
G2 = G.^2;
E = sum(G2,2);
rE = (G2 * U_ls) ./ (E*ones(1,3));
rE_mag = sqrt(sum(rE.^2,2));
spread_deg = 2*acos(2*rE_mag-1)*180/pi;

end

