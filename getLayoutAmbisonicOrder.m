function order = getLayoutAmbisonicOrder(ls_dirs)
%GETLAYOUTAMBIORDER Returns the equivalent HOA order of the speaker layout
% GETLAYOUTAMBIORDER computes an equivalent ambisonic decoding order of a 
% speaker layout, based on mean spread across the layout. The relation is
% given by Zotter & Frank in
%
%   Zotter, F., Frank, M. (2012). All-Round Ambisonic Panning and Decoding. 
%   Journal of the Audio Engineering Society, 60(10), 807:820.
%
% The code requires the VBAP library found in:
%
% <https://github.com/polarch/Vector-Base-Amplitude-Panning>
%
% Inputs:   
%   ls_dirs: speaker directions in [azi1 elev1; azi2 elev2;... ; aziL elevL]
%            convention, in degrees
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

ls_dirs_rad = ls_dirs*pi/180;
ls_groups = sphDelaunay(ls_dirs_rad);

% find center angles of each triplet
Nls = size(ls_dirs,1);
ls_vecs = zeros(Nls, 3);
[ls_vecs(:,1), ls_vecs(:,2), ls_vecs(:,3)] = sph2cart(ls_dirs_rad(:,1), ...
    ls_dirs_rad(:,2), 1);
Ntri = size(ls_groups,1);
center_dirs_rad = zeros(Ntri, 2);
for k=1:Ntri
    triplet_vecs = ls_vecs(ls_groups(k,:), :);
    center_vec = sum(triplet_vecs);
    [center_dirs_rad(k,1), center_dirs_rad(k,2)] = cart2sph(center_vec(1), ...
        center_vec(2), center_vec(3));
end
% VBAP gains for triplet centers
%ls_gain = vbap(center_dirs_rad, ls_dirs_rad, ls_groups);
center_dirs = center_dirs_rad*180/pi;
ls_gain = vbap(center_dirs, findLsTriplets(ls_dirs), invertLsMtx(ls_dirs, findLsTriplets(ls_dirs))).';

% compute energy vectors for each triplet center
rE = (ls_vecs' * ls_gain.^2);
% vbap spread for every triplet
rE_mag = sqrt(sum(rE.^2))';
spread = 2*acos(rE_mag)*180/pi;
% average spread
spread_avg = mean(spread);
% equivalent ambisonic order
order = round(2*(137.9/(spread_avg)) - 1.52);
if order==0, order=1; end

end

