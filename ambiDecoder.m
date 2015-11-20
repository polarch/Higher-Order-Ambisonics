function [D, order] = ambiDecoder(ls_dirs, method, rE_WEIGHT, order)
%AMIDECODER Summary of this function goes here
%   Detailed explanation goes here

% find speaker triplets
ls_dirs_rad = ls_dirs*pi/180;
ls_groups = sphDelaunay(ls_dirs_rad);
ls_num = size(ls_dirs,1);

if ~exist('order','var')
    % compute ambisonic equivalent order (from Zotter & Frank)
    order = getLayoutAmbisonicOrder(ls_dirs, ls_groups);
end

if ~exist('rE_WEIGHT','var')
    rE_WEIGHT = 1;
end

if rE_WEIGHT
    a_n = getMaxREweights(order);
end

% compute real SH matrix for layout
Y_ls = getRSH(order, ls_dirs);

switch method
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

function order = getLayoutAmbisonicOrder(ls_dirs, ls_groups)

ls_dirs_rad = ls_dirs*pi/180;

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
