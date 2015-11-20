function Dmdip = csad(ls_dirs, N)
%CSAD Implements the Constant angular spread HOA decoding of Epain, Jin & Zotter.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Archontis Politis, 15/11/2015
%   archontis.politis@aalto.fi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ls_dirs_rad = ls_dirs*pi/180;
[U_ls(:,1), U_ls(:,2), U_ls(:,3)] = sph2cart(ls_dirs_rad(:,1), ls_dirs_rad(:,2), 1);

[~, tdirs_rad] = getTdesign(4*N);
tdirs = tdirs_rad*180/pi;

ls_groups = findLsTriplets(ls_dirs);
Gvbip = getGainTable(ls_dirs, [5 5], 0, 'vbip');
E = sum(Gvbip.^2,2);
rE = (Gvbip.^2 * U_ls) ./ (E*ones(1,3));
rE_mag = sqrt(sum(rE.^2,2));
spread = 2*acos(2*rE_mag-1);
max_spread = max(spread);

ls_groups = findLsTriplets(ls_dirs);
ls_invMtx = invertLsMtx(ls_dirs,ls_groups);

ls_num = size(ls_dirs,1);
vs_num = size(tdirs,1);
Gmdip = zeros(ls_num, vs_num);
for nd = 1:size(tdirs,1);
%    disp(['Solving ' num2str(nd) '/' num2str(size(tdirs,1)) ' panning directions'])
    dir = tdirs(nd,:);
    function2minimise = @(x) findVbipGains(dir, ls_groups, ls_invMtx, x, U_ls) - max_spread*180/pi;
    
    opti_spread(nd) = fzero(function2minimise,max_spread*180/pi);
    Gmdip(:,nd) = vbip(dir, ls_groups, ls_invMtx, opti_spread(nd), 32, 4);
end

Y_td = getRSH(N, tdirs);
Dmdip = Gmdip * Y_td' * inv(Y_td*Y_td');

end

function spread_deg = findVbipGains(dir, ls_groups, ls_invMtx, target_spread_deg, U_ls)

G = vbip(dir, ls_groups, ls_invMtx, target_spread_deg, 32, 4);
G2 = G.^2;
E = sum(G2,2);
rE = (G2 * U_ls) ./ (E*ones(1,3));
rE_mag = sqrt(sum(rE.^2,2));
spread_deg = 2*acos(2*rE_mag-1)*180/pi;

end

