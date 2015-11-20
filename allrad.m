function D_allrad = allrad(ls_dirs, order)
%ALLRAD Summary of this function goes here
%   Detailed explanation goes here

        % t-value for the t-design
        t = 2*order + 1;
        % vbap gains for selected t-design
        [~, t_dirs_rad] = getTdesign(t);
        t_dirs = t_dirs_rad*180/pi;
        G_td = vbap(t_dirs , findLsTriplets(ls_dirs), invertLsMtx(ls_dirs, findLsTriplets(ls_dirs))).';
        
        % spherical harmonic matrix for t-design
        % convert to [azimuth zenith] for SH convention
        Y_td = getRSH(order, t_dirs).';
        
        % allrad decoder
        Ntd = size(t_dirs_rad,1);
        D_allrad = 4*pi/Ntd * G_td * Y_td;

end

