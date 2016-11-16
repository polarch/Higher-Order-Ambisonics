function D_allrad = allrad(ls_dirs, order)
%ALLRAD Implements the All-round Ambisonic Decoding of Zotter & Frank.
% ALLRAD implements the all-round ambisonic decoding published by Zotter &
% Frank in
%
%   Zotter, F., Frank, M. (2012). All-Round Ambisonic Panning and Decoding. 
%   Journal of the Audio Engineering Society, 60(10), 807:820.
%
% It is an ambisonic decoder which combines rendering to an ideal uniform
% virtual loudspeaker layout with energy preserving properties, then
% rendered to an arbitrary loudspeaker setup by means of vector-base
% amplitude panning (VBAP). The code requires the VBAP library and the
% Spherical Harmonic transform library (for VBAP and t-Designs accordingly)
% found in:
%
% <https://github.com/polarch/Vector-Base-Amplitude-Panning>
% <https://github.com/polarch/Spherical-Harmonic-Transform>
%
% Inputs:   
%   ls_dirs: speaker directions in [azi1 elev1; azi2 elev2;... ; aziL elevL]
%            convention, in degrees
%   order:   order of the HOA decoding matrix. For an irregular speaker
%            layout, a well behaved decoding order can be found by the
%            getLayoutAmbisonicOrder() function.
%
% Outputs:
%   D_allrad: [L x (order+1)^2] decoding matrix
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Archontis Politis, 15/11/2015
%   archontis.politis@aalto.fi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % t-value for the t-design
        %t = 2*order + 1;
        t = 20;
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

