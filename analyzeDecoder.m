function [E, A, rV_mag, rE_mag, rV_ang, rE_ang] = analyzeDecoder(Mtx, ls_dirs, type, angRes, PLOT_ON, INFO_ON)
%ANALYZEDECODER Analyzes energy and velocity/energy vectors for a decoder.
% ANALYZEDECODER evaluates total energy and amplitude gain for a grid of 
% directions, as well as energy and velocity vectors, and spread, for an 
% ambisonic decoder or for a panning function in general. If the type
% specified is 'decoder', an ambisonic decoding matrix should be passed to
% the function, and then the loudspeaker gains for a grid of directions,
% specified by angRes, are computed. If the type specified is 'panner', then 
% a matrix of loudspeaker gains should be passed, that correspond to the
% grid directions specified by angRes. In both cases, the results can be
% plotted or returned as outputs. If an additional option INFO_ON is
% enabled, then some extra information about the deviations of
% amplitude, energy, energy/velocity vector magnitudes, directional
% errors and perceptual spread, are printed in the console.
%
% Inputs:   
%   Mtx:     The [L x (order+1)^2] ambisonic decoding matrix, for
%            type='decoder', or the [(180/elevRes+1) x (360/aziRes+1) x L]
%            matrix of loudspeaker gains for a grid of 
%            angRes = [aziRes elevRes]. See the TEST_AMBI_SCRIPT.m for an
%            example on panning analysis using VBAP/VBIP.
%   ls_dirs: speaker directions in [azi1 elev1; azi2 elev2;... ; aziL elevL]
%            convention, in degrees, for L loudspeakers in the layout
%   type:    'decoder' for analysis of ambisonic decoding, or 'panner' for
%            analysis of panning functions evaluated on a grid, see also
%            the Mtx input variable. If not defined, the default is
%            'decoder'
%   angRes:  The resolution of the regular grid of direcitons to evaluate 
%            the metrics, with angRes = [azi_resolution elev_resolution] in
%            degrees. In case of type='panner', the resolution should match
%            the dimensions of Mtx, see the description of the Mtx variable
%            above. In the case of type = 'decoder', the resolution can be 
%            anything. If not defined, default is angRes = [5 5].
%   PLOT_ON: (0,1) Plot the various metrics. If not defined, default is 0.
%   INFO_ON: Displays ranges of the various metrics and spread values,
%            computed from the magnitude of the energy vector. In case of
%            type='decoder', some additional information is displayed, such
%            as equivalent ambisonic order, and theoretical spread. If not 
%            defined, default is 0.
%
% Outputs:
%   E:      Energy as the sum of squared gains for each grid point.
%   A:      Amplitude as sum of gains for each grid point.
%   rV_mag: Magnitude of the velocity vector, at each grid point.
%   rE_mag: Magnitude of the energy vector, at each grid point.
%   rV_ang: Directional error in degrees, between the direction of the grid 
%           point and the direction of the velocity vector for that grid point.
%   rE_ang: Directional error in degrees, between the direction of the grid 
%           point and the direction of the energy vector for that grid point.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Archontis Politis, 15/11/2015
%   archontis.politis@aalto.fi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<3, type = 'decoder'; angRes = [5 5]; PLOT_ON = 0; INFO_ON = 0;
elseif nargin<4, angRes = [5 5]; PLOT_ON = 0; INFO_ON = 0;
elseif nargin<5, PLOT_ON = 0; INFO_ON = 0;
elseif nargin<6, INFO_ON = 0;
end

ls_dirs_rad = ls_dirs*pi/180;
[U_ls(:,1), U_ls(:,2), U_ls(:,3)] = sph2cart(ls_dirs_rad(:,1), ls_dirs_rad(:,2), 1);

% regular grid of directions
aziRes = angRes(1);
polarRes = angRes(2);
azi = -180:aziRes:180;
elev = -90:polarRes:90;
[Azi, Elev] = meshgrid(azi, elev);
src_dirs_rad = [Azi(:) Elev(:)]*pi/180;
aziElev2aziIncl = @(dirs) [dirs(:,1) pi/2-dirs(:,2)];
[U_grid(:,1), U_grid(:,2), U_grid(:,3)] = sph2cart(src_dirs_rad(:,1), src_dirs_rad(:,2), 1);

switch type
    case 'decoder'
        M_dec = Mtx;
        
        order = sqrt(size(M_dec,2))-1;
        Y = getSH(order, aziElev2aziIncl(src_dirs_rad), 'real');
        
        G = Y * M_dec.';
        
    case 'panner'
        
        G = reshape(Mtx, [size(Mtx,1)*size(Mtx,2) size(Mtx,3)]);
end

A = sum(G,2);
rV = (G * U_ls) ./ (A*ones(1,3));
E = sum(G.^2,2);
rE = (G.^2 * U_ls) ./ (E*ones(1,3));

rV_mag = sqrt(sum(rV.^2,2));
rE_mag = sqrt(sum(rE.^2,2));
[rV_dirs_rad(:,1), rV_dirs_rad(:,2)] = cart2sph(rV(:,1), rV(:,2), rV(:,3));
[rE_dirs_rad(:,1), rE_dirs_rad(:,2)] = cart2sph(rE(:,1), rE(:,2), rE(:,3));
U_rV = rV./(sqrt(sum(rV.^2,2))*ones(1,3));
U_rE = rE./(sqrt(sum(rE.^2,2))*ones(1,3));
dot_rV = sum(U_rV.*U_grid,2);
dot_rV(dot_rV>1) = 1; dot_rV(dot_rV<-1) = -1;
dot_rE = sum(U_rE.*U_grid,2);
dot_rE(dot_rE>1) = 1; dot_rE(dot_rE<-1) = -1;
rV_ang = acos(dot_rV)*180/pi;
rE_ang = acos(dot_rE)*180/pi;

spread_rE_daniel = 2*acos(rE_mag)*180/pi;
spread_rE_frank = 186.4*(1-rE_mag)+10.7;
spread_rE_epain = 2*acos(2*rE_mag-1)*180/pi;

if INFO_ON
    disp(' ')
    disp(['Amplitude range: ' num2str(20*log10(min(A))) ' ~ ' num2str(20*log10(max(A))) ' dB'])
    disp(['Energy range: ' num2str(10*log10(min(E))) ' ~ ' num2str(10*log10(max(E))) ' dB'])
    disp(['Directional error range of rV: ' num2str(min(rV_ang)) ' ~ ' num2str(max(rV_ang)) ' deg'])
    disp(['Directional error range of rE: ' num2str(min(rE_ang)) ' ~ ' num2str(max(rE_ang)) ' deg'])
    disp(['Magnitude range of rV: ' num2str(min(rV_mag)) ' ~ ' num2str(max(rV_mag))])
    disp(['Magnitude range of rE: ' num2str(min(rE_mag)) ' ~ ' num2str(max(rE_mag))])
    disp(['Spread range of rE (Daniel): ' num2str(min(spread_rE_daniel)) ' ~ ' num2str(max(spread_rE_daniel)) ' deg'])
    disp(['Spread range of rE (Frank): ' num2str(min(spread_rE_frank)) ' ~ ' num2str(max(spread_rE_frank)) ' deg'])
    disp(['Spread range of rE (Epain et al.): ' num2str(min(spread_rE_epain)) ' ~ ' num2str(max(spread_rE_epain)) ' deg'])
    if isequal(type, 'decoder')
        Neq = getLayoutAmbisonicOrder(ls_dirs);
        disp(['Equivalent ambisonic order of layout: ' num2str(Neq)])
        disp(['Theoretical magnitude of rE for equivalent order: ' num2str(getTheoreticalEVmag(Neq))])
        disp(['Theoretical magnitude of rE for decoding order: ' num2str(getTheoreticalEVmag(order))])
        disp(['Theoretical spread of equivalent order (Daniel): ' num2str(2*acos(getTheoreticalEVmag(Neq))*180/pi)])
        disp(['Theoretical spread of decoding order (Daniel): ' num2str(2*acos(getTheoreticalEVmag(order))*180/pi)])
    end
    disp(' ')
end

A = reshape(A, [(180/polarRes+1) (360/aziRes+1)]);
E = reshape(E, [(180/polarRes+1) (360/aziRes+1)]);
rV_mag = reshape(rV_mag, [(180/polarRes+1) (360/aziRes+1)]);
rE_mag = reshape(rE_mag, [(180/polarRes+1) (360/aziRes+1)]);
rV_ang = reshape(rV_ang, [(180/polarRes+1) (360/aziRes+1)]);
rE_ang = reshape(rE_ang, [(180/polarRes+1) (360/aziRes+1)]);

if PLOT_ON
    
    figure
    subplot(231)
    plotSphericalGrid(A, angRes, ls_dirs, gca);
    title('Amplitude A')
    subplot(232)
    plotSphericalGrid(rV_mag, angRes, ls_dirs, gca);
    title('Velocity vector magnitude ||rV||')
    caxis([0 1])
    subplot(233)
    plotSphericalGrid(rE_mag, angRes, ls_dirs, gca);
    title('Energy vector magnitude ||rE||')
    caxis([0 1])
    
    subplot(234)
    plotSphericalGrid(E, angRes, ls_dirs, gca);
    title('Energy E')
    subplot(235)
    plotSphericalGrid(rV_ang, angRes, ls_dirs, gca);
    title('Velocity vector angular error (deg)')
    subplot(236)
    plotSphericalGrid(rE_ang, angRes, ls_dirs, gca);
    title('Energy vector angular error (deg)')
end

end
