function h_ax = plotSphericalGrid(F, ang_res, ls_dirs, h_ax)
%PLOTSPHERICALGRID Plots spherical data on a regular grid of directions.
% PLOTSPHERICALGRID plots data evaluated on a regular grid of azimuths and
% elevations, with a specific angular resolution.
%
% Inputs:   
%   F:       data on a matrix of size [(180/elev_res+1) x (360/azi_res+1)],
%            where each row is the data values for azimuths from
%            -180:azi_res:180deg, for a single elevation, and with elevations
%            from -90:elev_res:90deg.
%   ang_res: [azi_resolution elev_resolution] in degrees
%   ls_dirs: speaker directions in [azi1 elev1; azi2 elev2;... ; aziL elevL]
%            convention, to be plotted on top of the data (optional)
%   h_ax:    axis handle to plot, useful for e.g. subplots (optional)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Archontis Politis, 15/11/2015
%   archontis.politis@aalto.fi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if exist('ls_dirs','var') && ~isempty(ls_dirs), LS_PLOT = 1;
else LS_PLOT = 0; end

if nargin<4
    LS_PLOT = 1;
    figure
    h_ax = gca;
end

% regular grid of directions
aziRes = ang_res(1);
polarRes = ang_res(2);
azi = -180:aziRes:180;
elev = -90:polarRes:90;


linesX = cat(3, repmat([-180; 180], [1 19]), [-90:10:90; -90:10:90], ones(2,19));
linesY = cat(3, [-180:10:180; -180:10:180], repmat([-90; 90], [1 37]), ones(2,37));

surf(azi, elev, F), view(2), colorbar, axis([-180 180 -90 90])
shading interp
line([linesX(:,:,1) linesY(:,:,1)], [linesX(:,:,2) linesY(:,:,2)], 1.5*max(F(:))*[linesX(:,:,3) linesY(:,:,3)], 'color',[100,100,100]/255)

if LS_PLOT
    hold on, plot3(ls_dirs(:,1), ls_dirs(:,2), 2*max(F(:))*ones(length(ls_dirs),1), 'ro', 'MarkerFaceColor','r')
end
