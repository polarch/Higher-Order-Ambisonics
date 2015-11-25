function bfsig_rot = rotateBformat(bfsig, yaw, pitch, roll)
%ROTATEBFORMAT Rotate a B-format sound scene, by yaw-pitch-roll angles.
% ROTATEBFORMAT takes an input B-format signal set and according to a
% specified rotation of the coordinate system, returns the B-format signals
% on the rotated frame.
%
% The convention used is of a right-hand coordinate, active rotation, 
% yaw-pitch-roll convention, with yaw applied first, then pitch, then roll.
%
% The code requires the Spherical Harmonic transform library (for the 
% euler2rotationMatrix() functions) found in:
%
% <https://github.com/polarch/Spherical-Harmonic-Transform>
%
% Inputs:   
%   bfsig:  B-format signals
%   yaw:    yaw angle in degrees, around z (first rotation)
%   pitch:  pitch angle in degrees, around y' (second rotation)
%   roll:   roll angle in degrees, around x'' (third rotation)
% Outputs:
%   bfsig_rot:  rotated B-format signals
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Archontis Politis, 15/11/2015
%   archontis.politis@aalto.fi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get rotation matrix
Rzyx = euler2rotationMatrix(-yaw*pi/180, -pitch*pi/180, roll*pi/180, 'zyx');

% augment with zero order
Rbf = zeros(4,4);
Rbf(1) = 1;
Rbf(2:4,2:4) = Rzyx;

% apply to B-format signals
bfsig_rot = bfsig * Rbf.';

end
