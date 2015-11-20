function bfsig_rot = rotateBformat(bfsig, yaw, pitch, roll)
%ROTATEBFORMAT Rotate a B-format sound scene, by yaw-pitch-roll angles.
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
