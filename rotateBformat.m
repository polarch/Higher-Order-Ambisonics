function bfsig_rot = rotateBformat(bfsig, yaw, pitch, roll)
%ROTATEBFORMAT Summary of this function goes here
%   Detailed explanation goes here

% get rotation matrix
Rzyx = euler2rotationMatrix(-yaw*pi/180, -pitch*pi/180, roll*pi/180, 'zyx');

% augment with zero order
Rbf = zeros(4,4);
Rbf(1) = 1;
Rbf(2:4,2:4) = Rzyx;

% apply to B-format signals
bfsig_rot = bfsig * Rbf.';

end
