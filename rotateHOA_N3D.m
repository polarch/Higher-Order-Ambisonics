function hoasig_rot = rotateHOA_N3D(hoasig, yaw, pitch, roll)
%ROTATEHOA_N3D Summary of this function goes here
%   Detailed explanation goes here

% get order
Nchan = size(hoasig,2);
Nord = sqrt(Nchan)-1;

% get rotation matrix
Rzyx = euler2rotationMatrix(-yaw*pi/180, -pitch*pi/180, roll*pi/180, 'zyx');

% compute rotation matrix in the SH domain
Rshd = getSHrotMtx(Rzyx, Nord, 'real');

% apply to hoa signals
hoasig_rot = hoasig * Rshd.';

end
