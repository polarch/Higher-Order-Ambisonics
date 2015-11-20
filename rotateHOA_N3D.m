function hoasig_rot = rotateHOA_N3D(hoasig, yaw, pitch, roll)
%ROTATEHOA_N3D Rotate a HOA encoded sound scene, by yaw-pitch-roll angles.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Archontis Politis, 15/11/2015
%   archontis.politis@aalto.fi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
