function stereosig = decodeBformat2stereo(BFsig, dirCoeff, stereoAngle)
%DECODEBFORMAT2STEREO Returns stereo speaker signals from B-format signals.
% DECODEBFORMAT2STEREO decodes B-format signals to a stereo loudspeaker layout. 
%
% Inputs:
%   BFsig:      B-format signals
%   dirCoeff:   directivity coefficient of virtual microphone stereo pair
%               {1: omni, 0.5: cardioid, 0.366: supercardioid, 0.25: hypercardioid, 0:dipole}
%   stereoAngle:    stereo angle (angle span of virtual microphone pair)
%
% Outputs:
%   LSsig:  decoded loudspeaker signals
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Archontis Politis, 15/11/2015
%   archontis.politis@aalto.fi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

foasig_N3D = BFsig(:,1:4);
foasig_SN3D = convert_N3D_SN3D(foasig_N3D, 'n2sn');

vmic_dirs = [stereoAngle/2 0; -stereoAngle/2 0]*pi/180;
stereoDecMtx = VmicB(vmic_dirs, dirCoeff);
stereosig = foasig_SN3D * stereoDecMtx.';

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function vmicGains = VmicB(vmic_directions, a)
% VMICB Computes the gains on the B-format signals for virtual microphones

    % get the unit vectors of each vmic direction
    Nvmic = size(vmic_directions, 1);
    u_vmic = zeros(Nvmic, 3);
    [u_vmic(:,1), u_vmic(:,2), u_vmic(:,3)] = sph2cart(vmic_directions(:,1), ...
        vmic_directions(:,2), ones(Nvmic, 1));
    
    vmicGains = [a*ones(Nvmic, 1) (1-a)*u_vmic]';
    
end
