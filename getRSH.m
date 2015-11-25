function R_N = getRSH(N, dirs_deg)
%GETRSH Get vector of real orthonormal spherical harmonic values up to order N
%
% Inputs:
%   N:      maximum order of harmonics
%   dirs:   [azimuth_1 elevation_1; ...; azimuth_K elevation_K] angles
%           in degs for each evaluation point, where elevation is the
%           polar angle from the horizontal plane
%
% Outpus:
%   R_N:    [(N+1)^2 x K] matrix of SH values
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Archontis Politis, 10/10/2013
%   archontis.politis@aalto.fi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ndirs = size(dirs_deg, 1);
Nharm = (N+1)^2;

% convert to rads
dirs = dirs_deg*pi/180;
% initialize SH matrix
R_N = zeros(Nharm, Ndirs);
% zero order
R_N(1,:) = 1/sqrt(4*pi);
% higher orders
if N>0 
    idx_R = 1;
    for n=1:N
        
        m = (0:n)';
        % vector of unnormalised associated Legendre functions of current order
        Pnm = legendre(n, sin(dirs(:,2)'));
        % cancel the Condon-Shortley phase from the definition of
        % the Legendre functions to result in signless real SH
        uncondon = (-1).^[m(end:-1:2);m] * ones(1,Ndirs);
        Pnm = uncondon .* [Pnm(end:-1:2, :); Pnm];
        
        % normalisations
        norm_real = sqrt( (2*n+1)*factorial(n-m) ./ (4*pi*factorial(n+m)) );
        
        % convert to matrix, for direct matrix multiplication with the rest
        Nnm = norm_real * ones(1,Ndirs);
        Nnm = [Nnm(end:-1:2, :); Nnm];
        
        CosSin = zeros(2*n+1,Ndirs);
        % zero degree
        CosSin(n+1,:) = ones(1,size(dirs,1));
        % positive and negative degrees
        CosSin(m(2:end)+n+1,:) = sqrt(2)*cos(m(2:end)*dirs(:,1)');
        CosSin(-m(end:-1:2)+n+1,:) = sqrt(2)*sin(m(end:-1:2)*dirs(:,1)');
        Rnm = Nnm .* Pnm .* CosSin;
        
        R_N(idx_R + (1:2*n+1), :) = Rnm;
        idx_R = idx_R + 2*n+1;
    end
end

end
