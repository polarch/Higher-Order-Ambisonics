function rE_mag = getTheoreticalEVmag(order)
%GETTHEORETICALEVMAG Theoretical energy vector magnitude of an ideal decoder.
% GETTHEORETICALEVMAG returns the theoretical energy vector magnitude of a
% continuous source distribution (ideal loudspeaker setup), when max-rE
% weighting is considered. The relation is the one found in
%
%   Zotter, F., Frank, M. (2012). All-Round Ambisonic Panning and Decoding. 
%   Journal of the Audio Engineering Society, 60(10), 807:820.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Archontis Politis, 15/11/2015
%   archontis.politis@aalto.fi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n = (0:order)';

% if ~MAX_RE, rE_mag = 2*sum(n+1)/sum(2*n+1); end

a_n = zeros(order+2,1);
for l=0:order+1
    temp = legendre(l, cos(137.9*(pi/180)/(order+1.51)));
    a_n(l+1) = temp(1);
end
a_n2 = a_n(n+1).*a_n(n+1);
a_n_a_n1 = a_n(n+1).*a_n(n+2);

rE_mag = 2*sum( (n+1).*a_n_a_n1 )/sum( (2*n+1).*a_n2 );
