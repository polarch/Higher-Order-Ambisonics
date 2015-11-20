function rE_mag = getTheoreticalEVmag(order)
%GETTHEORETICALEVMAG Returns the theoretical energy vector magnitude of an ideal decoder.
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
