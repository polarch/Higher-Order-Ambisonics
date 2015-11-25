function outsig = convert_ACN_SID(insig, type)
%CONVERT_ACN_SID Converts from the ACN to SID convention for HOA channel indexing
% CONVERT_ACN_SID Converts between HOA signals in the "canonical" indexing
% of SH components, q=n^2+N+m+1, where q is the channel number, n is the
% order and m is the degree, and the SID channel indexing, as defined by
% Jerome Daniel in
%
%   Daniel, J. (2003). Spatial Sound Encoding Including Near Field Effect : 
%   Introducing Distance Coding Filters and a Viable , New Ambisonic Format. 
%   In 23rd International Conference of AES. Copenhagen, Denmark.
%
% The ACN indexing follows the order (0,0), (1,-1), (1,0), (1,1), ...
% The SID indexing follows the order (0,0), (1,1), (1,-1), (1,0), ...
%
% Inputs:
%   insig:  The ACN or SID HOA signals
%   type:   'acn2sid' to convert from ACN to SID, or 'sid2acn' to do the
%           opposite
%
% Outputs:
%   outsig:  the converted signals
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Archontis Politis, 15/11/2015
%   archontis.politis@aalto.fi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

order = sqrt(size(insig,2))-1;
idx_ACN = zeros((order+1)^2, 2);
idx_SID = zeros((order+1)^2, 2);
idx_ACN(1,:) = [0 0];
idx_SID(1,:) = [0 0];
idx_n = 1;
for n = 1:order
    idx_ACN(idx_n+(1:2*n+1),1) = n;
    idx_ACN(idx_n+(1:2*n+1),2) = -n:n;
    idx_SID(idx_n+(1:2*n+1),1) = n;
    for m = 1:n
        idx_SID(idx_n+2*m-1,2) = n-m+1;
        idx_SID(idx_n+2*m,2) = m-n-1;
    end
    idx_SID(idx_n+2*n+1,2) = 0;
    
    idx_n = idx_n + 2*n+1;
end

switch lower(type)
    case 'acn2sid'
        ACNsig = insig;
        
        exchange_idx = idx_SID(:,1).^2 + idx_SID(:,1) + idx_SID(:,2) + 1;
        SIDsig = ACNsig(:, exchange_idx);            

        outsig = SIDsig;
        
    case 'sid2acn'
        SIDsig = insig;
        
        [~, exchange_idx] = sort(idx_SID(:,1).^2 + idx_SID(:,1) + idx_SID(:,2) + 1);
        ACNsig = SIDsig(:, exchange_idx);

        outsig = ACNsig;
end


end

