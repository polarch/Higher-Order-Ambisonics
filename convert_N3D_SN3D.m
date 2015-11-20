function outsig = convert_N3D_SN3D(insig, type)
%CONVERT_N3D_TO_SN3D Summary of this function goes here
%   Detailed explanation goes here

order = sqrt(size(insig,2))-1;

switch lower(type)
    case 'sn2n'
        SN3Dsig = insig;
        N3Dsig = zeros(size(SN3Dsig));
        % transformation between N3D and SN3D is a simple scaling of same
        % order components
        for no=0:order
            idx_order = no^2 + (1:2*no+1);
            N3Dsig(:, idx_order) = sqrt(2*no+1) * SN3Dsig(:, idx_order);
        end
        outsig = N3Dsig;
        
    case 'n2sn'
        N3Dsig = insig;
        SN3Dsig = zeros(size(N3Dsig));
        % transformation between N3D and SN3D is a simple scaling of same
        % order components
        for no=0:order
            idx_order = no^2 + (1:2*no+1);
            SN3Dsig(:, idx_order) = N3Dsig(:, idx_order)/sqrt(2*no+1);
        end        

        outsig = SN3Dsig;
end

end

