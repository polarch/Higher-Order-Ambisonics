function LSsig = decodeHOA_N3D(hoasig, M_dec, cutoffs, fs)
%DECODEHOA_N3D Returns speaker signals from HOA signals and a decoding matrix.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Archontis Politis, 15/11/2015
%   archontis.politis@aalto.fi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<3, nBands = 1;
else nBands = length(cutoffs)+1; end

if size(M_dec,3) ~= nBands
    error('One decoding matrix should be provided for each frequency band')
end

if nBands == 1
    LSsig = hoasig * M_dec.';
else
    hoasig_band = applyFilterbank(hoasig, cutoffs, fs);
    
    nLS = size(M_dec,1);
    LSsig_band = zeros(size(hoasig,1), nLS, nBands);
    for nb=1:nBands
        LSsig_band(:,:,nb) = hoasig_band(:,:,nb) * M_dec(:,:,nb).';
    end
    LSsig = sum(LSsig_band,3);    
end

end


function sig_band = applyFilterbank(sig, cutoffs, fs)

lSig = size(sig,1);
nCh = size(sig,2);
nBand = length(cutoffs)+1;
% order of fir filters
lFilt = 100;

% zero pad
sig_pad = [sig; zeros(lFilt/2,nCh)];
sig_band = zeros([size(sig_pad) nBand]);


% create first and last lowpass and highpass in the filterbank
filters = zeros(lFilt+1, nBand);
filters(:,1) = fir1(lFilt, cutoffs(1)/(fs/2), 'low');
filters(:,nBand) = fir1(lFilt, cutoffs(nBand-1)/(fs/2), 'high');

if nBand > 2
    for i = 2:(nBand-1)
        filters(:,i) = fir1(lFilt, [cutoffs(i-1) cutoffs(i)]/(fs/2), 'bandpass');
    end
end

for i = 1:nBand
    sig_band(:,:,i) = fftfilt(filters(:,i), sig_pad);
end

sig_band = sig_band(lFilt/2+1:end,:,:);

end
