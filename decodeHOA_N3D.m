function LSsig = decodeHOA_N3D(hoasig, M_dec, cutoffs, fs)
%DECODEHOA_N3D Returns speaker signals from HOA signals and a decoding matrix.
% DECODEHOA_N3D decodes HOA signals to a specific loudspeaker layout. The 
% decoding matrix (or matrices), should be computed first (by ambiDecoder()
% function). Frequency-dependent decoding can be achieved by specifying the
% cutoff frequencies that switch between different decoding matrices, and
% by passing as many matrices as ranges.
%
% Inputs:   
%   hoasig: HOA signals (N3D normalization, ACN channel order)
%   M_dec:  [K x (order+1)^2 x (Ncutoff+1)] decoding matrix, where K is the
%           number of speakers, order is the HOA signals order, and Ncutoff
%           is the number of cutoff frequencies (number of decoding ranges
%           minus one). The third dimension contains the different decoder
%           matrices for the higher frequency ranges, see
%           TEST_AMBI_SCRIPT.m for an example.
%   cutoffs: vector of cutoff frequencies, for frequency-dependent decoding
%   fs:     sample rate
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
