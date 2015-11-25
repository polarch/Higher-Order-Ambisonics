function LSsig = decodeBformat(BFsig, M_dec, cutoffs, fs)
%DECODEBFORMAT Returns speaker signals from B-format signals and a decoding matrix.
% DECODEBFORMAT decodes B-format signals to a specific loudspeaker layout. 
% The decoding matrix (or matrices), should be computed first (by the 
% ambiDecoder() function). Frequency-dependent decoding can be achieved by 
% specifying the cutoff frequencies that switch between different decoding 
% matrices, and by passing as many matrices as ranges.
%
% Inputs:
%   BFsig:  B-format signals
%   M_dec:  [K x 4 x (Ncutoff+1)] decoding matrix, where K is the number of 
%           speakers, order is the HOA signals order, and Ncutoff is the 
%           number of cutoff frequencies (number of decoding ranges minus 
%           one). The third dimension contains the different decoder matrices 
%           for the higher frequency ranges, see TEST_AMBI_SCRIPT.m for an 
%           example.
%   cutoffs:    vector of cutoff frequencies, for frequency-dependent decoding
%   fs:         sample rate
%
% Outputs:
%   LSsig:  decoded loudspeaker signals
%
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


N3Dsig = convert_N3D_Bformat(BFsig, 'b2n');
if nBands == 1
    LSsig = N3Dsig * M_dec(:,1:4).';
else
    N3Dsig_band = applyFilterbank(N3Dsig, cutoffs, fs);
    
    nLS = size(M_dec,1);
    LSsig_band = zeros(size(N3Dsig,1), nLS, nBands);
    for nb=1:nBands
        LSsig_band(:,:,nb) = N3Dsig_band(:,:,nb) * M_dec(:,1:4,nb).';
    end
    LSsig = sum(LSsig_band,3);    
end

end


function sig_band = applyFilterbank(sig, cutoffs, fs)

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
