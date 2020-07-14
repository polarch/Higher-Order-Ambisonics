function rir_filt = synthesizeHOAreverb(order, fs, t60, fc, FLATTEN)
%SYNTHESIZEHOAREVERB Synthesizes a quick and dirty ambisonic late reverb tail
% SYNTHESIZEHOAREVERB generates expoenntially decaying white noise
% sequences, with a user prescribed energy decay curve, controlled through
% the reverberation time parameters. If more than one reverberation times
% are defined, then bandpassed noise tails are generated, for more
% realistic frequency-dependent decay.
%
% Inputs:
% order:    HOA order
% fs:       sample rate
% t60:      reverberation times in different bands
% fc:       center frequencies of reverberation time bands (octave bands)
% FLATTEN:  {0,1} force the late RIR to have a flat-spectrum
%
% Output:
% rir_filt: the later ambisonic RIR
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Archontis Politis, 15/11/2015
%   archontis.politis@aalto.fi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if nargin<5, FLATTEN = 0; end
    
    % number of HOA channels
    nSH = (order+1)^2;
    % number of frequency bands
    nBands = length(t60);
    % decay constants
    alpha = 3*log(10)./t60;
    % length of RIR
    lFilt = ceil(max(t60)*fs);
    t = (0:1/fs:t60-1/fs)';
    % generate envelopes
    env = exp(-t*alpha);
    % generate RIRs
    rir = randn(lFilt, nSH, nBands);
    for k = 1:nBands
        rir(:, :, k) = rir(:,:,k).*(env(:,k)*ones(1,nSH));
    end
    % get filterbank IRs for each band
    filterOrder = 200;
    h_filt = filterbank(fc, filterOrder, fs);
    % filter rirs
    rir_filt = zeros(lFilt+ceil(filterOrder/2), nSH);
    for n = 1:nSH
        h_temp = [squeeze(rir(:,n,:)); zeros(ceil(filterOrder/2), nBands)];
        rir_filt(:, n) = sum(fftfilt(h_filt, h_temp), 2);
    end
    
    if FLATTEN, rir_filt = equalizeMinphase(rir_filt); end
end


function h_filt = filterbank(fc, filterOrder, fs)
% fc:   the center frequencies of the bands
% Nord: order of hte FIR filter

    if length(fc) == 1
        h_filt = 1;

    elseif length(fc) == 2
        h_filt = zeros(filterOrder+1, 2);

        % lowpass
        f_ll = 2*fc(1)/sqrt(2);
        w_ll = f_ll/(fs/2);
        h_filt(:, 1) = fir1(filterOrder, w_ll);
        % highpass
        f_hh = fc(2)/sqrt(2);
        w_hh = f_hh/(fs/2);
        h_filt(:, 2) = fir1(filterOrder, w_hh, 'high');

    else
        Nbands = length(fc);
        h_filt = zeros(filterOrder+1, Nbands);

        % lowpass
        f_ll = 2*fc(1)/sqrt(2);
        w_ll = f_ll/(fs/2);
        h_filt(:, 1) = fir1(filterOrder, w_ll);
        % highpass
        f_hh = fc(end)/sqrt(2);
        w_hh = f_hh/(fs/2);
        h_filt(:, end) = fir1(filterOrder, w_hh, 'high');
        % bandpass
        for k = 2:Nbands-1
            fl = fc(k)/sqrt(2);
            fh = 2*fc(k)/sqrt(2);
            wl = fl/(fs/2);
            wh = fh/(fs/2);
            w = [wl wh];
            h_filt(:, k) = fir1(filterOrder, w, 'bandpass');
        end
    end

end


function rir_filt_flat = equalizeMinphase(rir_filt)
%MAKEFLATVERB Makes the decaying noise spectrally flat
%   Detailed explanation goes here

Nrir = size(rir_filt,2);
for n=1:Nrir
    % equalise TDI by its minimum phase form to unity magnitude response
    tdi_f = fft(rir_filt(:,n));
    tdi_min_f = exp(conj(hilbert(log(abs(tdi_f)))));
    tdi_eq = real(ifft(tdi_f./tdi_min_f));
    rir_filt_flat(:,n) = tdi_eq;
end

end

