function rats = GetModalRMSRatio(Sqq_lab,Sqq_fl,rms_inds,df)
% get RMS modal acceleration ratios
%
% [ratios] = GetModalRMSRatio(Sqq_lab,Sqq_fl,rms_inds,df)
%
% Sqq_lab = [Nmodes x Nmodes x Nf] spectral density matrix of modal accelerations (lab).
% Sqq_fl = [Nmodes x Nmodes x Nf] spectral density matrix of modal accelerations (flight).
% rms_inds = frequency indices over which to compute the RMS
% df = frequency spacing (Hz)
%
nmodes = size(Sqq_lab,1); % num modes
rats = zeros(nmodes,1); % initialize ratios vec
for ii = 1:nmodes
    numer = sqrt(sum(squeeze(Sqq_lab(ii,ii,rms_inds)))*df); % lab term = numerator
    denom = sqrt(sum(squeeze(Sqq_fl(ii,ii,rms_inds)))*df); % fl term = denom
    rats(ii) = numer/denom; % get ratio (lab/fl)
end
end