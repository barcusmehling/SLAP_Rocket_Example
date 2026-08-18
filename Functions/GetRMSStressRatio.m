function rat = GetRMSStressRatio(Ssigsig_lab,Ssigsig_fl,rms_inds,df)
% Calculate RMS stress ratio using stress PSDs
%
% rat = GetRMSStressRatio(Ssigsig_lab,Ssigsig_fl,rms_inds,df)
%
% where
%   Ssigsig_fl = matrix of stress PSDs
%   rms_inds = frequency indices over which to compute the RMS
%   df = frequency spacing in Hz
%
numer = sqrt(sum(Ssigsig_lab(rms_inds))*df); % calculate RMS lab stress (numerator)
denom = sqrt(sum(Ssigsig_fl(rms_inds))*df); % RMS fl stress (denominator)
rat = numer/denom; % divide lab / flight to get ratio
end
