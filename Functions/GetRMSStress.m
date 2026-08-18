function RMSStr = GetRMSStress(Ssigsig_fl,rms_inds,df)
% Calculate RMS stress from a matrix of stress PSDs
%
% RMSStr = GetRMSStress(Ssigsig_fl,rms_inds,df)
%
% where
%   Ssigsig_fl = matrix of stress PSDs
%   rms_inds = frequency indices over which to compute the RMS
%   df = frequency spacing in Hz
% 

RMSStr = sqrt(sum(Ssigsig_fl(rms_inds))*df); % RMS fl stress

end
