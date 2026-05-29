function [scaling,metrics] = SLAPfunc(Sxx_lab,Sxx_fl,phi,fs,fb_inds,rms_inds,bf,Ts,p,varargin)
% Apply Scaled Lab PSD method (SLAP) to create a specification for a
% vibration qualification test
% Marcus Behling | 10/15/2025
%
% [scaling,metrics] = SLAPfunc(Sxx_lab,Sxx_fl,phi,fs,fb_inds,rms_inds,bf,Ts,p);
%
%%%%% Inputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sxx_lab = (No x No) lab PSDs to be modal filtered
% Sxx_fl = (No x No) flight PSDs to be modal filtered
% phi = (No x N) mode shapes used in modal filtering - Should contain the
%    rigid body modes followed by the fixed-interface modes.
% fs = frequency vector for PSDs
% fb_inds = which columns of phi correspond to fixed-interface modes
% rms_inds = which freq. lines to include when calculating metrics
% bf = fatigue exponent (material property)
% Ts = durations of [lab flight]
% p = confidence level (between 0 and 1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%% Outputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% scaling = final scale factor applied to lab PSD to cause the metrics to
%    all be conservative (RMS is scaled by the square root of this value).
% metrics = 4x1 vector with [RMS Stress; Peak Stress; Fatigue, RMS FB modal resp.]
%    ratios between flight and test, respectively
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% To exclude a metric, use:
% [scaling,metrics] = SLAPfunc(Sxx_lab,Sxx_fl,phi,fs,fb_inds,rms_inds,bf,Ts,p,use_ind);
%   i.e. use_ind=[1,2,3] to use the first three metrics only.
%

if nargin>9
    use_ind=varargin{1};
else
    use_ind=[1:4]; % Use all metrics by default.
end

df = fs(2)-fs(1); % frequency spacing for calculating RMS
ws = 2*pi*fs; % angular frequency vector

% get fixed-base modal displacement RMS and PSDs
dflag = 1; % displacement flag = 1 -> calc modal displacement
Spp_lab = ModalFilterFunc(Sxx_lab,phi,fs,fb_inds,dflag); 
Spp_fl = ModalFilterFunc(Sxx_fl,phi,fs,fb_inds,dflag);

% Calculate stress-based damage metrics
[Ssigsig_lab, Ssigsig_fl] = GetStressPSDs(Spp_lab,Spp_fl); % approx VM stress PSDs
sig_rms_ratio = GetRMSStressRatio(Ssigsig_lab,Ssigsig_fl,rms_inds,df); % RMS stress ratio
sig_peak_ratio = GetPeakStressRatio(Ssigsig_lab,Ssigsig_fl,fs,rms_inds,Ts,p); % peak stress ratio
fatigue_ratio = GetFatigueRatio(Ssigsig_lab,Ssigsig_fl,bf,fs,rms_inds,Ts); % fatigue ratio

% Calculate acceleration-based damage metric
dflag = 0; % calculate acceleration
Sqq_lab = ModalFilterFunc(Sxx_lab,phi,fs,fb_inds,dflag);
Sqq_fl = ModalFilterFunc(Sxx_fl,phi,fs,fb_inds,dflag);
qrats = GetModalRMSRatio(Sqq_lab,Sqq_fl,rms_inds,df);

% Calculate candidate scaling factors (Eq. 26)
sc_factors(1) = sig_rms_ratio^-2; 
sc_factors(2) = sig_peak_ratio^-2;
sc_factors(3) = fatigue_ratio^(-2/bf);
sc_factors(4) = max(qrats.^-2);

% retain largest of the PSD scaling factors as the final scaling
scaling = abs(max(sc_factors(use_ind))); 

% scaled values of metrics
metrics = abs([sqrt(scaling)*sig_rms_ratio sqrt(scaling)*sig_peak_ratio scaling^(bf/2)*fatigue_ratio sqrt(scaling)*min(qrats)]); % get rid of residual (~0) imaginary parts, if any  

end