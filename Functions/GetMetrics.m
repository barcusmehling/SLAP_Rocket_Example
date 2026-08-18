function [metrics] = GetMetrics(Sxx_fl,phi,fs,fb_inds,rms_inds,bf,Tflight,p)
    % Apply Scaled Lab PSD method (SLAP) to create a specification for a
    % vibration qualification test
    % Marcus Behling | 10/15/2025
    
    %%%%% Inputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Sxx_lab = lab PSDs to be modal filtered
    % Sxx_fl = flight PSDs to be modal filtered
    % phi = mode shapes used in modal filtering
    % fs = frequency vector for PSDs
    % fb_inds = which columns of phi correspond to fixed-base modes
    % rms_inds = which freq. lines to include when calculating metrics
    % bf = fatigue exponent (material property)
    % Ts = durations of [lab flight]
    % p = confidence level (between 0 and 1)
    % varargin{1} = small component flag - use RB and FB modes as the
    % metric of interest
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%% Outputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % scaling = final scale factor applied to lab PSD (RMS is scaled by the
    % square root of this value)
    % metrics = 4x1 vector with [RMS Stress, Peak Stress, Fatigue, RMS FB Accel]
    % ratios, respectively
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    df = fs(2)-fs(1); % frequency spacing for calculating RMS
    ws = 2*pi*fs; % angular frequency vector
    
    % get fixed-base modal displacement RMS and PSDs
    dflag = 1; % displacement flag = 1 -> calc modal displacement
    Spp_fl = ModalFilterFunc(Sxx_fl,phi,fs,fb_inds,dflag);

    % Calculate stress-based damage metrics
    [Ssigsig_fl] = GetStressPSDs(Spp_fl); % approx VM stress PSDs
    sig_rms = GetRMSStress(Ssigsig_fl,rms_inds,df); % RMS stress
    sig_peak = GetPeakStress(Ssigsig_fl,fs,rms_inds,Tflight,p); % peak stress
    fatigue = GetFatigue(Ssigsig_fl,bf,fs,rms_inds,Tflight); % fatigue for flight
    % Calculate acceleration-based damage metric
    dflag = 0; % calculate acceleration
    Sqq_fl = ModalFilterFunc(Sxx_fl,phi,fs,fb_inds,dflag);
    qRMS = GetModalRMS(Sqq_fl,rms_inds,df);
    % qrats = 1e20;

    metrics = [sig_rms, sig_peak, fatigue, qRMS];

end