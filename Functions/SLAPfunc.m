function [scaling,metrics] = SLAPfunc(Sxx_lab,Sxx_fl,phi,fs,fb_inds,rms_inds,bf,Ts,p)
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
    % metrics = 3x1 vector with fb RMS, approx stress, and approx fatigue
    % ratios, respectively
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    % qrats = 1e20;

    % Calculate candidate scaling factors (Eq. 26)
    sc_factor1 = sig_rms_ratio^-2; 
    sc_factor2 = sig_peak_ratio^-2;
    sc_factor3 = fatigue_ratio^(-2/bf);
    sc_factor4 = max(qrats.^-2);

    % retain largest of the PSD scaling factors as the final scaling
    scaling = abs(max([sc_factor1 sc_factor2 sc_factor3 sc_factor4])); 

    % scaled values of metrics
    metrics = abs([sqrt(scaling)*sig_rms_ratio sqrt(scaling)*sig_peak_ratio scaling^(bf/2)*fatigue_ratio sqrt(scaling)*min(qrats)]); % get rid of residual (~0) imaginary parts, if any  
end