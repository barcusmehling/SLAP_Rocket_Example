function rd_est = GetFatigueRatio(Ssigsig_lab,Ssigsig_fl,bf,fs,freq_inds,Ts)
    % Calculate fatigue damage ratio (lab / flight) using VM stress PSDs

    Tlab = Ts(1); % duration of lab test
    Tflight = Ts(2); % duration of flight
    ws = 2*pi*fs; % convert to rad/s from Hz
    
    Ssigsig_lab = reshape(Ssigsig_lab,[length(Ssigsig_lab) 1]); % make VM stress PSDs column vectors (if not already)
    Ssigsig_fl = reshape(Ssigsig_fl,[length(Ssigsig_fl) 1]);
    ws = reshape(ws,[length(ws) 1]); % make ws column vector (for elementwise mult)

    scale_vec = ws.^(2/bf); % vector of weighted frequencies (higher freq = more load cycles = weighted heavier)

    numer = sum(scale_vec(freq_inds).*Ssigsig_lab(freq_inds)); % numerator / lab term weighted avg.

    denom = sum(scale_vec(freq_inds).*Ssigsig_fl(freq_inds)); % fl term

    rd_est = (numer/denom)^(bf/2)*Tlab/Tflight; % scale by relative time to get damage ratio
end