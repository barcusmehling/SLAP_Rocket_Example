function fatigue_damage = GetFatigue(Ssigsig_fl,bf,fs,freq_inds,Tflight)
% Calculate fatigue damage ratio (lab / flight) using VM stress PSDs
%
% fatigue_damage = GetFatigue(Ssigsig_fl,bf,fs,freq_inds,Tflight)
%
% Ssigsig_fl = nf x 1 matrix containing flight stress PSD
% fs = nf x 1 frequency vector
% freq_inds = frequencies over which to compute the metric
%
% Tflight = duration of flight
% bf - see paper

ws = 2*pi*fs; % convert to rad/s from Hz

Ssigsig_fl = reshape(Ssigsig_fl,[length(Ssigsig_fl) 1]);
ws = ws(:); % make ws column vector (for elementwise mult)

scale_vec = ws.^(2/bf); % vector of weighted frequencies (higher freq = more load cycles = weighted heavier)

denom = sum(scale_vec(freq_inds).*Ssigsig_fl(freq_inds)); % fl term

fatigue_damage = (denom)^(bf/2)*Tflight; % scale by relative time to get damage ratio
end