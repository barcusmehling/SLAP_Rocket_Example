function rat = GetPeakStressRatio(sig_lab,sig_fl,fs,rms_inds,Ts,p)
    % Calculate peak stress ratio

    %%%%%%%% Variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % sig_fl = flight von mises stress PSD
    % sig_lab = lab von mises stress PSD
    % fs = frequency vector
    % rms_inds = indices over which to calculate statistical moments m0 and
    % m2
    % Ts = [Tfl Tlab] durations of flight operation and lab test
    % p = confidence level (between 0 and 1)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    Tl = Ts(1); % lab test duration
    Tf = Ts(2); % flight duration

    df = fs(2)-fs(1); % freq spacing
    dw = 2*pi*df; % freq spacing in rad/s
    ws = 2*pi*fs; % freq vec rad/s

    sig_fl = reshape(sig_fl,[length(sig_fl) 1]); % make all column vectors for elementwise mult
    sig_lab = reshape(sig_lab,[length(sig_lab) 1]);
    ws = reshape(ws,[length(ws) 1]);

    % Calculate peak stress in lab and flight, divide to get ratio
    % See: Kabe Struct Dynamics Vol 2, p. 590
    lam = @(b,sig_psd) sum(1/2/pi*ws(rms_inds).^b.*sig_psd(rms_inds)*dw); % function handle that calculates bth order spectral moments 
    
    m0f = lam(0,sig_fl); % zeroth order spectral moment - flight
    m2f = lam(2,sig_fl); % second order spectral moment - flight
    v0f = 1/2/pi*sqrt(m2f/m0f); % zero crossing rate flight - p. 598 Kabe

    m0l = lam(0,sig_lab); % same but for lab env
    m2l = lam(2,sig_lab);
    v0l = 1/2/pi*sqrt(m2l/m0l);

    % flight and lab peak stresses - p. 598 Kabe
    peakf = sqrt(2*log(v0f*Tf/-log(p)))*sqrt(m0f); 
    peakl = sqrt(2*log(v0l*Tl/-log(p)))*sqrt(m0l);

    rat = peakl / peakf; % peak stress ratio
end