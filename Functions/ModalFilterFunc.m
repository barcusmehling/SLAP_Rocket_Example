function Sqq = ModalFilterFunc(Sxx,phi,fs,fb_inds,dispflag)
    % Sxx = flight or lab environment PSD (g^2/Hz)

    % phi = matrix containing DOF in Sxx, TS modes, and fixed-base modes

    % fs = frequency vector

    % fb_inds = indices in phi corresponding to fixed-base modes

    % dispflag == 1 --> convert acc to disp PSD
    % dispflag == 0 --> keep as acc PSD

    df = fs(2)-fs(1); % frequency resolution
    nf = length(fs); % num frequency lines
    nfb = length(fb_inds); % num fixed base modes

    ws = 2*pi*fs; % angular frequency

    nmodes = size(phi,2);
    Sqq = zeros(nmodes,nmodes,nf); % initialize modal disp (or acceleration) PSD matrix
    
    if dispflag == 1
        for ii = 1:nf % get disp PSD at each fline and modal filter
            Sxx_d = Sxx(:,:,ii) / (ws(ii)^4) * 9.81^2; % convert g^2/Hz to m^2/Hz
            Sqq(:,:,ii) = pinv(phi)*Sxx_d*pinv(phi)';
        end
        Sqq = Sqq(fb_inds,fb_inds,:); % retain fixed-base modes only
    else
        for ii = 1:nf
            Sqq(:,:,ii) = pinv(phi)*Sxx(:,:,ii)*pinv(phi)';
        end
    end

    
end