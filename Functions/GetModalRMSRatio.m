function rats = GetModalRMSRatio(Sqq_lab,Sqq_fl,rms_inds,df)
    % get RMS modal acceleration ratios

    nmodes = size(Sqq_lab,1); % num modes
    rats = zeros(nmodes,1); % initialize ratios vec
    for ii = 1:nmodes
        numer = sqrt(sum(squeeze(Sqq_lab(ii,ii,rms_inds)))*df); % lab term = numerator
        denom = sqrt(sum(squeeze(Sqq_fl(ii,ii,rms_inds)))*df); % fl term = denom
        rats(ii) = numer/denom; % get ratio (lab/fl)
    end
end