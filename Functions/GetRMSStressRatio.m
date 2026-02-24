function rat = GetRMSStressRatio(Ssigsig_lab,Ssigsig_fl,rms_inds,df)
    % Calculate RMS stress ratio using stress PSDs

    numer = sqrt(sum(Ssigsig_lab(rms_inds))*df); % calculate RMS lab stress (numerator)
    denom = sqrt(sum(Ssigsig_fl(rms_inds))*df); % RMS fl stress (denominator)
    rat = numer/denom; % divide lab / flight to get ratio
end
