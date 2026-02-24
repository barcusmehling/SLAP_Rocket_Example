function [sigrms,sigpsd,sigloc] = GetStressFunc(Hs,Sff,df,rms_inds,sh_inds)
    % Hs = nstresspoints x 1 cell, where each element contains a 6 x
    % nforces matrix relating force PSD to stress tensor PSD

    % Sff = nforces x nforces x nflines. Forces applied in environment or
    % in the lab

    % df = frequency spacing

    % rms_inds = frequency lines to include in RMS calculation. If you use
    % the 0 Hz frequency line, RMS stress usually gets really big, so this
    % is mainly here to exclude that line.

    A1 = [1 -0.5 -0.5;-0.5 1 -0.5;-0.5 -0.5 1]; % Assemble VM Stress Quadratic Form Matrix
    A2 = 3*eye(3);
    A = blkdiag(A1,A2);

    rmss = zeros(length(Hs),1);

    H_ex = Hs{1};
    nf = size(H_ex,3); % frequency lines

    sig_psds = zeros(length(Hs),nf);
    
    for ii = 1:length(rmss)
        Hii = Hs{ii};
        
        for jj = 1:nf
            Ssigsig = squeeze(Hii(:,sh_inds,jj))*Sff(:,:,jj)*squeeze(Hii(:,sh_inds,jj))';
            sig_psds(ii,jj) = sig_psds(ii,jj) + abs(sum(sum(Ssigsig.*A)));
        end
        rmss(ii) = sqrt(df*sum(sig_psds(ii,rms_inds)));
    end
    [sigrms,ind1] = max(rmss);

    sigpsd = sig_psds(ind1,:);
    sigloc = ind1;

end