function [Ssigsig_lab,Ssigsig_fl] = GetStressPSDs(Spp_lab,Spp_fl)
    % Get approximate stress PSDs to use in SLAP

    nf = size(Spp_lab,3); % num freq lines

    Ssigsig_lab = zeros(nf,1); % lab and flight (approx.) VM stress PSDs
    Ssigsig_fl = zeros(nf,1);

    for ii = 1:nf % estimate stress at each freq line (Eqs. 20 and 21, "method for conservative MIMO vibration testing")
        diag_lab = real(sum(diag(Spp_lab(:,:,ii)))); % sum FB PSDs for lab

        Ssigsig_lab(ii) = Ssigsig_lab(ii) + diag_lab; % add FB PSD terms to lab stress PSD

        Spp_off_lab = Spp_lab(:,:,ii)-diag(diag(Spp_lab(:,:,ii))); % retain FB CPSD terms

        off_lab = .5*sum(sum(abs(real(Spp_off_lab)))); % create off diagonal term for stress PSD 

        Ssigsig_lab(ii) = Ssigsig_lab(ii) - off_lab; % subtract off diag terms from lab VM stress PSD

        diag_fl = real(sum(diag(Spp_fl(:,:,ii)))); % repeat previous steps for flight environment

        Ssigsig_fl(ii) = Ssigsig_fl(ii) + diag_fl;

        Spp_off_fl = Spp_fl(:,:,ii)-diag(diag(Spp_fl(:,:,ii)));

        off_fl = .5*sum(sum(abs(real(Spp_off_fl))));

        Ssigsig_fl(ii) = Ssigsig_fl(ii) + off_fl; % add off diag terms to flight VM stress PSD 
    end

    inds = find(Ssigsig_lab < 0); % make all nonnegative terms zero

    Ssigsig_lab(inds) = 0;

end