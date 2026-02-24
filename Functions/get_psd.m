function psd = get_psd(Sxx) % retain diagonal PSD terms from CPSD matrix
    psd = zeros(size(Sxx,1),size(Sxx,3));
    for ii = 1:size(Sxx,1)
        for jj = 1:size(Sxx,3)
            psd(ii,jj) = real(Sxx(ii,ii,jj));
        end
    end
end