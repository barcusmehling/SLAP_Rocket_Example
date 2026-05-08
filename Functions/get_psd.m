function psd = get_psd(Sxx)
% Extract diagonal PSD terms from a CPSD matrix
% Sxx is a (No x No x Nf) autospectral density matrix
% psd is an (No x Nf) matrix containing only the diagonal PSDs Sxx(i,i,:)
%
    psd = zeros(size(Sxx,1),size(Sxx,3));
    for ii = 1:size(Sxx,1)
        for jj = 1:size(Sxx,3)
            psd(ii,jj) = real(Sxx(ii,ii,jj));
        end
    end
end