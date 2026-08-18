function [modal_RMS,varargout] = GetModalRMS(Sqq_fl,rms_inds,df)
% get RMS modal acceleration
%
% [modal_RMS] = GetModalRMS(Sqq_fl,rms_inds,df)
%
% Sqq_fl = [Nmodes x Nmodes x Nf] spectral density matrix of modal accelerations.
% rms_inds = frequency indices over which to compute the RMS
% df = frequency spacing (Hz)
%
% Output:
% modal_RMS = (scalar) RMS value for mode with the largest response.
%
nmodes = size(Sqq_fl,1); % num modes
modal_RMS = zeros(nmodes,1); % initialize ratios vec
for ii = 1:nmodes
    modal_RMS(ii) = sqrt(sum(squeeze(Sqq_fl(ii,ii,rms_inds)))*df); % fl term = denom
end

[modal_RMS,m_index] = max(abs(modal_RMS)); % imaginary parts are small artifacts of computation

if nargout>2
    varargout{1}=m_index;
end

end