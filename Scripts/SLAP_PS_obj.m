%% Objective function for SLAP-PS

function [f_obj,f_metrics,f_inp,metrics_ratio] = SLAP_PS_obj(F_in,metrics_fl,Winp,H_lab,fs,sh_inds,Sxx,phi_filt,filt_inds,fb_inds,rms_inds,bf,Tlab,p)
nsh = size(H_lab,2);
nf = length(fs); 

% create buzz test environment
Sff_PS = diag(F_in); % Force PSD matrix (assuming uncorrelated inputs)
Sxx_PS = zeros(size(Sxx(filt_inds,filt_inds,:)));
for ii = 1:nf
    Sxx_PS(:,:,ii) = H_lab(filt_inds,sh_inds,ii)*Sff_PS*H_lab(filt_inds,sh_inds,ii)';
end

% Compute metrics for this force
[metrics_lab] = GetMetrics(Sxx_PS,phi_filt,fs,fb_inds,rms_inds,bf,Tlab,p);

% Compute the SLAP damage metrics and the scaling required to make the
% test conservative.
metrics_ratio = metrics_lab./metrics_fl;
f_metrics = sum((log10(metrics_lab)-log10(metrics_fl)).^2);
f_inp = sum(F_in.^2);

f_obj =  f_metrics + Winp*f_inp;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end