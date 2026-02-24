% Modal Filter environment and expand out to reference DOF to check how
% accurately modal responses are estimated

clc;close all;clear all;
addpath ..\Functions; % for get_psd.m
load ..\ModeShapes\BARCAccelModes; % from GetModesAtAccels
load ..\Environment\RocketEnv;

Sxx_BARC = Sxx(28:end,28:end,:); % DUT channels only (first 27 baseplate)
filt_inds = 1:90;
ref_inds = 40:42;
filt_inds(ref_inds) = []; % modal filter using 87 DUT channels (exclude 3)

nmodes = size(phi,2);
nf = size(Sxx,3);
Sqq = zeros(nmodes,nmodes,nf); % initialize some variables
Sxx_est = zeros(size(Sxx_BARC));

for ii = 1:nf % estimate modal responses and expand back to physical coordinates
    Sqq(:,:,ii) = pinv(phi(filt_inds,:))*Sxx_BARC(filt_inds,filt_inds,ii)*pinv(phi(filt_inds,:))';
    Sxx_est(:,:,ii) = phi*Sqq(:,:,ii)*phi';
end 

refpsd_fl = get_psd(Sxx_BARC(ref_inds,ref_inds,:)); % retain diagonal terms of original and filter estimate for ref accel
refpsd_filt = get_psd(Sxx_est(ref_inds,ref_inds,:));

figure; % plot 
tlt = tiledlayout(3,1);
titles = {'Ref X','Ref Y','Ref Z'};

for ii = 1:3
    nexttile;
    semilogy(fs,abs(refpsd_fl(ii,:)),'k',fs,abs(refpsd_filt(ii,:)),'r--','Linewidth',2)
    grid on;
    xlim([10 2000])
    title(titles{ii})
    if ii == 1
        legend('Flight','Filtered')
    end
end

ylabel(tlt,'Acceleration PSD (g^2/Hz)','interpreter','tex')
xlabel(tlt,'Frequency (Hz)','interpreter','tex')
tlt.Padding = 'compact';
tlt.TileSpacing = 'tight';

%% plot modal responses if desired
% qpsds = get_psd(Sqq);
% qpsds_fb = qpsds(7:end,:);
% 
% figure;
% semilogy(fs,abs(qpsds_fb)./(ws.^4),'Linewidth',2)












