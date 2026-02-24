% Simulate a MIMO test controlling to flight environment
clc;close all;clear all;
addpath ..\Functions; 

load ..\Environment\RocketEnv; % flight environment
load ..\FRFs\LabFRF; % control FRF
load ..\FRFs\LabStressFRFs.mat; % lab stress FRFs to calculate stress in lab

ctrl_accs = 1:27; % control channels (baseplate accels, like a 6DOF test)
ref_accs = 67:69; % reference channels (a random triax on the DUT)

Sxx_est = zeros(size(Sxx)); % lab environment

nsh = size(H,2); % number of shakers
Sff = zeros(nsh,nsh,length(fs)); % lab force PSD matrix

for ii = 1:length(fs) % calculate shaker forces and lab env at each fline
    cnthresh = 0.01*max(svd(H(ctrl_accs,:,ii))); % condition number threshold = 0.01 * largest SV of FRF mat
    % following line calculates forces while implementing CN threshold. See
    % pinv.m documentation for more detail.
    Sff(:,:,ii) = pinv(H(ctrl_accs,:,ii),cnthresh)*Sxx(ctrl_accs,ctrl_accs,ii)*pinv(H(ctrl_accs,:,ii),cnthresh)';
    Sxx_est(:,:,ii) = H(:,:,ii)*Sff(:,:,ii)*H(:,:,ii)'; % calculate lab response at all DOF
end

envpsd = get_psd(Sxx(ref_accs,ref_accs,:)); % flight environment PSDs at DUT ref accels
labpsd = get_psd(Sxx_est(ref_accs,ref_accs,:)); % '' lab env

% Plot reference responses in flight vs. lab
figure('Units','normalized','Position',[0.1 0.1 0.8 0.4]);
tlt = tiledlayout(1,3);
titles = {'Ref X','Ref Y','Ref Z'};

for ii = 1:3
    nexttile;
    semilogy(fs,abs(envpsd(ii,:)),'k',fs,abs(labpsd(ii,:)),'b','Linewidth',2)
    xlabel('Frequency (Hz)','interpreter','tex')
    grid on;
    xlim([10 2000])
    if ii == 1
        legend('Flight','Test')
    end
    title(titles{ii})
end

ylabel(tlt,'Acceleration PSD (g^2/Hz)','interpreter','tex')
tlt.TileSpacing = 'tight';
tlt.Padding = 'tight';

% Stress analysis
ind1 = find(fs >= 95,1);
ind2 = find(fs >= 2000,1);
rms_inds = ind1:ind2; % calculate RMS stress only including freqs above 100 Hz
df = 5; % 5 Hz frequency spacing

[sigrms_lab,sigpsd_lab,sigloc_lab] = GetStressFunc(Hs,Sff,df,rms_inds,1:7); % calculate max lab VM stress PSD

load ..\Environment\FlightStressPSD; % max flight VM stress PSD

% compare flight and lab stress PSDs
figure;
semilogy(fs,abs(sigpsd_fl),'k',fs,abs(sigpsd_lab),'b','Linewidth',2)
grid on;
xlabel('Frequency (Hz)')
ylabel('Stress PSD (Pa^2/Hz)')
legend('Flight','Lab Test')
xlim([10 2000])







