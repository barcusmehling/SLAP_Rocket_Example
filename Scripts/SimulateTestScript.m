% Simulate a MIMO test controlling to flight environment
clc;close all;clear all;
addpath ..\Functions; 

load ..\LargeFiles\RocketEnv; % flight environment
    % Sxx - No x No x Nf matrix of spectral densities in flight.
    % fs - 1 x Nf frequency vector in Hz  
load ..\FRFs\LabFRF; % control FRF
    % H - No x Ni x Nf matrix of FRFs in the lab configuration
load ..\LargeFiles\LabStressFRFs.mat; % lab stress FRFs to calculate stress in lab
    % Hs - Nstr x 1 cell array of stress FRFs in the lab configuration
    % Each stress FRF is 6 x Ni x Nf where 6 is the number of stress
    % components.

% %%%%%%% Simulation Settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
shakers = 1:6; % which shakers in H to use in control
ctrl_accs = 1:27; % control channels (baseplate accels, like a 6DOF test)
ref_accs = 67:69; % reference channels (a random triax on the DUT)

% Beginning of computations
Sxx_est = zeros(size(Sxx)); % lab environment

nsh = length(shakers); % number of shakers
Sff = zeros(nsh,nsh,length(fs)); % lab force PSD matrix

for ii = 1:length(fs) % calculate shaker forces and lab env at each fline
    cnthresh = 0.01*max(svd(H(ctrl_accs,:,ii))); % condition number threshold = 0.01 * largest SV of FRF mat
    % following line calculates forces while implementing CN threshold. See
    % pinv.m documentation for more detail.
    Sff(:,:,ii) = pinv(H(ctrl_accs,shakers,ii),cnthresh)*Sxx(ctrl_accs,ctrl_accs,ii)*pinv(H(ctrl_accs,shakers,ii),cnthresh)';
    Sxx_est(:,:,ii) = H(:,shakers,ii)*Sff(:,:,ii)*H(:,shakers,ii)'; % calculate lab response at all DOF
end

% Create plot comparing spectra:
plotTriaxPSD(fs,Sxx,Sxx_est,ref_accs)

% Stress analysis
ind1 = find(fs >= 95,1);
ind2 = find(fs >= 2000,1);
rms_inds = ind1:ind2; % calculate RMS stress only including freqs above 100 Hz
df = 5; % 5 Hz frequency spacing

[sigrms_lab,sigpsd_lab,sigloc_lab] = GetStressFunc(Hs,Sff,df,rms_inds,shakers); % calculate max lab VM stress PSD

load ..\Environment\FlightStressPSD; % max flight VM stress PSD

% compare flight and lab stress PSDs
figure;
semilogy(fs,abs(sigpsd_fl),'k',fs,abs(sigpsd_lab),'b','Linewidth',2)
grid on;
xlabel('Frequency (Hz)')
ylabel('Stress PSD (Pa^2/Hz)')
legend('Flight','Lab Test')
xlim([10 2000])







