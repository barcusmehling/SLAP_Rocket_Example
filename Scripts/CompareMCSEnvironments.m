% Compare the flight environments used in the SLAPScript for the Monte
% Carlo Simulations by overlaying response spectrum at the center of the
% plate for each.
% M. Allen, May 2026
clc;close all;clear all;
%% Load Files (Takes ~30 seconds because files are big)
addpath ..\Functions\;

disp('Loading Flight FRF...')
load ..\FRFs\FlightFRF; % load and save flight acceleration/force FRF
%   H - No x Ni x Nf - (# outputs) x (# inputs) x (# frequency lines)
%   fs - frequency vector in Hz (Nf x 1)
% MSA: Need to also load something that tells us what these points are!
% %%%%%%%%%
H_fl = H; clear H;
nf = length(fs);
ws = 2*pi*fs;
df = fs(2)-fs(1);

load ..\Environment\FlightForces; % nominal force spectra for the 8 forcing vectors on rocket - starting point for making environments
%   fmat - Ni x Nf matrix of nominal force spectra

load ..\LargeFiles\RocketEnv; % flight environment
    % Sxx - No x No x Nf matrix of spectral densities in flight.
    % fs - 1 x Nf frequency vector in Hz  
nacc = size(H_fl,1); % number of accel channels in H_fl
% Assume that the 5th accelerometer is the one in the center, plot that
ref_accs = 13:15;

%% Simulate Tests
nsims = 100; % number of simulations to perform (each takes 2-3 seconds on my computer...)

scales = zeros(nsims,1); % save SLAP scaling factors for each test
metric_vals = zeros(nsims,4); % SLAP metric values after scaling
actual_metrics = zeros(nsims,3); % actual elastic damage metric values after scaling

bf = 7.3; % fatigue exponent for Aluminum (acc. to Larsen and Irvine "Review of spectral fatigue methods"...)
Ts = [1;5]; % relative durations of test and flight, respectively. 
p = 0.99; % 99 percent confidence level - "we have 99 percent confidence that the peak value will be below the value we calculate"

nflforces = size(fmat,1); % number of flight forces
fl_force_inds = 1:nflforces;

Sxx_fl = zeros(nacc,nacc,length(fs),nsims); % initialize flight and lab environments and force PSD
Sff_fl = zeros(nflforces,nflforces,length(fs),nsims);

for ii = 1:nsims

    % Create the environment for this simulation.
    fmatrand = 1/3.04*(10.^(2*rand(size(fmat))-1)+1i*10.^(2*rand(size(fmat))-1)); % Make new forces for new env
    fmatfin = fmat.*fmatrand; % "final" force vectors to apply to make flight env

    for jj = 1:length(fs)
        Sff_fl(:,:,jj,ii) = (fmatfin(:,jj)*fmatfin(:,jj)'); % calculate flight force PSD
        Sxx_fl(:,:,jj,ii) = squeeze(H_fl(:,:,jj))*Sff_fl(:,:,jj,ii)*squeeze(H_fl(:,:,jj)');
        
    end
   
    % Plot reference responses in flight vs. lab
    if ii==1
        figure(1);
        directions = {'X','Y','Z'};
    end
    if ii==nsims; lcol='r'; else lcol='k'; end
    for n = 1:3
        subplot(3,1,n)
        semilogy(fs,abs(squeeze(Sxx_fl(ref_accs(n),ref_accs(n),:,ii))),lcol); hold on;
        if ii==1
            grid on;
            xlim([10 2000])
            ylabel([directions{n},'-PSD (g^2/Hz)'],'interpreter','tex')
            if n==3; xlabel('Frequency (Hz)','interpreter','tex'); end
        end
    end

end
% Add nominal on top
for n=1:3
    subplot(3,1,n)
    semilogy(fs,abs(squeeze(Sxx(ref_accs(n),ref_accs(n),:))),'b','Linewidth',2);
end

% % Find statistics of random environments
% Sxx_fl_mean=mean(abs(Sxx_fl),4);
% Sxx_fl_std=std(Sxx_fl,[],4); % This doesn't work well, would need to use
% "quantile" or something to get a better upper/lower bound.
% 
% % Plot comparison of PSDs
% figure(10);
% tlt = tiledlayout(3,1);
% %titles = {'Ref X','Ref Y','Ref Z'};
% 
% for ii = 1:3
%     nexttile;
%     semilogy(fs,abs(squeeze(Sxx_fl_mean(ref_accs(ii),ref_accs(ii),:))),'k',...
%         fs,abs(squeeze(Sxx_fl_mean(ref_accs(ii),ref_accs(ii),:))+2*squeeze(Sxx_fl_std(ref_accs(ii),ref_accs(ii),:))),'--k',...
%         fs,abs(squeeze(Sxx_fl_mean(ref_accs(ii),ref_accs(ii),:))-2*squeeze(Sxx_fl_std(ref_accs(ii),ref_accs(ii),:))),'--k',...
%         fs,abs(squeeze(Sxx(ref_accs(ii),ref_accs(ii),:))),'b','Linewidth',2)
%     xlabel('Frequency (Hz)','interpreter','tex')
%     grid on;
%     xlim([10 2000])
%     if ii == 1
%         legend('mean(Flight)','+2\sigma','-2\sigma','Nominal')
%     end
% end
% ylabel(tlt,'Acceleration PSD (g^2/Hz)','interpreter','tex')
% tlt.TileSpacing = 'tight';
% tlt.Padding = 'tight';