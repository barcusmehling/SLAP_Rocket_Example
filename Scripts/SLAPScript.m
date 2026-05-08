% Simulate many SLAP tests and make plots for UNSGC 2026 paper
% Marcus Behling | 2/18/2026
clc;close all;clear all;
%% Load Files (Takes ~30 seconds because files big)
addpath ..\Functions\;

disp('Loading Flight FRF...')
load ..\FRFs\Flight_FRF; % load and save flight acceleration/force FRF
%   H - No x Ni x Nf - (# outputs) x (# inputs) x (# frequency lines)
%   fs - frequency vector in Hz (Nf x 1)
% MSA: Need to also load something that tells us what these points are!
% %%%%%%%%%
H_fl = H; clear H;
nf = length(fs);
ws = 2*pi*fs;
df = fs(2)-fs(1);

load ..\Environment\Flight_Forces; % nominal force spectra for the 8 forcing vectors on rocket - starting point for making environments
%   fmat - Ni x Nf matrix of nominal force spectra

disp('Loading Lab FRF...')
load ..\FRFs\Lab_FRF; % control FRF
H_lab = H; clear H;
sh_inds = 1:7; % Select which of the potential shaker locations to use
nsh = length(sh_inds); % number of shakers

disp('Loading Flight Stress FRFs...')
load ..\LargeFiles\Flight_Stress_FRFs.mat; % FRF to calculate stress in flight env
Hs_fl = Hs; clear Hs;

disp('Loading Lab Stress FRFs...')
load ..\LargeFiles\Lab_Stress_FRFs.mat; % FRF to calculate stress in lab
Hs_lab = Hs; clear Hs;

load ..\ModeShapes\BARC_Accel_Modes; % modes for modal filtering (needed to do SLAP)
%   phi - (Nacc*3) x (6 + Nfb)
%       First six columns are rigid body modes,
%       The remaining columns are fixed-interface modes
%       ADD DOCUMENTATION - ACCEL SET??
phi_filt = phi; % use 6 rigid body + 6 fixed-base modes in modal filter
    fb_inds = 7:12; % these columns are fixed-base modes in phi_filt

% Define indicies to match between fixed-interface modes and flight
nacc = size(H_fl,1); % number of accel channels in H_fl
    ctrl_inds = 1:27; % Indicies in H_fl that correspond to control channels
    filt_inds = 28:117; % Indicies in H_fl that correspond to the DUT accel channels.
        % These are the same as the (Nacc*3) channels in phi from BARCAccelModes
    rms_inds = 51:401; % frequency indicies (in fs) to use when computing RMS stress
    % 51:401 uses everything 250 and 2000 Hz.  (Motions below 250 Hz can be
    % large in magnitude but don't cause stress because the DUT is basically rigid)

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

ctrl = 0; % SLAP-Control if 1, SLAP-Buzz if 0

%%%%%%%%%%%%%%%% SLAP-Buzz params %%%%%%%%%%%%%%%%%%%%%%%%%%%
Sff_lab = eye(nsh,nsh); % all diagonal terms 1 N^2/Hz (can change this to scale shakers preferentially! This is just a starting point)

% create off-diagonal terms 
% set coherence to be small : 0.05, as inputs uncorrelated in buzz test
% randomize phase between shakers uniformly between -pi and pi
% assumes same force at each fline - could change this as well...
coh = 0.05;
for ii = 1:nsh-1 % increment through all top nondiagonal terms (bottom ones are complex conjugates)
    for jj = ii+1:nsh
        phase_angle = 2*pi*(rand-0.5);
        Sfxfx = Sff_lab(ii,ii);
        Sfyfy = Sff_lab(jj,jj);
        re = sqrt(coh*Sfxfx*Sfyfy/(1+tan(phase_angle)^2)); % real part
        im = re*tan(phase_angle); % imaginary part
        Sfxfy = re+1i*im;
        Sff_lab(ii,jj) = Sfxfy;
        Sff_lab(jj,ii) = Sfxfy';
    end
end

% create buzz test environment
Sxx_lab = zeros(nacc,nacc,nf);
Sff_lab_temp = zeros(nsh,nsh,nf); % define at all flines
for ii = 1:nf
    Sxx_lab(:,:,ii) = H_lab(:,sh_inds,ii)*Sff_lab*H_lab(:,sh_inds,ii)';
    Sff_lab_temp(:,:,ii) = Sff_lab; 
end
Sff_lab = Sff_lab_temp; % swap with original to define at all flines
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for ii = 1:nsims
    disp(['Run ' num2str(ii) ' of ' num2str(nsims)])
    Sxx_fl = zeros(nacc,nacc,length(fs)); % initialize flight and lab environments and force PSD
    Sff_fl = zeros(nflforces,nflforces,length(fs));

    if ctrl==1
        Sxx_lab = zeros(nacc,nacc,length(fs));
        Sff_lab = zeros(nsh,nsh,length(fs));
    end

    % Create the environment for this simulation.
    fmatrand = 1/3.04*(10.^(2*rand(size(fmat))-1)+1i*10.^(2*rand(size(fmat))-1)); % Make new forces for new env
    fmatfin = fmat.*fmatrand; % "final" force vectors to apply to make flight env

    for jj = 1:length(fs)
        Sff_fl(:,:,jj) = (fmatfin(:,jj)*fmatfin(:,jj)'); % calculate flight force PSD
        Sxx_fl(:,:,jj) = squeeze(H_fl(:,:,jj))*Sff_fl(:,:,jj)*squeeze(H_fl(:,:,jj)');
        
        if ctrl == 1
            cnthresh = 0.01*max(svd(H_lab(ctrl_inds,sh_inds,jj))); % condition number threshold
            Sff_lab(:,:,jj) = pinv(H_lab(ctrl_inds,sh_inds,jj),cnthresh)*Sxx_fl(ctrl_inds,ctrl_inds,jj)*pinv(H_lab(ctrl_inds,sh_inds,jj)',cnthresh); % calculate lab shaker forces
            Sxx_lab(:,:,jj) = H_lab(:,sh_inds,jj)*Sff_lab(:,:,jj)*H_lab(:,sh_inds,jj)'; % calculate lab environment
        end
    end

    [sigrms_fl,sigpsd_fl,sigloc_fl] = GetStressFunc(Hs_fl,Sff_fl,df,rms_inds,fl_force_inds); % RMS stress and stress PSD in flight env
    [sigrms_lab,sigpsd_lab,sigloc_lab] = GetStressFunc(Hs_lab,Sff_lab,df,rms_inds,sh_inds); % RMS stress and stress PSD in flight env

    [scaling,metric_vals(ii,:)] = SLAPfunc(Sxx_lab(filt_inds,filt_inds,:),Sxx_fl(filt_inds,filt_inds,:),phi_filt,fs,fb_inds,rms_inds,bf,Ts,p); % Apply SLAP

    % Create a plot comparing the spectra at a point:
        % {
        ref_accs = 67:69; % reference channels (a random triax on the DUT)
        plotTriaxPSD(fs,Sxx_fl,Sxx_lab*scaling,ref_accs)
        pause
        %}

    sig_lab_rms = sigrms_lab*sqrt(scaling); % control
    sigpsd_lab = sigpsd_lab*scaling;
    stress_ratio = sig_lab_rms / sigrms_fl; % RMS stress ratio
    fatigue_ratio = GetFatigueRatio(sigpsd_lab,sigpsd_fl,bf,fs,rms_inds,Ts); % fatigue ratio
    peak_ratio = GetPeakStressRatio(sigpsd_lab,sigpsd_fl,fs,rms_inds,Ts,p);

    actual_metrics(ii,:) = [peak_ratio stress_ratio fatigue_ratio];
    scales(ii) = scaling;
    
end

readme = 'SLAP-Control, shakers 1-7, checking if I get consistent results.';
save('..\Results\Results_SLAP_Control','actual_metrics','scales','metric_vals','readme')
%% Load and plot results (need to have variables before big for loop above loaded in)
load ..\Results\Results_SLAP_Control.mat;

titles = {'Peak Stress Ratio','RMS Stress Ratio','Fatigue Damage Ratio'};
nbins = 16;

figure('Units','normalized','Position',[0.1 0.1 0.8 0.4]);   % [left bottom width height]
tlt = tiledlayout(1,3);
for jj = 1:3
    nexttile;
    h = histogram(rmoutliers(actual_metrics(:,jj)),'NumBins',nbins);
    h.FaceColor = 'r';
    h.FaceAlpha = 1;  % Fully opaque
    h.EdgeColor = 'k';  % Black edges for contrast
    xlabel(titles{jj})
    if jj == 1
        ylabel('Count')
    end
    grid on;
end

tlt.Padding = 'compact';
tlt.TileSpacing = 'tight';
sgtitle('Actual Damage Metrics: SLAP-Control, 4 Shakers')

figure('Units','normalized','Position',[0.1 0.1 0.8 0.4]);   % [left bottom width height]
tlt = tiledlayout(1,3);
for jj = 1:3
    nexttile;
    if jj < 3
        h = histogram(rmoutliers(actual_metrics(:,jj)./sqrt(scales)),'NumBins',nbins);
    else
        h = histogram(rmoutliers(actual_metrics(:,jj)./(scales.^(bf/2))),'NumBins',nbins);
    end
    h.FaceColor = 'r';
    h.FaceAlpha = 1;  % Fully opaque
    h.EdgeColor = 'k';  % Black edges for contrast
    xlabel(titles{jj})
    if jj == 1
        ylabel('Count')
    end
    grid on;
end

tlt.Padding = 'compact';
tlt.TileSpacing = 'tight';
sgtitle('Actual Damage Metrics: Unscaled 4 Shaker Control')
