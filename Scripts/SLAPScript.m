% Simulate many SLAP tests and make plots for UNSGC 2026 paper
% Marcus Behling | 2/18/2026
clc;close all;clear all;
%% Load Files (Takes ~30 seconds because files big)
addpath ..\Functions\;

disp('Loading Flight FRF...')
load ..\FRFs\FlightFRF; % load and save flight force - acceleration FRF for making new environments
H_fl = H;
nf = length(fs);
ws = 2*pi*fs;
df = fs(2)-fs(1);
clear H;

load ..\Environment\FlightForces; % nominal forces for the 8 forcing vectors on rocket - starting point for making environments

disp('Loading Lab FRF...')
load ..\FRFs\LabFRF; % control FRF
H_lab = H;
sh_inds = 1:7;
nsh = length(sh_inds); % number of shakers
clear H;

disp('Loading Flight Stress FRFs...')
load ..\FRFs\FlightStressFRFs.mat; % FRF to calculate stress in flight env
Hs_fl = Hs;
clear Hs;

disp('Loading Lab Stress FRFs...')
load ..\FRFs\LabStressFRFs.mat; % FRF to calculate stress in lab
Hs_lab = Hs;
clear Hs;

load ..\ModeShapes\BARCAccelModes; % modes for modal filtering (needed to do SLAP)
nacc = size(H_fl,1); % number of accels
ctrl_inds = 1:27;
filt_inds = 28:117; % only use DUT accel channels in filtering
phi_filt = phi; % use 6 rigid body + 6 fixed-base modes in modal filter
fb_inds = 7:12; % these columns are fixed-base modes in phi mat
rms_inds = 51:401; % everything between 250 and 2000 Hz include in RMS stress (stuff below 250 Hz can be large in magnitude but doesn't cause stress because basically rigid)
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
clc;

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

    sig_lab_rms = sigrms_lab*sqrt(scaling); % control
    sigpsd_lab = sigpsd_lab*scaling;
    stress_ratio = sig_lab_rms / sigrms_fl; % RMS stress ratio
    fatigue_ratio = GetFatigueRatio(sigpsd_lab,sigpsd_fl,bf,fs,rms_inds,Ts); % fatigue ratio
    peak_ratio = GetPeakStressRatio(sigpsd_lab,sigpsd_fl,fs,rms_inds,Ts,p);

    actual_metrics(ii,:) = [peak_ratio stress_ratio fatigue_ratio];
    scales(ii) = scaling;
    
end

% readme = 'SLAP-Buzz, shakers 1-8, checking if I get consistent results.';
% save('Results_SLAPBuzz_2','actual_metrics','scales','metric_vals','readme')
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





