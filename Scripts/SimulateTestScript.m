% Simulate a MIMO test controlling to a flight environment using a defined
% set of shakers and control accelerometers. Computes and compares SLAP
% metrics for Control, Buzz, and Per-Shaker (PS) test configurations.
%
close all; clear all;
addpath ..\Functions;

load ..\LargeFiles\Rocket_Env; % flight environment
    % Sxx = [No x No x Nf] spectral density matrix of the flight environment on the DUT
    % fs  = [Nf x 1] frequency vector (Hz)
    % Node definitions: ../FRFs/Flight_Accel_Nodes.mat, ../FRFs/Flight_Force_Nodes.mat

load ..\FRFs\Lab_FRF; % control FRF
    % H   = [No x Ni x Nf] FRF matrix measured in the lab
    % fs  = [Nf x 1] frequency vector (Hz) — must match above
    H_lab = H; clear H;
    % DOFH encodes channel IDs as node.direction (e.g., 3.2 = node 3, Y-direction)
    DOFH = kron([1:39].',ones(3,1))+kron(ones(39,1),[1:3].'/10);
    % Accelerometer and shaker node locations:
    %   load('../FRFs/Lab_Accel_Nodes.mat')  → accel_nodes
    %   load('../FRFs/Lab_Shaker_Nodes.mat') → shaker_nodes
    %   load('..\ModeShapes\BARC_Baseplate_Modes.mat') → nodes

load ..\LargeFiles\Lab_Stress_FRFs.mat; % lab stress FRFs
    % Hs = {nstress x 1} cell array; each element is a [6 x nforces] FRF
    %   mapping force spectrum to stress tensor [s_xx, s_yy, s_zz, s_xy, s_xz, s_yz]

load ..\ModeShapes\BARC_Accel_Modes; % mode shapes for SLAP modal filter
    % phi = [(Nacc*3) x (6+Nfb)] — columns 1:6 are rigid-body, 7:end are fixed-base modes
    phi_filt  = phi;
    fb_inds   = 7:12;    % fixed-base mode columns in phi_filt
    filt_inds = 28:117;  % rows in H_lab corresponding to DUT accel channels (matches phi)

ctrl_accs = [1:3,7:9,19:21,25:27]; % control accelerometer channels
    % Accels 1:9 (three triax) capture 6DOF base-plate motion for a minimal test;
    % adding more channels up to all 27DOF improves conditioning. Check observability:
    disp('Condition number of rigid-body modes at control accels:');
    cond(phi_filt(ctrl_accs,1:6))

ref_accs = 67:69; % reference triax on the DUT (for monitoring, not control)
    %{
    % To visualize this accelerometer's location:
    node_ind = ref_accs(3)/3;
    FEM = load('..\ModeShapes\BARC_Baseplate_Modes.mat');
    AL  = load('../FRFs/Lab_Accel_Nodes.mat');
    figure(5);
    plot3(FEM.nodes(:,2),FEM.nodes(:,3),FEM.nodes(:,4),'.',...
        FEM.nodes(AL.accel_nodes(node_ind,2),2),FEM.nodes(AL.accel_nodes(node_ind,2),3),...
        FEM.nodes(AL.accel_nodes(node_ind,2),4),'ro'); axis equal
    %}

%% Compute Lab Response for 6DOF Control Test
Sxx_est = zeros(size(Sxx));

nsh     = 3; %size(H_lab,2); % number of shakers
sh_inds = 1:nsh;
Sff_lab = zeros(nsh,nsh,length(fs));

for ii = 1:length(fs)
    cnthresh = 0.01*max(svd(H_lab(ctrl_accs,1:nsh,ii))); % regularization at 1% of largest SV
    Sff_lab(:,:,ii) = pinv(H_lab(ctrl_accs,1:nsh,ii),cnthresh)*Sxx(ctrl_accs,ctrl_accs,ii)*pinv(H_lab(ctrl_accs,1:nsh,ii),cnthresh)';
    Sxx_est(:,:,ii) = H_lab(:,1:nsh,ii)*Sff_lab(:,:,ii)*H_lab(:,1:nsh,ii)';
end

envpsd_ctrl = get_psd(Sxx(ctrl_accs,ctrl_accs,:));     % flight PSD at control accels
labpsd_ctrl = get_psd(Sxx_est(ctrl_accs,ctrl_accs,:)); % lab PSD at control accels

envpsd = get_psd(Sxx(ref_accs,ref_accs,:));     % flight PSD at DUT ref accels
labpsd = get_psd(Sxx_est(ref_accs,ref_accs,:)); % lab PSD at DUT ref accels

%% Stress Analysis
ind1     = find(fs >= 95, 1);
ind2     = find(fs >= 2000, 1);
rms_inds = ind1:ind2; % RMS integration range (above 95 Hz to exclude low-freq noise)
df       = fs(2)-fs(1);

[sigrms_lab,sigpsd_lab,sigloc_lab] = GetStressFunc(Hs,Sff_lab,df,rms_inds,1:nsh);

load ..\Environment\Flight_Stress_PSD; % sigrms_fl, sigpsd_fl, sigloc_fl
disp(['Max RMS stress — Flight: ',num2str(sigrms_fl/1e3),' kPa,  Lab: ',num2str(sigrms_lab/1e3),' kPa']);
disp(['Lab/Flight stress ratio: ',num2str(sigrms_lab/sigrms_fl)]);

figure(2);
semilogy(fs,abs(sigpsd_fl),'k',fs,abs(sigpsd_lab),'b','Linewidth',2)
grid on;
xlabel('Frequency (Hz)')
ylabel('Stress PSD (Pa^2/Hz)')
legend('Flight','Lab Test')
title('\bfStress in Flight and Lab 6DOF Test');
xlim([20 2000])

%% SLAP-Control Test
bf = 7.3;    % fatigue exponent for aluminum (Larsen & Irvine, "Review of spectral fatigue methods")
Ts = [1;5];  % test and flight durations (relative)
p  = 0.99;   % confidence level for peak stress estimate

[scaling_control,metric_vals] = SLAPfunc(Sxx_est(filt_inds,filt_inds,:),Sxx(filt_inds,filt_inds,:),...
    phi_filt,fs,fb_inds,rms_inds,bf,Ts,p,[1:3]);
disp('SLAP Metrics-Control: [RMS Stress; Peak Stress; Fatigue; RMS FB modal resp.]');
scaling_control
metric_vals

figure(2); hold on;
semilogy(fs,abs(sigpsd_lab)*scaling_control,'--','Linewidth',2)
hold off;
legend('Flight','Lab Test','SLAP-Control Test')

%% SLAP-Buzz Test

% Initial forcing: unit autospectrum per shaker with low inter-shaker coherence
Sff_buzz = eye(nsh,nsh);
nf = length(fs);
coh = 0.05; % target inter-shaker coherence

for ii = 1:nsh-1
    for jj = ii+1:nsh
        phase_angle = 2*pi*(rand-0.5);
        Sfxfx = Sff_buzz(ii,ii);
        Sfyfy = Sff_buzz(jj,jj);
        re = sqrt(coh*Sfxfx*Sfyfy/(1+tan(phase_angle)^2));
        im = re*tan(phase_angle);
        Sfxfy = re+1i*im;
        Sff_buzz(ii,jj) = Sfxfy;
        Sff_buzz(jj,ii) = Sfxfy';
    end
end

Sxx_buzz = zeros(size(Sxx));
Sff_buzz_temp = zeros(nsh,nsh,nf);
for ii = 1:nf
    Sxx_buzz(:,:,ii)      = H_lab(:,sh_inds,ii)*Sff_buzz*H_lab(:,sh_inds,ii)';
    Sff_buzz_temp(:,:,ii) = Sff_buzz;
end
Sff_buzz = Sff_buzz_temp; clear Sff_buzz_temp;

buzzpsd     = get_psd(Sxx_buzz(ref_accs,ref_accs,:));
buzzctrlacc = get_psd(Sxx_buzz(ctrl_accs,ctrl_accs,:));

[scaling_buzz,buzz_metric_vals] = SLAPfunc(Sxx_buzz(filt_inds,filt_inds,:),Sxx(filt_inds,filt_inds,:),...
    phi_filt,fs,fb_inds,rms_inds,bf,Ts,p,[1:3]);
disp('SLAP-Buzz Metrics: [RMS Stress; Peak Stress; Fatigue; RMS FB modal resp.]');
disp(buzz_metric_vals)

[sigrms_buzz,sigpsd_buzz,sigloc_buzz] = GetStressFunc(Hs,Sff_buzz,df,rms_inds,1:nsh);

figure(2); hold on;
semilogy(fs,abs(sigpsd_buzz)*scaling_buzz,'-.','Linewidth',2)
hold off;
legend('Flight',[num2str(nsh),' DOF'],'SLAP-Control','SLAP-Buzz')

% Reference accelerometer responses
figure(3); set(gcf,'Units','normalized','Position',[0.1 0.23 0.4 0.25]);
tlt    = tiledlayout(1,3);
titles = {'Ref X','Ref Y','Ref Z'};
for ii = 1:3
    nexttile;
    semilogy(fs,abs(envpsd(ii,:)),'k',fs,abs(labpsd(ii,:)),'b',...
        fs,scaling_control*abs(labpsd(ii,:)),'--',fs,scaling_buzz*abs(buzzpsd(ii,:)),'Linewidth',2)
    xlabel('Frequency (Hz)','interpreter','tex')
    grid on; xlim([10 2000])
    if ii == 1
        legend('Flight',[num2str(nsh),' DOF'],'SLAP-Control','SLAP-Buzz')
    end
    title(titles{ii})
end
ylabel(tlt,'Acceleration PSD (g^2/Hz)','interpreter','tex')
tlt.TileSpacing = 'tight'; tlt.Padding = 'tight';

% Control accelerometer responses
figure(4); set(gcf,'Units','normalized','Position',[0.1 0.55 0.4 0.25]);
tlt    = tiledlayout(1,3);
titles = {'Control 1X','Control 1Y','Control 1Z'};
for ii = 1:3
    nexttile;
    semilogy(fs,abs(envpsd_ctrl(ii,:)),'k',fs,abs(labpsd_ctrl(ii,:)),'b',...
        fs,scaling_control*abs(labpsd_ctrl(ii,:)),'--',...
        fs,scaling_buzz*abs(buzzctrlacc(ii,:))); set(get(gca,'Children'),'LineWidth',2)
    xlabel('Frequency (Hz)','interpreter','tex')
    grid on; xlim([10 2000])
    if ii == 1
        legend('Flight',[num2str(nsh),' DOF'],'SLAP-Control','SLAP-Buzz')
    end
    title(titles{ii})
end
ylabel(tlt,'Acceleration PSD (g^2/Hz)','interpreter','tex')
tlt.TileSpacing = 'tight'; tlt.Padding = 'tight';

%% SLAP-PS Test (Per Shaker) — IN PROGRESS
% Adjusts shaker-by-shaker forcing to minimize shaker voltage while meeting the environment.

[metrics_fl] = GetMetrics(Sxx(filt_inds,filt_inds,:),phi,fs,fb_inds,rms_inds,bf,Ts(2),p);

return % ---- remainder is incomplete placeholder ----

    F_in_0 = diag(Sff_buzz(:,:,1)); % initial force autospectrum estimate

[sigrms_PS,sigpsd_PS,sigloc_PS] = GetStressFunc(Hs,Sff_PS,df,rms_inds,1:nsh);
