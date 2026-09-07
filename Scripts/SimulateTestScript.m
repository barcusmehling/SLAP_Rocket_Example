% Simulate a MIMO test controlling to a flight environment with a defined
% set of shakers, using a set of control accelerometers.
%
close all;clear all;
addpath ..\Functions; 

load ..\LargeFiles\Rocket_Env; % flight environment
    % Sxx = [No x No x Nf] Spectral Density Matrix defining the flight
    %   environment on the DUT. 
    % fs = [Nf x 1] Frequency vector (Hz)
    % The FEM nodes included in No outputs are defined in the following
    % files, in the same format as explained below for the lab nodes.
    % ../FRFs/Flight_Accel_Nodes.mat and ../FRFs/Flight_Force_Nodes.mat
load ..\FRFs\Lab_FRF; % control FRF
    % H = [No x Ni x Nf] Frequency Response Function Matrix measured in the
    %   Lab, which is used to generate the control signals for the [Ni]
    %   shakers. [(# outputs) x (# inputs) x (# frequency lines)]
    % fs = [Nf x 1] Frequency vector (Hz) (must match that above)
    H_lab = H; clear H;
    % The accelerometers are placed in the FEM nodes in "accel_nodes" in:
    %   load('../FRFs/Lab_Accel_Nodes.mat')
    % Create DOF vector to keep track of these.
    DOFH = kron([1:39].',ones(3,1))+kron(ones(39,1),[1:3].'/10);
    % The shaker nodes and directions are defined in "shaker_nodes"
    %   load('../FRFs/Lab_Shaker_Nodes.mat')
    % The node locations are defined in "nodes" in:
    %   load('..\ModeShapes\BARC_Baseplate_Modes.mat');
load ..\LargeFiles\Lab_Stress_FRFs.mat; % lab stress FRFs to calculate stress in lab
    % Hs = (nstresspoints x 1) cell array, where each element contains a [6 x nforces]
    %   frequency response matrix relating force spectrum to the spectrum of
    %   the [6x1] stress tensor: [s_xx, s_yy, s_zz, s_xy, s_xz, s_yz] where s
    %   denotes the stress.
% Load the fixed base modes of the DUT for the SLAP metrics:
load ..\ModeShapes\BARC_Accel_Modes; % modes for modal filtering (needed to do SLAP)
    % phi - (Nacc*3) x (6 + Nfb)
    %   First six columns are rigid body modes,
    %   The remaining columns are fixed-interface modes
    %   See above for how to correlate the accel locations to the FEM nodes
phi_filt = phi; % use 6 rigid body + 6 fixed-base modes in modal filter
    fb_inds = 7:12; % these columns are fixed-base modes in phi_filt
    filt_inds = 28:117; % Indicies in H_fl that correspond to the DUT accel channels.
        % These are the same as the (Nacc*3) channels in phi from BARCAccelModes

ctrl_accs = [1:3,7:9,19:21,25:27];%1:27; % control channels
    % Using accels 1:9 simulates a minimal 6DOF test, with three triax accels
    % that capture the 6DOF motion of the base plate.  One can add more
    % channels, up to all 27DOF (9 accels) on the base plate. The 6DOF test
    % will not be successful unless the accels chosen capture the rigid
    % body motion of the base plate; this is checked below.
    disp('Condition Number of Rigid Body Modes with Control Accels:');
    cond(phi_filt(ctrl_accs,1:6))
ref_accs = 67:69; % reference channels (a triax on the DUT)
    % To see the location of this accelerometer, one could use:
    %{
    node_ind=ref_accs(3)/3; % 3 DOF per node
    FEM=load('..\ModeShapes\BARC_Baseplate_Modes.mat')
    AL=load('../FRFs/Lab_Accel_Nodes.mat')
    
    figure(5);
    plot3(FEM.nodes(:,2),FEM.nodes(:,3),FEM.nodes(:,4),'.',...
        FEM.nodes(AL.accel_nodes(node_ind,2),2),FEM.nodes(AL.accel_nodes(node_ind,2),3),...
        FEM.nodes(AL.accel_nodes(node_ind,2),4),'ro'); axis equal
    %}

%% Compute the Response of the DUT when controlled in the 6DOF Lab Test
Sxx_est = zeros(size(Sxx)); % environment obtained in the lab

nsh = 3; %size(H_lab,2); % number of shakers - Use all by default.
sh_inds = 1:nsh; % Select which of the potential shaker locations to use
Sff_lab = zeros(nsh,nsh,length(fs)); % lab force PSD matrix

for ii = 1:length(fs) % calculate shaker forces and lab env at each frequency
    cnthresh = 0.01*max(svd(H_lab(ctrl_accs,1:nsh,ii))); % condition number threshold = 0.01 * largest SV of FRF mat
    % following line calculates forces while implementing CN threshold. See
    % pinv.m documentation for more detail.
    Sff_lab(:,:,ii) = pinv(H_lab(ctrl_accs,1:nsh,ii),cnthresh)*Sxx(ctrl_accs,ctrl_accs,ii)*pinv(H_lab(ctrl_accs,1:nsh,ii),cnthresh)';
    Sxx_est(:,:,ii) = H_lab(:,1:nsh,ii)*Sff_lab(:,:,ii)*H_lab(:,1:nsh,ii)'; % calculate lab response at all DOF
end

% Plot response at a control accelerometer
envpsd_ctrl = get_psd(Sxx(ctrl_accs,ctrl_accs,:)); % flight environment PSDs at control accels
labpsd_ctrl = get_psd(Sxx_est(ctrl_accs,ctrl_accs,:)); % '' lab env

% Extract diagonal terms from SDMs to plot only those PSDs below.
envpsd = get_psd(Sxx(ref_accs,ref_accs,:)); % flight environment PSDs at DUT ref accels
labpsd = get_psd(Sxx_est(ref_accs,ref_accs,:)); % '' lab env

%% Stress analysis
ind1 = find(fs >= 95,1);
ind2 = find(fs >= 2000,1);
rms_inds = ind1:ind2; % calculate RMS stress only including freqs above 100 Hz
df = fs(2)-fs(1); % frequency spacing

% Find the true peak RMS stress in the lab, as well as the PSD at the
% location where it occurs:
[sigrms_lab,sigpsd_lab,sigloc_lab] = GetStressFunc(Hs,Sff_lab,df,rms_inds,1:nsh); % calculate max lab VM stress PSD

% Load the true peak RMS stress from the flight.
load ..\Environment\Flight_Stress_PSD; % max flight VM stress PSD
    % Contains: sigrms_fl,sigpsd_fl,sigloc_fl
disp(['Max RMS stress in Flight: ',num2str(sigrms_fl/1e3),' kPa, and Lab: ',num2str(sigrms_lab/1e3), ' kPa']);
disp(['Ratio of Max RMS stresses: (lab/flight): ',num2str(sigrms_lab/sigrms_fl)]);

% Plot to compare flight and lab stress PSDs
figure(2);
semilogy(fs,abs(sigpsd_fl),'k',fs,abs(sigpsd_lab),'b','Linewidth',2)
grid on;
xlabel('Frequency (Hz)')
ylabel('Stress PSD (Pa^2/Hz)')
legend('Flight','Lab Test')
title('\bfStress in Flight and Lab 6DOF Test');
xlim([20 2000])

%% Simulate a SLAP-Control test
bf = 7.3; % fatigue exponent for Aluminum (acc. to Larsen and Irvine "Review of spectral fatigue methods"...)
Ts = [1;5]; % relative durations of test and flight, respectively. 
p = 0.99; % 99 percent confidence level - "we have 99 percent confidence that the peak value will be below the value we calculate"

% Compute the SLAP damage metrics and the scaling required to make the
% test conservative.
[scaling_control,metric_vals] = SLAPfunc(Sxx_est(filt_inds,filt_inds,:),Sxx(filt_inds,filt_inds,:),...
    phi_filt,fs,fb_inds,rms_inds,bf,Ts,p,[1:3]); % Apply SLAP
disp('SLAP Metrics-Control: [RMS Stress; Peak Stress; Fatigue, RMS FB modal resp.]');
scaling_control
metric_vals

% Add line for SLAP lab stress PSD - in practice one wouldn't know this,
% only the response, such as those shown in figure(1).
figure(2);
hold on;
semilogy(fs,abs(sigpsd_lab)*scaling_control,'--','Linewidth',2)
hold off;
legend('Flight','Lab Test','SLAP-Control Test')

%% Simulate a SLAP-Buzz test

    %%%%%%%%%%%%%%%% SLAP-Buzz Environment %%%%%%%%%%%%%%%%%%%%%%%%%%%
    Sff_buzz = eye(nsh,nsh); % all diagonal terms 1 N^2/Hz (can change this to scale shakers preferentially! This is just a starting point)
    %%%%%%%%% MSA Hack in here and update forcing per SLAP-PS %%%%%%%%%%%%
        % Sff_buzz = diag([0.005e-5,5e-5,5e-5,5e-5,5e-5,5e-5]);
    nf = length(fs);

    % Create off-diagonal terms for forcing spectral density matrix.
    % Set coherence to be small : 0.05, as inputs are uncorrelated in buzz 
    % test, randomize phase between shakers uniformly between -pi and pi,
    % assumes same force at each frequency.
    coh = 0.05;
    for ii = 1:nsh-1 % increment through all top nondiagonal terms (bottom ones are complex conjugates)
        for jj = ii+1:nsh
            phase_angle = 2*pi*(rand-0.5);
            Sfxfx = Sff_buzz(ii,ii);
            Sfyfy = Sff_buzz(jj,jj);
            re = sqrt(coh*Sfxfx*Sfyfy/(1+tan(phase_angle)^2)); % real part
            im = re*tan(phase_angle); % imaginary part
            Sfxfy = re+1i*im;
            Sff_buzz(ii,jj) = Sfxfy;
            Sff_buzz(jj,ii) = Sfxfy';
        end
    end
    
    % create buzz test environment
    Sxx_buzz = zeros(size(Sxx));
    Sff_buzz_temp = zeros(nsh,nsh,nf); % define at all flines
    for ii = 1:nf
        Sxx_buzz(:,:,ii) = H_lab(:,sh_inds,ii)*Sff_buzz*H_lab(:,sh_inds,ii)';
        Sff_buzz_temp(:,:,ii) = Sff_buzz; 
    end
    Sff_buzz = Sff_buzz_temp; % replace original with this version, which has values at all frequency lines.
        clear Sff_buzz_temp
    
    % Pull out Buzz PSD for plotting
    buzzpsd = get_psd(Sxx_buzz(ref_accs,ref_accs,:)); % '' lab env
    buzzctrlacc = get_psd(Sxx_buzz(ctrl_accs,ctrl_accs,:)); % '' lab env
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute the SLAP damage metrics and the scaling required to make the
% test conservative.
[scaling_buzz,buzz_metric_vals] = SLAPfunc(Sxx_buzz(filt_inds,filt_inds,:),Sxx(filt_inds,filt_inds,:),...
    phi_filt,fs,fb_inds,rms_inds,bf,Ts,p,[1:3]); % Apply SLAP
disp('SLAP-Buzz Metrics: [RMS Stress; Peak Stress; Fatigue, RMS FB modal resp.]');
disp(buzz_metric_vals)

% Find the true peak RMS stress in the lab, as well as the PSD at the
% location where it occurs:
[sigrms_buzz,sigpsd_buzz,sigloc_buzz] = GetStressFunc(Hs,Sff_buzz,df,rms_inds,1:nsh); % calculate max lab VM stress PSD

% Add line for SLAP lab stress PSD - in practice one wouldn't know this,
% only the response, such as those shown in figure(1).
figure(2);
hold on;
semilogy(fs,abs(sigpsd_buzz)*scaling_buzz,'-.','Linewidth',2)
hold off;
legend('Flight',[num2str(nsh), ' DOF'],'SLAP-Control','SLAP-Buzz')

% Plot reference responses in flight vs. lab
figure(3); set(gcf,'Units','normalized','Position',[0.1     0.23          0.4         0.25]);
tlt = tiledlayout(1,3);
titles = {'Ref X','Ref Y','Ref Z'};
for ii = 1:3
    nexttile;
    semilogy(fs,abs(envpsd(ii,:)),'k',fs,abs(labpsd(ii,:)),'b',...
        fs,scaling_control*abs(labpsd(ii,:)),'--',fs,scaling_buzz*abs(buzzpsd(ii,:)),'Linewidth',2)
    xlabel('Frequency (Hz)','interpreter','tex')
    grid on;
    xlim([10 2000])
    if ii == 1
        legend('Flight',[num2str(nsh), ' DOF'],'SLAP-Control','SLAP-Buzz')
    end
    title(titles{ii})
end
ylabel(tlt,'Acceleration PSD (g^2/Hz)','interpreter','tex')
tlt.TileSpacing = 'tight';
tlt.Padding = 'tight';

% Plot reference responses in flight vs. lab
figure(4); set(gcf,'Units','normalized','Position',[0.1  0.55  0.4  0.25]);
tlt = tiledlayout(1,3);
titles = {'Control 1X','Control 1Y','Control 1Z'};
for ii = 1:3
    nexttile;
    semilogy(fs,abs(envpsd_ctrl(ii,:)),'k',fs,abs(labpsd_ctrl(ii,:)),'b',...
        fs,scaling_control*abs(labpsd_ctrl(ii,:)),'--',...
        fs,scaling_buzz*abs(buzzctrlacc(ii,:))); set(get(gca,'Children'),'LineWidth',2)
    xlabel('Frequency (Hz)','interpreter','tex')
    grid on;
    xlim([10 2000])
    if ii == 1
        legend('Flight',[num2str(nsh), ' DOF'],'SLAP-Control','SLAP-Buzz')
    end
    title(titles{ii})
end
ylabel(tlt,'Acceleration PSD (g^2/Hz)','interpreter','tex')
tlt.TileSpacing = 'tight';
tlt.Padding = 'tight';    

%% Simulate a SLAP-PS Test "Per Shaker"
% Adjust shaker-by-shakerk to minimize shaker voltage while also meeting
% environment.

% Compute the metrics for the flight environment
[metrics_fl] = GetMetrics(Sxx(filt_inds,filt_inds,:),phi,fs,fb_inds,rms_inds,bf,Ts(2),p);

return
    % Initial estimate for force autospectra
    F_in_0 = diag(Sff_buzz(:,:,1)); % all frequency lines are identical

    


% Find the true peak RMS stress in the lab, as well as the PSD at the
% location where it occurs:
[sigrms_PS,sigpsd_PS,sigloc_PS] = GetStressFunc(Hs,Sff_PS,df,rms_inds,1:nsh); % calculate max lab VM stress PSD



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

