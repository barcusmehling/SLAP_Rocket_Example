% Calculate FRFs relating forces to stress in flight and lab 
% Flight
clc;close all;clear all;
addpath ..\Functions\;
load ..\ModeShapes\FullRocketModes; % load flight disp. modes, force locs, and stress modes
load ..\FRFs\FlightForceNodes.mat;
load ..\ModeShapes\FlightStressModes.mat;

nmodes = length(fn); % 91 below 3000 Hz
mode_inds = 1:nmodes;
fb_inds = [15 40 46 47 52 80]; % DUT deformation modes

zts = 0.05*ones(nmodes,1); % rocket elastic modes 5 percent damping
zts(1:6) = 0; % no rigid body damping
zts(fb_inds) = 0.01; % BARC modes 1 percent damping

wns = 2*pi*fn; % natural frequencies and freq vector for FRF
fs = 0:5:3000;
ws = 2*pi*fs;

% get modes at desired locations
% 1 - flight force ("shaker") locations
phi_sh = zeros(size(shaker_nodes,1),nmodes);
for ii = 1:size(shaker_nodes,1)
    if shaker_nodes(ii,3)~=0
        phi_sh(ii,:) = sign(shaker_nodes(ii,3))*phi(shaker_nodes(ii,2)*6-5,1:nmodes);
    elseif shaker_nodes(ii,4)~=0
        phi_sh(ii,:) = sign(shaker_nodes(ii,4))*phi(shaker_nodes(ii,2)*6-4,1:nmodes);
    else
        phi_sh(ii,:) = sign(shaker_nodes(ii,5))*phi(shaker_nodes(ii,2)*6-3,1:nmodes);
    end
end

n_stress_dof = size(psi,1); % 459 stress els
Hs = cell(n_stress_dof,1); % store FRFs for each element in a cell array
FRF_type = 1; % displacement

for ii = 1:length(Hs) % calculate FRF for each element
    Acl = {squeeze(psi(ii,:,:)),phi_sh};
    H = permute(cl_model(wns, zts, Acl,ws,FRF_type),[2 3 1]); 
    Hs{ii} = H; % save to cell
end

% save('FlightStressFRFs','Hs')

% Calculate stress PSDs in flight for the fun of it

fvec = 0.5*logspace(0,-2,601);

fmat = [fvec;fvec;fvec;fvec;fvec;fvec;fvec;fvec];

fmat(1:2,:) = 0.2;

% create Sff
Sff = zeros(size(fmat,1),size(fmat,1),size(fmat,2));
df = 5; % frequency spacing (Hz)

for ii = 1:length(fs)
    Sff(:,:,ii) = (fmat(:,ii)*fmat(:,ii)')/df; % force PSD
end

ind1 = find(fs >= 95,1); % indices to use in RMS
ind2 = find(fs >= 2000,1);
rms_inds_fl = ind1:ind2;

[sigrms_fl,sigpsd_fl,sigloc_fl] = GetStressFunc(Hs,Sff,df,rms_inds_fl,1:8);

% save('FlightStressPSD','sigrms_fl','sigpsd_fl','sigloc_fl','rms_inds_fl')

% Lab
clc;close all;clear all;
load ..\ModeShapes\BARC_BaseplateModes;
load ..\FRFs\LabShakerNodes
load ..\ModeShapes\LabStressModes.mat;

% get FRF matrix
clc;close all;
nmodes = 12; % 12 modes below 3,000 Hz
wns = 2*pi*fn(1:nmodes);
zts = 0.01*ones(nmodes,1); % assume 1 percent damping 
zts(1:6) = 0; % no damping for rigid body modes

fs = 0:5:3000;
ws = 2*pi*fs;

% get modes at desired locations
% 1 - shaker locations
phi_sh = zeros(size(shaker_nodes,1),nmodes);
for ii = 1:size(shaker_nodes,1)
    if shaker_nodes(ii,3)~=0
        phi_sh(ii,:) = sign(shaker_nodes(ii,3))*phi(shaker_nodes(ii,2)*6-5,1:nmodes);
    elseif shaker_nodes(ii,4)~=0
        phi_sh(ii,:) = sign(shaker_nodes(ii,4))*phi(shaker_nodes(ii,2)*6-4,1:nmodes);
    else
        phi_sh(ii,:) = sign(shaker_nodes(ii,5))*phi(shaker_nodes(ii,2)*6-3,1:nmodes);
    end
end

n_stress_dof = size(psi,1);
Hs = cell(n_stress_dof,1);
FRF_type = 1; % displacement

for ii = 1:length(Hs)
    Acl = {squeeze(psi(ii,:,1:nmodes)),phi_sh};
    H = permute(cl_model(wns, zts, Acl,ws,FRF_type),[2 3 1]); 
    Hs{ii} = H;
end

% save('LabStressFRFs','Hs')








