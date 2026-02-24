% Create FRF relating lab shaker force to acceleration
% requires FEM modes and locations of accel and shaker nodes from
% GetLabLocationsScript.m
clc;close all;clear all;
addpath ..\Functions\; % for cl_model function

load ..\ModeShapes\BARC_BaseplateModes; % load lab modes and accel and shaker locs
load ..\FRFs\LabAccelNodes.mat;
load ..\FRFs\LabShakerNodes;

% get lab FRF matrix
clc;
nmodes = 12; % 12 modes below 3,000 Hz (include residual effects for 2000 Hz bandwidth)
wns = 2*pi*fn(1:nmodes);
zts = 0.01*ones(nmodes,1); % assume 1 percent damping for elastic modes
zts(1:6) = 0; % no damping for rigid body modes

% get modes at desired locations (retains only selected degrees of freedom)
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

% 3 - accel locs (xyz for each  - triaxes)
phi_r = zeros(3*size(accel_nodes,1),nmodes);
c = 1;
for ii = 1:size(accel_nodes,1)
    phi_r(c:c+2,:) = phi(accel_nodes(ii,2)*6-5:accel_nodes(ii,2)*6-3,1:nmodes);
    c = c + 3;
end

Acl= {phi_r,phi_sh}; % residue matrix cell

fs = 0:5:3000; % frequency vector of interest
ws = 2*pi*fs;
FRF_type = 2; % acceleration div. by force (see cl_model doc.)

H = permute(cl_model(wns, zts,Acl,ws,FRF_type),[2 3 1])/9.81; % convert to g

% save('LabFRF','H','fs')








