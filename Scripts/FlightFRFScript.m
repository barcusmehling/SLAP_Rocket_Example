% get FRF relating flight forces to accelerations 
clc; close all; clear;
addpath ..\Functions\; % for cl_model function

% Load rocket modes and accelerometer and force locations from GetRocketLocations.m
load ..\ModeShapes\FullRocketModes;
load ..\FRFs\FlightAccelNodes.mat;
load ..\FRFs\FlightForceNodes.mat; % note - these are labeled as "shaker_nodes"

% get flight FRF matrix
clc;close all;
nmodes = length(fn); % 91 modes below 3000 Hz
mode_inds = 1:nmodes;
fb_inds = [15 40 46 47 52 80]; % fixed-base DUT modes 

zts = 0.05*ones(nmodes,1); % rocket elastic modes 5 percent damping because of bolts and stuff in rocket
zts(1:6) = 0; % no damping for rigid body modes
zts(fb_inds) = 0.01; % BARC modes 1 percent damping

wns = 2*pi*fn; % natural frequencies in rad/s
fs = 0:5:3000; % frequency vector at which to define FRF (Hz)
ws = 2*pi*fs; % rad/s

FRF_type = 2; % acceleration (see cl_model.m)

% get modes at desired locations
% 1 - flight force ("shaker") locations
phi_sh = zeros(size(shaker_nodes,1),nmodes); % initialize matrix
for ii = 1:size(shaker_nodes,1) % go through each selected node
    if shaker_nodes(ii,3)~=0 % if X force
        phi_sh(ii,:) = sign(shaker_nodes(ii,3))*phi(shaker_nodes(ii,2)*6-5,1:nmodes); % save value from FEM to new phi matrix
    elseif shaker_nodes(ii,4)~=0 % if Y force
        phi_sh(ii,:) = sign(shaker_nodes(ii,4))*phi(shaker_nodes(ii,2)*6-4,1:nmodes);
    else % if Z force
        phi_sh(ii,:) = sign(shaker_nodes(ii,5))*phi(shaker_nodes(ii,2)*6-3,1:nmodes);
    end
end

% 2 - get modes at accel locations
phi_acc = zeros(3*size(accel_nodes,1),nmodes);
c = 1;
for ii = 1:size(accel_nodes,1)
    phi_acc(c:c+2,:) = phi(accel_nodes(ii,2)*6-5:accel_nodes(ii,2)*6-3,1:nmodes); % save X, Y, Z value for each accel (triaxes)
    c = c + 3;
end

Acl = {phi_acc,phi_sh}; % see cl_model.m documentation

H = permute(cl_model(wns, zts, Acl,ws,FRF_type),[2 3 1])/9.81; % get FRF and convert from m/s^2 to g

% save('FlightFRF','H','fs') 















