% Get Rigid Body + Fixed-Base Modes at BARC Accel Locs
clc;close all;clear all;
load ..\FRFs\LabAccelNodes;
barc_accelnodes = accel_nodes(10:end,2);

load ..\ModeShapes\BARC_BaseplateModes;

phi_acc_rb = zeros(3*size(barc_accelnodes,1),6);
c = 1;
for ii = 1:size(barc_accelnodes,1)
    phi_acc_rb(c:c+2,:) = phi(barc_accelnodes(ii)*6-5:barc_accelnodes(ii)*6-3,1:6);
    c = c + 3;
end

clear phi;
load ..\ModeShapes\BARC_FixedBaseModes;

phi_acc_el = zeros(3*size(barc_accelnodes,1),size(phi,2));
c = 1;
for ii = 1:size(barc_accelnodes,1)
    phi_acc_el(c:c+2,:) = phi(barc_accelnodes(ii)*6-5:barc_accelnodes(ii)*6-3,:);
    c = c + 3;
end

phi_acc = [phi_acc_rb phi_acc_el];
clear phi;

phi = phi_acc;
readme = '6 Rigid Body + Fixed-Base Modes';

% save('BARCAccelModes','phi','readme')