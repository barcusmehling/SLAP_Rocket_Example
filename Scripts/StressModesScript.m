% Extract stress modes from abaqus .csv files
% I exported these as csvs from abaqus because I was afraid to break the
% MATLAB functions
% 2/19/2026
% BARC Fixed-Base Stress Modes (I don't actually use these for anything at
% the moment...)
clc;close all;clear all;
data = readtable('BARC_FixedBaseStressModes.csv');
ndof = 459; % num. stress elements to calculate stress on. can check this in csv file.
nmodes = 6; % 6 fixed-base modes of the DUT 
psi = zeros(ndof,6,nmodes+1); % stress mode matrix - 6 (num stress tensor components) x nmodes for each stress el

c = 1;
for ii = 1:nmodes+1 % save stress tensor values to psi
    psi(:,:,ii) = [data.S_S11(c:c+ndof-1) data.S_S22(c:c+ndof-1) data.S_S33(c:c+ndof-1) data.S_S12(c:c+ndof-1) data.S_S13(c:c+ndof-1) data.S_S23(c:c+ndof-1)];
    c = c + ndof;
end

psi(:,:,1) = []; % get rid of fake 1st mode - doesn't actually correspond to a mode

% save('FixedBaseStressModes','psi')

%% Flight stress modes
clc;close all;clear all;
data = readtable('FlightStressModes.csv');
ndof = 459;
nmodes = 91;
psi = zeros(ndof,6,nmodes+1);

c = 1;
for ii = 1:nmodes+1
    psi(:,:,ii) = [data.S_S11(c:c+ndof-1) data.S_S22(c:c+ndof-1) data.S_S33(c:c+ndof-1) data.S_S12(c:c+ndof-1) data.S_S13(c:c+ndof-1) data.S_S23(c:c+ndof-1)];
    c = c + ndof;
end

psi(:,:,1) = [];

% save('FlightStressModes','psi')

%% Lab setup stress modes
clc;close all;clear all;
data = readtable('BARC_baseplate_StressModes.csv');
ndof = 459;
nmodes = 15;
psi = zeros(ndof,6,nmodes+1);

c = 1;
for ii = 1:nmodes+1
    psi(:,:,ii) = [data.S_S11(c:c+ndof-1) data.S_S22(c:c+ndof-1) data.S_S33(c:c+ndof-1) data.S_S12(c:c+ndof-1) data.S_S13(c:c+ndof-1) data.S_S23(c:c+ndof-1)];
    c = c + ndof;
end

psi(:,:,1) = [];

% save('BARCBaseplateStressModes','psi')










