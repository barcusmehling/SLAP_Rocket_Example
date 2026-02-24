% Load modes from abaqus to MATLAB using .inp

clc;close all;clear all;
model = 'BARC_FixedBase';
strABQDir = cd;
controls.batchFile='abaqus.bat'; % Name or path/name to Abaqus *.bat file
controls.silent=true; 
controls.numCores=1; % Specify if using more than 1 core.
[phi,fn,dof,nodes,elestruct] = getModesAbq(model,strABQDir,controls);
uimplot5({nodes,elestruct.eles,fn,phi})
% save('BARC_FixedBaseModes','phi','fn','dof','nodes','elestruct')




