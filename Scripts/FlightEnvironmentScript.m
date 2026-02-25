% Make flight environment using FRF from FlightFRFScript.m
clc;close all;clear all;
load ..\FRFs\FlightFRF;

% Define flight forcing PSDs
fvec = 0.5*logspace(0,-2,601); % make force PSDs linearly decay on loglog plot

fmat = [fvec;fvec;fvec;fvec;fvec;fvec;fvec;fvec]; % 8 times for 8 flight forces

fmat(1:2,:) = 0.2; % replace two launch forces with constant values rather than decay

% save('FlightForces','fmat')

% % Plot force PSDs
% figure;
% t = tiledlayout(4,2);
% chans = {'Force 1','Force 2','Force 3','Force 4','Force 5','Force 6','Force 7','Force 8'};
% 
% for ii = 1:size(fmat,1)
%     nexttile;
%     semilogy(fs,abs(fmat(ii,:)),'k','Linewidth',2)
%     title(chans{ii},'interpreter','tex')
%     grid on;
%     xlim([10 3000])
% %     xlabel('Frequency (Hz)','interpreter','tex')
% %     ylabel('Force PSD (N^2/Hz)','interpreter','tex')
% end
% t.Padding = 'compact';
% t.TileSpacing = 'tight';
% 
% xlabel(t,'Frequency (Hz)','interpreter','tex')
% ylabel(t, 'Force PSD (N^2/Hz)','interpreter','tex')
% title(t, 'Nominal Flight Force PSDs','interpreter','tex')

% Create flight environment PSD matrix
Sxx = zeros(size(H,1),size(H,1),length(fs)); % initialize
df = 5; % frequency spacing (Hz)

for ii = 1:length(fs)
    Sxx(:,:,ii) = H(:,:,ii)*(fmat(:,ii)*fmat(:,ii)')*H(:,:,ii)'/df; % x = H*f in PSD form. Divide by df to get /Hz units.
end

% save('RocketEnv.mat','Sxx','fs')
