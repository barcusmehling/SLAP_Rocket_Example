function plotTriaxPSD(fs,Sxx,Sxx_est,ref_accs)
% Create a plot comparing the flight PSD with that acheived in the Lab at a
% single triax.
%
% plotTriaxPSD(fs,Sxx,Sxx_est,ref_accs)
%   fs = frequency vector 1 x Nf
%   Sxx = No x No x Nf SDM from flight
%   Sxx_est = No x No x Nf SDM from lab test

% Extract the diagonal terms of the spectral density matrix for only the
% reference channels (defined above).
envpsd = get_psd(Sxx(ref_accs,ref_accs,:)); % flight environment PSDs at DUT ref accels
labpsd = get_psd(Sxx_est(ref_accs,ref_accs,:)); % '' lab env

% Plot reference responses in flight vs. lab
figure('Units','normalized','Position',[0.1 0.1 0.8 0.4]);
tlt = tiledlayout(1,3);
titles = {'Ref X','Ref Y','Ref Z'};

for ii = 1:3
    nexttile;
    semilogy(fs,abs(envpsd(ii,:)),'k',fs,abs(labpsd(ii,:)),'b','Linewidth',2)
    xlabel('Frequency (Hz)','interpreter','tex')
    grid on;
    xlim([10 2000])
    if ii == 1
        legend('Flight','Test')
    end
    title(titles{ii})
end

ylabel(tlt,'Acceleration PSD (g^2/Hz)','interpreter','tex')
tlt.TileSpacing = 'tight';
tlt.Padding = 'tight';