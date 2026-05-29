function [sigrms,sigpsd,sigloc,varargout] = GetStressFunc(Hs,Sff,df,rms_inds,sh_inds)
% Function that calculates the stress metrics given an input force spectrum
% Sff and matrix of stress transfer functions Hs.
%
% [sigrms,sigpsd,sigloc] = GetStressFunc(Hs,Sff,df,rms_inds,sh_inds)
%
% INPUTS:
% Hs = nstresspoints x 1 cell, where each element contains a [6 x nforces]
%   frequency response matrix relating force spectrum to the spectrum of
%   the [6x1] stress tensor: [s_xx, s_yy, s_zz, s_xy, s_xz, s_yz] where s
%   denotes the stress. 
%
% Sff = nforces x nforces x nflines. Spectral Density Matrix of the forces
%   applied to the system
%
% df = frequency spacing (Hz)
%
% rms_inds = frequency lines to include in RMS calculation. If you use
%   the 0 Hz frequency line, RMS stress usually gets really big, so this
%   is mainly here to exclude that line.
% sh_inds = (nforces x 1) vector of shakers to use.  If Hs has more shakers
%   than there are applied force, this vector tells the function which
%   shaker indicies to apply the forces to.
%
% OUTPUTS:
% sigrms = RMS Von Mises stress at the strain gauge with the largest RMS value
% sigpsd = PSD of Von Mises stress for the strain gauge with the largest
%   RMS value
% sigloc = index of strain gauge with the largest RMS Von Mises stress.
%
% To obtain the Von Mises stress PSD at all gauges:
%
% [sigrms,sigpsd,sigloc,sig_psds_all] = GetStressFunc(Hs,Sff,df,rms_inds,sh_inds)
%
% Marcus Behling, 2026, BYU Mechanical Engineering
% Edits by M.S. Allen & B.M. Bahr, May 2026
%

    A1 = [1 -0.5 -0.5;-0.5 1 -0.5;-0.5 -0.5 1]; % Assemble VM Stress Quadratic Form Matrix
    A2 = 3*eye(3);
    A = blkdiag(A1,A2);

    rmss = zeros(length(Hs),1);

    H_ex = Hs{1};
    nf = size(H_ex,3); % frequency lines

    sig_psds = zeros(length(Hs),nf);
    
    for ii = 1:length(rmss)
        Hii = Hs{ii};
        
        for jj = 1:nf
            Ssigsig = squeeze(Hii(:,sh_inds,jj))*Sff(:,:,jj)*squeeze(Hii(:,sh_inds,jj))';
            sig_psds(ii,jj) = sig_psds(ii,jj) + abs(sum(sum(Ssigsig.*A)));
                % This implements Eq. (4) in Behling, et al, AIAA-JSCR, 2026
        end
        rmss(ii) = sqrt(df*sum(sig_psds(ii,rms_inds)));
    end
    % Find maximum RMS stress and strain gauge where it occurs:
    [sigrms,ind1] = max(rmss);

    sigpsd = sig_psds(ind1,:);
    sigloc = ind1;

    % This only exports the PSD at the strain gauge with the largest RMS
    % strain. To plot a few of the maximum ones, use:
    %{
    [mvs,sind]=sort(rmss);
    figure; semilogy([0:1:size(sig_psds,2)-1]*df,sig_psds(sind(1:10),:));
    %}
    if nargout>3
        varargout{1}=sig_psds;
    end

end