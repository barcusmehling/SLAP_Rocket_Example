function [H_fit] = cl_model(wns, zts, Acl,ws,FRF_type,varargin);
% [H_fit] = cl_model(wn,zt,Acl,ws,FRF_type,sum_flag,LRUR)
%
% Function which creates an FRF model of the modes described by 'lam' and 'A_fit'.
% Optimized for fast processing
%
% 'wn' is a vector of natural frequencies
% 'zt' is a vector of damping ratios
% 'Acl' is an array of residues of size [N,No,Ni] where N is the number
%       of modes, No the number of outputs and Ni the number of inputs.
%       For classical, mass normalized modes these are:
%           Acl(r,:,:) = phi_r_R*phi_r_DP
%       where 'R' denotes the response points and 'DP' the drive points.
% OR 'Acl' = {phi_r,phi_d} where these are the mode shapes at the drive and
%       response points of interest respectively.
%       phi_r => [No,N], phi_d => [Ni,N]
% The rest of the arguments are passed to ss_model:
% 'ws' is the frequency vector for the FRF
% 'FRF_type' is a flag indicating the type of FRF
%       1 - Displacement/Force
%       2 - Acceleration/Force
%       3 - Velocity/Force
% 'LRUR' contains the lower and upper residuals arranged like A_fit above.
%       H(w) (displacement) = Modes + LR/(i*w)^2 + UR;
%   Optional paramter, sum_flag tells whether to return the sum of modal FRFs, 
%   or individual FRFs.  sum_flag = 
%       's' - Returns the sum of modal FRFs, (H_fit is up to three dimensional)
%       'i' - Returns individual modal FRFs, (H_fit is four dimensional)
%           so that:  H_fit(:,:,:,r) is the rth modal FRF.
%
% If lam_fit and A_fit are empty [] then the algorithm reconstructs the
% residual terms LR,UR only.
%
% By Matt Allen - Aug. 2006
%

if isempty(zts); zts = zeros(sizes(wns)); end

if length(wns)~= length(zts)
    error('Length of wn and zt vectors are not equal');
end

if iscell(Acl)
    phir=Acl{1}; phid=Acl{2};
    if size(phir,2)~=size(phid,2);
        disp('Mode vectors supplied rather than residue matrix');
        error('Mode shape matrices for response and drive points must have the same number of modes');
    end
    % Form residue matrices for each mode
    % Note, this overwrites the Acl that was passed in
    Acl=zeros(size(phir,1),size(phid,1),length(wns));
        % The residue matrix for the rth mode is Acl(:,:,r)
    for k=1:length(wns);
        Acl(:,:,k)=phir(:,k)*phid(:,k).';
    end
    Acl=permute(Acl,[3,1,2]); % Permute to size requested by AMI:
end

% if any(abs(wns) < 100*eps);
%     disp('Rigid Body Modes Detected');
%     rgbm_flg = (abs(wns) < 100*eps);
%     Pex = sum(Acl(rgbm_flg,:,:),1);
%     wns = wns(~rgbm_flg); zts = zts(~rgbm_flg); Acl = Acl(~rgbm_flg,:,:);
%     if length(varargin) > 3;
%         if size(varargin{4}(1,:,:)) ~= size(Pex);
%             error('Pex supplied by user not the same size as modal Residues');
%         end
%         varargin{4}(1,:,:) = varargin{4}(1,:,:)+Pex;
%     else
%         varargin{4} = Pex;
%     end
% end
% 
% lam = -zts.*wns+i*wns.*sqrt(1-zts.^2);
% 
% % Find state space residues
% Ni = size(Acl,3);
% Acl = mimo2simo_rs(Acl);
% if any(imag(lam) == 0); warning('Imaginary Natural Frequencies - Unstable Modes'); end
% A = i*diag((2*imag(lam)).^-1)*Acl;
% Acl = simo2mimo_rs(Acl,Ni);
% 
% % Pass to ss_model
% Hfit = ss_model(lam,A,varargin{:});

%% From ss_model
[a,b] = size(ws);
    if b>a;
        ws = ws.';
        [a,b] = size(ws);
    end % Assure that ws is a column
    if b>1; warning('ws is not a column vector'); end

[Na,No,Ni] = size(Acl);
N = length(wns);
	if Na ~= N
        warning('Sizes of Acl and wns are not the same, minimum used');
        N = min([N,Na]);
	end

if ~isempty(varargin)
    sum_flag = lower(varargin{1});
else
    sum_flag = 's'; % Default is to return the sum of modal FRFs
end

if nargin > 6;
    Pex = varargin{2};
else
    Pex = [];
end

% Allow for FRF_type to be a character or number
if ischar(FRF_type);
    if strcmp(upper(FRF_type),'D'); FRF_type = 1;
    elseif strcmp(upper(FRF_type),'A'); FRF_type = 2;
    elseif strcmp(upper(FRF_type),'V'); FRF_type = 3;
    else
        error('Unrecognized FRF_type String');
    end
end

if ~isempty(wns) && ~isempty(Acl); % find out of band terms only.
    Acl = mimo2simo_rs(Acl); % Reshape for easy processing
    if sum_flag == 's';
        if N > 1;
            H_fit = zeros(length(ws),No*Ni);
            for k = 1:N
                if FRF_type == 1; % Displacement
                    H_fit = H_fit + (wns(k).^2-ws.^2+2*1i*zts(k)*wns(k)*ws).^-1*Acl(k,:);
                elseif FRF_type == 2; % Acceleration Transfer Function
                    H_fit = H_fit + ((-ws.^2).*(wns(k).^2-ws.^2+2*1i*zts(k)*wns(k)*ws).^-1)*Acl(k,:);
                elseif FRF_type == 3; % Velocity Transfer Function
                    H_fit = H_fit + ((1i*ws).*(wns(k).^2-ws.^2+2*1i*zts(k)*wns(k)*ws).^-1)*Acl(k,:);
                else
                    error('Unrecognized FRF type'); FRF_type
                end
            end
        else % one mode only
            % This extra code probably contributes little to speeding
            % things up.
            if FRF_type == 1; % Displacement
                H_fit = (wns.^2-ws.^2+2*1i*zts*wns*ws).^-1*Acl;
            elseif FRF_type == 2; % Acceleration Transfer Function
                H_fit = ((-ws.^2).*(wns.^2-ws.^2+2*1i*zts*wns*ws).^-1)*Acl;
            elseif FRF_type == 3; % Velocity Transfer Function
                H_fit = ((1i*ws).*(wns.^2-ws.^2+2*1i*zts*wns*ws).^-1)*Acl;
            else
                error('Unrecognized FRF type'); FRF_type
            end
        end
    else
        H_fit = zeros(length(ws),No*Ni,N);
        for k = 1:N
            if FRF_type == 1; % Displacement
                H_fit(:,:,k) = (wns(k).^2-ws.^2+2*1i*zts(k)*wns(k)*ws).^-1*Acl(k,:);
            elseif FRF_type == 2; % Acceleration Transfer Function
                H_fit(:,:,k) = ((-ws.^2).*(wns(k).^2-ws.^2+2*1i*zts(k)*wns(k)*ws).^-1)*Acl(k,:);
            elseif FRF_type == 3; % Velocity Transfer Function
                H_fit(:,:,k) = ((1i*ws).*(wns(k).^2-ws.^2+2*1i*zts(k)*wns(k)*ws).^-1)*Acl(k,:);
            else
                error('Unrecognized FRF type'); FRF_type
            end
        end
    end
    
else
    H_fit = zeros(length(ws),No*Ni);
end

% Find residual terms
if ~isempty(Pex);
    Pex = mimo2simo_rs(Pex);
    if size(Pex,1) ~= 2
        Pex(2,:,:) = 0; % No UR included
    end
    P_fit = zeros(length(ws),No*Ni);
    if FRF_type == 1; % Displacement
        P_fit = P_fit + (1i*ws).^-2*Pex(1,:) + Pex(2*ones(size(ws)),:);
    elseif FRF_type == 2; % Acceleration Transfer Function
        P_fit = P_fit + Pex(ones(size(ws)),:) + (1i*ws).^2*Pex(2,:);
    elseif FRF_type == 3; % Velocity Transfer Function
        P_fit = P_fit + (1i*ws).^-1*Pex(1,:) + (1i*ws)*Pex(2,:);
    else
        error('Unrecognized FRF type'); FRF_type
    end
    
    H_fit = H_fit + P_fit;
    
end

% Reshape
H_fit = simo2mimo_rs(H_fit,Ni);
