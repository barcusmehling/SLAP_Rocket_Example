function [Xout] = simo2mimo_rs(Xin,Ni);
% [Xout] = simo2mimo_rs(Xin,Ni);
%
% Reshapes the array Xin from SIMO to MIMO convention
% For example, if Xin is 100 x 27, and Ni = 3 Xout is 100 x 9 x 3
%
% This also reshapes x_model, so if size(Xin) = [Nf,NiNo,Nmodes],
% the function loops over the modes so that size(X_out) = [Nf,Ni,No,Nmodes]
%
[a b c] = size(Xin);
    b = b/Ni;
    if round(b) ~= b; error('Matrix not compatible with Number of Inputs'); end

Xout = zeros(a,b,Ni,c);
for m = 1:c
    for ii = 1:1:Ni
        Xout(:,:,ii,m) = Xin(:,[((ii-1)*b+1):1:ii*b],m);
	end
end
