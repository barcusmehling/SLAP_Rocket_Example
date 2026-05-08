function [Xout] = mimo2simo_rs(Xin);
% Subfunction of AMI:  Matt Allen 11/20/2003
% MIMO AMI data is stored in an array H(a,b,c) where
%   a = frequency point
%   b = Output index (size(H,2) = No)
%   c = Input inces (size(H,3) = Ni)
%
% This function stacks the FRF matrices for each drive k, H(:,:,k)
% into a single matrix such that size(H) = [w, Ni*No].

[a b c] = size(Xin);

Xout = zeros(a,b*c);
for ii = 1:1:c
    Xout(:,[((ii-1)*b+1):1:ii*b]) = Xin(:,:,ii);
end