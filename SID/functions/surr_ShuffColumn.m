function [B] = surr_ShuffColumn(B,ii)
%Performs joint random permutation of [ii] columns of 2D matrix

if ~exist('ii','var'), ii=1; end

p=randperm(size(B,1));
B(p,ii) = B(:,ii);

end

