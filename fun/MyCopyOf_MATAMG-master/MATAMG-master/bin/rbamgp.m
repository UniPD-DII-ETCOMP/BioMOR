function P = rbamgp(ia,ja,aa,C,x)
% RBAMGP rBAMG interpolation.
%   P = RBAMGP(IA,JA,AA,C,X)
%
%   See also amgp, aamgp, ibamgp, pdata.

%   Written by Minho Park
%   Department of Applied Mathematics
%   University of Colorado at Boulder
%   parkmh@colorado.edu
%   Date : 2008-11-24
if nargin < 5
    error('MATLAB:rbamgp:NotEnoughInputs','Not enough input arguments.');
end

m = length(ia)-1;           % Size of A

nc = sum(C);                % Number of C point

carray = zeros(1,m);    % carray maps index in fine gris to that in coarse grid

[nx lengthx] = size(x);

% Calculate residuals
r = zeros(nx,lengthx);

for i = 1 : nx
    r(i,:) = mm(ia,ja,aa,x(i,:));
end


% maxindex is used to allocate memory for P .
% Note that this is just an approximation.
maxindex = intpsetup(C,ia,carray);

ip = zeros(maxindex,1);
jp = zeros(maxindex,1);
ap = zeros(maxindex,1);

% Call mex function
% pindex is nummber of nonzeros in P
pindex = rbamgp_mex(C,ia,ja,aa,ip,jp,ap,carray,x,r);

% Build the interpolation matrix P
P = sparse(ip(1:pindex-1)',jp(1:pindex-1)',ap(1:pindex-1).',m,nc);


