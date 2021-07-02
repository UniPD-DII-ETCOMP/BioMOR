function P = aamgp(ia,ja,aa,C,x)
% AAMGP Adaptive AMG interpolation.
% AAMGP use one smooth vector x.
%   P = AAMGP(IA,JA,AA,C,X)
%
%   See also strong_connection, amg_coarsening, amgp, bamgp, ibamgp, rbamgp, pdata.

%   Written by Minho Park
%   Department of Applied Mathematics
%   University of Colorado at Boulder
%   parkmh@colorado.edu
%   $Date: 2008/11/27 $


if nargin < 5
    error('MATLAB:aamgp:NotEnoughInputs','Not enough input arguments.');
end



m = length(ia)-1;           % Size of A

nc = sum(C);                % Number of C point

carray = zeros(1,m);    % carray maps index in fine gris to that in coarse grid

% maxindex is used to allocate memory for P .
% Note that this is just an approximation.
maxindex = intpsetup(C,ia,carray);

% Memory allocation
ip = zeros(maxindex,1);
jp = zeros(maxindex,1);
ap = zeros(maxindex,1);


% Call mex function
% pindex is nummber of nonzeros in P
pindex = aamgp_mex(C,ia,ja,aa,ip,jp,ap,carray,x);

% Build the interpolation matrix P
P = sparse(ip(1:pindex-1)',jp(1:pindex-1)',ap(1:pindex-1).',m,nc);

