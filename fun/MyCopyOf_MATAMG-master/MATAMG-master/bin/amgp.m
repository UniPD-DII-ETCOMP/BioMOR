function P = amgp(ia,ja,aa,is,js,C)
% AMGP Classical AMG interpolation.
%   P = AMGP(IA,JA,AA,S,C)
%
%   See also strong_connection, amg_coarsening, aamgp, bamgp, ibamgp

%   Written by Minho Park
%   Department of Applied Mathematics
%   University of Colorado at Boulder
%   parkmh@colorado.edu
%   $Date: 2008/11/27 $
%   Update : 2010/8/17 changed to access a_k_j once.

if nargin < 6
    error('MATLAB:amgp:NotEnoughInputs','Not enough input arguments.');
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
pindex = amgp_mex(C,ia,ja,aa,ip,jp,ap,is,js,carray);


% Build the interpolation matrix P
P = sparse(ip(1:pindex-1)',jp(1:pindex-1)',ap(1:pindex-1).',m,nc);

