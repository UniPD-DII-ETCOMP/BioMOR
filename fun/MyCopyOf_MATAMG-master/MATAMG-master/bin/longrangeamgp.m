function P = longrangeamgp(ia,ja,aa,inbd,jnbd,C2)
% AMGP Classical AMG interpolation.
%   P = LONGRANGEAMGP(IA,JA,AA,S,C)
%
%   See also length2strong, aggcoarsening, amgp, aamgp, bamgp, ibamgp

%   Written by Minho Park
%   min.park@nottingham.ac.uk
%   $Date: 2011/07/07 $

% Changelog
% v.0.1.2 - Added long range AMG interpolation (longrangeamgp.m)

if nargin < 6
    error('MATLAB:longrangeamgp:NotEnoughInputs','Not enough input arguments.');
end

m = length(ia)-1;           % Size of A

nc = sum(C2);                % Number of C point

carray = zeros(1,m);    % carray maps index on the fine gris into that ont the coarse grid

% maxindex is used to allocate memory for P .
% Note that this is just an approximation.
maxindex = intpsetup(C2,inbd,carray);

% Memory allocation
ip = zeros(maxindex,1);
jp = zeros(maxindex,1);
ap = zeros(maxindex,1);

% Call mex function
% pindex is nummber of nonzeros in P
pindex = longrangeamgp_mex(C2,ia,ja,aa,ip,jp,ap,inbd,jnbd,carray);

% [ip jp ap]
% Build the interpolation matrix P
P = sparse(ip(1:pindex-1)',jp(1:pindex-1)',ap(1:pindex-1).',m,nc);