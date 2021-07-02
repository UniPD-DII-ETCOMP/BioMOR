function C = amgcoarsening(ia,ja,aa,is,js,it,jt,nofc)
% AMG_COARSENING Determine coarse grid points from matrix A(ia,ja,aa)
%   C = AMG_COARSENING(IA,JA,AA,S)
%   will generate vector C which has 1 at coarse grid point
%   and 0 at fine grid point. S is strong connection matrix.
%
%   See also strong_connection, cr_coarsening, cdata.

%   Written by Minho Park
%   Department of Applied Mathematics
%   University of Colorado at Boulder
%   parkmh@colorado.edu
%   $Date: 2008/11/22 $
%   Update : 2010/01/22 Double linked list
%            2010/01/25 MEX for First coloring and second coloring
if nargin < 8
    error('MATLAB:amgcoarsening:NotEnoughInputs','Not enough input arguments.');
end


m = length(ia)-1;       % Size of matrix A

% Memory allocation
C = zeros(1,m);

% First Coloring Pass
firstcoloring(C,it,jt,is,js,nofc);

% Second Coloring Pass
secondcoloring(C,ia,ja,aa,is,js);

