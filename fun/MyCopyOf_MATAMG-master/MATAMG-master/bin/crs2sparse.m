function A = crs2sparse(ia,ja,aa,m,n)  
% CRS2SPARSE Convert compressed row stroage into sparse matrix.
%   A = CRS2SPARSE(IA,JA,AA) convert row index IA, column index JA,
%   and values AA INTO matrix a.
%   [ia ja aa] = sparse2crs(A)  
%
%   REMARKS:
%   1. Every diagonal element must be nonzero.
% 
%   See also mm, sparse2crs.

%   Written by Minho Park
%   Department of Applied Mathematics
%   University of Colorado at Boulder
%   parkmh@colorado.edu
%   $Date: 2009/08/18 $
if nargin < 4
m = length(ia)-1;
end
if nargin < 5
n = length(ia)-1;
end

iia = zeros(size(ja));
for i = 1:m
    for j = ia(i) : ia(i+1)-1
        iia(j) = i;
    end
end

A = sparse(iia,ja,aa);
        
