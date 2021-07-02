function [ia ja aa] = sparse2crs(A)  
% SPARSE2CRS Convert sparse matrix to compressed row stroage.
%   [IA JA AA] = SPARSE2CRS(A) convert matrix A into row index IA,
%   column index JA, and values AA.
%   [ia ja aa] = sparse2crs(A)  
%
%   REMARKS:
%   1. Every diagonal element must be nonzero.
% 
%   See also mm, crs2sparse.

%   Written by Minho Park
%   Department of Applied Mathematics
%   University of Colorado at Boulder
%   parkmh@colorado.edu
%   $Date: 2009/08/18 $
%   $Update : 2010/01/22 Make algorithm fast
m = length(A);

[ja i aa] = find(A.');       % get sparsity information of A;

ia = zeros(m+1,1);          % ia is index where first element in 
                            % each row starts

diff = zeros(m,1);                            

% nnzrow(i,diff);
% 
% ia(1) = 1;
% for k = 1:m
%     ia(k+1) = ia(k) + diff(k);
% end
      
crsrow(i,diff,ia);


