function a_i_j = getA(ia,ja,aa,row,col)
% GETA Return the matrix value, A(row,col). 
%   A_I_J = GETA(IA,JA,AA,ROW,COL)
%
% See also sparse2crs, crs2sparse, mm

%   Written by Minho Park
%   Department of Applied Mathematics
%   University of Colorado at Boulder
%   parkmh@colorado.edu
%   $Date: 2009/08/19 $

a_i_j = 0;

for j = ia(row):ia(row+1)-1    
    if ja(j) == col
        a_i_j = aa(j);
        break
    end
end