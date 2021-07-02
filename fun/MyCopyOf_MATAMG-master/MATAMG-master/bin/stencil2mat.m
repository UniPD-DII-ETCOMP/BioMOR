function A = stencil2mat(stencil,n)
% STENCIL2MAT Generate the matrix from stencil. 
%   A = STENCIL2MAT(STENCIL, N)
%   
%   length(stencil)  : 3(1D) [1]C [2]W [3]E
%                    : 5(2D) [1]C [2]W [3]E [4]S [5]N
%                    : 9(2D) [1]C [2]W [3]E [4]SW [5]S [6]SE [7]NW [8]N [9]NE
%                    : 7(3D) [1]C [2]W [3]E [4]S [5]N [6]D [7]U
%   See also --------

%   Minho Park
%   Department of Applied Mathematics
%   University of Colorado at Boulder
%   parkmh@colorado.edu
%   Date : 2010.02.06
%   Update : 
% generate matrix from stencil 


[ia ja aa] = stencil2mat_mex(stencil,n);
A = sparse(ia,ja,aa);