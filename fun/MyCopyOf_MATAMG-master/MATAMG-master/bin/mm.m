function b = mm(ia,ja,aa,x)
% MM Matrix-vector multiplication.
%   B = MM(IA, JA, AA, X) calculate multiplication between A and X
%   , and IA, JA and AA are crs form of A.
%
% See also sparse2crs

%   Written by Minho Park
%   Department of Applied Mathematics
%   University of Colorado at Boulder
%   parkmh@colorado.edu
%   $Date: 2009/08/18 $

% memory allocation for the result
b = zeros(length(x),1);

% Matrix vector multiplication
for i = 1:length(ia)-1
    for j = ia(i):ia(i+1)-1
        b(i) = b(i) + aa(j)*x(ja(j));        
    end
end
    