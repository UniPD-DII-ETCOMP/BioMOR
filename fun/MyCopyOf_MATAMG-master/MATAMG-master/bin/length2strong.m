function [S2 inbd jnbd] = length2strong(ia,ja,C,path)
% LENGTH2STRONG determines the point that are strongly n-connected C points
% w.r.t the number paths, p and path length, l. In this code we only
% consider length 2 and it is enough for most cases. 
%   S2 = LENGTH2STRONG(IA,JA,S,C,path)
%
%   See also strongconnection, aggresivecoarsening.

%   Wirtten by Minho Park
%   $Date: 2011/5/17 $

% Changelog
% v.0.1.2 - Added length 2 strong connection (length2strong.m)
% v.0.1.2 - Fixed length2strong_mex to return neighborhood points (inbd and jnbd)(length2strong.m)



if nargin ~= 4
    error('MATLAB:length2strong:notEnoughInputs','Not enough input arguments.');
end

[is2 js2 as2 index inbd jnbd] = length2strong_mex(ia,ja,C,path);

S2 = sparse(is2(1:index),js2(1:index),as2(1:index),length(C),length(C));