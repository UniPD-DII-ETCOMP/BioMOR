function C2 = aggcoarsening(S2,C,ia,ja,inbd,jnbd)
% AGGCOARSENING generate the C/F splitting by more aggresive way. 
%
% 
%   See also length2strong

%   Written by Minho Park
%   email : parkmh78@gmail.com
%   $Date : 2011/05/18

% Changelog
% v.0.1.2 - Added aggresive first coloring scheme (aggcoarsening.m)
% v.0.1.2 - Added aggresive second coloring scheme (aggcoarsening.m)

m = length(S2);
% save S2 S2
% Memory allocation
C2 = zeros(1,m);

[it jt] = sparse2crs(S2');
[is js] =  sparse2crs(S2);

% First coloring
aggfirstcoloring(C,C2,it,jt,is,js);
% save C2 C2
% Second coloring
% aggsecondcoloring(C2,ia,ja,inbd,jnbd,is,js);
% save C3 C2