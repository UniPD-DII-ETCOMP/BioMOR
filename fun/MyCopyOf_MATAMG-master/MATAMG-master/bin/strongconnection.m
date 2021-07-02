 function [nofc is js it jt] = strongconnection(ia,ja,aa,theta,offdiag)
% STRONG_CONNECTION (MEX)Determine the strong dependency and influence of
% matrix A(ia,ja,aa) with respect to the thresold theta.
%   S = STRONG_CONNECTION(IA,JA,AA,THETA)
%   will generate matirx S and its transpost T.
%   IA, JA, and AA is crs of A and THETA is thereshold.
%
%   See also amgcoarsening.

%   Written by Minho Park
%   parkmh78@gmail.com
%   $Date: 2008/11/28 $
%   Update : 2010/02/06 can use both MEX-File and M-File.
%          : 2010/08/16 added dynamic threshold for strong connection.
% Changelog :
% v.0.1.2 - Removed unnecessary computation caused by using MATLAB index in mex file.

% Set default
if nargin < 4
    error('MATLAB:strongconnection:NotEnoughInputs','Not enough input arguments.');
end

% % Consider only negative off-diagonal values of matrix A
% strongness = 'negative';

% Calculate strong connections of A
if strcmp(offdiag,'negative')
    [ns nst is js it jt] = sc(ia,ja,aa,theta);
    % Build sparse matrix from the matrix in CRS form
elseif strcmp(offdiag,'absolute')
    [ns nst is js it jt] = sc_positive(ia,ja,aa,theta);
    % Build sparse matrix from the matrix in CRS form
end

nofc = ns + nst;

