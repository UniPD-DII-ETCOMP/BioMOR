% This M-file shows how to use the MATAMG toolbox.
clear
close all
clc
dad=pwd;
cd ..;
cd ..;
addpath(genpath(pwd));
cd(dad)
% Build a matrix A
stencil = [4 -1 -1 -1 -1];
n = 2000;
A = stencil2mat(stencil,n);

figure
spy(A)

b = A * ones(length(A),1);

% Generate MATAMG option file
opt = amgset;                               % default option file
opt = amgset(opt,'coarsest',20);            % set the number of levels        
opt = amgset(opt,'PreCond','pcg');          % set the Krylov method
opt = amgset(opt,'PrintOnScreen','on');    % turn off the option to print the log on the screen 
opt = amgset(opt,'SaveCsn','on');           % save the set of coarse-grid points
opt = amgset(opt,'CsnType','amg');          % choose the coarsening method
opt = amgset(opt,'LogFile','on');  
opt = amgset(opt,'Log','on');  
opt = amgset(opt,'LogNum','on');  
opt = amgset(opt,'TolAMG',1e-8);  
% Initial vector
v = rand(length(A),1);

% Solve a linear system
global load_amgdata
load_amgdata ='on';

[amgdata] = amginitsetup(A,opt);

tic
[v1, amgout, amgtime] = amg(A,v,b,opt);
toc

err=norm(v1-ones(length(A),1))/norm(ones(length(A),1))


% tic
% [v2, amgout, amgtime] = myRTamg(A,v,b,opt);
% toc

