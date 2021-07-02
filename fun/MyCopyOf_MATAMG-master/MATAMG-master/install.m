% For more information, see <a href="matlab: 
% web('http://www.grandmaster.colorado.edu/~parkmh')">Minho Park's Web site</a>.
cd bin
Ver = amgversion;
user = 'Minho Park';
email = 'min.park@nottingham.ac.uk';
cd ..

clc
fprintf(' ********************************************\n')
fprintf('\n')
fprintf(' Matlab Algebraic Multigrid Toolbox  %s \n',Ver) 
fprintf('\n')
fprintf('%20s %s\n','Written by', user);
fprintf('%13s %s\n','email :',email);
fprintf(' ********************************************\n')

cwd = pwd;
matamgroot = cwd;

% Add path
% Generate amgpath.m file
fid = fopen([pwd filesep 'bin' filesep 'amgpath.m'],'w');
fprintf(fid,'function matamgpath = amgpath\n');

addpath(matamgroot);
addpath(fullfile(matamgroot,'bin'));

fprintf(fid,'matamgpath = ''%s'';\n',matamgroot);

fclose(fid);
savepath;


% Compile mex files
cd(fullfile(cwd,'src'));
copyfile('*.c',fullfile(cwd,'bin'))

cd(fullfile(cwd,'bin'));

fprintf('\n2. Compile mex -O Files\n')
cd([cwd filesep 'bin'])

try
    fprintf('Compiling crsrow.c\n')
    mex -O crsrow.c
    fprintf('Compiling stencil2mat_mex.c\n')
    mex -O stencil2mat_mex.c
    fprintf('Compiling gs.c\n');
    mex -O gs.c
    fprintf('Compiling jacobi.c\n');
    mex -O jacobi.c
    fprintf('Compiling bgs.c\n');
    mex -O bgs.c
    fprintf('Compiling cr.c\n');
    mex -O cr.c
    fprintf('Compiling sc.c\n')
    mex -O sc.c
    fprintf('Compiling sc_positive.c\n')
    mex -O sc_positive.c
    fprintf('Compiling firstcoloring.c\n')
    mex -O firstcoloring.c
    fprintf('Compiling secondcoloring.c\n')
    mex -O secondcoloring.c
    fprintf('Compiling intpsetup.c\n')
    mex -O intpsetup.c
    fprintf('Compiling amgp_mex.c\n')
    mex -O amgp_mex.c
    fprintf('Compiling aamgp_mex.c\n')
    mex -O aamgp_mex.c
    fprintf('Compiling rbamgp_mex.c\n')
    mex -O rbamgp_mex.c
    fprintf('Compiling length2strong_mex.c\n')
    mex -O length2strong_mex.c
    fprintf('Compiling aggfirstcoloring.c\n')
    mex -O aggfirstcoloring.c
    fprintf('Compiling aggsecondcoloring.c\n')
    mex -O aggsecondcoloring.c
    fprintf('Compiling longrangeamgp_mex.c\n')
    mex -O longrangeamgp_mex.c
catch err
    cd ..
    error('mex -O compile error')
end

% fprintf('\n3. Delete Source Files\n')
% delete *.c
% delete *.h


cd ..
clear all

