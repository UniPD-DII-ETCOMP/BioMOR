function [v, amgout, amgtime] = myRTamg(A,v,b,opt,returnMtx)
%AMG Run amg algorithm based on AMG OPTIONS parameters.
%   V = AMG(A,V,B,OPTIONS) use algebraic multigrid algorithm to solve
%   the linear system A*x = B based on algebraic multigrid algorithm.
%
%  See amgset.

% todo - adaptive process
% todo - update contents.m
% todo - update info.xml
% todo - GUI
% todo - parallel coarsening
% todo - plot coarsenging
% todo - compute overall efficiency
% todo - html
% v.0.1.2 - Fixed bug in rbamgp.c
% v.0.1.2 - Changed OutFile to LogFile, and added Log and LogNum;
amgout.residual=[];
amgout.errors=[];
if nargin < 5
    returnMtx = 'off';
end

if nargin < 4
    error('MATLAB:amg:NotEnoughInputs','Not enough input arguments.');
end

if ~isa(opt,'struct')
    error('MATLAB:amg:Arg4NotStruct',...
        'Fouth argument must be an options structure by AMGSET.');
end

[m n] =  size(A);
if m~= n
    error('MATLAB:amg:NotSquared','Matrix must be squared.')
end

if (n~= length(v)) ||  (n~= length(b))
    error('MATLAB:amg:DiffDim','Inner matrix dimensions must agree')
end

logflag = amgget(opt,'Log','off');

if strcmp(logflag,'on')
    % Get current matamg version
    ver = amgversion;
    
    lognum = amgget(opt,'LogNum','off');
    
    chapt = 1; % Chater number
    % Load index number from amg.idx
    
    if strcmp(lognum,'on')
        fid1 = fopen('amg.idx','r');
        if fid1 == -1
            fid1 = fopen('amg.idx','w');
            fprintf(fid1,'%d',1);
            
        else
            idx = fscanf(fid1,'%d',1) + 1;  % increase index by 1
            fclose(fid1);
            fid1 = fopen('amg.idx','w');
            fprintf(fid1,'%d',idx);
        end
        fclose(fid1);
    end
    
    % Open file for output
    logfile     = amgget(opt,'LogFile','out');
    if strcmp(lognum,'on')
        filename    = sprintf('%s%2.3d.rst',logfile,idx);
    else
        filename = sprintf('%s.rst',logfile);
    end
    fid2        = fopen(filename,'w');
    
    % Get current date and time
    c1 = clock;
end
maxCycle = amgget(opt,'MaxCycle',20);
tol     = amgget(opt,'TolAMG',1e-10);

residual    = zeros(maxCycle+1,1);       % stroage for L2 norm of residual || b - A*x ||
errors      = zeros(maxCycle+1,1);         % stroage for L2 norm of error || v ||. This is meaningful when b is a zero vector.
residual(1) = norm(b-A*v);
errors(1)   = norm(v);

nb = norm(b);
if nb == 0
    nb = 1;
end

if strcmp(logflag,'on')
    fprintf(fid2,'\t>>>> MATAMG %s <<<<\n',ver);
    fprintf(fid2,'%10s : %s\n','File',filename);
    fprintf(fid2,'%10s : %d/%d/%d\n%10s : %d:%d:%2.3f\n','Date', c1(2), c1(3), c1(1),'Time', c1(4), c1(5), c1(6));
end
% if intpType is 'off', then only relaxation algorithm of relType
intpType = amgget(opt,'IntpType','amg');
switch intpType
    case {'off'}
        [IA, JA, AA] = sparse2crs(A);
        rType   = amgget(opt,'RelType','gs');
        rPara   = amgget(opt,'RelPara',1);
        
        if strcmp(logflag,'on')
            % Display amg OPTION
            fprintf(fid2,'%d. Relaxation Options\n',chapt);
            chapt = chapt + 1;
            
            fprintf(fid2,'%10s : ''%s''\n','RelType',rType);
            fprintf(fid2,'%10s : %1.1f\n','RelPara',rPara);
            fprintf(fid2,'%10s : %d\n','MaxCycle',maxCycle);
            fprintf(fid2,'%10s : %2.2e\n','Tol', tol);
            
            % Convergence Factor
            fprintf(fid2,'\n%d. Convergence History\n', chapt);
            chapt = chapt + 1;
            fprintf(fid2,'=====================================================\n');
            fprintf(fid2,'%5s ||%10s%6s%6s |%9s%6s%6s\n','CYCLE','RESIDUAL', 'CF', 'AVG','ERROR','CF','AVG');
            fprintf(fid2,'=====================================================\n');
            fprintf(fid2,'%3s   ||%10.2e%6s%6s |%9.2e%6s%6s\n','0',residual(1), '', '',errors(1),'','');
        end
        if residual(1)/nb < tol     % initial is very close to solution
            return
        end
        
        cycleTic = tic;
        switch rType
            case {'gs'} % Gauss Seidel iterative method
                for i = 1:maxCycle
                    gs(IA,JA,AA,v,b,rPara,1)
                    residual(i+1)   = norm(b-A*v);
                    errors(i+1)     = norm(v);
                    if strcmp(logflag,'on')
                        fprintf(fid2,'%3d   ||%10.2e%6.2f%6.2f |%9.2e%6.2f%6.2f\n',i,residual(i+1), residual(i+1)/residual(i),...
                            (residual(i+1)/residual(1))^(1/i),errors(i+1),errors(i+1)/errors(i),...
                            (errors(i+1)/errors(1))^(1/i));
                    end
                    if (residual(i+1)/residual(1)) < tol %|| r(i+1)/r(i) > .95
                        break
                    end
                end
            case {'jacobi'}     % Jacobi iterative method
                for i = 1:maxCycle
                    jacobi(IA,JA,AA,v,b,rPara,1,2)
                    residual(i+1)   = norm(b-A*v);
                    errors(i+1)     = norm(v);
                    if strcmp(logflag,'on')
                        fprintf(fid2,'%3d   ||%10.2e%6.2f%6.2f |%9.2e%6.2f%6.2f\n',i,residual(i+1), residual(i+1)/residual(i),...
                            (residual(i+1)/residual(1))^(1/i),errors(i+1),errors(i+1)/errors(i),...
                            (errors(i+1)/errors(1))^(1/i));
                    end
                    if (residual(i+1)/residual(1)) < tol %|| r(i+1)/r(i) > .95
                        break
                    end
                end
            otherwise
                errid   = 'MATLAB:funfun:amg:relType:notARelType';
                errmsg  = sprintf('Invalid value for OPTIONS parameter %s: must be ''gs'' or ''jacobi''.',options.relType);
                error(errid, errmsg);
        end
        amgout.cycle = toc(cycleTic);
    case {'amg','aamg','rbamg','lramg'}
        global load_amgdata
        if strcmp(load_amgdata,'off')
             [amgdata] = amginitsetup(A,opt);
        elseif strcmp(load_amgdata,'on')
            if exist('CsnData.mat','file')
                load CsnData.mat amgdata
            else
                [amgdata] = amginitsetup(A,opt);
            end
        end
        
%         [amgdata amgtime] = amginitsetup(A,opt);
        
        % Update coarsest level if changed
        vdata.cst     = amgdata.cst;
        
        
        % Build vdata structure
        vdata.Iter1   = amgget(opt,'Iter1',1);
        vdata.Iter2   = amgget(opt,'Iter2',1);
        vdata.RelType = amgget(opt,'RelType','gs');
        vdata.RelPara = amgget(opt,'RelPara',1);
        vdata.Mu      = amgget(opt,'Mu',1);
        
        precnd = amgget(opt,'PreCond','off');
        
        if strcmp(logflag,'on')
            % Display amg OPTION
            fprintf(fid2,'_____________________________________________________\n\n');
            fprintf(fid2,'%d. AMG OPTIONS\n',chapt);chapt = chapt + 1;
            fprintf(fid2,'%10s : %d\n','Coarsest',vdata.cst);
            fprintf(fid2,'%10s : %d\n','MaxCycle',maxCycle);
            fprintf(fid2,'%10s : %2.2e\n','Tol', tol);
            fprintf(fid2,'%10s : %d\n','Iter1', vdata.Iter1);
            fprintf(fid2,'%10s : %d\n','Iter2', vdata.Iter2);
            fprintf(fid2,'%10s : %d\n','MuCycle', vdata.Mu);
            fprintf(fid2,'%10s : ''%s''\n','CsnType',amgget(opt,'CsnType','amg'));
            fprintf(fid2,'%10s : %2.2f\n','Theta',amgget(opt,'Theta',.25));
            fprintf(fid2,'%10s : %d\n','Nu',amgget(opt,'Nu',5));
            fprintf(fid2,'%10s : %2.2f\n','Alpha',amgget(opt,'Alpha',.7));
            
            fprintf(fid2,'%10s : ''%s''\n','IntpType',amgget(opt,'IntpType','amg'));
            fprintf(fid2,'%10s : %d\n','NumVec', amgget(opt,'NumVec',5));
            fprintf(fid2,'%10s : %d\n','NumRel', amgget(opt,'NumRel',5));
            fprintf(fid2,'%10s : ''%s''\n','RelType',vdata.RelType);
            fprintf(fid2,'%10s : %1.1f\n','RelPara',vdata.RelPara);
            fprintf(fid2,'%10s : %d\n','MaxCycle',maxCycle);
            fprintf(fid2,'%10s : ''%s''\n','PreCond',precnd);
            
            % Display Complexity
            fprintf(fid2,'\n_____________________________________________________\n\n');
            fprintf(fid2,'%d. Complexity\n',chapt);
            allnnzs = 0; alldofs = 0;
            fprintf(fid2,'\n%d.1. Operator Complexity\n',chapt);
            fprintf(fid2,'=================\n');
            fprintf(fid2,'LEVEL ||     NNZs\n');
            fprintf(fid2,'=================\n');
            for i = 1 : vdata.cst
                fprintf(fid2,'%5d || %8d\n',i,amgdata.nnz(i));
                allnnzs = allnnzs + amgdata.nnz(i);
            end
            fprintf(fid2,'*** Operator Complexity : %2.2f [(nnz(1)+...+nnz(%d))/nnz(1)]\n',allnnzs/amgdata.nnz(1),vdata.cst);
            
            fprintf(fid2,'\n%d.2. Grid Complexity\n',chapt);
            fprintf(fid2,'=================\n');
            fprintf(fid2,'LEVEL ||     DOFs\n');
            fprintf(fid2,'=================\n');
            for i = 1 : vdata.cst
                fprintf(fid2,'%5d || %8d\n',i,amgdata.nofdof(i));
                alldofs = alldofs + amgdata.nofdof(i);
            end
            fprintf(fid2,'*** Grid Complexity     : %2.2f [(dof(1)+...+dof(%d))/dof(1)]\n',alldofs/amgdata.nofdof(1),vdata.cst);
            
            fprintf(fid2,'\n%d.3. Density\n',chapt);chapt = chapt + 1;
            fprintf(fid2,'=================\n');
            fprintf(fid2,'LEVEL ||  Density\n');
            fprintf(fid2,'=================\n');
            for i = 1 : vdata.cst
                fprintf(fid2,'%5d || %6.1f %% \n',i,amgdata.nnz(i)/amgdata.nofdof(i)^2*100);
            end
            
            
            % Display Convergence Factor
            fprintf(fid2,'_____________________________________________________\n\n');
            fprintf(fid2,'%d. Convergence History\n',chapt); chapt = chapt+1;
            fprintf(fid2,'=====================================================\n');
            fprintf(fid2,'%5s ||%10s%6s%6s |%9s%6s%6s\n','CYCLE','RESIDUAL', 'CF', 'AVG','ERROR','CF','AVG');
            fprintf(fid2,'=====================================================\n');
            fprintf(fid2,'%3s   ||%10.2e%6s%6s |%9.2e%6s%6s\n','0',residual(1), '', '',errors(1),'','');
        end
        
        if residual(1)/nb < tol     % initial is very close to solution
            return
        end
        
        cycleTic = tic;
        switch precnd
            case {'off'}
                for i = 1 : maxCycle
                    v  = vcycle(A,v,b,1,vdata,amgdata);
                    
                    residual(i+1)   = norm(b-A*v);
                    errors(i+1)     = norm(v);
                    
                    if strcmp(logflag,'on')
                        fprintf(fid2,'%3d   ||%10.2e%6.2f%6.2f |%9.2e%6.2f%6.2f\n',i,residual(i+1), residual(i+1)/residual(i),...
                            (residual(i+1)/residual(1))^(1/i),errors(i+1),errors(i+1)/errors(i),...
                            (errors(i+1)/errors(1))^(1/i));
                    end
                    if (residual(i+1)/residual(1)) < tol %|| r(i+1)/r(i) > .95
                        %                 lastindex = i+1;
                        break
                    end
                    
                end
                amgout.amgcycle = i;
                amgout.error_cf = errors(i+1)/errors(i);
                amgout.res_cf = residual(i+1)/residual(i);
            case {'pcg'}
                r = b - A*v;
                z = zeros(size(b));
                z = vcycle2(A,z,r,1,vdata,amgdata);
                p = z;
                for i = 1 : maxCycle
                    nu = r'*z;
                    Ap = A*p;
                    alpha = nu/(p'*Ap);
                    v = v + alpha * p;
                    r = r - alpha * Ap;
                    residual(i+1) = norm(b - A*v);
                    errors(i+1) = norm(v);
                    if strcmp(logflag,'on')
                        fprintf(fid2,'%3d   ||%10.2e%6.2f%6.2f |%9.2e%6.2f%6.2f\n',i,residual(i+1), residual(i+1)/residual(i),...
                            (residual(i+1)/residual(1))^(1/i),errors(i+1),errors(i+1)/errors(i),...
                            (errors(i+1)/errors(1))^(1/i));
                    end
                    if (residual(i+1)/residual(1)) < tol %|| r(i+1)/r(i) > .95
                        %                 lastindex = i+1;
                        break
                    end
                    z = zeros(size(b));
                    z = vcycle2(A,z,r,1,vdata,amgdata);
                    beta = r'*z/nu;
                    p = z + beta * p;
                end
                amgout.amgcycle = i;
                amgout.error_cf = errors(i+1)/errors(i);
                amgout.res_cf = residual(i+1)/residual(i);
            otherwise
        end
        amgout.residual = residual;
        amgout.errors = errors;
        amgtime.cycle = toc(cycleTic);
    otherwise
        error('MATLAB:funfun:amg:intpType:notAIntpType',...
            'Invalid value for OPTIONS parameter %s: must be ''amg'', ''aamg'', ''rbamg'' or ''off''.',intpType);
end
if strcmp(logflag,'on')
    fclose(fid2);
end
amgout.nofdof = amgdata.nofdof;
amgout.nnz = amgdata.nnz;
amgout.cst = amgdata.cst;

if strcmp(returnMtx,'on')
    amgout.P     = amgdata.P;
    amgout.PT    = amgdata.PT;
    amgout.Mtx   = amgdata.Mtx;
    amgout.IA    = amgdata.IA;
    amgout.JA    = amgdata.JA;
    amgout.AA    = amgdata.AA;
end
