% function [amgdata amgtime] = amginitsetup(A,opt)
function [amgdata]= amginitsetup(A,opt)

%AMGINITSETUP Generate AMGDATA which includes interpolation matrices,
%  coarser grid matrices and target vector(s) if necessary.
%  AMGDATA = AMGINITSETUP(A,OPTIONS) computes all matrices from
%  AMG initial setup process.
%
%  See also amg, amgget.

% changelog
% v.0.1.2 - Added transpose of interpolation to amgdata.
% v.0.1.2 - Added aggressive coarsening and interpolation

% tTotal = tic;
% Get values from amg OPTIONS file
cst     = amgget(opt,'Coarsest',2);
csnType = amgget(opt,'CsnType','amg');
intpType= amgget(opt,'IntpType','amg');
relType = amgget(opt,'RelType','gs');
relPara = amgget(opt,'RelPara',1);
nVec    = amgget(opt,'NumVec',5);
nRel    = amgget(opt,'NumRel',5);
theta   = amgget(opt,'Theta',.25);
nu      = amgget(opt,'Nu',5);
alpha   = amgget(opt,'Alpha',.7);
saveCsn = amgget(opt,'SaveCsn','off');

% Allocate all cells
amgdata.P   = cell(cst-1,1); % Interpolation
amgdata.PT  = cell(cst-1,1); % Transpose of Interpolation
amgdata.IA  = cell(cst-1,1); % Row index
amgdata.JA  = cell(cst-1,1); % Column index
amgdata.AA  = cell(cst-1,1); % Matric values
amgdata.Mtx = cell(cst ,1);  % Matrix in matlab format
amgdata.X   = cell(cst-1,1); % target vectors
amgdata.C   = cell(cst-1,1); % coarsening data

% % Set all times zero
% amgtime.P = zeros(cst-1,1);
% amgtime.RAP = zeros(cst-1,1);
% amgtime.C = zeros(cst-1,1);

if strcmp(csnType,'amg')
    amgdata.IS = cell(cst-1,1);
    amgdata.JS = cell(cst-1,1);
    amgdata.IT = cell(cst-1,1);
    amgdata.JT = cell(cst-1,1);
end

if strcmp(csnType,'a21') || strcmp(csnType,'a22')
    amgdata.S = cell(cst-1,1);
    amgdata.S2 = cell(cst-1,1);
    amgtime.C2 = zeros(cst-1,1);
end

amgdata.Mtx{1} = A;

amgdata.nofdof = zeros(1,cst);  % number of degree of freedoms
amgdata.nnz    = zeros(1,cst);  % number of nonzeros

% Load the saved coarsening data
if strcmp(csnType,'load')
    load CsnData;   
end

% Initail set up
for i = 1 : cst - 1     
    % save the fine grid matrix data into amgdata
    amgdata.nnz(i) = nnz(amgdata.Mtx{i});
    [amgdata.IA{i} amgdata.JA{i} amgdata.AA{i}] = sparse2crs(amgdata.Mtx{i});
    amgdata.nofdof(i) = length(amgdata.IA{i})-1;
    sTime = 0; sTime2 = 0; cTime = 0; cTime2 = 0;
    % Coarsening process
    switch csnType
        case {'amg'}
%             sTic = tic;

            [nofc amgdata.IS{i} amgdata.JS{i} amgdata.IT{i} amgdata.JT{i}] ...
                = strongconnection(amgdata.IA{i},amgdata.JA{i},amgdata.AA{i},theta,'negative');
%             sTime = toc;
%             sTime = toc(sTic);
%             cTic = tic;
            
            amgdata.C{i} = amgcoarsening(amgdata.IA{i},amgdata.JA{i},amgdata.AA{i},amgdata.IS{i},amgdata.JS{i},...
                                         amgdata.IT{i},amgdata.JT{i},nofc);
%             cTime = toc;
%             cTime = toc(cTic);
        case {'amg(abs)'}
%             sTic = tic;
            [nofc amgdata.IS{i} amgdata.JS{i} amgdata.IT{i} amgdata.JT{i}] ...
                = strongconnection(amgdata.IA{i},amgdata.JA{i},amgdata.AA{i},theta,'absolute');
%             sTime = toc;
%             sTime = toc(sTic);
%             cTic = tic;
            amgdata.C{i} = amgcoarsening(amgdata.IA{i},amgdata.JA{i},amgdata.AA{i},amgdata.IS{i},amgdata.JS{i},...
                                         amgdata.IT{i},amgdata.JT{i},nofc);
%             cTime = toc;
%             cTime = toc(cTic);
        case {'cr'}
%             cTic = tic;
            amgdata.C{i} = cr_coarsening(amgdata.IA{i},amgdata.JA{i},amgdata.AA{i},nu,alpha);
%             cTime = toc;
%             cTime = toc(cTic);
        case {'load'}
            
        case {'a21'}
            if i == 1
%                 sTic = tic;
                [nofc IS JS IT JT] ...
                    = strongconnection(amgdata.IA{i},amgdata.JA{i},amgdata.AA{i},theta,'negative');
%                 sTime = toc;
%                 sTime = toc(sTic);
%                 cTic = tic;
                C1 = amgcoarsening(amgdata.IA{i},amgdata.JA{i},amgdata.AA{i},IS,JS,...
                    IT,JT,nofc);
%                 save C1 C1
%                 cTime = toc;
%                 cTime = toc(cTic);
%                 sTic2 = tic;
                [S inbd jnbd] = length2strong(amgdata.IA{i},amgdata.JA{i},C1,1);
%                 sTime2 = toc;
%                 sTime2 = toc(sTic2);
%                 cTic2 = tic;
                amgdata.C{i} = aggcoarsening(S,C1,amgdata.IA{i},amgdata.JA{i},inbd,jnbd);
%                 cTime2 = toc;
%                 cTime2 = toc(cTic2);
            else
%                 sTic = tic;
                [nofc amgdata.IS{i} amgdata.JS{i} amgdata.IT{i} amgdata.JT{i}] ...
                    = strongconnection(amgdata.IA{i},amgdata.JA{i},amgdata.AA{i},theta,'negative');
%                 sTime = toc;
%                 sTime = toc(sTic);
                
%                 cTic = tic;
                amgdata.C{i} = amgcoarsening(amgdata.IA{i},amgdata.JA{i},amgdata.AA{i},amgdata.IS{i},amgdata.JS{i},...
                    amgdata.IT{i},amgdata.JT{i},nofc);
%                 cTime = toc;
%                 cTime = toc(cTic);
            end
        case {'a22'}
            if i == 1
%                 sTic = tic;
                [nofc IS JS IT JT] ...
                    = strongconnection(amgdata.IA{i},amgdata.JA{i},amgdata.AA{i},theta,'negative');
%                 sTime = toc;
%                 sTime = toc(sTic);
%                 cTic = tic;
                C1 = amgcoarsening(amgdata.IA{i},amgdata.JA{i},amgdata.AA{i},IS,JS,...
                    IT,JT,nofc);
%                 cTime = toc;
%                 cTime = toc(cTic);
%                 sTic2 = tic;
                [amgdata.S{i} inbd jnbd] = length2strong(amgdata.IA{i},amgdata.JA{i},C1,2);
%                 sTime2 = toc;
%                 sTime2 = toc(sTic2);
%                 cTic2 = tic;
                amgdata.C{i} = aggcoarsening(amgdata.S{i},C1,amgdata.IA{i},amgdata.JA{i},inbd,jnbd);
%                 cTime2 = toc;
%                 cTime2 = toc(cTic2);
            else
%                  sTic = tic;
                [nofc amgdata.IS{i} amgdata.JS{i} amgdata.IT{i} amgdata.JT{i}] ...
                    = strongconnection(amgdata.IA{i},amgdata.JA{i},amgdata.AA{i},theta,'negative');
%                 sTime = toc;
%                 sTime = toc(sTic);
%                 cTic = tic;
                amgdata.C{i} = amgcoarsening(amgdata.IA{i},amgdata.JA{i},amgdata.AA{i},amgdata.IS{i},amgdata.JS{i},...
                    amgdata.IT{i},amgdata.JT{i},nofc);
%                 cTime = toc;
%                 cTime = toc(cTic);
            end
        otherwise
            error('MATLAB:funfun:amginitget:csnType:notACsnType',...
                'Invalid value for OPTIONS parameter %s: must be ''amg'', ''cr'', ''load'', ''a1'' or ''a2''.',csnType);
    end
%     cTime
%     amgtime.C(i) = cTime;
%     amgtime.S(i) = sTime;
%     amgtime.C2(i) = cTime2;
%     amgtime.S2(i) = sTime2;
%     pTime = tic;
    
    % Compute interpolation weight    
    switch intpType
        case {'amg'}
            amgdata.P{i} = amgp(amgdata.IA{i},amgdata.JA{i},amgdata.AA{i},...
                amgdata.IS{i},amgdata.JS{i},...
                amgdata.C{i});            
        case {'lramg'}
            if i == 1
                amgdata.P{i} = longrangeamgp(amgdata.IA{i},amgdata.JA{i},amgdata.AA{i},inbd,jnbd,amgdata.C{i});                
%                 P = amgdata.P{i};
%                 save P P
            else
                amgdata.P{i} = amgp(amgdata.IA{i},amgdata.JA{i},amgdata.AA{i},...
                amgdata.IS{i},amgdata.JS{i},...
                amgdata.C{i});
            end
%             P = amgdata.P{i};
        case {'aamg'}
            if i == 1   % if finest level
                %TODO zero-crossing & random vector
                amgdata.X{i} = rand(amgdata.nofdof(i),1);
            else        % if not finest level
                amgdata.X{i} = (amgdata.X{i-1}(amgdata.C{i-1}==1));
            end
            
            % Apply relaxation to target vector
            switch relType
                case {'gs'}
                    gs(amgdata.IA{i},amgdata.JA{i},amgdata.AA{i},...
                        amgdata.X{i},zeros(length(IA)-1),relPara,nRel);
                case {'jacobi'}
                    jacobi(amgdata.IA{i},amgdata.JA{i},amgdata.AA{i},...
                        amgdata.X{i},zeros(length(IA)-1),relPara,nRel);
                otherwise
                    error('MATLAB:funfun:amginitsetup:relType:notARelType',...
                        'Invalid value for OPTIONS parameter %s: must be ''gs'' or ''jacobi''.',relType);
            end
            
            % Build interpolation operator
            amgdata.P{i} = aamgp(amgdata.IA{i},amgdata.JA{i},amgdata.AA{i},amgdata.C{i},amgdata.X{i});
        case {'rbamg'}
            amgdata.X{i} = zeros(nVec,length(amgdata.IA{i})-1);
            for j = 1 : nVec
                if i == 1    % if finest level
                    x = rand(length(amgdata.IA{i})-1,1);
                    x = x/norm(x);   % normalied in L2 norm
                else         % if not finest level
                    xf = amgdata.X{i-1}(j,:).';
                    x = xf(amgdata.C{i-1} == 1); % Injection
                end
                
                % Apply relaxation to target vector
                switch relType
                    case {'gs'}
                        gs(amgdata.IA{i},amgdata.JA{i},amgdata.AA{i},x,zeros(length(amgdata.IA{i})-1),...
                            relPara,nRel);
                    case {'jacobi'}
                        jacobi(amgdata.IA{i},amgdata.JA{i},amgdata.AA{i},x,zeros(length(amgdata.IA{i})-1),...
                            relPara,nRel);
                    otherwise
                        error('MATLAB:funfun:amginitsetup:relType:notARelType',...
                            'Invalid value for OPTIONS parameter %s: must be ''gs'' or ''jacobi''.',relType);
                end
                amgdata.X{i}(j,:) = x.'/sqrt(x'*amgdata.Mtx{i}*x);
            end
            amgdata.P{i} = rbamgp(amgdata.IA{i},amgdata.JA{i},amgdata.AA{i},amgdata.C{i},amgdata.X{i});
            
            for j = 1 : nVec
            subplot(2,1,1)
            plot(amgdata.X{i}(j,:)' - amgdata.P{i} *amgdata.X{i}(j,amgdata.C{i}==1)','*-' )
            pause
            end
        otherwise
            error('MATLAB:funfun:amginitsetup:intpType:notAIntpType',...
                'Invalid value for OPTIONS parameter %s: must be ''amgp'', ''aamgp'', ''rbamgp'' or ''off''.',intpType);
    end
    amgdata.PT{i} = amgdata.P{i}';
%     amgtime.P(i) = toc;
%     amgtime.P(i) = toc(pTime);
    
%     rapTime = tic;
    % Build coarse grid operator
    % TODO HPD
    amgdata.Mtx{i+1} = amgdata.PT{i}*amgdata.Mtx{i}*amgdata.P{i};
%     amgtime.RAP(i) = toc;
%     amgtime.RAP(i) = toc(rapTime);
    
    if length(amgdata.Mtx{i+1}) < 101
%             warning('[NOTE] Since the size of matrix on level %d is less than or equal to 100, \nwe stop the coarsening at this point\n',i+1);
        break
    end
    if nnz(amgdata.Mtx{i+1})/length(amgdata.Mtx{i+1})^2 > .7
        break
    end
end

amgdata.nnz(i+1) = nnz(amgdata.Mtx{i+1});

amgdata.nofdof(i+1) = length(amgdata.Mtx{i+1});
amgdata.cst = i+1;

if strcmp(saveCsn,'on')
    C = amgdata.C;
    IS = amgdata.IS;
    JS = amgdata.JS;
    IT = amgdata.IT;
    JT = amgdata.JT;
    if strcmp(csnType,'a21') || strcmp(csnType,'a22')
    save CsnData.mat C IS JS IT JT inbd jnbd -v7.3    
    else
    save CsnData.mat amgdata -v7.3
    end
end

% amgtime.total = toc(tTotal);