function v  = vcycle2(A,v,f,level,vdata,amgdata)
% VCYCLE   Run V(or mu) cycle.
%   V =
%   VCYCLE(A,V,F,ITER1,ITER2,LEVEL,RDATA,COARSEST,PS,MS,IA,JA,AA,MUCYCLE)
%
%   See also amginitsetup, rdata.

%   Minho Park
%   Department of Applied Mathematics
%   University of Colorado at Boulder
%   parkmh@colorado.edu
%   Date : 2008-11-24
%   Update : 

if level == vdata.cst
%     size(A)
    v =  A\f;
elseif level < vdata.cst
    % Pre Relaxation
    switch vdata.RelType
        case {'gs'}  % Gauss_Seidel method
            gs(amgdata.IA{level},amgdata.JA{level},amgdata.AA{level},...
                v,f,vdata.RelPara,vdata.Iter1);
        case {'jacobi'}  % Jacobi method
            jacobi(amgdata.IA{level},amgdata.JA{level},amgdata.AA{level},...
                v,f,vdata.RelPara,vdata.Iter1);
        otherwise
    end
    
    % Restrict the fine grid residual00
    f_2h = amgdata.PT{level}*(f-A*v);
    
    v_2h = zeros(size(f_2h));
    
    if level < vdata.cst -1
        for i = 1 : vdata.Mu
            v_2h = vcycle2(amgdata.Mtx{level+1},v_2h,f_2h,level+1,vdata,amgdata);            
        end
    else
        v_2h = vcycle2(amgdata.Mtx{level+1},v_2h,f_2h,level+1,vdata,amgdata);
    end
    
    % Correction
    v = v + amgdata.P{level}*v_2h;
    
    % Post Relaxation
    switch vdata.RelType
        case {'gs'}  % Gauss_Seidel method
            bgs(amgdata.IA{level},amgdata.JA{level},amgdata.AA{level},...
                v,f,vdata.RelPara,vdata.Iter2);
        case {'jacobi'}  % Jacobi method
            jacobi(amgdata.IA{level},amgdata.JA{level},amgdata.AA{level},...
                v,f,vdata.RelPara,vdata.Iter2);
        otherwise
    end    
else
    
end
