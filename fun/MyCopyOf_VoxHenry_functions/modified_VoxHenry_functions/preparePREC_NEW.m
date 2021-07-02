%%
% This function has been inspired by 
% lse_sparse_precon_prepare.m which is available at
% https://github.com/acyucel/VoxHenry and in the directory VoxHenry_functions
%%
function [A_inv,LL,UU,PP,QQ,RR,Sch_comp] =  preparePREC_NEW(d,...
    z_realF,idxFx,idxFy,idxFz,...
    st_sparse_preconL,AeeR,Aee,Kt,prec,opt_amg_MATMAG,opt_amg_AGMG)
A_inv=[];
LL=[];
UU=[];
PP=[];
QQ=[];
RR=[];
Sch_comp=[];
dx = d(1); dy = d(2); dz = d(3);
num_curr = size(Aee,2);
num_currx=length(idxFx);
num_curry=length(idxFy);
num_currz=length(idxFz);
diag_pulse=zeros(num_curr,1);
diag_pulse(1:num_currx,1)                                        =1./(z_realF(idxFx)*dx/(dy*dz) + st_sparse_preconL(1));
diag_pulse(num_currx+1:num_currx+num_curry,1)                    =1./(z_realF(Kt+idxFy)*dy/(dz*dx) + st_sparse_preconL(2));
diag_pulse(num_currx+num_curry+1:num_currx+num_curry+num_currz,1)=1./(z_realF(Kt+Kt+idxFz)*dz/(dx*dy) + st_sparse_preconL(3));
inds=zeros(num_curr,3);
inds(1:num_curr,1)=[1:1:num_curr];
inds(1:num_curr,2)=inds(1:num_curr,1);
inds(:,3)=(diag_pulse);
if strcmp(prec,'lu')
    A_inv=sparse(inds(:,1),inds(:,2),inds(:,3));
    Sch_comp=  - (AeeR*A_inv*AeeR.');  
    [LL,UU,PP,QQ,RR] = lu(Sch_comp);    
elseif strcmp(prec,'amg_MATMAG')
    A_inv=sparse(inds(:,1),inds(:,2),abs(inds(:,3)));
    Sch_comp=  - (AeeR*A_inv*AeeR.');        
    global amgdata
    if isempty(amgdata)
        [amgdata] = amginitsetup(-Sch_comp,opt_amg_MATMAG);
    end
elseif strcmp(prec,'amg_AGMG')
    A_inv=sparse(inds(:,1),inds(:,2),abs(inds(:,3)));
    Sch_comp=  - (AeeR*A_inv*AeeR.');  
    verbose=1;
    x0=zeros(size(Sch_comp,1),1);
    ijob=1;
    [x,flag,relres,iter,resvec]=agmg(-Sch_comp,x0,opt_amg_AGMG.restart,opt_amg_AGMG.tol,opt_amg_AGMG.maxit,verbose,x0,ijob);        
end
end
