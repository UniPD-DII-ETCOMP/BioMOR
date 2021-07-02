%%
% This function has been inspired by 
% lse_sparse_precon_multiply.m which is available at
% https://github.com/acyucel/VoxHenry
%%
function [JOut_full_out] = multiplyPREC_CAP_NEW(JOut_full_in,AeR,A_inv,...
    LL,UU,PP,QQ,RR,Sch_comp,prec,opt_amg_MATMAG,opt_amg_AGMG)
    num_nodeR=size(AeR,1); 
    num_curr=size(AeR,2);
    JOut_full_out=zeros(num_nodeR+num_curr,1);
if strcmp(prec,'lu')
    warning off
    JOut_full_out(num_curr+1:num_curr+num_nodeR) = ...
        QQ * (UU \ (LL \ (PP * (RR \ (JOut_full_in(num_curr+1:num_curr+num_nodeR) - (AeR * (A_inv * JOut_full_in(1:num_curr))) ) ))));
    warning on
elseif strcmp(prec,'amg_MATMAG')
        v=zeros(size(Sch_comp,1),1);
        JOut_full_out(num_curr+1:num_curr+num_nodeR) = (JOut_full_in(num_curr+1:num_curr+num_nodeR) - ( AeR * ( A_inv * JOut_full_in(1:num_curr)) ) );        
        bb=JOut_full_out(num_curr+1:num_curr+num_nodeR);
        [xx_re, amgout] = myRTamg2(-Sch_comp,v,real(-bb),opt_amg_MATMAG);
        [xx_im, amgout] = myRTamg2(-Sch_comp,v,imag(-bb),opt_amg_MATMAG);
         xx=xx_re+1j*xx_im;
         JOut_full_out(num_curr+1:num_curr+num_nodeR)=xx;
elseif strcmp(prec,'amg_AGMG')
        v=zeros(size(Sch_comp,1),1);
        JOut_full_out(num_curr+1:num_curr+num_nodeR) = (JOut_full_in(num_curr+1:num_curr+num_nodeR) - ( AeR * ( A_inv * JOut_full_in(1:num_curr)) ) );        
        bb=JOut_full_out(num_curr+1:num_curr+num_nodeR);
        [xx_re,flag,relres_amg,iter,resvec]=agmg(-Sch_comp,real(-bb),opt_amg_AGMG.restart,opt_amg_AGMG.tol,opt_amg_AGMG.maxit,0,v,2);
        [xx_im,flag,relres_amg,iter,resvec]=agmg(-Sch_comp,imag(-bb),opt_amg_AGMG.restart,opt_amg_AGMG.tol,opt_amg_AGMG.maxit,0,v,2);
        xx=xx_re+1j*xx_im;
        JOut_full_out(num_curr+1:num_curr+num_nodeR)=xx;
end
    JOut_full_out(1:num_curr) = A_inv*(JOut_full_in(1:num_curr) - ((AeR.')*JOut_full_out(num_curr+1:num_curr+num_nodeR)));

end