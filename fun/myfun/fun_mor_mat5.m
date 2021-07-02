function [yR,yA,yL,V,b_hat]=fun_mor_mat5(N_check,MOR_tol,MOR_maxdim,rhs_vectR,...
                        R,Aee,AeeR,opt_amg_MATMAG,...
                        opt_amg_AGMG,L,M,N,...
                        dx,dy,dz,Ind,idxF,...
                        C,idxFx,idxFy,idxFz,st_sparse_preconM,opCirculantM_all,...
                        Kt,prec,...
                        inner_it,outer_it,tol)
% 
N_dofs=length(rhs_vectR);
%%
V = []; % 
Nmat=size(Ind,2);
point_rho=zeros(Nmat,1);
mat_rho=zeros(Nmat,2);
for ii = 1:Nmat
    yR(ii).R=[];
    if strcmp(Ind(ii).tag,'cond')
       mat_rho(ii,:)=Ind(ii).rho;  
    end
end
%%
%%
count=1;
disp('...start mor...')
check_stop=0;
residual_stor=zeros(MOR_maxdim,1);
brhs=rhs_vectR;
Nf=length(idxF);
Np=N_dofs-Nf;
e=brhs(1:Nf);
miwb=C.'*e;
clear e
% R fixed
% z_realx_loc=R(1).Rx*0;
% z_realy_loc=R(1).Ry*0;
% z_realz_loc=R(1).Rz*0;
% mat_fixed=abs(1-[Ind(:).mor]);
% rhopoint=0.5*(mat_rho(:,2)+mat_rho(:,1))+point_rho.*0.5.*(mat_rho(:,2)-mat_rho(:,1));
% for ii =1:Nmat
%     z_realx_loc=z_realx_loc+mat_fixed(ii)*rhopoint(ii)*R(ii).Rx;
%     z_realy_loc=z_realy_loc+mat_fixed(ii)*rhopoint(ii)*R(ii).Ry;
%     z_realz_loc=z_realz_loc+mat_fixed(ii)*rhopoint(ii)*R(ii).Rz;
% end
% rloc_fixed=[(dx/(dy*dz))*z_realx_loc;(dy/(dx*dz))*z_realy_loc;(dz/(dy*dx))*z_realz_loc];
fMVM_L = @(J) multiplyMATVECT_EDDY(J,opCirculantM_all,R(1).Rx*0,...
R(1).Ry*0,R(1).Rz*0,idxFx,idxFy,idxFz,[dx dy dz],0*AeeR,L,M,N);
    
while check_stop==0%r > in_mor_tol
    %% 
    disp('...create R...')
    tic
    rhopoint=0.5*(mat_rho(:,2)+mat_rho(:,1))+point_rho.*0.5.*(mat_rho(:,2)-mat_rho(:,1));
    % 
    z_realx_loc=R(1).Rx*0;
    z_realy_loc=R(1).Ry*0;
    z_realz_loc=R(1).Rz*0;
    for ii =1:Nmat
        z_realx_loc=z_realx_loc+rhopoint(ii)*R(ii).Rx;
        z_realy_loc=z_realy_loc+rhopoint(ii)*R(ii).Ry;
        z_realz_loc=z_realz_loc+rhopoint(ii)*R(ii).Rz;
    end
    rloc=[(dx/(dy*dz))*z_realx_loc;...
          (dy/(dx*dz))*z_realy_loc;...
          (dz/(dy*dx))*z_realz_loc];
    toc
    %% 
    disp('...solve...')
    tic
    [A_inv,LL,UU,PP,QQ,RR,Sch_comp] = ...
        preparePREC_NEW2(rloc,idxFx,idxFy,...
        idxFz,st_sparse_preconM,AeeR,Aee,Kt,prec,opt_amg_MATMAG,opt_amg_AGMG);
    fPMV = @(xx)multiplyPREC_CAP_NEW(xx,AeeR,A_inv,LL,UU,PP,QQ,RR,...
        Sch_comp,prec,opt_amg_MATMAG,opt_amg_AGMG);
    fMVM = @(J) multiplyMATVECT_EDDY(J,opCirculantM_all,z_realx_loc,...
    z_realy_loc,z_realz_loc,idxFx,idxFy,idxFz,[dx dy dz],AeeR,L,M,N);
    [v, flag, relres, iter, resvec] = pgmres_mod_no_plot(@(J)fMVM(J),brhs, inner_it, tol, outer_it, @(JOut_full_in)fPMV(JOut_full_in) );  
    toc
    disp(' ')
    %% 
    disp('...gso...')
    tic
    [V, ~]=mgson_update([V,v]);  % Gram-Schmidt
    toc
    disp(' ')
    %% 
    % R
    disp('...VT*R*V...')  
    tic
    for ii = 1:Nmat
        z_realx_loc=R(ii).Rx;
        z_realy_loc=R(ii).Ry;
        z_realz_loc=R(ii).Rz;
        rloc=[(dx/(dy*dz))*z_realx_loc;(dy/(dx*dz))*z_realy_loc;(dz/(dy*dx))*z_realz_loc];  
        clear z_realx_loc z_realy_loc z_realz_loc
        %         tmpR(ii).R(:,count)=[rloc.*V(1:Nf,count);zeros(Np,1)];
        tmpRv=[rloc.*V(1:Nf,count);zeros(Np,1)];
%         tmpvR=[V(1:Nf,count).*rloc;zeros(Np,1)];
%         clear rloc
        tmp1=V'*tmpRv;
        tmp2=tmpRv'*V;
%         clear tmpRv tmpRv
        yR(ii).R=[yR(ii).R,tmp1(1:end-1);
            tmp2];
        
    end
    toc
% L and A
   disp('...VT*A*V...') 
   tic
   tmpL(:,count)=[fMVM_L(V(:,count))];
   tmpA(:,count)=[AeeR.'*V(Nf+1:end,count);...
         AeeR*V(1:Nf,count)];
    yA=V' * tmpA;
    yL=V' * tmpL;
    toc
    %% rhs poblema ridotto
    disp('...rhs...') 
    tic
    b_hat = V' * brhs;
    toc
    %% random mode
    disp('...check...')
    tic
    residual = 0;
    N_check2=N_check;%max(N_check,count);
    res_sto=zeros(N_check2,1);
    Irho_store=zeros(Nmat,N_check2);
    for jj = 1:N_check2
        Irho=(rand(Nmat,1)-0.5)*2;
        Irho_store(:,jj)=Irho; 
        rhopoint=0.5*(mat_rho(:,2)+mat_rho(:,1))+Irho.*0.5.*(mat_rho(:,2)-mat_rho(:,1));
        z_realx_loc=R(1).Rx*0;
        z_realy_loc=R(1).Ry*0;
        z_realz_loc=R(1).Rz*0;
        yRloc=yR(1).R*0;
        for ii = 1:Nmat
            yRloc=yRloc+rhopoint(ii)*yR(ii).R;
            z_realx_loc=z_realx_loc+rhopoint(ii)*R(ii).Rx;
            z_realy_loc=z_realy_loc+rhopoint(ii)*R(ii).Ry;
            z_realz_loc=z_realz_loc+rhopoint(ii)*R(ii).Rz;            
        end
        rloc=[(dx/(dy*dz))*z_realx_loc;(dy/(dx*dz))*z_realy_loc;(dz/(dy*dx))*z_realz_loc];
        A_hat = yRloc+yA+yL;
        %
        x_hat = A_hat \ (b_hat);
        xx=V*x_hat;
        
        lxx=fMVM_L(xx);
        z=rloc.*xx(1:Nf)+lxx(1:Nf);
        rhs=C.'*z(1:Nf);
        residual = norm(rhs-miwb) / norm(miwb);
        res_sto(jj)=residual;
    end
    [res_max,ins_res_max]=max(res_sto);
    residual_stor(count,1)=res_max;
    NcheckOK=nnz(find(res_sto<MOR_tol));
    if res_max<MOR_tol || count>MOR_maxdim
       disp(['NcheckOK = ',num2str(NcheckOK),' ',num2str(res_max)]) 
       check_stop=1;
    else
        point_rho=Irho_store(:,ins_res_max);
        disp(['NcheckOK = ',num2str(NcheckOK),' ',num2str(res_max)])         
    end
    disp('')
        figure(667)
        semilogy(count,(res_max),'o')
        hold on
        semilogy(count*ones(N_check2,1),(res_sto),'.')
		text(count,(res_max),num2str(NcheckOK))
        title('mor')
        ylim([sort([MOR_tol/10 max(residual_stor)])])
        drawnow
        count=count+1;
        toc
        %%
end
%%

end