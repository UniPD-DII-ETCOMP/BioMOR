function [yR,yL,V,b_hat]=fun_mor_freq(N_check,MOR_tol,MOR_maxdim,w_min,w_max,rhs_w1,...
                        d,z_realx_loc,z_realy_loc,z_realz_loc,...
                        idxFx,idxFy,idxFz,AeeR,Aee,Kt,inner_it,tol,...
                        outer_it,opCirculantM_all,st_sparse_preconM,...
                        prec,L,M,N,z_realF,opt_amg_MATMAG,opt_amg_AGMG)
% 
N_dofs=length(rhs_w1);
%%
V = []; % 
point_freq=0.0;
%%
count=1;
disp('...start mor...')
check_stop=0;
residual_stor=zeros(MOR_maxdim,1);
figure(667)
while check_stop==0%
    %% 
    wpoint=0.5*(w_max+w_min)+point_freq*0.5*(w_max-w_min);
    %% 
    % PRECONDITIONER
    disp('----COMPUTING PRECONDITIONER--------------------------------')
    mytic_prec=tic;
    [A_inv,LL,UU,PP,QQ,RR,Sch_comp] = ...
        preparePREC_NEW(d,z_realF,idxFx,idxFy,...
        idxFz,st_sparse_preconM*wpoint,AeeR,Aee,Kt,prec,opt_amg_MATMAG,opt_amg_AGMG);
    fPMV = @(xx)multiplyPREC_CAP_NEW(xx,AeeR,A_inv,LL,UU,PP,QQ,RR,...
        Sch_comp,prec,opt_amg_MATMAG,opt_amg_AGMG);
    disp([' Total time for computing preconditioner ::: ' ,num2str(toc(mytic_prec))]);
    disp(' ')
    % SYSTEM: A*X
    fMVM = @(J) multiplyMATVECT_EDDY(J,opCirculantM_all*wpoint,z_realx_loc,...
        z_realy_loc,z_realz_loc,idxFx,idxFy,idxFz,d,AeeR,L,M,N);
    %% 
    % RHS
    brhs=(wpoint*rhs_w1);
    %   
    disp('solving...')
    [v, flag, relres, iter, resvec] =  pgmres_mod(@(J)fMVM(J),brhs, inner_it, tol, outer_it, @(JOut_full_in)fPMV(JOut_full_in) );    
    %
    bre=fMVM(v);
    err_rhs=norm(bre-brhs)/norm(brhs);
    disp(['error on rhs = ',num2str(err_rhs)])
    if err_rhs>MOR_tol
        warning('MOR tollerance has been changed...')
        disp(['from ',num2str(MOR_tol)]);
        MOR_tol=err_rhs;
        disp(['to ',num2str(MOR_tol)]);
    end
    %% 
    [V, ~]=mgson([V,v]);  % Gram-Schmidt
    %% 
    YR = zeros(N_dofs, size(V, 2));
    yR = zeros(size(V, 2), size(V, 2));
    YL = zeros(N_dofs, size(V, 2));
    yL = zeros(size(V, 2), size(V, 2));
    % L
    fMVM_L = @(J) (1/1j)*multiplyMATVECT_EDDY(J,opCirculantM_all,...
        0*z_realx_loc,0*z_realy_loc,0*z_realz_loc,idxFx,idxFy,idxFz,d,0*AeeR,L,M,N);
    for cc=1:size(V,2)
        YL(:, cc) = fMVM_L(V(:,cc));%HL*V;
    end
    yL(:, :) = V' * YL(:, :);
    % R
    fMVM_R = @(J) multiplyMATVECT_EDDY(J,0*opCirculantM_all,...
        z_realx_loc,z_realy_loc,z_realz_loc,idxFx,idxFy,idxFz,d,AeeR,L,M,N);
    for cc=1:size(V,2)
        YR(:, cc) = fMVM_R(V(:,cc));%HR*V;
    end
    yR(:, :) = V' * YR(:, :);
    %% 
    b_hat = V' * rhs_w1;
    %% 
    N_check2=max(N_check,count);
    res_sto=zeros(N_check2,1);
    Ifreq_store=zeros(N_check2,1);
    for jj = 1:N_check2
        Ifreq=(rand(1,1)-0.5)*2;
        Ifreq_store(jj)=Ifreq;
        wpoint=0.5*(w_max+w_min)+Ifreq*0.5*(w_max-w_min);
        A_hat = yR(:, :)+1j*wpoint*yL(:, :);
        %
        x_hat = A_hat \ (wpoint*b_hat);
        z = (YR(:,:)+1j*wpoint*YL(:, :)) * x_hat;
        zz=z-(wpoint*rhs_w1);
        residual = norm(zz) / norm(wpoint*rhs_w1);
        res_sto(jj)=residual;
    end
    [res_max,ins_res_max]=max(res_sto);
    residual_stor(count,1)=res_max;
    NcheckOK=nnz(find(res_sto<MOR_tol));
    if res_max<MOR_tol || count>MOR_maxdim
       disp(['NcheckOK = ',num2str(NcheckOK),' ',num2str(res_max)]) 
       check_stop=1;
    else
        point_freq=Ifreq_store(ins_res_max);
        disp(['NcheckOK = ',num2str(NcheckOK),' ',num2str(res_max)])         
    end
        figure(667)
        semilogy(count,(res_max),'o')
        hold on
        semilogy(count*ones(N_check2,1),(res_sto),'.')
		text(count,(res_max),num2str(NcheckOK))
        title('mor')
        ylim([MOR_tol/10 max(residual_stor)])
        drawnow
        count=count+1;
        disp(' ')
end
end