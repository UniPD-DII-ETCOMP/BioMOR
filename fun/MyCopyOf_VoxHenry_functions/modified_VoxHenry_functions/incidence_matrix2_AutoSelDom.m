%% 
% This function has been inspired by lse_compute_Ae_matrix.m available at
% https://github.com/acyucel/VoxHenry/ 
function [Ae,Aee,AeeR,idxF,idxFx,idxFy,...
    idxFz,Ae1x,Ae1y,Ae1z,id_sel,dim_domain,bins2,Ndom] = ...
    incidence_matrix2_AutoSelDom(Kt,LMN,idxV)
L = LMN(1); M = LMN(2); N = LMN(3);
Ae=sparse(Kt,3*Kt);
Ae=Ae+sparse(1:Kt,1:Kt,ones(Kt,1),Kt,3*Kt); % 
Ae=Ae+sparse(1:Kt,Kt+1:2*Kt,ones(Kt,1),Kt,3*Kt); % 
Ae=Ae+sparse(1:Kt,2*Kt+1:3*Kt,ones(Kt,1),Kt,3*Kt); % 
idfx_fv=zeros(Kt,2); % 
idfy_fv=zeros(Kt,2); % 
idfz_fv=zeros(Kt,2); % 
for ii = 1:L-1
    for jj = 1:M
        for kk = 1:N
         idfx_fv(ii+(jj-1)*L+(kk-1)*M*L,1)=ii  +(jj  -1)*L+(kk  -1)*M*L;
         idfx_fv(ii+(jj-1)*L+(kk-1)*M*L,2)=ii+1+(jj  -1)*L+(kk  -1)*M*L;
        end
    end
end
Ae1x=zeros(2,Kt); 
Ae1x(1,:)=1:Kt;
Ae1x(2,:)=idfx_fv(:,2);
for ii = 1:L
    for jj = 1:M-1
        for kk = 1:N
         idfy_fv(ii+(jj-1)*L+(kk-1)*M*L,1)=ii  +(jj  -1)*L+(kk  -1)*M*L;
         idfy_fv(ii+(jj-1)*L+(kk-1)*M*L,2)=ii  +(jj+1-1)*L+(kk  -1)*M*L;
        end
    end
end
Ae1y=zeros(2,Kt);
Ae1y(1,:)=1:Kt;
Ae1y(2,:)=idfy_fv(:,2);
for ii = 1:L
    for jj = 1:M
        for kk = 1:N-1
         idfz_fv(ii+(jj-1)*L+(kk-1)*M*L,1)=ii  +(jj  -1)*L+(kk  -1)*M*L;
         idfz_fv(ii+(jj-1)*L+(kk-1)*M*L,2)=ii  +(jj  -1)*L+(kk+1-1)*M*L;
        end
    end
end
Ae1z=zeros(2,Kt);
Ae1z(1,:)=1:Kt;
Ae1z(2,:)=idfz_fv(:,2);
sel=find(idfx_fv(:,2));
Ae=Ae+sparse(idfx_fv(sel,2),idfx_fv(sel,1),-ones(length(sel),1),Kt,3*Kt);
sel=find(idfy_fv(:,2));
Ae=Ae+sparse(idfy_fv(sel,2),Kt+idfy_fv(sel,1),-ones(length(sel),1),Kt,3*Kt);
sel=find(idfz_fv(:,2));
Ae=Ae+sparse(idfz_fv(sel,2),2*Kt+idfz_fv(sel,1),-ones(length(sel),1),Kt,3*Kt);
clear idfx_fv idfx_fv idfz_fv
clear sel
Aee = Ae(idxV,:); % 
idxFx = find(sum(abs(Aee(:,1:Kt)),1)==2).'; % 
idxFy = 0*Kt+find(sum(abs(Aee(:,Kt+1:2*Kt)),1)==2).'; % 
idxFz = 0*2*Kt+find(sum(abs(Aee(:,2*Kt+1:3*Kt)),1)==2).'; % 
idxF = find(sum(abs(Aee),1)==2).'; % 
Aee = Aee(:,idxF); 
id_sel=[];
dim_domain=[];
bins2=[];
Ndom=[];
%%
[rowp,colp] = find(Aee==1);
[rown,coln] = find(Aee==-1);
G1=zeros(2,size(Aee,2));
G1(1,colp)=rowp;
G1(2,coln)=rown;
% disp('Start computing Graph');
mytic_graph = tic; 
Gmatlab2 = graph(G1(1,:).',G1(2,:).',[],size(Aee,1));
toc(mytic_graph);
% disp('Computing number of components');
% mytic_comp = tic;
try
  [bins2,~] = conncomp(Gmatlab2);    
catch
    warning('conncomp MATLAB functions was not detected, (slow) function graph_connected_components.m is used instead')
  CC=abs(Aee)*abs(Aee).'-diag(diag(abs(Aee)*abs(Aee).'));
  [bins2,~] = graph_connected_components(CC);    
end
% toc(mytic_comp);
unique_id=unique(bins2);
Ndom=length(unique_id);
disp([' Number of detected domains is ',num2str(Ndom)])
%%
%     warning('no imposed potential found !')
id_sel = zeros(Ndom,1);
dim_domain = zeros(Ndom,1);
for ii = 1:Ndom
    wwho=find(bins2==unique_id(ii));
    id_sel(ii)=wwho(1);
    dim_domain(ii)=length(wwho);
end
AeeR=Aee;
AeeR(id_sel,:)=[];
disp([' Size of full incidence matrix: ', num2str(size(Ae))])
disp([' Size of incidence matrix relative to non-empty voxels: ', num2str(size(Aee))])
disp([' Size of reduced incidence matrix: ', num2str(size(AeeR))])
end