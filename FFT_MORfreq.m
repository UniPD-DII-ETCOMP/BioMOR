close all
clear global
clear
clc
restoredefaultpath
warning on
dad = pwd;
cd('fun'); addpath(genpath(pwd)); cd(dad)
cd('fortran'); addpath(pwd); cd(dad)
%%
%% BEGIN USER SETTINGS
%%
%% Directory
name_dir='test1';
%% Frequency
freq_min = 50; %[Hz]
freq_max = 10000; %[Hz]
%% Selections
plot_vectorsJ_flag = 1; %quiver plot of real and imag of J
paraview_export_flag = 1; % export to paraviw
Integration_flag = 'NumNum'; %'NumAn'; 'NumNum' (Integration: NumericalNumerical or AnalyticalNumerical)
ext_field_flag = 1; % exernal field
% below you can write the external magnetic vector potential as a function of x,y,z
% and omega. Active only if ext_field_flag=1
Ax_ext = @(x,y,z) -1j*0.5*y; 
Ay_ext = @(x,y,z)  1j*0.5*x; 
Az_ext = @(x,y,z)  0*z; 
%% Solver parameters
prec =  'amg_MATMAG'; % 'amg_MATMAG' 'amg_AGMG' 'lu';
tol = 1e-6;
inner_it = 40;
outer_it = 5;
% amg_MATMAG
opt_amg_MATMAG = amgset;                               % default option file
opt_amg_MATMAG = amgset(opt_amg_MATMAG,'coarsest',10);            % set the number of levels
opt_amg_MATMAG = amgset(opt_amg_MATMAG,'PreCond','pcg');          % set the Krylov method
opt_amg_MATMAG = amgset(opt_amg_MATMAG,'PrintOnScreen','on');    % turn off the option to print the log on the screen 
opt_amg_MATMAG = amgset(opt_amg_MATMAG,'SaveCsn','off');           % save the set of coarse-grid points
opt_amg_MATMAG = amgset(opt_amg_MATMAG,'CsnType','amg');          % choose the coarsening method
global amgdata
amgdata=[];
%
% amg_AGMG
opt_amg_AGMG.restart=1;
opt_amg_AGMG.tol=1e-8;
opt_amg_AGMG.maxit=200;
%% MOR
mor.N_check=10;
mor.MOR_tol=1e-3;
mor.MOR_maxdim=30;
%%
%% END USER SETTINGS
%%
%% Add Path
if strcmp(prec,'amg_AGMG')
cd fun
if ~exist('AGMG', 'dir')
   warning('You are trying to use "AGMG" but the directory is not detected.')
   disp(' ')
   warning('You have to download AGMG toolbox from "http://www.agmg.eu/" and place it inside  "/fun/AGMG" directory')
   disp(' ')
   warning('prec is set to amg_MATMAG')
   prec =  'amg_MATMAG';
end
cd ..
end
cd('data'); cd(name_dir); load('data.mat'); 
fileList = dir('*.stl');
figure
hold on
xmin=[];xmax=[];ymin=[];ymax=[];zmin=[];zmax=[];ccolor=distinguishable_colors(size(fileList,1));
for ii = 1:size(fileList,1)
    [stlcoords] = READ_stl(fileList(ii).name);
    xco = squeeze( stlcoords(:,1,:) )';
    yco = squeeze( stlcoords(:,2,:) )';
    zco = squeeze( stlcoords(:,3,:) )';
    [hpat] = patch(xco,yco,zco,ccolor(ii,:),'edgecolor','none');
    alpha(0.1) % for trasparency
    axis equal
    xlabel('x');
    ylabel('y');
    zlabel('z');
    view(3)
    title('stl (original, not scaled)')
    drawnow
end
cd(dad)
modelname = name_dir;
%% EM constants
mu = 4*pi*1e-7;
co = 299792458;
eo = 1/co^2/mu;
omega = 2*pi*[freq_min freq_max];
%% extract data information 
nVoxesl=L*M*N;
rhoVoxel=zeros(nVoxel,1);
idxV=[]; rhomin=Inf; hhi=1; hho=0;
for ii = 1:Nmat
    Ind(ii).ind=reshape(Ind(ii).ind,length(Ind(ii).ind),1);
    if strcmp(Ind(ii).tag,'air') || strcmp(Ind(ii).tag,'mag') || strcmp(Ind(ii).tag,'diel') || strcmp(Ind(ii).tag,'pot')
        % nothing to do here (?)
    elseif strcmp(Ind(ii).tag,'cond')
        idxV=[idxV;Ind(ii).ind];  
        rhoVoxel(Ind(ii).ind,1)=Ind(ii).rho;  
        rhomin=min([rhomin,Ind(ii).rho]);
    end
end
idxV=unique(idxV);
idxVR=idxV;
%% Grid Definition
disp('----DOMAIN--------------------------------')
%%Grid resolution
disp([' Number of voxels in x direction: ', num2str(L)])
disp([' Number of voxels in y direction: ', num2str(M)])
disp([' Number of voxels in z direction: ', num2str(N)])
disp(' Resolution:')
dx = smeshx; dy = smeshy; dz = smeshz;
disp([' dx = ',num2str(dx),' m']); disp([' dy = ',num2str(dy),' m']); disp([' dz = ',num2str(dz),' m'])
d = [dx dy dz]; 
Kt = nVoxel; %total number of voxels
K = length(idxV); %number of non-empty voxels
%% Set Material Properties
rho_eV=reshape(rhoVoxel,L,M,N); % 
clear rhoVoxel
%%
disp([' Total number of voxels: ', num2str(Kt)])
disp([' Number of non-empty voxels: ', num2str(K)])
disp(' ')
%% Incidence Matix A
disp('----COMPUTING INCIDENCE--------------------------------')
mytic=tic;
[Ae,Aee,AeeR,idxF,idxFx,idxFy,...
    idxFz,Ae1x,Ae1y,Ae1z,indPotv,dim_dom,id_dom,N_dom] = ...
    incidence_matrix2_AutoSelDom(Kt,[L M N],idxV);
disp([' Number of DoFs: ', num2str(size(AeeR,1)+size(AeeR,2))])
disp([' Time for computing incidence ::: ' ,num2str(toc(mytic))]);
disp(' ')
%% Forcing Term: Incident E field
% NOTE: Electric field with components (-iwy/2,iwx/2,0)
% since component x (y) does not depend on x(y), field is calculated at voxel
% barycenter; but in general must be calculated in barycenters of faces!!!
% Thus, we are introducing an approximation here. 
if ext_field_flag
    Ax = Ax_ext(xyz(:,:,:,1),xyz(:,:,:,2),xyz(:,:,:,3));
    Ay = Ay_ext(xyz(:,:,:,1),xyz(:,:,:,2),xyz(:,:,:,3));
    Az = Az_ext(xyz(:,:,:,1),xyz(:,:,:,2),xyz(:,:,:,3));
    %%RHS array: <Einc,f_a>_V = int_V dot(Einc,f_a) dV
    Gram = dx*dy*dz; %volume of cubic element
    Vx = (Gram.*Ax)./(dy*dz);
    Vy = (Gram.*Ay)./(dx*dz);
    Vz = (Gram.*Az)./(dx*dy);
    clear Ax Ay Az
else
    Vx=zeros(L*M*N,1);
    Vy=zeros(L*M*N,1);
    Vz=zeros(L*M*N,1);
end
%% Matrices Z_real and Z_imag
rho_eF=0.5*(abs(Ae(:,:)).'*rho_eV(:)); clear Ae rho_eV
z_realF=rho_eF;
indFneq=setdiff([1:3*Kt].',[idxFx;idxFy+Kt;idxFz+2*Kt]);
z_realF(indFneq,:)=0;z_realx=zeros(L,M,N);z_realx(idxFx)=z_realF(idxFx);
z_realx_loc = z_realx(idxFx);clear z_realx
z_realy=zeros(L,M,N);z_realy(idxFy)=z_realF(Kt+idxFy);
z_realy_loc = z_realy(idxFy);clear z_realy
z_realz=zeros(L,M,N);z_realz(idxFz)=z_realF(2*Kt+idxFz);
z_realz_loc = z_realz(idxFz);clear z_realz
%% Compute Green Tensor 
disp('----COMPUTING GREEN TENSOR--------------------------------')
mytic_G=tic;
[Gmn] = computeGREEN(d,L,M,N,Integration_flag);
disp([' Total time for getting Green tensor ::: ' ,num2str(toc(mytic_G))]);
disp(' ')
%% Compute Circulant Tensors
disp('----COMPUTING CIRCULANT TENSOR--------------------------------')
disp(' Circulant Tensors related to L matrix')
mytic_cir=tic;
[opCirculantM_all,st_sparse_preconM] = computeCIRCULANT(Gmn,d,'L');
%%Add constants to Circulants
opCirculantM_all = (mu)*opCirculantM_all;
st_sparse_preconM = (1j*mu)*st_sparse_preconM;
disp([' Total time for getting circulant tensors ::: ' ,num2str(toc(mytic_cir))])
clear Gmn %Green tensor is not used anymore
disp(' ')
%% Generating RHS vector
num_node = size(Aee,1); %all potential nodes in non-empty voxels 
num_nodeR = size(AeeR,1); %all potential nodes in non-empty voxels excluding ones with given potential
num_curr = size(Aee,2); %all currens in non-empty voxels 
%%Define RHS: current + potentials
rhs_vect = [Vx(idxFx);Vy(idxFy);Vz(idxFz);zeros(num_node,1)]; 
clear Vx Vy Vz
%%Reduced RHS: 
rhs_vectR = rhs_vect;
rhs_vectR(num_curr+indPotv,:) = [];
clear Vx Vy Vz
%%
disp('----GENERATING MOR-------------------------------------------------')
[yR,yL,V,b_hat]=fun_mor_freq(mor.N_check,mor.MOR_tol,mor.MOR_maxdim,...
                        omega(1),omega(2),rhs_vectR,...
                        d,z_realx_loc,z_realy_loc,z_realz_loc,...
                        idxFx,idxFy,idxFz,AeeR,Aee,Kt,inner_it,tol,...
                        outer_it,opCirculantM_all,st_sparse_preconM,...
                        prec,L,M,N,z_realF,opt_amg_MATMAG,opt_amg_AGMG);
disp(' ')                    
%% use mor
disp('----USE ROM--------------------------------')
disp(' ')
disp(['ROM dimension: ',num2str(length(b_hat))])
wpoint=0.5*(omega(2)+omega(1))/2; % select the desired frequency
tic
A_hat = yR(:, :)+1j*wpoint*yL(:, :); % create MOR for the desired frequency 
x_hat = A_hat \ (wpoint*b_hat); % solving ROM for the desired frequency 
vsol=V*x_hat; % from ROM to FOM (solution)
disp([' Total time for solving ROM ::: ' ,num2str(toc(tic)),' s'])
disp(' ')
%% extract solution
Jout = zeros(L,M,N,3);
Jout(idxF) = vsol(1:num_curr) ; % return to global variables
%%
%% POST PROCESSING
%%
%% Post Processing J
disp('----POST PROCESSING J------------------------------')
mytic_prec=tic;
[J,XYZ] = fun_my_postRT2(Jout,Kt,Ae1x,Ae1y,Ae1z,xyz,L,M,N,d);
potval=zeros(Kt,1);
indLocPotDofs=setdiff(idxV,idxV(indPotv));
potval(indLocPotDofs,1)=vsol(num_curr+1:end);
disp([' Total time for post processing J ::: ' ,num2str(toc(mytic_prec))]);
disp(' ')
%% Plot Vectors
if plot_vectorsJ_flag
jjR = real(J);%reshape(real(Jout),L*M*N,3)/(l^2);
figure
subplot(1,2,1)
normJR=sqrt(jjR(:,1).^2+jjR(:,2).^2+jjR(:,3).^2);
quiver3_c_scal(XYZ(:,1),XYZ(:,2),XYZ(:,3),jjR(:,1),jjR(:,2),jjR(:,3),...
          normJR,4);
axis equal
c1=colorbar;
caxis([min(normJR) max(normJR)]);
xlabel('x')
ylabel('y')
zlabel('z')
title('Current Density Vector \Re Part')
c1.Location = 'southoutside';
xlim([min(XYZ(:,1))-dx max(XYZ(:,1))+dx])
ylim([min(XYZ(:,2))-dy max(XYZ(:,2))+dy])
zlim([min(XYZ(:,3))-dz max(XYZ(:,3))+dz])
%
jjI = imag(J); %reshape(imag(Jout),L*M*N,3)/(l^2);
subplot(1,2,2)
normJI=sqrt(jjI(:,1).^2+jjI(:,2).^2+jjI(:,3).^2);
quiver3_c_scal(XYZ(:,1),XYZ(:,2),XYZ(:,3),jjI(:,1),jjI(:,2),jjI(:,3),...
          normJI,4);
axis equal
c1=colorbar;
caxis([min(normJI) max(normJI)]);
xlabel('x')
ylabel('y')
zlabel('z')
title('Current Density Vector \Im Part')
c1.Location = 'southoutside';
xlim([min(XYZ(:,1))-dx max(XYZ(:,1))+dx])
ylim([min(XYZ(:,2))-dy max(XYZ(:,2))+dy])
zlim([min(XYZ(:,3))-dz max(XYZ(:,3))+dz])
end
%% paraview
if paraview_export_flag
disp('----EXPORT TO PARAVIEW------------------------------')
xd=xyz(:,:,:,1);
yd=xyz(:,:,:,2);
zd=xyz(:,:,:,3);
xidx=xd(idxV);
yidx=yd(idxV);
zidx=zd(idxV);
P0=[...
    [xidx-dx/2,yidx-dy/2,zidx+dz/2];...
    [xidx-dx/2,yidx-dy/2,zidx-dz/2];...
    [xidx+dx/2,yidx-dy/2,zidx-dz/2];...
    [xidx+dx/2,yidx-dy/2,zidx+dz/2];...
    [xidx-dx/2,yidx+dy/2,zidx+dz/2];...
    [xidx-dx/2,yidx+dy/2,zidx-dz/2];...
    [xidx+dx/2,yidx+dy/2,zidx-dz/2];...
    [xidx+dx/2,yidx+dy/2,zidx+dz/2]];
VP=[1:K;...
    K+1:2*K;...
    2*K+1:3*K;...
    3*K+1:4*K;...
    4*K+1:5*K;...
    5*K+1:6*K;...
    6*K+1:7*K;...
    7*K+1:8*K];
warning off
[~] = ...
    fun_for_ParaView_vec_HEXA(...
    jjR(idxV,:),jjI(idxV,:),P0,VP,dad,[modelname,'J']);
warning on
end
%% to nii
Jnorm3D=reshape(sqrt(J(:,1).*conj(J(:,1))+J(:,2).*conj(J(:,2))+J(:,3).*conj(J(:,3))),L,M,N);
Jnii=make_nii(Jnorm3D);
view_nii(Jnii);
