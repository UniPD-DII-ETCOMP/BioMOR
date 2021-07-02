clear 
clear global
close all
clc
global meshXmin meshXmax meshYmin meshYmax meshZmin meshZmax
%% BEGIN USER SETTINGS
paraview_export_flag = 1;
x_ray_flag = 1;
model_name='SheppLoganHeadModel';
%
stl_files(10).name = 'el1.stl'; 
stl_files(10).tag = 'cond';
stl_files(10).cur=[];
stl_files(10).rho=8.26e-2;
stl_files(10).rhomin=8.26e-2*0.8; % only used if pMOR on parameters is performed
stl_files(10).rhomax=8.26e-2*1.2; % only used if pMOR on parameters is performed
%
stl_files(9).name = 'el2.stl'; 
stl_files(9).tag = 'cond';
stl_files(9).cur=[];
stl_files(9).rho=1.15e-1;
stl_files(9).rhomin=1.15e-1*0.8;
stl_files(9).rhomax=1.15e-1*1.2;
%
stl_files(8).name = 'el3.stl'; 
stl_files(8).tag = 'cond';
stl_files(8).cur=[];
stl_files(8).rho=2;
stl_files(8).rhomin=2*0.8;
stl_files(8).rhomax=2*1.2;
%
stl_files(7).name = 'el4.stl'; 
stl_files(7).tag = 'cond';
stl_files(7).cur=[];
stl_files(7).rho=2;
stl_files(7).rhomin=2*0.8;
stl_files(7).rhomax=2*1.2;
%
stl_files(6).name = 'el5.stl'; 
stl_files(6).tag = 'cond';
stl_files(6).cur=[];
stl_files(6).rho=6.95e-2;
stl_files(6).rhomin=6.95e-2*0.8;
stl_files(6).rhomax=6.95e-2*1.2;
%
stl_files(5).name = 'el6.stl'; 
stl_files(5).tag = 'cond';
stl_files(5).cur=[];
stl_files(5).rho=3e-1;
stl_files(5).rhomin=3e-1*0.8;
stl_files(5).rhomax=3e-1*1.2;
%
stl_files(4).name = 'el7.stl'; 
stl_files(4).tag = 'cond';
stl_files(4).cur=[];
stl_files(4).rho=3e-1;
stl_files(4).rhomin=3e-1*0.8;
stl_files(4).rhomax=3e-1*1.2;
%
stl_files(3).name = 'el8.stl'; 
stl_files(3).tag = 'cond';
stl_files(3).cur=[];
stl_files(3).rho=3e-1;
stl_files(3).rhomin=3e-1*0.8;
stl_files(3).rhomax=3e-1*1.2;
%
stl_files(2).name = 'el9.stl'; 
stl_files(2).tag = 'cond';
stl_files(2).cur=[];
stl_files(2).rho=3e-1;
stl_files(2).rhomin=3e-1*0.8;
stl_files(2).rhomax=3e-1*1.2;
%
stl_files(1).name = 'el10.stl'; 
stl_files(1).tag = 'cond';
stl_files(1).cur=[];
stl_files(1).rho=3e-1;
stl_files(1).rhomin=3e-1*0.8;
stl_files(1).rhomax=3e-1*1.2;
% to scale a stl file from any unit to meters
scal_geomery.x=1; scal_geomery.y=1; scal_geomery.z=1;
% Box 
% number of voxels in the x y z directions
Nx=40;%110;
Ny=50;%150;
Nz=45;%148;
% corners
flag_auto=1; % if 1, user_data below are ignored
% user_data
meshXmin = -0.06; % (m)
meshXmax = +0.06;  % (m)
meshYmin = -0.06;% (m)
meshYmax = +0.06;   % (m)
meshZmin = -0.06;% (m)
meshZmax = +0.06;  % (m)
%% END USER SETTINGS
how_many_stl=size(stl_files,2);
%%
dad=pwd;
cd ..; cd ..; cd('fun'); addpath(genpath(pwd)); cd(dad)
%% Plot the original STL mesh
ccolor=distinguishable_colors(how_many_stl);
figure
hold on
xmin=[];
xmax=[];
ymin=[];
ymax=[];
zmin=[];
zmax=[];
for ii = 1:how_many_stl
[stlcoords] = READ_stl(stl_files(ii).name);
xco = squeeze( stlcoords(:,1,:) )';
yco = squeeze( stlcoords(:,2,:) )';
zco = squeeze( stlcoords(:,3,:) )';
[hpat] = patch(xco*scal_geomery.x,yco*scal_geomery.y,zco*scal_geomery.z,ccolor(ii,:),'edgecolor','none');
alpha(0.1) % for trasparency
axis equal
xlabel('x');
ylabel('y');
zlabel('z');
view(3)
title('stl')
drawnow
pause(1)
% 
xmin=min([xmin;xco(:)]);
xmax=max([xmax;xco(:)]);
ymin=min([ymin;yco(:)]);
ymax=max([ymax;yco(:)]);
zmin=min([zmin;zco(:)]);
zmax=max([zmax;zco(:)]);
end
%%
if flag_auto
    meshXmin=xmin;
    meshXmax=xmax;
    meshYmin=ymin;
    meshYmax=ymax;
    meshZmin=zmin;
    meshZmax=zmax;
else
    meshXmin = meshXmin/scal_geomery.x;
    meshXmax = meshXmax/scal_geomery.x; 
    meshYmin = meshYmin/scal_geomery.y; 
    meshYmax = meshYmax/scal_geomery.y;
    meshZmin = meshZmin/scal_geomery.z;
    meshZmax = meshZmax/scal_geomery.z;    
end
%% 
disp('===================================================================')
disp('voxelize...')
for ii = 1:how_many_stl
    [o(ii).OUTPUTgrid,...
    o(ii).gridCOx,...
    o(ii).gridCOy,...
    o(ii).gridCOz] = VOXELISE_mod(Nx,Ny,Nz,stl_files(ii).name,'xyz');
    o(ii).gridCOx=o(ii).gridCOx*scal_geomery.x;
    o(ii).gridCOy=o(ii).gridCOy*scal_geomery.y;
    o(ii).gridCOz=o(ii).gridCOz*scal_geomery.z;
end
xyz= grid3dRT2(o(1).gridCOx,o(1).gridCOy,o(1).gridCOz);
xd=xyz(:,:,:,1);
yd=xyz(:,:,:,2);
zd=xyz(:,:,:,3);
%
idx=zeros(Nx,Ny,Nz);
for ii = 1:how_many_stl
idx=idx+(o(ii).OUTPUTgrid);
end
idx=find(idx);
% plot
xidx=xd(idx);
yidx=yd(idx);
zidx=zd(idx);
disp(' ')
%%
disp('===================================================================')
disp('x_ray...')
if x_ray_flag
    for ii = 1:how_many_stl
        figure
        subplot(1,3,1);
        title(['xray object',num2str(ii),' ZY'])
        hold on
        imagesc(squeeze(sum(o(ii).OUTPUTgrid,1)));
        colormap(gray(256));
        xlabel('Z-direction');
        ylabel('Y-direction');
        axis equal 
        subplot(1,3,2);
        title(['xray object',num2str(ii),' ZX'])    
        hold on
        imagesc(squeeze(sum(o(ii).OUTPUTgrid,2)));
        colormap(gray(256));
        xlabel('Z-direction');
        ylabel('X-direction');
        axis equal 
        subplot(1,3,3);
        title(['xray object',num2str(ii),' YX'])            
        hold on        
        imagesc(squeeze(sum(o(ii).OUTPUTgrid,3)));
        colormap(gray(256));
        xlabel('Y-direction');
        ylabel('X-direction');
        axis equal 
    end
end
drawnow
%%
if Nx>1
dx=abs(o(1).gridCOx(2)-o(1).gridCOx(1));
else
dx=(max(xco(:))-min(xco(:)))*scal_geomery.x;    
end
if Ny>1
dy=abs(o(1).gridCOy(2)-o(1).gridCOy(1));
else
dy=(max(yco(:))-min(yco(:)))*scal_geomery.y;    
end
if Nz>1
dz=abs(o(1).gridCOz(2)-o(1).gridCOz(1));
else
dz=(max(zco(:))-min(zco(:)))*scal_geomery.z;
end
xidx=xd([idx]);
yidx=yd([idx]);
zidx=zd([idx]);
disp(' ')
%%  paraview
disp('===================================================================')
disp('paraview...')
if paraview_export_flag
P0=[...
    [xidx-dx/2,yidx-dy/2,zidx+dz/2];...
    [xidx-dx/2,yidx-dy/2,zidx-dz/2];...
    [xidx+dx/2,yidx-dy/2,zidx-dz/2];...
    [xidx+dx/2,yidx-dy/2,zidx+dz/2];...
    [xidx-dx/2,yidx+dy/2,zidx+dz/2];...
    [xidx-dx/2,yidx+dy/2,zidx-dz/2];...
    [xidx+dx/2,yidx+dy/2,zidx-dz/2];...
    [xidx+dx/2,yidx+dy/2,zidx+dz/2]];
nVox=length([idx]);
VP=[1:nVox;...
    nVox+1:2*nVox;...
    2*nVox+1:3*nVox;...
    3*nVox+1:4*nVox;...
    4*nVox+1:5*nVox;...
    5*nVox+1:6*nVox;...
    6*nVox+1:7*nVox;...
    7*nVox+1:8*nVox];
val=zeros(length(o(1).gridCOx)*length(o(1).gridCOy)*length(o(1).gridCOz),1);
for ii = 1:how_many_stl
    val(find(o(ii).OUTPUTgrid))=stl_files(ii).rho;
end
val3D=reshape(val,length(o(1).gridCOx),length(o(1).gridCOy),length(o(1).gridCOz));
val=val(idx);
dad=pwd;
[Ricc] = fun_for_ParaView_sca_HEXA(...
    val,val,P0,VP,dad,model_name);
end
%%
L=length(o(1).gridCOx);
M=length(o(1).gridCOy);
N=length(o(1).gridCOz);
nVoxel=L*M*N;
smeshx=dx;smeshy=dy;smeshz=dz;
%%  
for ii = 1:how_many_stl
    Ind(ii).ind= find(o(ii).OUTPUTgrid);
    Ind(ii).tag=stl_files(ii).tag;
    Ind(ii).rho=stl_files(ii).rho;
    Ind(ii).rhomin=stl_files(ii).rhomin;
    Ind(ii).rhomax=stl_files(ii).rhomax;
end
%%
Nmat = how_many_stl;
save data.mat Ind L M N Nmat smeshx smeshy smeshz xyz -v7.3
%% to nii
RES=make_nii(val3D);
view_nii(RES);
