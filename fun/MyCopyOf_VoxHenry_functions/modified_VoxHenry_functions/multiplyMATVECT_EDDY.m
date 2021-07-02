%%
% This function has been inspired by 
% lse_matvect_mult.m which is available at
% https://github.com/acyucel/VoxHenry
%%
function [sol] = multiplyMATVECT_EDDY(JIn0,fN,z_realx_loc,z_realy_loc,...
    z_realz_loc,idxFx,idxFy,idxFz,d,AeeR,L,M,N)
dx = d(1); dy = d(2); dz = d(3);
num_node=size(AeeR,1); %num potential nodes
num_curr=size(AeeR,2); %num current faces
sol=zeros(num_curr+num_node,1);
if dx == dy && dy == dz 
    indloc = 1:length(idxFx);
    fJ=zeros(L, M, N);
    fJ(idxFx) = JIn0(indloc);
    fJ = fftn(fJ,[2*L, 2*M, 2*N]);
    Jout = (fN(:,:,:,1)) .* fJ; % 
    clear fJ
    Jout = ifftn(Jout);
    Jout =Jout(1:L,1:M,1:N);
    sol(indloc)=Jout(idxFx);
    clear Jout
    sol(indloc) = (dx/(dy*dz)) * (z_realx_loc .* JIn0(indloc)) + 1j*sol(indloc);
  
    indloc = length(idxFx)+1:length(idxFx)+length(idxFy);
    fJ=zeros(L, M, N);
    fJ(idxFy) = JIn0(indloc);
    fJ = fftn(fJ,[2*L, 2*M, 2*N]);
    Jout = (fN(:,:,:,1)) .* fJ; % 
    clear fJ
    Jout = ifftn(Jout);
    Jout =Jout(1:L,1:M,1:N);
    sol(indloc)=Jout(idxFy);
    clear Jout
    sol(indloc) = (dy/(dx*dz)) * (z_realy_loc .* JIn0(indloc)) + 1j*sol(indloc);

    indloc = length(idxFx)+length(idxFy)+1:length(idxFx)+length(idxFy)+length(idxFz);
    fJ=zeros(L, M, N);
    fJ(idxFz) = JIn0(indloc);
    fJ = fftn(fJ,[2*L, 2*M, 2*N]);
    Jout = (fN(:,:,:,1)) .* fJ; % 
    clear fJ
    Jout = ifftn(Jout);
    Jout =Jout(1:L,1:M,1:N);
    sol(indloc)=Jout(idxFz);
    clear Jout
    sol(indloc) = (dz/(dx*dy)) * (z_realz_loc .* JIn0(indloc)) + 1j*sol(indloc);
else
    indloc = 1:length(idxFx);
    fJ=zeros(L, M, N);
    fJ(idxFx) = JIn0(indloc);
    fJ = fftn(fJ,[2*L, 2*M, 2*N]);
    Jout = (fN(:,:,:,1)) .* fJ; % 
    clear fJ
    Jout = ifftn(Jout);
    Jout =Jout(1:L,1:M,1:N);
    sol(indloc)=Jout(idxFx);
    clear Jout
    sol(indloc) = (dx/(dy*dz)) * (z_realx_loc .* JIn0(indloc)) + 1j*sol(indloc);
  
    indloc = length(idxFx)+1:length(idxFx)+length(idxFy);
    fJ=zeros(L, M, N);
    fJ(idxFy) = JIn0(indloc);
    fJ = fftn(fJ,[2*L, 2*M, 2*N]);
    Jout = (fN(:,:,:,2)) .* fJ; % 
    clear fJ
    Jout = ifftn(Jout);
    Jout =Jout(1:L,1:M,1:N);
    sol(indloc)=Jout(idxFy);
    clear Jout
    sol(indloc) = (dy/(dx*dz)) * (z_realy_loc .* JIn0(indloc)) + 1j*sol(indloc);

    indloc = length(idxFx)+length(idxFy)+1:length(idxFx)+length(idxFy)+length(idxFz);
    fJ=zeros(L, M, N);
    fJ(idxFz) = JIn0(indloc);
    fJ = fftn(fJ,[2*L, 2*M, 2*N]);
    Jout = (fN(:,:,:,3)) .* fJ; % Gxx*Jx
    clear fJ
    Jout = ifftn(Jout);
    Jout =Jout(1:L,1:M,1:N);
    sol(indloc)=Jout(idxFz);
    clear Jout
    sol(indloc) = (dz/(dx*dy)) * (z_realz_loc .* JIn0(indloc)) + 1j*sol(indloc);    
end
sol(1:num_curr) = sol(1:num_curr) + (AeeR.'*JIn0(num_curr+1:num_curr+num_node)) ;
sol(num_curr+1:num_curr+num_node) = AeeR*JIn0(1:num_curr);
end

