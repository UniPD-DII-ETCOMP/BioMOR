function [C] = fun_make_C(L,M,N,K,idxV,Kt,idxF)
VP=zeros(8,K);
Lp=L+1;
Mp=M+1;
Np=N+1;
for ii = 1:K
    lmn=idx2triplet(idxV(ii),L,M,N);
    l=lmn(1);
    m=lmn(2);
    n=lmn(3);
    VP(1,ii)=l  +Lp*(m  -1)+Lp*Mp*(n+1-1);
    VP(2,ii)=l  +Lp*(m  -1)+Lp*Mp*(n-1  );
    VP(3,ii)=l+1+Lp*(m  -1)+Lp*Mp*(n-1  );
    VP(4,ii)=l+1+Lp*(m  -1)+Lp*Mp*(n+1-1);
    VP(5,ii)=l  +Lp*(m+1-1)+Lp*Mp*(n+1-1);
    VP(6,ii)=l  +Lp*(m+1-1)+Lp*Mp*(n-1  );
    VP(7,ii)=l+1+Lp*(m+1-1)+Lp*Mp*(n-1  );
    VP(8,ii)=l+1+Lp*(m+1-1)+Lp*Mp*(n+1-1);
end
[G1,C1,D1,F1]=gcd_mexed_HEXA(VP);
nf=size(C1,2);
ne=size(G1,2);
Cpos = (C1+abs(C1))/2;
Cneg = (C1-abs(C1))/2;
[r_pos_c,c_pos_c,val_pos_c] = find(Cpos);
[r_neg_c,c_neg_c,val_neg_c] = find(Cneg);
Matrix_C = sparse(c_pos_c,val_pos_c,ones(size(c_pos_c,1),1),nf,ne);
Matrix_C = Matrix_C+sparse(c_neg_c,abs(val_neg_c),-ones(size(c_neg_c,1),1),nf,ne);
ind_face_ext=find(F1(6,:)==0).';
n_face_ext=length(ind_face_ext);
ind_face_int=find(F1(6,:)~=0).';
corr_Fgcd_Fvoxel=zeros(length(ind_face_int),1);
corr_Fgcd_Fvoxel_flagXYZ=zeros(length(ind_face_int),1);
for ii = 1:length(ind_face_int)
    v1=F1(5,ind_face_int(ii));
    v2=F1(6,ind_face_int(ii));
    lmn1=idx2triplet(idxV(v1),L,M,N);
    lmn2=idx2triplet(idxV(v2),L,M,N);
    if lmn2(1)-lmn1(1)==1
        corr_Fgcd_Fvoxel(ii)=     lmn1(1)+L*(lmn1(2)-1)+L*M*(lmn1(3)-1);
        corr_Fgcd_Fvoxel_flagXYZ(ii)=1;
    elseif lmn2(2)-lmn1(2)==1
        corr_Fgcd_Fvoxel(ii)=Kt  +lmn1(1)+L*(lmn1(2)-1)+L*M*(lmn1(3)-1);
        corr_Fgcd_Fvoxel_flagXYZ(ii)=2;
    elseif lmn2(3)-lmn1(3)==1
        corr_Fgcd_Fvoxel(ii)=2*Kt+lmn1(1)+L*(lmn1(2)-1)+L*M*(lmn1(3)-1);
        corr_Fgcd_Fvoxel_flagXYZ(ii)=3;        
    end
end
%
ind_edge_ext=unique(abs(reshape(C1(:,ind_face_ext),1,4*n_face_ext))).';
ind_edge_int=setdiff(1:ne,ind_edge_ext).';                                  
[~,IA,IB]=intersect(idxF,corr_Fgcd_Fvoxel,'stable');
flag_edgeXYZ=zeros(length(ind_edge_int),1);
for ii = 1:length(ind_edge_int)
    p1=G1(1,ind_edge_int(ii));
    p2=G1(2,ind_edge_int(ii));
    lmn1=idx2triplet((p1),L+1,M+1,N+1);
    lmn2=idx2triplet((p2),L+1,M+1,N+1);
    if lmn2(1)-lmn1(1)==1
        flag_edgeXYZ(ii)=1;
    elseif lmn2(2)-lmn1(2)==1
        flag_edgeXYZ(ii)=2;
    elseif lmn2(3)-lmn1(3)==1
        flag_edgeXYZ(ii)=3;        
    end    
end
loc_edge_x=(find(flag_edgeXYZ==1));
loc_edge_y=(find(flag_edgeXYZ==2));
loc_edge_z=(find(flag_edgeXYZ==3));
Matrix_C=Matrix_C(ind_face_int,ind_edge_int);
C=Matrix_C(IB,:);
C=C(:,[loc_edge_x;loc_edge_y;loc_edge_z]);
end

