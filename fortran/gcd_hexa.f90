subroutine gcd_hexa(vp,gshort,c,d,nn,ne,nf,nv,f)
use commonstuff_hexa
implicit none
integer*8 vp(8,nf)
integer*8 nn,ne,nf,nv
integer*8,allocatable,dimension(:,:) :: gshort,c,d
integer*8 i,j,esiz,fsiz,ehtpsiz,elocsiz,fhtpsiz,flocsiz,quattro
integer*8 nnd2(2),nnd2ori(2),nnd4(4),nnd4rot(4),ind,htsrcedge,htsrcfacehexa,ied
double precision log2
integer*8,allocatable,dimension(:) :: eht,elinks,fht,flinks
integer*8,allocatable,dimension(:,:) :: f,g
external htsrcedge,htsrcfacehexa
robust=.true.
quattro = int(4,8)
nn=maxval(vp(1:8,1:nv))
esiz=max(nn*10,10000)
ehtpsiz=int(log2(esiz))+1
elocsiz=2**ehtpsiz
allocate(eht(elocsiz),elinks(elocsiz))
allocate(g(2,elocsiz))
!G matrix
ne=0
eht(:)=0
do i=1,nv  
  do j=1,12
     nnd2(:)=vp(locedg_hexa(j,:),i)
     nnd2ori(:)=nnd2(:)
     call ordsw2(nnd2(1),nnd2(2))
     ind=htsrcedge(nnd2,ehtpsiz,g,eht,elinks)
     ied=ind
     if(ind.eq.0)then
       ne=ne+1
       if(ne.gt.elocsiz)then
         write(6,*) 'ne > locsiz in gcd, aborting... nv,locsiz=',nv,elocsiz
         stop
       endif
       ied=ne
       call htinsedge (ne,nnd2,ehtpsiz,g,eht,elinks )
!recheck
       if(robust)then
         ind=htsrcedge (nnd2,ehtpsiz,g,eht,elinks )
         if(ind.ne.ne)then
           write(6,*) 'impossible error in gcd_quad'
           stop
         endif
       endif      
     endif
  enddo
enddo
!F and D matrices
fsiz=max(nv*6,10000)
fhtpsiz=int(log2(fsiz))+1
flocsiz=2**fhtpsiz
allocate(fht(flocsiz),flinks(flocsiz))
allocate(f(8,flocsiz)) !4 nodes defining face, 2 attached elements, 2 local faces
allocate(d(6,nv))
nf=0
fht(:)=0
do i=1,nv
  do j=1,6
    nnd4(1:4)=vp(locfac_hexa(j,1:4),i) ! prendo i 4 nodi della faccia locale j 
    ! nnd4rot(1:3)=nnd4(1:3)
    ! call rotate3(nnd3rot(1),nnd3rot(2),nnd3rot(3))
    call orderk(quattro,nnd4) ! faccio un sort dei nodi della faccia 
	ind=htsrcfacehexa ( nnd4, fhtpsiz,f, fht, flinks )
	if(ind.eq.0)then
      nf=nf+1
       if(nf.gt.flocsiz)then
         write(6,*) 'nf > locsiz in vol2fac, aborting... nv,locsiz=',nv,flocsiz
         stop
       endif
!make sure that inside facs we have real existing faces with smallest index first
      call htinsfacehexa ( nf, nnd4, fhtpsiz,f, fht, flinks )
      f(1:4,nf)=vp(locfac_hexa(j,:),i)
	  f(5,nf)=i
	  f(7,nf)=j
      f(6,nf)=0
      f(8,nf)=0
      !regions(1,nf)=reg(i)
      d(j,i)=nf
    else
      !regions(2,ind)=reg(i)
      d(j,i)=-ind
	  f(6,ind)=i
      f(8,ind)=j
    endif
  enddo
enddo
!what follows is a consistency check
if(robust)then
do i=1,nv
  do j=1,6
    nnd4(1:4)=vp(locfac_hexa(j,1:4),i)
    !nnd3rot(1:3)=nnd3(1:3)
    !call rotate3(nnd3rot(1),nnd3rot(2),nnd3rot(3))
    call orderk(quattro, nnd4)
    ind=htsrcfacehexa ( nnd4, fhtpsiz,f, fht, flinks )
    if(ind.eq.0)then
	  write(66,*) 'nnd4'
	  write(66,*) nnd4
	  write(66,*) 'ind'
	  write(66,*) ind
      write(66,*) 'Impossible error in gcd...'
      stop
    endif
enddo  
enddo
endif
!C matrix
allocate(c(4,nf))
do i=1,nf
  do j=1,4
     nnd2(1:2)=f(locfacedg_hexa(j,1:2),i)
     nnd2ori(:)=nnd2(:)
     call ordsw2(nnd2(1),nnd2(2))
     ind=htsrcedge ( nnd2, ehtpsiz,g, eht, elinks )
     ied=ind
     if(ind.eq.0)then
       write(6,*) 'Immpossible error XXX in gcd'
       stop
     endif
     if(nnd2ori(1).eq.nnd2(1).and.nnd2ori(2).eq.nnd2(2))then
       c(j,i)=ied
     else
       c(j,i)=-ied
     endif
      enddo
enddo
deallocate(elinks,eht)
deallocate(fht,flinks)
allocate(gshort(2,ne))
gshort(1:2,1:ne)=g(1:2,1:ne)
deallocate(g)
end subroutine gcd_hexa