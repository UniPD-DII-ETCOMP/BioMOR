double precision function log2(n)
implicit none
integer*8 n
log2=log(dble(n))/log(dble(2))
end function

subroutine ordsw2(aa,bb)
implicit none
integer*8 aa,bb,kk1
  if(aa.gt.bb)then
    kk1=aa
    aa=bb
    bb=kk1
  endif
end subroutine

subroutine ordsw3(aa,bb,cc) ! sort del vettore
implicit none
integer*8 aa,bb,cc,kk1,kk2
  if(aa.gt.bb)then
    kk1=aa
    aa=bb
    bb=kk1
  endif
  if(cc.ge.bb)then
    return
  endif
  if(cc.le.aa)then
     kk1=aa
     kk2=bb
     aa=cc
     bb=kk1
     cc=kk2
     return
  endif
  kk1=bb
  bb=cc
  cc=kk1
end subroutine

subroutine rotate3(aa,bb,cc)
implicit none
! rotate aa,bb,cc so that smallest node is in position 1
integer*8 imin,aa,bb,cc,kk1,kk2,kk3
imin=min(aa,bb,cc)
if(aa.eq.imin) return
if(bb.eq.imin)then
  kk1=aa
  kk2=bb
  kk3=cc
  aa=kk2
  bb=kk3
  cc=kk1
  return
endif
  kk1=aa
  kk2=bb
  kk3=cc
  aa=kk3
  bb=kk1
  cc=kk2
end subroutine


integer*8 function fhash(arr,siz,seed,cut)
implicit none
integer*8 arr(*),siz,seed,cut
integer*8 i,shft
!FNV hash
fhash = 1966136261
do i=1,siz
   fhash=ieor(fhash*16777619,arr(i))
enddo
shft=ishft(1,cut)-1
fhash=iand(fhash,shft)+1
!write(6,*) 'fhash:',arr(1:siz),fhash,cut
end function

subroutine htinsedge ( ind, ab, htedgpowersiz,edg, edght, edglinks )
  implicit none
  integer*8 ab(2)
  integer*8 ind,htedgpowersiz
  integer*8 k
  integer*8 leni
  integer*8 edg(2,*),edght(*),edglinks(*)
  integer*8 fhash
  external fhash
  leni=2
  k=fhash(ab, leni, htedgpowersiz, htedgpowersiz)
  edg(1:2,ind)=ab(1:2)
  edglinks(ind) = edght(k)
  edght(k) = ind
!  write(6,*) 'ins:',ind,ab,k,htedgpowersiz
end

function htsrcedge ( ab, htedgpowersiz,edg, edght, edglinks )
  implicit none
  integer*8 ab(2),htedgpowersiz
  integer*8 htsrcedge
  integer*8 ind
  integer*8 k
  integer*8 leni
  integer*8 edg(2,*),edght(*),edglinks(*)
  integer*8 fhash
  external fhash
  leni=2
  k=fhash(ab, leni, htedgpowersiz, htedgpowersiz)
!  write(6,*) 'src:',ab,k,htedgpowersiz
  ind = edght(k)
  do
    if ( ind == 0 ) then
      exit
    end if
    if ( edg(1,ind) == ab(1) .and. edg(2,ind) == ab(2) )then
            exit
    endif
    ind = edglinks(ind)
  end do
  htsrcedge = ind
end

function htsrcface ( abc, htfacpowersiz,fc, fcht, fclinks )
  implicit none
  integer*8 abc(3),htfacpowersiz
  integer*8 htsrcface
  integer*8 ind
  integer*8 k
  integer*8 cod(3),leni
  integer*8 fc(5,*),fcht(*),fclinks(*)
  integer*8 abcloc(3)
  integer*8 fhash
  external fhash
  leni=3
  k=fhash(abc, leni, htfacpowersiz, htfacpowersiz)
!  write(6,*) 'bef k=',k,htfacpowersiz
  ind = fcht(k)
!  write(6,*) 'aft ind=',ind
  do
    if ( ind == 0 ) then
      exit
    end if
    abcloc(1:3)=fc(1:3,ind)
    call ordsw3 ( abcloc(1), abcloc(2),abcloc(3) )    
    if ( abcloc(1) == abc(1) .and. abcloc(2) == abc(2) .and. abcloc(3) == abc(3))then
            exit
    endif
    ind = fclinks(ind)
  end do
  htsrcface = ind
end

function htsrcfacehexa ( abcd, htfacpowersiz,fc, fcht, fclinks )
  implicit none
  integer*8 abcd(4),htfacpowersiz
  integer*8 htsrcfacehexa
  integer*8 ind
  integer*8 k
  integer*8 cod(4),leni
  integer*8 fc(8,*),fcht(*),fclinks(*)
  integer*8 abcdloc(4)
  integer*8 fhash
  external fhash
  leni=4
  k=fhash(abcd, leni, htfacpowersiz, htfacpowersiz)
!  write(6,*) 'bef k=',k,htfacpowersiz
  ind = fcht(k)
!  write(6,*) 'aft ind=',ind
  do
    if ( ind == 0 ) then
      exit
    end if
    abcdloc(1:4)=fc(1:4,ind)
    call orderk ( 4, abcdloc )
    if ( abcdloc(1) == abcd(1) .and. abcdloc(2) == abcd(2) .and. abcdloc(3) == abcd(3) .and. abcdloc(4) == abcd(4))then
            exit
    endif
    ind = fclinks(ind)
  end do
  htsrcfacehexa = ind
end

subroutine htinsface ( ind, abc, htfacpowersiz,fc, fcht, fclinks )
  implicit none
  integer*8 abc(3)
  integer*8 ind,htfacpowersiz
  integer*8 k
  integer*8 leni
  integer*8 cod(3)
  integer*8 fc(5,*),fcht(*),fclinks(*)
  integer*8 fhash
  external fhash
  leni=3
  k=fhash(abc, leni, htfacpowersiz, htfacpowersiz)
  fc(1:3,ind)=abc(1:3)
  fclinks(ind) = fcht(k)
  fcht(k) = ind
!  write(6,*) 'htinsfac:',ind,fc(1:3,ind)
end

subroutine htinsfacehexa ( ind, abcd, htfacpowersiz,fc, fcht, fclinks )
  implicit none
  integer*8 abcd(4)
  integer*8 ind,htfacpowersiz
  integer*8 k
  integer*8 leni
  integer*8 cod(4)
  integer*8 fc(8,*),fcht(*),fclinks(*)
  integer*8 fhash
  external fhash
  leni=4
  k=fhash(abcd, leni, htfacpowersiz, htfacpowersiz)
  fc(1:4,ind)=abcd(1:4)
  fclinks(ind) = fcht(k)
  fcht(k) = ind
!  write(6,*) 'htinsfac:',ind,fc(1:3,ind)
end


subroutine orderk ( k, ind )
!*****************************************************************************80
!
!! ORDERK reorders K elements of an array in nondecreasing order.
!
!  Discussion: 
!
!    This routine reorders K elements of array IND in nondecreasing order.
!
!    It is assumed that K is small, say <= 15, so that insertion sort
!    is used. If K is larger, a faster sort such as heapsort should
!    be used.
!
!  Author:
!
!    Barry Joe, 
!    Department of Computing Science, 
!    University of Alberta,
!    Edmonton, Alberta, Canada  T6G 2H1
!    Email: barry@cs.ualberta.ca
!
!  Parameters:
!
!    Input, integer K, the size of array IND.
!
!    Input/output, integer IND(1:K), an array, which is sorted on output.
!
  implicit none
  integer*8 k
  integer*8 i
  integer*8 ind(k)
  integer*8 j
  integer*8 s
  integer*8 t
  do i = 2, k
    t = ind(i)
    j = i
    do
      s = ind(j-1)
      if ( s <= t ) then
        exit
      end if
      ind(j) = s
      j = j - 1
      if ( j <= 1 ) then
        exit
      end if
    end do
    ind(j) = t
  end do
  return
end



