program test

implicit none
 
interface
 subroutine gcd_hexa(vp,g,c,d,nn,ne,nf,nv,f)
    integer*8 vp(8,*)
    integer*8 nn,ne,nf,nv
    integer*8,allocatable,dimension(:,:) :: g,c,d,f
  end subroutine  
end interface

integer*8 vp(8,2),i
    
integer*8 nn,ne,nf,nv
integer*8,allocatable,dimension(:,:) :: g,c,d,f

vp(1:8,1)=(/1,2,3,4,5,6,7,8/)
vp(1:8,2)=(/5,6,7,8,9,10,11,12/)
nv=2
call gcd_hexa(vp,g,c,d,nn,ne,nf,nv,f)

write(66,*) 'G'
!write(66,*) g
!write(66,*) 'ne'
!write(66,*) ne
do i=1,ne
  write(66,*) i,g(1:2,i)
enddo

write(66,*) 'D'
!write(66,*) d
!write(66,*) 'nv'
!write(66,*) nv
do i=1,nv
  write(66,*) i,d(1:6,i)
enddo

write(66,*) 'C'
!write(66,*) c
!write(66,*) 'nf'
!write(66,*) nf
do i=1,nf
  write(66,*) i,c(1:4,i)
enddo

write(66,*) 'F'
do i=1,nf
  write(66,*) i
  write(66,*) f(1:8,i)
enddo

if(allocated(g)) deallocate(g)
if(allocated(c)) deallocate(c)
if(allocated(d)) deallocate(d)

end program