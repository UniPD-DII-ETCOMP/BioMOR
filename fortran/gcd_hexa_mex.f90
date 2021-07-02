!---- never change ------

#include "fintrf.h"      

subroutine mexFunction(nlhs, plhs, nrhs, prhs)

   implicit none
   
interface
 subroutine gcd_hexa(vp,gshort,c,d,nn,ne,nf,nv,f)
    integer*8 vp(8,*)
    integer*8 nn,ne,nf,nv
    integer*8,allocatable,dimension(:,:) :: gshort,c,d,f
  end subroutine  
end interface

! mwPointer mexFunction arguments:
      mwPointer plhs(*), prhs(*)
      integer nlhs, nrhs
!---- end never change ------
	  
	  
! mwSize stuff for mexing
      mwSize mo,no,siz
	  
! mwPointer stuff for mexing
      mwPointer mxCreateDoubleMatrix
      mwPointer mxGetPr
	  mwPointer vp_pr,g_pr,c_pr,d_pr,nn_pr,ne_pr,nf_pr,nv_pr, f_pr
      mwPointer m, n
      mwPointer mxGetM, mxGetN

!integer (normal not 4/8)
      integer flag,i,j,ii,jj,icount
      integer mxIsNumeric 
      integer*4 ComplexFlag
	  
! fortran subroutine arguments
	  real*8,allocatable,dimension(:,:) :: vp_r
      integer*8,allocatable,dimension(:,:) :: vp,g,c,d,f
	  real*8,allocatable,dimension(:,:) :: greal,creal,ddreal,freal
	  integer*8 nn, ne, nf, nv
!	  real*8 
	  character*80 msg
      logical debu
       
      debu = .true. ! .true. o .false. per attivare o disattivare il debug
	  if(debu) open(unit=66,file='log.txt',status='unknown')

!    Check to see input is numeric.
!	  do ii = 1,1
        if (mxIsNumeric(prhs(1)) .ne. 1) then
          call mexErrMsgIdAndTxt ('MATLAB:gcd:NonNumeric', &
                                'Inputs must be numeric.')
        endif
!	  enddo
	  if(debu) write(66,*) 'sono arrivato al check numerico'
	  
	  
!     Check that input #1 is scalar matrix and fetch it
      m = mxGetM(prhs(1))
      n = mxGetN(prhs(1))
	  if(debu) write(66,*) m
	  if(debu) write(66,*) n
     ! if(m .ne. 8 .or. n .ne. N_hexa*3) then
      if(m .ne. 8) then !il check può essere fatto solo sulle righe (??)            TAKE A LOOK HERE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         call mexErrMsgIdAndTxt ('MATLAB:convec:NonRowVector', &
                                'Input 1 must 8 x NV')
      endif	  
	  nv=n
      siz = m*n
      vp_pr = mxGetPr(prhs(1))
	  write(66,*) 'qui1'
	  allocate(vp_r(8,n))
	  write(66,*) 'qui2'
      allocate(vp(8,n))
	  write(66,*) 'qui3'
      call mxCopyPtrToReal8(vp_pr, vp_r, siz) ! da double precision a reale
	  vp = int(vp_r,8) !da reale a intero
	  if(debu) write(66,*) 'input 1'
	  deallocate(vp_r)

	  	  

! call the computational subroutine.

       call gcd_hexa(vp,g,c,d,nn,ne,nf,nv,f)
       if(debu) write(66,*) 'ho chiamato la subroutine e creato le matrici'

	  deallocate(vp)
! Create a matrix for the return argument 1 (g)
      mo=2
      no=ne
	  ComplexFlag = 0
      plhs(1) = mxCreateDoubleMatrix(mo, no, ComplexFlag)
! Load the output 1 into a MATLAB array.
      g_pr = mxGetPr(plhs(1))
      siz=mo*no
	  
	  allocate(greal(2,ne))
	  greal=real(g,8)
	  deallocate(g)		  
      
	  call mxCopyReal8ToPtr(greal, g_pr, siz)
      if(debu) write(66,*) 'fatta G'
      deallocate(greal)	
	  
! Create a matrix for the return argument 2 (c)
      mo=4
      no=nf
	  ComplexFlag = 0
	  if(debu) write(66,*) 'prima di mxCreateDoubleMatrix'
      plhs(2) = mxCreateDoubleMatrix(mo, no, ComplexFlag)
! Load the output2 into a MATLAB array.
      c_pr = mxGetPr(plhs(2))
      siz=mo*no
	  
	  allocate(creal(4,nf))
	  creal=real(c,8)
	  deallocate(c)	
	  
	  if(debu) write(66,*) 'prima di mxCopyReal8ToPtr, sizes C and nf:',size(c,1),size(c,2),nf
      call mxCopyReal8ToPtr(creal, c_pr, siz)
      if(debu) write(66,*) 'fatta C'
	  deallocate(creal)
	  
! Create a matrix for the return argument 3 (d)
      mo=6
      no=nv
	  ComplexFlag = 0
      plhs(3) = mxCreateDoubleMatrix(mo, no, ComplexFlag)
! Load the output 3 into a MATLAB array.
      d_pr = mxGetPr(plhs(3))
      siz=mo*no
	  
	  allocate(ddreal(6,nv))
	  ddreal=real(d,8)
	  deallocate(d)		  
	  
      call mxCopyReal8ToPtr(ddreal, d_pr, siz)
      if(debu) write(66,*) 'fatta D'
	  deallocate(ddreal)
	  
! Create a matrix for the return argument 4 (f)
      mo=8
      no=nf
	  ComplexFlag = 0
	  if(debu) write(66,*) 'prima di mxCreateDoubleMatrix'
      plhs(4) = mxCreateDoubleMatrix(mo, no, ComplexFlag)
! Load the output2 into a MATLAB array.
      f_pr = mxGetPr(plhs(4))
      siz=mo*no
	  
	  allocate(freal(8,nf))
	  freal=real(f(1:8,1:nf),8)
	  deallocate(f)		  
	  
	  if(debu) write(66,*) 'prima di mxCopyReal8ToPtr, sizes F and nf:',size(f,1),size(f,2),nf
      call mxCopyReal8ToPtr(freal, f_pr, siz)
      if(debu) write(66,*) 'fatta F'	  
      deallocate(freal)

      if(debu) write(66,*) 'ho convertito le matrici per matlab'
      
      if(debu) write(66,*) 'ho deallocato la matrice e ora chiudo il log.txt'
	  
      
      if(debu) close(66)
      return
      end