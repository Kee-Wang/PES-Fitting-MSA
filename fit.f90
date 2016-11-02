program fit
use basis
implicit none

  external dgelss
  integer, parameter:: dp=kind(0.d0)
  real(dp) :: rmse, mue, mse, a0
  integer :: data_size
  real(dp),allocatable::xyz(:,:,:),x(:)
  real(dp),allocatable::v(:),b(:),p(:)
  real(dp),allocatable::yij(:,:), A(:,:)
  real(dp),allocatable::coeff(:),v_out(:),s(:)
  real(dp) :: work(150000), dr(3)
  integer :: ncoeff, natm, ndis
  integer :: i, j, k, m, info, rank
  character(len=32) :: data_file
  character :: symb

  natm=XXX         ! change to the number of atoms
  ncoeff=YYY       ! change to the number of coeff. (size of c in basis.f90)
  data_size=ZZZ    ! change to the number of data points in pts.dat
  a0=2.0_dp

  ndis=natm*(natm-1)/2

  open(10,file='coeff.dat',status='unknown')
  open(11,FILE='points.eng',status='unknown')
  open(12,file='points.dat',status='old')
  
  allocate(x(ndis))
  allocate(xyz(data_size,natm,3))
  allocate(v(data_size),v_out(data_size),b(data_size),coeff(ncoeff),s(ncoeff))
  
  do i=1,data_size
     read(12,*)
     read(12,*) v(i)
     do j=1,natm
        read(12,*) symb,xyz(i,j,:)
     end do
  end do
    
  allocate(yij(data_size,ndis))

  do m=1,data_size
     k = 1
     do i=1,natm-1
        do j=i+1,natm

           yij(m,k)=0.0_dp
           dr=xyz(m,i,:)-xyz(m,j,:)
           yij(m,k)=sqrt(dot_product(dr,dr))
           yij(m,k)=yij(m,k)/0.5291772083_dp
           yij(m,k)=exp(-yij(m,k)/a0)

           k=k+1
        end do
     end do
  end do
  
  deallocate(xyz)
  allocate(p(ncoeff))
  allocate(A(data_size,ncoeff))
  
  do m=1,data_size
     x=yij(m,:)
     call bemsav(x,p)
     A(m,:)=p
  end do
  b=v
  
  call dgelss(data_size,ncoeff,1,A,data_size,b,data_size,s,1.0d-8,rank,work,150000,info)
  
  coeff(:)=b(1:ncoeff)
  
  do i=1,ncoeff
     write(10,*) coeff(i)
  end do
  
  mse=0.0_dp
  rmse=0.0_dp
  mue=0.0_dp
  do i=1,data_size
     v_out(i)=emsav(yij(i,:),coeff)
     write (11,*) v(i),v_out(i),abs(v(i)-v_out(i))*219474.63
     mse=mse+abs(v(i)-v_out(i))
     rmse=rmse+(v(i)-v_out(i))**2
     mue=mue+sqrt((v(i)-v_out(i))**2)
  end do
  
  mse=mse/dble(data_size)
  rmse=sqrt(rmse/dble(data_size))
  mue=mue/dble(data_size)
  
  print*, 'MSE  = ', mse , ' Hartree'
  print*, 'RMSE = ', rmse , ' Hartree'
  print*, 'MUE  = ', mue , ' Hartree'
  
  close (10)
  close (11)
  close (12)

end program
