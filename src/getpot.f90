program main
use pes_shell
implicit none

  character (len=1) :: symb
  integer :: i, j, natm, numero
  real :: v
  real, dimension (:,:), allocatable :: xyz
  real, dimension (:), allocatable :: grad
  character(len=32)::fname, bname
  real,parameter::aucm = 219474.63
  real,parameter::auang= 0.5291772083d0

  call getarg(1,fname)
  i=index(fname,'.',.true.)
  bname=fname(1:i-1)

  open(12,file=trim(fname),status='old')
  open(13,file=trim(bname)//'.out',status='unknown')


  call pes_init()
  read (12,*) natm
  read (12,*)
  
  allocate(xyz(3,natm),grad(3*natm))

  do i=1,natm
     read (12,*)symb, xyz(:,i)
  end do

  xyz=xyz/auang
  v=f(xyz)
  write(13,'(A)') "Energy (in a.u.) from the fit:"
  write(13,'(F20.10)') v
  grad=g(xyz)

  write(13,'(A)') "Gradient (in a.u.) from the fit:"
  do i=1,size(grad)
     write(13,'(F20.10)') grad(i)
  end do

  close (12)
  close (13)

end program
