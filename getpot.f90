program main
use pes_shell
implicit none

  character (len=1) :: symb
  integer :: i, j, natm, numero
  real :: v
  real, dimension (:,:), allocatable :: xyz
  real, dimension (:), allocatable :: g1
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
  
  allocate(xyz(3,natm),g1(3*natm))

  do i=1,natm
     read (12,*)symb, xyz(:,i)
  end do

  xyz=xyz/auang
  v=f(xyz)
  write (13,'(3F20.8)') v

  g1=g(xyz)

  do i=1,size(g1)
     write(13,'(3F20.10)') g1(i)
  end do

  close (12)
  close (13)

end program
