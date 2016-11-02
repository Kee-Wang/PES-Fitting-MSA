module pes_shell
use basis
use gradient
implicit none

  real::coeff(1:XXX) ! change to number of coefficients
                     ! (size of c in bemsa.f90)
  save coeff

contains
  !==================================!
  ! read the coefficients of the PES !
  !==================================!
  subroutine pes_init()
    !::::::::::::::::::
    integer::i

    open(10,file='coeff.dat',status='old')

    do i=1,size(coeff) 
       read (10,*) coeff(i)
    end do

    return
    close (10)
  end subroutine pes_init

  !====================================!
  ! Function to evaluate the potential !
  !====================================!
  function f(xyz)
    real,dimension(:,:),intent(in)::xyz
    real::f
    !::::::::::::::::::::::::::::::
    real,dimension(size(xyz,2)*(size(xyz,2)-1)/2)::x
    real,dimension(3)::dr
    real::a0  ! the same as the fitting code
    integer::i,j,k

    a0 = 2.d0

    k = 1
    do i=1,size(xyz,2)-1
       do j=i+1,size(xyz,2)
          dr = xyz(:,i) - xyz(:,j)
          x(k) = sqrt(dot_product(dr,dr))
          k = k+1
       end do
    end do

    do i=1,size(x)
       x(i)=exp(-x(i)/a0)
    end do

    f=emsav(x,coeff)

    return
  end function f

  !===========================!
  ! function to calculate the !
  ! analytical gradient       !
  !===========================!
  function g(xyz)
    real,dimension(:,:),intent(in)::xyz
    real,dimension(size(xyz,2)*3)::g
    !::::::::::::::::::::
    real,dimension(size(xyz,2)*(size(xyz,2)-1)/2)::r,x
    real,dimension(3,size(xyz,2)*(size(xyz,2)-1)/2)::dr
    real,dimension(size(xyz,2)*3,size(xyz,2)*(size(xyz,2)-1)/2)::drdx
    real,dimension(1:YYY)::p   ! change to number of popynomials
                               ! (size of p in bemsa.f90)
    real,dimension(1:ZZZ)::m  ! change to number of monomials
                              ! (size of m in bemsa.f90
    real::a0
    integer::i,j,k

    a0 = 2.d0  ! the same as in fitting code

    k = 1
    drdx = 0.d0
    do i=1,size(xyz,2)-1
       do j=i+1,size(xyz,2)
          dr(:,k) = xyz(:,i) - xyz(:,j)
          r(k) = sqrt(dot_product(dr(:,k),dr(:,k)))

          drdx(3*i-2:3*i,k) = dr(:,k)/r(k)
          drdx(3*j-2:3*j,k) = -drdx(3*i-2:3*i,k)
          k = k+1
       end do
    end do

    do i=1,size(x)
       x(i)=exp(-r(i)/a0)
    end do

    call evmono(x,m)
    call evpoly(m,p)

    do i=1,3*size(xyz,2)
       g(i) = demsav(drdx,coeff,m,p,i)
    end do

    return
  end function g

end module pes_shell
