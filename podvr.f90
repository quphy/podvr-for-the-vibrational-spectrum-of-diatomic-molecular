module podvr_module
 use, intrinsic:: iso_fortran_env, only : d8 => real64
 implicit none
 real(kind=d8),parameter ::pi = dacos(-1.0_d8)

 !1-d potential
 abstract interface
  subroutine potential_ab(x,V)
      import :: d8
      real(kind=d8), intent(in)::x 
      real (kind=d8), intent(inout)::V 
  end subroutine
 end interface

 type podvr

 !dvr
  integer ::n_dvr 
  real(kind=d8) :: weight
  real(kind=d8) :: mass
  real(kind=d8),allocatable,dimension(:) :: dvr_grid
  real(kind=d8),allocatable,dimension(:,:) :: dvr_coefficient
  real(kind=d8),allocatable,dimension(:) ::dvr_eigenvalue

 !podvr
  integer:: n_podvr
  real(kind=d8),allocatable,dimension(:) :: podvr_grid
  real(kind=d8),allocatable,dimension(:,:,:) :: podvr_coefficient
  real(kind=d8),allocatable,dimension(:) ::podvr_eigenvalue
  real(kind=d8),allocatable,dimension(:,:) :: energy_rovib
  procedure(potential_ab),pointer,nopass::potential=>null()

  contains

  procedure,public ::dvr_calculation
  procedure,public ::podvr_calculation
 end type podvr

 contains

 subroutine dvr_calculation(this,n_dvrgrid,n_podvrgrid,x_start,x_end,mass_au,potential)
 implicit none
 class(podvr) ::this
 integer,intent(in)::n_dvrgrid
 integer,intent(in)::n_podvrgrid
 real(kind=d8),intent(in)::x_start,x_end,mass_au
 procedure(potential_ab)::potential
 real(kind=d8),allocatable ::work(:)
 real(kind=d8) :: L, V
 integer :: i, j, lwork, info

 this%n_dvr=n_dvrgrid
 this%n_podvr=n_podvrgrid
 this%mass=mass_au
 this%potential=>potential

 allocate(this%dvr_grid(this%n_dvr))
 allocate(this%dvr_coefficient(this%n_dvr,this%n_dvr))
 allocate(this%dvr_eigenvalue(this%n_dvr))

 allocate(this%podvr_grid(this%n_podvr))
 allocate(this%podvr_eigenvalue(this%n_podvr))

 L=x_end-x_start

 associate(n=>this%n_dvr,C=>this%dvr_coefficient,m=>this%mass)

 do i =1,n
     this%dvr_grid(i)=x_start+i*L/(n+1)
 end do
 
 this%weight=L/(n+1)

 do i = 1,n
  do j = 1,i-1
    C(i,j) = pi**2 / (4.0_d8 * m * L**2) * (dsin(pi * (i-j) / (2.0_d8 * (n+1)))**(-2) - &
     dsin(pi * (i+j) / (2.0_d8 * (n+1)))**(-2)) * (-1.0_d8) ** (i-j)
    C(j,i)=C(i,j)
  end do
  call this%potential(this%dvr_grid(i), V)
  C(i,i) = pi**2 / (4.0_d8 * m * L**2) * ((2.0_d8*(n+1)**2 + 1.0_d8) / 3.0_d8 - &
  dsin(pi * i / (n+1))**(-2)) + V
 end do 

 allocate(work(1))
 call dsyev('V','U',n,C,n,this%dvr_eigenvalue,work,-1,info)
 lwork=int(work(1))
 deallocate(work)

 allocate(work(lwork))
 call dsyev('V','U',n,C,n,this%dvr_eigenvalue,work,lwork,info)
 deallocate(work)
 end associate
 end subroutine 

 subroutine podvr_calculation(this,vmax,jmax)
  implicit none
  class(podvr) this
  integer,intent(in)::vmax
  integer,intent(in)::jmax
  real(kind=d8),dimension(this%n_dvr,this%n_podvr) :: podvr_trans
  real(kind=d8),dimension(this%n_podvr,this%n_podvr) :: x_matrix
  real(kind=d8),dimension(this%n_podvr,this%n_podvr) :: H_ref_matrix
  real(kind=d8),allocatable,dimension(:,:) :: energy_matrix
  real(kind=d8), allocatable :: work(:)
  integer :: i, j, l, lwork, info

  allocate(this%energy_rovib(0:vmax, 0:jmax))
  allocate(this%podvr_coefficient(this%n_podvr, 0:vmax, 0:jmax))
  allocate(energy_matrix(this%n_podvr,this%n_podvr))

  associate(np=>this%n_podvr,nd=>this%n_dvr)
  do i=1,nd
   do j=1,np
   podvr_trans(i,j)=this%dvr_coefficient(i,j)
   end do
  end do
 
 this%energy_rovib=0.0
 do i=1,np
  call phase(nd,podvr_trans(:,i))
 end do
 
 block 
  real(kind=d8)::temp
  do i=1,np
   do j=i,np
    temp=0.0
    do l=1,nd
    temp=temp+podvr_trans(l,i)*this%dvr_grid(l)* podvr_trans(l,j)
    end do
    x_matrix(i,j)=temp
    x_matrix(j,i)=temp
   end do    
  end do
 end block
 
 allocate(work(1))
 call dsyev('V','U',np,x_matrix,np,this%podvr_grid,work,-1,info)
 lwork=int(work(1))
 deallocate(work)

 allocate(work(lwork))
 call dsyev('V','U',np,x_matrix,np,this%podvr_grid,work,lwork,info)
 deallocate(work) 
 
 do i =1,np
  call phase(np,x_matrix(:,i))
 end do
 
 block 
 real(kind=d8)::temp
 integer::j_qn,v_qn
 H_ref_matrix=0.0_d8
 do i=1,np
  do j=1,i
   temp=0.0
   do l=1,np
    temp=temp+x_matrix(l,i)*this%dvr_eigenvalue(l)*x_matrix(l,j)
   end do
   H_ref_matrix(i,j)=temp
   H_ref_matrix(j,i)=temp
  end do
 end do
 
do j_qn = 0, jmax
  energy_matrix(:,:) = H_ref_matrix(:,:)
  do i = 1,np
    energy_matrix(i,i) = energy_matrix(i,i) + (j_qn*(j_qn+1))/(2.0_d8*this%mass*this%podvr_grid(i)**2)
      end do 

    allocate(work(1))
    call dsyev('V', 'U', np, energy_matrix, np, this%podvr_eigenvalue, work, -1, info)
    lwork = int(work(1))
    deallocate(work)
                
    allocate(work(lwork))
    call dsyev('V', 'U', np, energy_matrix, np, this%podvr_eigenvalue, work, lwork, info)
    deallocate(work)

  do i = 1, np
   call phase(np, energy_matrix(:,i))
  end do 

  do v_qn = 0, vmax
   this%energy_rovib(v_qn, j_qn) = this%podvr_eigenvalue(v_qn+1)
   end do 
 do i = 1, np
   do v_qn = 0, vmax 
    this%podvr_coefficient(i,v_qn,j_qn) = energy_matrix(i,v_qn+1)
   end do 
 end do 
end do 
 end block
 end associate
 deallocate(energy_matrix)
 end subroutine



 subroutine phase(n,vector)
  implicit none
  integer, intent(in)::n
  real(kind=d8),intent(inout)::vector(n)
  integer ::i
  logical :: nonzero=.false.
  do i =1,n 
    if (abs(vector(i))>1.0e-4)then
     nonzero=.true.
     exit
    end if
  end do

  if(.not. nonzero)then
   error stop "elements of eigenvector of DVR are all near zero"
  end if

  if(vector(i)<0.0)then
  vector=-vector
  end if

 end subroutine
end module
