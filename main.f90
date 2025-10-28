program H2_podvr
    use podvr_module
    use, intrinsic:: iso_fortran_env, only : d8 => real64
    implicit none

    type(podvr):: H2podvr

    integer :: dvrpoints, podvrpoints
    real(kind=d8) :: x_start, x_end, mass

    dvrpoints= 1000          
    podvrpoints=500         
    x_start=0.5_d8   
    x_end=5.0_d8     
    mass=0.5*1836.0_d8 
    call H2podvr%dvr_calculation(dvrpoints, podvrpoints, x_start, x_end, mass, potential_morse)
    call H2podvr%podvr_calculation(5, 10)


    print *, "v=2,j=3 energy (hartree)"
    print '(5F12.6)', H2podvr%energy_rovib(0,:)

contains

    subroutine potential_morse(x, V)
        use, intrinsic :: iso_fortran_env, only : d8 => real64
        real(kind=d8), intent(in)  :: x
        real(kind=d8), intent(inout) :: V
        real(kind=d8) :: De, a, re
        De = 0.1745_d8     
        a  = 1.028_d8    
        re = 1.4_d8      
        V  = De * (1.0_d8 - exp(-a * (x - re)))**2
    end subroutine 

end program 