module user
  use neko
  implicit none

  real(kind=rp) :: Ra = 0
  real(kind=rp) :: Pr = 0
  real(kind=rp) :: ta2 = 0
  type(coef_t), pointer :: c_Xh

contains
  ! Register user defined functions (see user_intf.f90)
  subroutine user_setup(u)
    type(user_t), intent(inout) :: u
    u%user_init_modules => set_Pr
    u%fluid_user_ic => set_ic
    u%fluid_user_f_vector => forcing
    u%user_mesh_setup => deform
  end subroutine user_setup
  
 
  subroutine deform(msh)
    type(mesh_t), intent(inout) :: msh
    integer :: i, p, nvert
    real(kind=rp) :: d, z

    nvert = size(msh%points)
    do i = 1, nvert
       z = msh%points(i)%x(3)
       if (z .lt. 0.5_rp) then
          msh%points(i)%x(3) = 2.0_rp*z**2.0_rp
       else
          msh%points(i)%x(3) = 1.0_rp - 2.0_rp*(z-1.0_rp)**2.0_rp
       end if
    end do
 
  end subroutine deform
 

  !> Dummy user initial condition
  subroutine set_ic(u, v, w, p, params)
    type(field_t), intent(inout) :: u
    type(field_t), intent(inout) :: v
    type(field_t), intent(inout) :: w
    type(field_t), intent(inout) :: p
    type(param_t), intent(inout) :: params
    type(field_t), pointer :: s
    integer :: i, j, k, e
    real(kind=rp) :: rand
    s => neko_field_registry%get_field('s')

    call rzero(u%x,u%dof%size())
    call rzero(v%x,v%dof%size())
    call rzero(w%x,w%dof%size())
    call rzero(s%x,w%dof%size())
    do i = 1, s%dof%size()
       s%x(i,1,1,1) = 1-s%dof%z(i,1,1,1) 
      ! This seems like enough to trigger into turbulence
      ! s%x(i,1,1,1) = 1-s%dof%z(i,1,1,1) + 0.01*sin(40*pi*s%dof%x(i,1,1,1)) &
      ! * sin(40*pi*s%dof%y(i,1,1,1))*sqrt((0.5**2-(0.5_rp-s%dof%z(i,1,1,1))**2))
    end do
    ! perturb not on element boundaries
    do e = 1, s%msh%nelv
       do k = 2,s%Xh%lx-1
          do j = 2,s%Xh%lx-1
             do i = 2,s%Xh%lx-1

                call random_number(rand)
                s%x(i,j,k,e) = 1-s%dof%z(i,j,k,e) + 0.01*rand*sin(40*pi*s%dof%x(i,j,k,e)) &
                * sin(40*pi*s%dof%y(i,j,k,e))*sqrt((0.5**2-(0.5_rp-s%dof%z(i,j,k,e))**2))
             end do
          end do
       end do
    end do
    if ((NEKO_BCKND_CUDA .eq. 1) .or. (NEKO_BCKND_HIP .eq. 1) &
       .or. (NEKO_BCKND_OPENCL .eq. 1)) then
       call device_memcpy(s%x,s%x_d,s%dof%size(),HOST_TO_DEVICE)
    end if

  end subroutine set_ic

  subroutine set_Pr(t, u, v, w, p, coef, params)
    real(kind=rp) :: t
    type(field_t), intent(inout) :: u
    type(field_t), intent(inout) :: v
    type(field_t), intent(inout) :: w
    type(field_t), intent(inout) :: p
    type(coef_t), intent(inout) :: coef
    type(param_t), intent(inout) :: params
    ! Reset the relevant nondimensional parameters
    ! Pr = input Pr
    ! Ra = input Re
    ! Re = 1/Pr
    Pr = params%Pr
    Ra = params%Re
    params%Re = sqrt(Ra / Pr)
    call save_coef(coef)
  end subroutine set_Pr


  subroutine save_coef(coef)
    type(coef_t), target :: coef
    c_Xh => coef
  end subroutine save_coef

  !> Forcing
  subroutine forcing(f)
    class(source_t) :: f
    integer :: i
    type(field_t), pointer :: u, v, w, s
    real(kind=rp) :: rapr, ta2pr
    u => neko_field_registry%get_field('u')
    v => neko_field_registry%get_field('v')
    w => neko_field_registry%get_field('w')
    s => neko_field_registry%get_field('s')
    ta2pr = ta2*Pr

    if ((NEKO_BCKND_CUDA .eq. 1) .or. (NEKO_BCKND_HIP .eq. 1) &
       .or. (NEKO_BCKND_OPENCL .eq. 1)) then
       call device_rzero(f%u_d,f%dm%size())
       call device_rzero(f%v_d,f%dm%size())
       call device_copy(f%w_d,s%x_d,f%dm%size())
    else
       call cmult2(f%u,v%x,Ta2Pr,f%dm%size())
       call cmult2(f%v,u%x,Ta2Pr,f%dm%size())
       call copy(f%w,s%x,f%dm%size())
    end if
  end subroutine forcing
end module user
