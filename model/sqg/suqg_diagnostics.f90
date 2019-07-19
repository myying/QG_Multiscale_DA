module qg_diagnostics

  !************************************************************************
  ! Energetics and other diagnostics
  ! 
  ! Routines:  Get_energetics, Get_spectra, energy, enstrophy
  !
  ! Dependencies: qg_arrays, qg_params, io_tools, op_rules, transform_tools
  !
  !************************************************************************

  implicit none
  private
  save

  public :: Get_energetics, Get_spectra, energy, enstrophy

contains

  !*************************************************************************

  real function energy(psi)

    ! Get total energy
    
    use op_rules,  only: operator(*)
    use qg_arrays, only: ksqd_

    complex,dimension(:,:,:),intent(in) :: psi

    energy = sum(sqrt(ksqd_)*psi*conjg(psi))

  end function energy

  !*************************************************************************

  real function enstrophy(q)

    ! Get total enstrophy
    
    use op_rules,  only: operator(*)
    
    complex,dimension(:,:,:),intent(in) :: q

    enstrophy = sum( q*conjg(q) )

  end function enstrophy

  !*************************************************************************

  function Get_energetics(framein) result(dframe)

    !************************************************************************
    ! Calculate KE, APE, and ENS, as well as all rates of generation
    ! and dissipation.  Also calculates current eddy rotation period.
    ! All energetics factors multiplied by 2 since we are only using
    ! upper-half plane spectral fields
    !************************************************************************
    
    use op_rules,  only: operator(+), operator(-), operator(*)
    use io_tools,  only: Message, Write_field
    use qg_arrays, only: ksqd_,kx_,dz,psi,q,force_o,drho, &
                         tracer,qdrag,filter,filter_t,force_ot
    use qg_params, only: dt,kmax,therm_drag,i,pi,F,&
                         time,cntr,uscale,use_forcing,use_tracer,&
                         quad_drag,filter_type,filter_type_t, &
                         use_forcing_t,use_mean_grad_t

    integer,intent(in) :: framein
    integer            :: dframe
    real,dimension(:,:),allocatable  :: temp
    real               :: ke=0., ens=0.
    real               :: gen_rmf_rate=0.
    real               :: thd_rate=0.,qd_rate=0.,filter_rate=0.
    real               :: tvar=0.,tpsi=0.,filter_rate_t=0.
    real               :: gen_rmf_t=0., gen_tg_t=0.

    dframe = framein + 1        ! Update diagnostics frame counter
    call Write_field(time,'diag1_time',dframe) ! Track diagnostic-writes
    
    ke =  -sum(q*conjg(psi))
    call Write_field(ke,'energy',dframe)

    ens = sum(q*conjg(q)) 
    call Write_field(ens,'enstrophy',dframe) 

    if (trim(filter_type)/='none') then
       allocate(temp(-kmax:kmax,0:kmax)); temp = 0.
       where (filter/=0) temp = (filter**(-1)-1)/(2*dt)
       filter_rate = 2*sum(real( temp*(conjg(psi)*q)))
       deallocate(temp)
       call Write_field(filter_rate,'filter_rate',dframe)
    endif
    if (therm_drag/=0) then
       thd_rate = 2*sum(therm_drag*conjg(psi)*q)
       call Write_field(thd_rate,'thd_rate',dframe)
    endif
    if (quad_drag/=0) then
       qd_rate = 2*sum(conjg(psi)*qdrag)
       call Write_field(qd_rate,'qd_rate',dframe)
    endif
    if (use_forcing) then
       gen_rmf_rate = -2*sum(real(conjg(psi)*force_o))
       call Write_field(gen_rmf_rate,'gen_rmf_rate',dframe)
    endif
    if (use_tracer) then
       tvar = sum(dz*tracer*conjg(tracer))
       call Write_field(tvar,'tvar',dframe)
       tpsi = sum(dz*tracer*conjg(psi))
       call Write_field(tpsi,'tpsi',dframe)
       if (trim(filter_type_t)/='none') then
          allocate(temp(-kmax:kmax,0:kmax)); temp = 0.
          where (filter_t/=0) temp = (filter_t**(-1)-1)/(2*dt)
          filter_rate_t = -2*sum(real( temp &
               *(dz*(conjg(tracer)*tracer)) ))
          call Write_field(filter_rate_t,'filter_rate_t',dframe)
          deallocate(temp)
       endif
       if (use_mean_grad_t) then
          gen_tg_t = - 2*sum(real(dz*tracer*conjg(i*kx_*psi)))
          call Write_field(gen_tg_t,'gen_tg_t',dframe)
       endif
       if (use_forcing_t) then
          gen_rmf_t = - 2*sum(real(conjg(tracer)*force_o))
          call Write_field(gen_rmf_t,'gen_rmf_t',dframe)
       endif
    endif

    if (.not.ieee_is_finite(ke)) call Message('INF or NAN in ke - quitting!',fatal='y')

    call Message('')          
    call Message('time step     =',tag=cntr)
    call Message('energy        =',r_tag=ke)
    call Message('enstrophy     =',r_tag=ens)
    if (use_forcing)  call Message('gen_rmf_rate  =',r_tag=gen_rmf_rate)
    if (therm_drag/=0)call Message('thermdrag_rate=',r_tag=thd_rate)
    if (quad_drag/=0) call Message('quaddrag_rate =',r_tag=qd_rate)
    if (trim(filter_type)/='none') &
                      call Message('filter_rate   =',r_tag=filter_rate)  ! temp
    if (use_tracer) then
       call Message('tvariance     =',r_tag=tvar)
       call Message('tvar_gen      =',r_tag=gen_rmf_t+gen_tg_t)
       call Message('tvar_dissip   =',r_tag=filter_rate_t)
    endif

  end function Get_energetics

  !*************************************************************************

  function Get_spectra(framein) result(dframe)

    !************************************************************************
    ! Calculate the isotropic horizontal wavenumber vs vertical
    ! wavenumber spectra of modal and layered energetics
    !************************************************************************

    use op_rules,        only: operator(+), operator(-), operator(*)
    use io_tools,        only: Message, Write_field
    use transform_tools, only: Spec2grid, Jacob
    use numerics_lib,    only: Ring_integral, sub2ind
    use strat_tools,     only: Layer2mode
    use qg_arrays,       only: ksqd_,kx_,ky_,kxv,kyv,dz,psi,q, &
                               filter,force_o,drho,tracer,qdrag,filter_t
    use qg_params,       only: time,dt,nx,ny,kmax,nz,therm_drag, &
                               i,pi,uscale,cntr,quad_drag,&
                               use_forcing,use_tracer,do_xfer_spectra,&
                               filter_type,&
                               filter_type_t,use_mean_grad_t,do_aniso_spectra

    integer,intent(in)                     :: framein
    integer                                :: dframe
    real,dimension(:,:,:),allocatable      :: field
    real,dimension(:,:),allocatable        :: spec, temp

    dframe = framein + 1
    call Write_field(time,'diag2_time',dframe) ! Track diagnostic-writes
    
    allocate(spec(kmax,nz),field(-kmax:kmax,0:kmax,nz))

    ! KE spectra

    field = -q*conjg(psi)
    spec = Ring_integral(real(field),kxv,kyv,kmax)
    call Write_field(spec,'kes',dframe)

    ! Spectra of energy along kx and ky axes if anisotropy expected
    
    if (do_aniso_spectra) then
       spec = real(field(1:kmax,0,:))
       call Write_field(spec,'kesx',dframe)
       spec = real(field(0,1:kmax,:))
       call Write_field(spec,'kesy',dframe)
    endif    ! Calculate eddy generation spectra
    
    ! Calculate eddy generation spectra
    
    if (use_forcing) then        ! From random Markovian forcing
       field(:,:,1) = -conjg(psi(:,:,1))*force_o
       spec = Ring_integral(2*real(field),kxv,kyv,kmax)
       call Write_field(spec,'gens_rmf',dframe)
    endif

    if (do_xfer_spectra) then
       field = real(conjg(psi)*jacob(psi,q))
       spec = Ring_integral(2*real(field),kxv,kyv,kmax)
       call Write_field(spec,'xfers',dframe)
    endif

    if (therm_drag/=0) then     ! Thermal drag dissipation 
       field = therm_drag*conjg(psi) 
       spec = Ring_integral(2*real(field),kxv,kyv,kmax)
       call Write_field(spec,'thds',dframe)
    endif

    if (trim(filter_type)/='none') then    ! Filter dissipation
       allocate(temp(-kmax:kmax,0:kmax)); temp = 0.
       where (filter/=0) temp = (filter**(-1)-1)/(2*dt)
       field = temp*(q*conjg(psi))
       spec = Ring_integral(2*real(field),kxv,kyv,kmax)
       deallocate(temp)
       call Write_field(spec,'filters',dframe)        
    endif
    
    if (quad_drag/=0) then      ! Quadratic drag dissipation
       field = conjg(psi)*qdrag
       spec = Ring_integral(2*real(field),kxv,kyv,kmax)
       call Write_field(spec,'qds',dframe)
    endif

    ! Write tracer spectra

    if (use_tracer) then
       field = real(dz*tracer*conjg(tracer))        ! tracer variance <t't'>
       spec = Ring_integral(field,kxv,kyv,kmax)
       call Write_field(spec,'tvars',dframe)
       field = -real(dz*tracer*conjg(psi(:,:,1:nz))) ! tracer psi cor <t'psi'>
       spec = Ring_integral(field,kxv,kyv,kmax)
       call Write_field(spec,'tpsis',dframe)
       if (use_mean_grad_t) then
          field = -real(dz*tracer*conjg(i*kx_*psi(:,:,1:nz))) ! <t'v'> vs. k
          spec = Ring_integral(2*field,kxv,kyv,kmax)
          call Write_field(spec,'tfluxs',dframe)
       endif
    endif

    deallocate(field,spec)

  end function Get_spectra

  !*********************************************************************

   elemental logical function ieee_is_nan(x)
     real,intent(in):: x
     ieee_is_nan = isnan(x)
   end function ieee_is_nan
   elemental logical function ieee_is_finite(x)
     real,intent(in):: x
     ieee_is_finite = .not. (isnan(x) .or. abs(x) > huge(x))
   end function ieee_is_finite

  !*********************************************************************

end module qg_diagnostics
