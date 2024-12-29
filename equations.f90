    ! Equations module for background and ! To avoid circular module issues, some things are not part of module

    ! Background evolution, return d tau/ d a, where tau is the conformal time
    function dtauda(this,a)
    use results
    implicit none
    class(CAMBdata) :: this
    real(dl), intent(in) :: a
    real(dl) :: dtauda

    dtauda = 1 / this%CP%DarkField%splines%dadtau%GetValue(a)

    end function dtauda

    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    !Gauge-dependent perturbation equations

    module GaugeInterface
    use precision
    use classes
    use results
    use MassiveNu
    implicit none
    public

    ! Equations and relation to synchronous gauge variables documented in the notes:
    ! https://cosmologist.info/notes/CAMB.pdf

    !Description of this file. Change if you make modifications.
    character(LEN=*), parameter :: Eqns_name = 'cdm_gauge'

    logical, parameter :: plot_evolve = .false. !for outputing time evolution

    integer, parameter :: ix_etak=1, ix_clxc=2, ix_df=3, ix_clxb=12, ix_vb=13 !Scalar array indices for each quantity
    integer, parameter :: basic_num_eqns = ix_vb

    logical :: DoTensorNeutrinos = .true.

    logical, parameter :: second_order_tightcoupling = .true.

    real(dl) :: Magnetic = 0._dl
    !Vector mode anisotropic stress in units of rho_gamma
    real(dl) :: vec_sig0 = 1._dl
    !Vector mode shear
    integer, parameter :: max_l_evolve = 256 !Maximum l we are ever likely to propagate
    !Note higher values increase size of Evolution vars, hence memory

    type, extends(TCambComponent) :: EvolutionVars
        real(dl) q, q2
        real(dl) k_buf,k2_buf ! set in initial

        integer Tg_ix !index of matter temerature perturbation
        integer reion_line_ix !index of matter temerature perturbation

        integer xe_ix !index of x_e perturbation
        integer Ts_ix !index of Delta_{T_s}

        integer r_ix !Index of the massless neutrino hierarchy
        integer g_ix !Index of the photon neutrino hierarchy

        integer q_ix !index into q_evolve array that gives the value q
        logical TransferOnly

        !       nvar  - number of scalar (tensor) equations for this k
        integer nvar,nvart, nvarv

        !Max_l for the various hierarchies
        integer lmaxg,lmaxnr,lmaxnu,lmaxgpol,MaxlNeeded
        integer lmaxnrt, lmaxnut, lmaxt, lmaxpolt, MaxlNeededt
        logical EvolveTensorMassiveNu(max_nu)
        integer lmaxnrv, lmaxv, lmaxpolv
        integer lmaxline !21cm multipoles for getting reionization effect

        integer polind  !index into scalar array of polarization hierarchy

        !array indices for massive neutrino equations
        integer nu_ix(max_nu), nu_pert_ix
        integer nq(max_nu), lmaxnu_pert
        logical has_nu_relativistic

        !Initial values for massive neutrino v*3 variables calculated when switching
        !to non-relativistic approx
        real(dl) G11(max_nu),G30(max_nu)
        !True when using non-relativistic approximation
        logical MassiveNuApprox(max_nu)
        real(dl) MassiveNuApproxTime(max_nu)

        !True when truncating at l=2,3 when k*tau>>1 (see arXiv:1104.2933)
        logical high_ktau_neutrino_approx

        !Massive neutrino scheme being used at the moment
        integer NuMethod

        !True when using tight-coupling approximation (required for stability at early times)
        logical TightCoupling, TensTightCoupling
        real(dl) TightSwitchoffTime

        !Numer of scalar equations we are propagating
        integer ScalEqsToPropagate
        integer TensEqsToPropagate
        !beta > l for _closed models
        integer FirstZerolForBeta
        !Tensor vars
        real(dl) aux_buf

        real(dl) pig, pigdot
        real(dl) poltruncfac

        logical no_nu_multpoles, no_phot_multpoles
        integer lmaxnu_tau(max_nu)  !lmax for massive neutinos at time being integrated
        logical nu_nonrelativistic(max_nu)

        real(dl) denlk(max_l_evolve),denlk2(max_l_evolve), polfack(max_l_evolve)
        real(dl) Kf(max_l_evolve)

        integer E_ix, B_ix !tensor polarizatisdon indices
        real(dl) denlkt(4,max_l_evolve),Kft(max_l_evolve)

        logical :: saha !still high x_e
        logical :: evolve_TM !\delta T_g evolved separately

        !Workaround for ifort, gives class pointer to avoid creating temps and huge slow down
        class(TThermoData), pointer :: ThermoData => null()

        real(dl), pointer :: OutputSources(:) => null()
        integer :: OutputStep = 0

    end type EvolutionVars

    class(CAMBdata), pointer :: State !Current state.
    !(Due to ifort bug State needs to be class pointer to avoid making temp class pointers when calling functions)
    type(CAMBParams), pointer :: CP   !Current parameters (part of state)

    !precalculated arrays
    real(dl) polfac(max_l_evolve),denl(max_l_evolve),vecfac(max_l_evolve),vecfacpol(max_l_evolve)

    real(dl), parameter :: ep0=1.0d-2
    integer, parameter :: lmaxnu_high_ktau=4 !Jan2015, increased from 3 to fix mpk for mnu~6eV

    real(dl) epsw
    real(dl), allocatable :: nu_tau_notmassless(:,:)
    real(dl) nu_tau_nonrelativistic(max_nu), nu_tau_massive(max_nu)

    contains

    subroutine SetActiveState(P)
    class(CAMBdata), target :: P

    State => P
    CP => P%CP

    end subroutine SetActiveState


    subroutine GaugeInterface_ScalEv(EV,y,tau,tauend,tol1,ind,c,w)
    type(EvolutionVars) EV
    real(dl) c(24),w(EV%nvar,9), y(EV%nvar), tol1, tau, tauend
    integer ind
    procedure(TClassDverk) :: dverk

    call dverk(EV,EV%ScalEqsToPropagate,derivs,tau,y,tauend,tol1,ind,c,EV%nvar,w)
    if (ind==-3) then
        call GlobalError('Dverk error -3: the subroutine was unable  to  satisfy  the  error ' &
            //'requirement  with a particular step-size that is less than or * ' &
            //'equal to hmin, which may mean that tol is too small' &
            //'--- but most likely you''ve messed up the y array indexing; ' &
            //'compiling with bounds checking may (or may not) help find the problem.',error_evolution)
    end if
    end subroutine GaugeInterface_ScalEv

    function f_K(x)
    real(dl) :: f_K
    real(dl), intent(in) :: x

    f_K = State%rofChi(x)

    end function f_K

    function next_nu_nq(nq) result (next_nq)
    integer, intent(in) :: nq
    integer q, next_nq

    if (nq==0) then
        next_nq=1
    else
        q = int(State%NuPerturbations%nu_q(nq))
        if (q>=10) then
            next_nq = State%NuPerturbations%nqmax
        else
            next_nq = nq+1
        end if
    end if

    end function next_nu_nq

    recursive subroutine GaugeInterface_EvolveScal(EV,tau,y,tauend,tol1,ind,c,w)
    use Recombination, only : CB1
    type(EvolutionVars) EV, EVout
    real(dl) c(24),w(EV%nvar,9), y(EV%nvar), yout(EV%nvar), tol1, tau, tauend
    integer ind, nu_i
    real(dl) cs2, opacity, dopacity
    real(dl) tau_switch_ktau, tau_switch_nu_massless, tau_switch_nu_massive, next_switch
    real(dl) tau_switch_no_nu_multpoles, tau_switch_no_phot_multpoles,tau_switch_nu_nonrel
    real(dl) noSwitch, smallTime
    !Sources
    real(dl) tau_switch_saha, Delta_TM, xe,a,tau_switch_evolve_TM

    noSwitch= State%tau0+1
    smallTime =  min(tau, 1/EV%k_buf)/100

    tau_switch_ktau = noSwitch
    tau_switch_no_nu_multpoles= noSwitch
    tau_switch_no_phot_multpoles= noSwitch

    !Massive neutrino switches
    tau_switch_nu_massless = noSwitch
    tau_switch_nu_nonrel = noSwitch
    tau_switch_nu_massive= noSwitch

    !Sources
    tau_switch_saha=noSwitch
    tau_switch_evolve_TM=noSwitch

    !Evolve equations from tau to tauend, performing switches in equations if necessary.

    if (.not. EV%high_ktau_neutrino_approx .and. .not. EV%no_nu_multpoles ) then
        tau_switch_ktau=  max(20, EV%lmaxnr-4)/EV%k_buf
    end if

    if (CP%Num_Nu_massive /= 0) then
        do nu_i = 1, CP%Nu_mass_eigenstates
            if (EV%nq(nu_i) /= State%NuPerturbations%nqmax) then
                tau_switch_nu_massless = min(tau_switch_nu_massless,nu_tau_notmassless(next_nu_nq(EV%nq(nu_i)),nu_i))
            else if (.not. EV%nu_nonrelativistic(nu_i)) then
                tau_switch_nu_nonrel = min(nu_tau_nonrelativistic(nu_i),tau_switch_nu_nonrel)
            else if (EV%NuMethod==Nu_trunc .and..not. EV%MassiveNuApprox(nu_i)) then
                tau_switch_nu_massive = min(tau_switch_nu_massive,EV%MassiveNuApproxTime(nu_i))
            end if
        end do
    end if

    if (.not. EV%no_nu_multpoles) & !!.and. .not. EV%has_nu_relativistic .and. tau_switch_nu_massless ==noSwitch)  &
        tau_switch_no_nu_multpoles= &
        max(15/EV%k_buf*CP%Accuracy%AccuracyBoost,min(State%taurend,EV%ThermoData%matter_verydom_tau))

    if (.not. EV%no_phot_multpoles .and. (.not.CP%WantCls .or. EV%k_buf>0.03*CP%Accuracy%AccuracyBoost)) &
        tau_switch_no_phot_multpoles =max(15/EV%k_buf,State%taurend)*CP%Accuracy%AccuracyBoost

    next_switch = min(tau_switch_ktau, tau_switch_nu_massless,EV%TightSwitchoffTime, tau_switch_nu_massive, &
        tau_switch_no_nu_multpoles, tau_switch_no_phot_multpoles, tau_switch_nu_nonrel, noSwitch, &
        tau_switch_saha, tau_switch_evolve_TM)

    if (next_switch < tauend) then
        if (next_switch > tau+smallTime) then
            call GaugeInterface_ScalEv(EV, y, tau,next_switch,tol1,ind,c,w)
            if (global_error_flag/=0) return
        end if

        EVout=EV

        if (next_switch == EV%TightSwitchoffTime) then
            !TightCoupling
            EVout%TightCoupling=.false.
            EVout%TightSwitchoffTime = noSwitch
            call SetupScalarArrayIndices(EVout)
            call CopyScalarVariableArray(y,yout, EV, EVout)
            EV=EVout
            y=yout
            ind=1
            !Set up variables with their tight coupling values
            y(EV%g_ix+2) = EV%pig
            call EV%ThermoData%Values(tau,a, cs2,opacity,dopacity)

            if (second_order_tightcoupling) then
                ! Francis-Yan Cyr-Racine November 2010

                y(EV%g_ix+3) = (3._dl/7._dl)*y(EV%g_ix+2)*(EV%k_buf/opacity)*(1._dl+dopacity/opacity**2) + &
                    (3._dl/7._dl)*EV%pigdot*(EV%k_buf/opacity**2)*(-1._dl)

                y(EV%polind+2) = EV%pig/4 + EV%pigdot*(1._dl/opacity)*(-5._dl/8._dl- &
                    (25._dl/16._dl)*dopacity/opacity**2) + &
                    EV%pig*(EV%k_buf/opacity)**2*(-5._dl/56._dl)
                y(EV%polind+3) = (3._dl/7._dl)*(EV%k_buf/opacity)*y(EV%polind+2)*(1._dl + &
                    dopacity/opacity**2) + (3._dl/7._dl)*(EV%k_buf/opacity**2)*((EV%pigdot/4._dl)* &
                    (1._dl+(5._dl/2._dl)*dopacity/opacity**2))*(-1._dl)
            else
                y(EV%g_ix+3) = 3./7*y(EV%g_ix+2)*EV%k_buf/opacity
                y(EV%polind+2) = EV%pig/4
                y(EV%polind+3) =y(EV%g_ix+3)/4
            end if
        else if (next_switch==tau_switch_ktau) then
            !k tau >> 1, evolve massless neutrino effective fluid up to l=2
            EVout%high_ktau_neutrino_approx=.true.
            EVout%nq(1:CP%Nu_mass_eigenstates) = State%NuPerturbations%nqmax
            call SetupScalarArrayIndices(EVout)
            call CopyScalarVariableArray(y,yout, EV, EVout)
            y=yout
            EV=EVout
        else if (next_switch == tau_switch_nu_massless) then
            !Mass starts to become important, start evolving next momentum mode
            do nu_i = 1, CP%Nu_mass_eigenstates
                if (EV%nq(nu_i) /= State%NuPerturbations%nqmax) then
                    if (next_switch == nu_tau_notmassless(next_nu_nq(EV%nq(nu_i)),nu_i)) then
                        EVOut%nq(nu_i) = next_nu_nq(EV%nq(nu_i))
                        call SetupScalarArrayIndices(EVout)
                        call CopyScalarVariableArray(y,yout, EV, EVout)
                        EV=EVout
                        y=yout
                        exit
                    endif
                end if
            end do
        else if (next_switch == tau_switch_nu_nonrel) then
            !Neutrino becomes non-relativistic, don't need high L
            do nu_i = 1, CP%Nu_mass_eigenstates
                if (.not. EV%nu_nonrelativistic(nu_i) .and.  next_switch==nu_tau_nonrelativistic(nu_i) ) then
                    EVout%nu_nonrelativistic(nu_i)=.true.
                    call SetupScalarArrayIndices(EVout)
                    call CopyScalarVariableArray(y,yout, EV, EVout)
                    EV=EVout
                    y=yout
                    exit
                end if
            end do
        else if (next_switch == tau_switch_nu_massive) then
            !Very non-relativistic neutrinos, switch to truncated velocity-weight hierarchy
            call EV%ThermoData%Values(tau,a, cs2,opacity)
            do nu_i = 1, CP%Nu_mass_eigenstates
                if (.not. EV%MassiveNuApprox(nu_i) .and.  next_switch== EV%MassiveNuApproxTime(nu_i) ) then
                    call SwitchToMassiveNuApprox(EV, a, y, nu_i)
                    exit
                end if
            end do
        else if (next_switch==tau_switch_no_nu_multpoles) then
            !Turn off neutrino hierarchies at late time where slow and not needed.
            ind=1
            EVout%no_nu_multpoles=.true.
            EVOut%nq(1:CP%Nu_mass_eigenstates ) = State%NuPerturbations%nqmax
            call SetupScalarArrayIndices(EVout)
            call CopyScalarVariableArray(y,yout, EV, EVout)
            y=yout
            EV=EVout
        else if (next_switch==tau_switch_no_phot_multpoles) then
            !Turn off photon hierarchies at late time where slow and not needed.
            ind=1
            EVout%no_phot_multpoles=.true.
            call SetupScalarArrayIndices(EVout)
            call CopyScalarVariableArray(y,yout, EV, EVout)
            y=yout
            EV=EVout
        else if (next_switch==tau_switch_saha) then
            !Sources
            ind=1
            EVout%saha = .false.
            call SetupScalarArrayIndices(EVout)
            call CopyScalarVariableArray(y,yout, EV, EVout)
            call EV%ThermoData%Values(tau,a, cs2,opacity)
            y=yout
            EV=EVout
            Delta_Tm = y(EV%g_ix)/4 ! assume delta_TM = delta_T_gamma
            xe= CP%Recomb%x_e(a)
            y(EV%xe_ix) = (1-xe)/(2-xe)*(-y(ix_clxb) + (3./2+  CB1/(CP%TCMB/a))*Delta_TM)
        else if (next_switch==tau_switch_evolve_TM) then
            !Sources
            ind=1
            EVout%evolve_TM = .true.
            call SetupScalarArrayIndices(EVout)
            call CopyScalarVariableArray(y,yout, EV, EVout)
            y=yout
            EV=EVout
            y(EV%Tg_ix) =y(EV%g_ix)/4 ! assume delta_TM = delta_T_gamma
        end if

        call GaugeInterface_EvolveScal(EV,tau,y,tauend,tol1,ind,c,w)
        return
    end if

    call GaugeInterface_ScalEv(EV,y,tau,tauend,tol1,ind,c,w)

    end subroutine GaugeInterface_EvolveScal

    function DeltaTimeMaxed(a1,a2, tol) result(t)
    real(dl) a1,a2,t
    real(dl), optional :: tol
    if (a1>1._dl) then
        t=0
    elseif (a2 > 1._dl) then
        t = State%DeltaTime(a1,1.01_dl, tol)
    else
        t = State%DeltaTime(a1,a2, tol)
    end if
    end function DeltaTimeMaxed

    subroutine GaugeInterface_Init
    !Precompute various arrays and other things independent of wavenumber
    integer j, nu_i
    real(dl) a_nonrel, a_mass,a_massive, time, nu_mass

    epsw = 100/State%tau0

    if (CP%WantScalars) then
        do j=2,max_l_evolve
            polfac(j)=real((j+3)*(j-1),dl)/(j+1)
        end do
    end if

    do j=1,max_l_evolve
        denl(j)=1._dl/(2*j+1)
    end do

    if (CP%Nu_Mass_eigenstates>0) then
        associate(nqmax => State%NuPerturbations%nqmax, nu_q => State%NuPerturbations%nu_q)
            if (allocated(nu_tau_notmassless)) deallocate(nu_tau_notmassless)
            allocate(nu_tau_notmassless(nqmax,max_nu))
            do nu_i=1, CP%Nu_Mass_eigenstates
                nu_mass = max(0.1_dl,State%nu_masses(nu_i))
                a_mass =  1.e-1_dl/nu_mass/CP%Accuracy%lAccuracyBoost
                time=State%DeltaTime(0._dl,State%NuPerturbations%nu_q(1)*a_mass)
                nu_tau_notmassless(1, nu_i) = time
                do j=2,nqmax
                    !times when each momentum mode becomes signficantly nonrelativistic
                    time= time + DeltaTimeMaxed(nu_q(j-1)*a_mass,nu_q(j)*a_mass, 0.01_dl)
                    nu_tau_notmassless(j, nu_i) = time
                end do

                a_nonrel =  2.5d0/nu_mass*CP%Accuracy%AccuracyBoost
                nu_tau_nonrelativistic(nu_i) =DeltaTimeMaxed(0._dl,a_nonrel)
                a_massive =  17.d0/nu_mass*CP%Accuracy%AccuracyBoost
                nu_tau_massive(nu_i) =nu_tau_nonrelativistic(nu_i) + DeltaTimeMaxed(a_nonrel,a_massive)
            end do
        end associate
    end if

    end subroutine GaugeInterface_Init


    subroutine SetupScalarArrayIndices(EV, max_num_eqns)
    !Set up array indices after the lmax have been decided
    use MassiveNu
    !Set the numer of equations in each hierarchy, and get total number of equations for this k
    type(EvolutionVars) EV
    integer, intent(out), optional :: max_num_eqns
    integer neq, maxeq, nu_i

    neq=basic_num_eqns
    maxeq=neq
    if (.not. EV%no_phot_multpoles) then
        !Photon multipoles
        EV%g_ix=basic_num_eqns+1
        if (EV%TightCoupling) then
            neq=neq+2
        else
            neq = neq+ (EV%lmaxg+1)
            !Polarization multipoles
            EV%polind = neq -1 !polind+2 is L=2, for polarizationthe first calculated
            neq=neq + EV%lmaxgpol-1
        end if
    end if
    if (.not. EV%no_nu_multpoles) then
        !Massless neutrino multipoles
        EV%r_ix= neq+1
        if (EV%high_ktau_neutrino_approx) then
            neq=neq + 3
        else
            neq=neq + (EV%lmaxnr+1)
        end if
    end if
    maxeq = maxeq +  (EV%lmaxg+1)+(EV%lmaxnr+1)+EV%lmaxgpol-1

    !Massive neutrinos
    if (CP%Num_Nu_massive /= 0) then
        EV%has_nu_relativistic = any(EV%nq(1:CP%Nu_Mass_eigenstates)/=State%NuPerturbations%nqmax)
        if (EV%has_nu_relativistic) then
            EV%lmaxnu_pert=EV%lmaxnu
            EV%nu_pert_ix=neq+1
            neq = neq+ EV%lmaxnu_pert+1
            maxeq=maxeq+ EV%lmaxnu_pert+1
        else
            EV%lmaxnu_pert=0
        end if

        do nu_i=1, CP%Nu_Mass_eigenstates
            if (EV%high_ktau_neutrino_approx) then
                EV%lmaxnu_tau(nu_i) = int(lmaxnu_high_ktau *CP%Accuracy%lAccuracyBoost)
            else
                EV%lmaxnu_tau(nu_i) =max(min(nint(0.8_dl*EV%q*nu_tau_nonrelativistic(nu_i) &
                    *CP%Accuracy%lAccuracyBoost),EV%lmaxnu),3)
                !!!Feb13tweak
                if (EV%nu_nonrelativistic(nu_i)) EV%lmaxnu_tau(nu_i)= &
                    min(EV%lmaxnu_tau(nu_i),nint(4*CP%Accuracy%lAccuracyBoost))
            end if
            if (State%nu_masses(nu_i) > 5000) &
                EV%lmaxnu_tau(nu_i) = EV%lmaxnu_tau(nu_i)*2 !megadamping
            EV%lmaxnu_tau(nu_i)=min(EV%lmaxnu,EV%lmaxnu_tau(nu_i))

            EV%nu_ix(nu_i)=neq+1
            if (EV%MassiveNuApprox(nu_i)) then
                neq = neq+4
            else
                neq = neq+ EV%nq(nu_i)*(EV%lmaxnu_tau(nu_i)+1)
            endif
            maxeq = maxeq + State%NuPerturbations%nqmax*(EV%lmaxnu+1)
        end do
    else
        EV%has_nu_relativistic = .false.
    end if

    EV%ScalEqsToPropagate = neq
    if (present(max_num_eqns)) then
        max_num_eqns=maxeq
    end if

    end subroutine SetupScalarArrayIndices

    subroutine CopyScalarVariableArray(y,yout, EV, EVout)
    type(EvolutionVars) EV, EVOut
    real(dl), intent(in) :: y(EV%nvar)
    real(dl), intent(out) :: yout(EVout%nvar)
    integer lmax,i, nq
    integer nnueq,nu_i, ix_off, ix_off2, ind, ind2
    real(dl) q, pert_scale

    yout=0
    yout(1:basic_num_eqns) = y(1:basic_num_eqns)

    if (.not. EV%no_phot_multpoles .and. .not. EVout%no_phot_multpoles) then
        if (EV%TightCoupling .or. EVOut%TightCoupling) then
            lmax=1
        else
            lmax = min(EV%lmaxg,EVout%lmaxg)
        end if
        yout(EVout%g_ix:EVout%g_ix+lmax)=y(EV%g_ix:EV%g_ix+lmax)
        if (.not. EV%TightCoupling .and. .not. EVOut%TightCoupling) then
            lmax = min(EV%lmaxgpol,EVout%lmaxgpol)
            yout(EVout%polind+2:EVout%polind+lmax)=y(EV%polind+2:EV%polind+lmax)
        end if
    end if

    if (.not. EV%no_nu_multpoles .and. .not. EVout%no_nu_multpoles) then
        if (EV%high_ktau_neutrino_approx .or. EVout%high_ktau_neutrino_approx) then
            lmax=2
        else
            lmax = min(EV%lmaxnr,EVout%lmaxnr)
        end if
        yout(EVout%r_ix:EVout%r_ix+lmax)=y(EV%r_ix:EV%r_ix+lmax)
    end if

    if (CP%Num_Nu_massive /= 0) then
        do nu_i=1,CP%Nu_mass_eigenstates
            ix_off=EV%nu_ix(nu_i)
            ix_off2=EVOut%nu_ix(nu_i)
            if (EV%MassiveNuApprox(nu_i) .and. EVout%MassiveNuApprox(nu_i)) then
                nnueq=4
                yout(ix_off2:ix_off2+nnueq-1)=y(ix_off:ix_off+nnueq-1)
            else if (.not. EV%MassiveNuApprox(nu_i) .and. .not. EVout%MassiveNuApprox(nu_i)) then
                lmax=min(EV%lmaxnu_tau(nu_i),EVOut%lmaxnu_tau(nu_i))
                nq = min(EV%nq(nu_i), EVOut%nq(nu_i))
                do i=1,nq
                    ind= ix_off + (i-1)*(EV%lmaxnu_tau(nu_i)+1)
                    ind2=ix_off2+ (i-1)*(EVOut%lmaxnu_tau(nu_i)+1)
                    yout(ind2:ind2+lmax) = y(ind:ind+lmax)
                end do
                do i=nq+1, EVOut%nq(nu_i)
                    lmax = min(EVOut%lmaxnu_tau(nu_i), EV%lmaxnr)
                    ind2=ix_off2+ (i-1)*(EVOut%lmaxnu_tau(nu_i)+1)
                    yout(ind2:ind2+lmax) = y(EV%r_ix:EV%r_ix+lmax)

                    !Add leading correction for the mass
                    q=State%NuPerturbations%nu_q(i)
                    pert_scale=(State%nu_masses(nu_i)/q)**2/2
                    lmax = min(lmax,EV%lmaxnu_pert)
                    yout(ind2:ind2+lmax) = yout(ind2:ind2+lmax) &
                        + y(EV%nu_pert_ix:EV%nu_pert_ix+lmax)*pert_scale
                end do
            end if
        end do

        if (EVOut%has_nu_relativistic .and. EV%has_nu_relativistic) then
            lmax = min(EVOut%lmaxnu_pert, EV%lmaxnu_pert)
            yout(EVout%nu_pert_ix:EVout%nu_pert_ix+lmax)=  y(EV%nu_pert_ix:EV%nu_pert_ix+lmax)
        end if
    end if
    !Sources
    if (.not. EV%saha .and. .not. EVOut%saha) then
        yout(EVOut%xe_ix) =y(EV%xe_ix)
    end if

    end subroutine CopyScalarVariableArray


    subroutine GetNumEqns(EV)
    use MassiveNu
    !Set the numer of equations in each hierarchy, and get total number of equations for this k
    type(EvolutionVars) EV
    real(dl) scal, max_nu_mass
    integer nu_i,q_rel,j

    if (CP%Num_Nu_massive == 0) then
        EV%lmaxnu=0
        max_nu_mass=0
    else
        max_nu_mass = maxval(State%nu_masses(1:CP%Nu_mass_eigenstates))
        do nu_i = 1, CP%Nu_mass_eigenstates
            !Start with momentum modes for which t_k ~ time at which mode becomes non-relativistic
            q_rel=0
            do j=1, State%NuPerturbations%nqmax
                !two different q's here EV%q ~k
                if (State%NuPerturbations%nu_q(j) > State%nu_masses(nu_i)*State%adotrad/EV%q) exit
                q_rel = q_rel + 1
            end do

            if (q_rel>= State%NuPerturbations%nqmax-2) then
                EV%nq(nu_i)=State%NuPerturbations%nqmax
            else
                EV%nq(nu_i)=q_rel
            end if
            !q_rel = nint(nu_masses(nu_i)*adotrad/EV%q) !two dffierent q's here EV%q ~k
            !EV%nq(nu_i)=max(0,min(nqmax0,q_rel)) !number of momentum modes to evolve intitially
            EV%nu_nonrelativistic(nu_i) = .false.
        end do

        EV%NuMethod = CP%MassiveNuMethod
        if (EV%NuMethod == Nu_Best) EV%NuMethod = Nu_Trunc
        !l_max for massive neutrinos
        EV%lmaxnu=nint(25*CP%Accuracy%lAccuracyBoost)
    end if

    EV%FirstZerolForBeta= 100000 !a large number

    EV%high_ktau_neutrino_approx = .false.
    if (CP%WantScalars) then
        EV%TightCoupling=.true.
        EV%no_phot_multpoles =.false.
        EV%no_nu_multpoles =.false.
        EV%MassiveNuApprox=.false.
        !Sources
        EV%saha = .true.
        EV%Evolve_TM = .false.

        EV%lmaxg  = max(nint(11*CP%Accuracy%lAccuracyBoost),3)
        EV%lmaxnr = max(nint(14*CP%Accuracy%lAccuracyBoost),3)
        if (max_nu_mass>700) EV%lmaxnr = max(nint(32*CP%Accuracy%lAccuracyBoost),3) !Feb13 tweak

        EV%lmaxgpol = EV%lmaxg

        if (EV%q < 0.05) then
            !Large scales need fewer equations
            scal  = 1
            scal = 4  !But need more to get polarization right
            EV%lmaxgpol=max(3,nint(min(8,nint(scal* 150* EV%q))*CP%Accuracy%lAccuracyBoost))
            EV%lmaxnr=max(3,nint(min(7,nint(sqrt(scal)* 150 * EV%q))*CP%Accuracy%lAccuracyBoost))
            if (EV%lmaxnr < EV%lmaxnu) then
                ! Nov 2020 change following Pavel Motloch report
                EV%lmaxnr = EV%lmaxnu
                !EV%lmaxnu = min(EV%lmaxnu, EV%lmaxnr) ! may be better but have not tested and makes small result changes
            endif
            EV%lmaxg=max(3,nint(min(8,nint(sqrt(scal) *300 * EV%q))*CP%Accuracy%lAccuracyBoost))
            !Sources
            EV%lmaxg=EV%lmaxg*4
            EV%lmaxgpol=EV%lmaxgpol*2
        end if

        if (EV%TransferOnly) then
            EV%lmaxgpol = min(EV%lmaxgpol,nint(5*CP%Accuracy%lAccuracyBoost))
            EV%lmaxg = min(EV%lmaxg,nint(6*CP%Accuracy%lAccuracyBoost))
        end if
        EV%lmaxnr=max(nint(45*CP%Accuracy%lAccuracyBoost),3)
        if (EV%q > 0.04 .and. EV%q < 0.5) then !baryon oscillation scales
            EV%lmaxg=max(EV%lmaxg,10)
        end if

        EV%poltruncfac=real(EV%lmaxgpol,dl)/max(1,(EV%lmaxgpol-2))
        EV%MaxlNeeded=max(EV%lmaxg,EV%lmaxnr,EV%lmaxgpol,EV%lmaxnu)
        if (EV%MaxlNeeded > max_l_evolve) call MpiStop('Need to increase max_l_evolve')
        call SetupScalarArrayIndices(EV,EV%nvar)
        EV%lmaxt=0
    else
        EV%nvar=0
    end if

    EV%nvarv=0

    end subroutine GetNumEqns

    subroutine SwitchToMassiveNuApprox(EV,a, y, nu_i)
    !When the neutrinos are no longer highly relativistic we use a truncated
    !energy-integrated hierarchy going up to third order in velocity dispersion
    type(EvolutionVars) EV, EVout
    integer, intent(in) :: nu_i
    real(dl) a,a2,pnu,clxnu,dpnu,pinu,rhonu
    real(dl) qnu
    real(dl) y(EV%nvar), yout(EV%nvar)

    a2=a*a
    EVout=EV
    EVout%MassiveNuApprox(nu_i)=.true.
    call SetupScalarArrayIndices(EVout)
    call CopyScalarVariableArray(y,yout, EV, EVout)

    !Get density and pressure as ratio to massles by interpolation from table
    call ThermalNuBack%rho_P(a*State%nu_masses(nu_i),rhonu,pnu)

    !Integrate over q
    call Nu_Integrate_L012(EV, y, a, nu_i, clxnu,qnu,dpnu,pinu)
    !clxnu_here  = rhonu*clxnu, qnu_here = qnu*rhonu
    dpnu=dpnu/rhonu
    qnu=qnu/rhonu
    clxnu = clxnu/rhonu
    pinu=pinu/rhonu

    yout(EVout%nu_ix(nu_i))=clxnu
    yout(EVout%nu_ix(nu_i)+1)=dpnu
    yout(EVout%nu_ix(nu_i)+2)=qnu
    yout(EVout%nu_ix(nu_i)+3)=pinu

    call Nu_Intvsq(EV,y, a, nu_i, EVout%G11(nu_i),EVout%G30(nu_i))
    !Analytic solution for higher moments, proportional to a^{-3}
    EVout%G11(nu_i)=EVout%G11(nu_i)*a2*a/rhonu
    EVout%G30(nu_i)=EVout%G30(nu_i)*a2*a/rhonu

    EV=EVout
    y=yout

    end subroutine SwitchToMassiveNuApprox

    subroutine MassiveNuVarsOut(EV,y,yprime,a,adotoa,grho,gpres,dgrho,dgq,dgpi, dgpi_diff,pidot_sum,clxnu_all)
    implicit none
    type(EvolutionVars) EV
    real(dl) :: y(EV%nvar), yprime(EV%nvar),a, adotoa
    real(dl), optional :: grho,gpres,dgrho,dgq,dgpi, dgpi_diff,pidot_sum,clxnu_all
    !grho = a^2 kappa rho
    !gpres = a^2 kappa p
    !dgrho = a^2 kappa \delta\rho
    !dgp =  a^2 kappa \delta p
    !dgq = a^2 kappa q (heat flux)
    !dgpi = a^2 kappa pi (anisotropic stress)
    !dgpi_diff = a^2 kappa (3*p -rho)*pi

    integer nu_i
    real(dl) pinudot,grhormass_t, rhonu, pnu,  rhonudot
    real(dl) grhonu_t,gpnu_t
    real(dl) clxnu, qnu, pinu, dpnu, grhonu, dgrhonu

    grhonu=0
    dgrhonu=0
    do nu_i = 1, CP%Nu_mass_eigenstates
        grhormass_t=State%grhormass(nu_i)/a**2

        !Get density and pressure as ratio to massless by interpolation from table
        call ThermalNuBack%rho_P(a*State%nu_masses(nu_i),rhonu,pnu)

        if (EV%MassiveNuApprox(nu_i)) then
            clxnu=y(EV%nu_ix(nu_i))
            !dpnu = y(EV%iq0+1+off_ix)
            qnu=y(EV%nu_ix(nu_i)+2)
            pinu=y(EV%nu_ix(nu_i)+3)
            pinudot=yprime(EV%nu_ix(nu_i)+3)
        else
            !Integrate over q
            call Nu_Integrate_L012(EV, y, a, nu_i, clxnu,qnu,dpnu,pinu)
            !clxnu_here  = rhonu*clxnu, qnu_here = qnu*rhonu
            !dpnu=dpnu/rhonu
            qnu=qnu/rhonu
            clxnu = clxnu/rhonu
            pinu=pinu/rhonu
            rhonudot = ThermalNuBack%drho(a*State%nu_masses(nu_i),adotoa)

            call Nu_pinudot(EV,y, yprime, a,adotoa, nu_i,pinudot)
            pinudot=pinudot/rhonu - rhonudot/rhonu*pinu
        endif

        grhonu_t=grhormass_t*rhonu
        gpnu_t=grhormass_t*pnu

        grhonu = grhonu  + grhonu_t
        if (present(gpres)) gpres= gpres + gpnu_t

        dgrhonu= dgrhonu + grhonu_t*clxnu
        if (present(dgq)) dgq  = dgq   + grhonu_t*qnu
        if (present(dgpi)) dgpi = dgpi  + grhonu_t*pinu
        if (present(dgpi_diff)) dgpi_diff = dgpi_diff + pinu*(3*gpnu_t-grhonu_t)
        if (present(pidot_sum)) pidot_sum = pidot_sum + grhonu_t*pinudot
    end do
    if (present(grho)) grho = grho  + grhonu
    if (present(dgrho)) dgrho= dgrho + dgrhonu
    if (present(clxnu_all)) clxnu_all = dgrhonu/grhonu

    end subroutine MassiveNuVarsOut

    subroutine Nu_Integrate_L012(EV,y,a,nu_i,drhonu,fnu,dpnu,pinu)
    type(EvolutionVars) EV
    !  Compute the perturbations of density and energy flux
    !  of one eigenstate of massive neutrinos, in units of the mean
    !  density of one eigenstate of massless neutrinos, by integrating over
    !  momentum.
    integer, intent(in) :: nu_i
    real(dl), intent(in) :: a, y(EV%nvar)
    real(dl), intent(OUT) ::  drhonu,fnu
    real(dl), optional, intent(OUT) :: dpnu,pinu
    real(dl) tmp, am, aq,v, pert_scale
    integer iq, ind

    !  q is the comoving momentum in units of k_B*T_nu0/c.

    drhonu=0
    fnu=0
    if (present(dpnu)) then
        dpnu=0
        pinu=0
    end if
    am=a*State%nu_masses(nu_i)
    ind=EV%nu_ix(nu_i)
    associate(nu_q=>State%NuPerturbations%nu_q, nu_int_kernel=>State%NuPerturbations%nu_int_kernel)
        do iq=1,EV%nq(nu_i)
            aq=am/nu_q(iq)
            v=1._dl/sqrt(1._dl+aq*aq)
            drhonu=drhonu+ nu_int_kernel(iq)* y(ind)/v
            fnu=fnu+nu_int_kernel(iq)* y(ind+1)
            if (present(dpnu)) then
                dpnu=dpnu+  nu_int_kernel(iq)* y(ind)*v
                pinu=pinu+ nu_int_kernel(iq)*y(ind+2)*v
            end if
            ind=ind+EV%lmaxnu_tau(nu_i)+1
        end do
        ind = EV%nu_pert_ix
        do iq=EV%nq(nu_i)+1,State%NuPerturbations%nqmax
            !Get the rest from perturbatively relativistic expansion
            aq=am/nu_q(iq)
            v=1._dl/sqrt(1._dl+aq*aq)
            pert_scale=(State%nu_masses(nu_i)/nu_q(iq))**2/2
            tmp = nu_int_kernel(iq)*(y(EV%r_ix)  + pert_scale*y(ind))
            drhonu=drhonu+ tmp/v
            fnu=fnu+nu_int_kernel(iq)*(y(EV%r_ix+1)+ pert_scale*y(ind+1))
            if (present(dpnu)) then
                dpnu=dpnu+ tmp*v
                pinu = pinu+ nu_int_kernel(iq)*(y(EV%r_ix+2)+ pert_scale*y(ind+2))*v
            end if
        end do
    end associate
    if (present(dpnu)) then
        dpnu = dpnu/3
    end if

    end subroutine Nu_Integrate_L012

    subroutine Nu_pinudot(EV,y, ydot, a,adotoa, nu_i,pinudot)
    type(EvolutionVars) EV
    integer, intent(in) :: nu_i
    real(dl), intent(in) :: a,adotoa, y(EV%nvar), ydot(EV%nvar)

    !  Compute the time derivative of the mean density in massive neutrinos
    !  and the shear perturbation.
    real(dl) pinudot
    real(dl) aq,q,v,aqdot,vdot
    real(dl) psi2,psi2dot
    real(dl) am, pert_scale
    integer iq,ind

    !  q is the comoving momentum in units of k_B*T_nu0/c.
    pinudot=0._dl
    ind=EV%nu_ix(nu_i)+2
    am=a*State%nu_masses(nu_i)
    do iq=1,EV%nq(nu_i)
        q=State%NuPerturbations%nu_q(iq)
        aq=am/q
        aqdot=aq*adotoa
        v=1._dl/sqrt(1._dl+aq*aq)
        vdot=-aq*aqdot/(1._dl+aq*aq)**1.5d0
        pinudot=pinudot+State%NuPerturbations%nu_int_kernel(iq)*(ydot(ind)*v+y(ind)*vdot)
        ind=ind+EV%lmaxnu_tau(nu_i)+1
    end do
    ind = EV%nu_pert_ix+2
    do iq=EV%nq(nu_i)+1,State%NuPerturbations%nqmax
        q=State%NuPerturbations%nu_q(iq)
        aq=am/q
        aqdot=aq*adotoa
        pert_scale=(State%nu_masses(nu_i)/q)**2/2
        v=1._dl/sqrt(1._dl+aq*aq)
        vdot=-aq*aqdot/(1._dl+aq*aq)**1.5d0
        psi2dot=ydot(EV%r_ix+2)  + pert_scale*ydot(ind)
        psi2=y(EV%r_ix+2)  + pert_scale*y(ind)
        pinudot=pinudot+State%NuPerturbations%nu_int_kernel(iq)*(psi2dot*v+psi2*vdot)
    end do

    end subroutine Nu_pinudot

    !cccccccccccccccccccccccccccccccccccccccccccccc
    subroutine Nu_Intvsq(EV,y, a, nu_i, G11,G30)
    type(EvolutionVars) EV
    integer, intent(in) :: nu_i
    real(dl), intent(in) :: a, y(EV%nvar)
    real(dl), intent(OUT) ::  G11,G30

    !  Compute the third order variables (in velocity dispersion)
    !by integrating over momentum.
    real(dl) aq,q,v, am
    integer iq, ind

    !  q is the comoving momentum in units of k_B*T_nu0/c.
    am=a*State%nu_masses(nu_i)
    ind=EV%nu_ix(nu_i)
    G11=0._dl
    G30=0._dl
    if (EV%nq(nu_i)/=State%NuPerturbations%nqmax) call MpiStop('Nu_Intvsq nq/=nqmax0')
    do iq=1, EV%nq(nu_i)
        q=State%NuPerturbations%nu_q(iq)
        aq=am/q
        v=1._dl/sqrt(1._dl+aq*aq)
        G11=G11+State%NuPerturbations%nu_int_kernel(iq)*y(ind+1)*v**2
        if (EV%lmaxnu_tau(nu_i)>2) then
            G30=G30+State%NuPerturbations%nu_int_kernel(iq)*y(ind+3)*v**2
        end if
        ind = ind+EV%lmaxnu_tau(nu_i)+1
    end do

    end subroutine Nu_Intvsq


    subroutine MassiveNuVars(EV,y,a,gpres,dgrho,dgq, wnu_arr)
    implicit none
    type(EvolutionVars) EV
    real(dl) :: y(EV%nvar), a, gpres,dgrho,dgq
    real(dl), intent(out), optional :: wnu_arr(max_nu)
    !gpres = a^2 kappa p
    !dgrho = a^2 kappa \delta\rho
    !dgp =  a^2 kappa \delta p
    !dgq = a^2 kappa q (heat flux)
    integer nu_i
    real(dl) grhormass_t, rhonu, qnu, clxnu, grhonu_t, gpnu_t, pnu

    do nu_i = 1, CP%Nu_mass_eigenstates
        grhormass_t=State%grhormass(nu_i)/a**2

        !Get density and pressure as ratio to massless by interpolation from table
        call ThermalNuBack%rho_P(a*State%nu_masses(nu_i),rhonu,pnu)

        if (EV%MassiveNuApprox(nu_i)) then
            clxnu=y(EV%nu_ix(nu_i))
            qnu=y(EV%nu_ix(nu_i)+2)
        else
            !Integrate over q
            call Nu_Integrate_L012(EV, y, a, nu_i, clxnu,qnu)
            !clxnu_here  = rhonu*clxnu, qnu_here = qnu*rhonu
            qnu=qnu/rhonu
            clxnu = clxnu/rhonu
        endif

        grhonu_t=grhormass_t*rhonu
        gpnu_t=grhormass_t*pnu

        gpres= gpres + gpnu_t
        dgrho= dgrho + grhonu_t*clxnu
        dgq  = dgq   + grhonu_t*qnu

        if (present(wnu_arr)) then
            wnu_arr(nu_i) =pnu/rhonu
        end if
    end do

    end subroutine MassiveNuVars


    subroutine output(EV, y, j, tau,sources)
    type(EvolutionVars) EV
    real(dl) y(EV%nvar), yprime(EV%nvar)
    integer, intent(in) :: j
    real(dl) tau
    real(dl), target :: sources(:)

    yprime = 0
    EV%OutputSources => Sources
    EV%OutputStep = j
    call derivs(EV,EV%ScalEqsToPropagate,tau,y,yprime)
    nullify(EV%OutputSources)

    end subroutine output


    subroutine initial(EV,y, tau)
    !  Scalar initial conditions.
    implicit none

    type(EvolutionVars) EV
    real(dl) y(EV%nvar)
    real(dl) Rp15,tau,x,x2,x3,om,omtau, &
        Rv,grhonu,chi
    real(dl) k,k2
    real(dl) a,a2, iqg, rhomass,a_massive, ep
    integer l,i, nu_i, j, ind
    integer, parameter :: i_clxg=1,i_clxr=2,i_clxb=3, &
        i_qg=4,i_qr=5,i_vb=6,i_pir=7,i_eta=8,i_aj3r=9
    integer, parameter :: i_max = i_aj3r
    real(dl) initv(1:i_max), initvec(1:i_max)

    nullify(EV%OutputSources)

    EV%k_buf=EV%q
    EV%k2_buf=EV%q2
    EV%Kf(1:EV%MaxlNeeded)=1._dl

    k=EV%k_buf
    k2=EV%k2_buf

    do j=1,EV%MaxlNeeded
        EV%denlk(j)=denl(j)*k*j
        EV%denlk2(j)=denl(j)*k*EV%Kf(j)*(j+1)
        EV%polfack(j)=polfac(j)*k*EV%Kf(j)*denl(j)
    end do

    !Get time to switch off tight coupling
    !The numbers here are a bit of guesswork
    !The high k increase saves time for very small loss of accuracy
    !The lower k ones are more delicate. Nead to avoid instabilities at same time
    !as ensuring tight coupling is accurate enough
    if (EV%k_buf > epsw) then
        if (EV%k_buf > epsw*5) then
            ep=ep0*5/CP%Accuracy%AccuracyBoost*0.65
        else
            ep=ep0
        end if
    else
        ep=ep0
    end if
    if (second_order_tightcoupling) ep=ep*2
    EV%TightSwitchoffTime = min(EV%ThermoData%tight_tau, EV%ThermoData%OpacityToTime(EV%k_buf/ep))

    y=0

    !  k*tau, (k*tau)**2, (k*tau)**3
    x=k*tau
    x2=x*x
    x3=x2*x
    rhomass =  sum(State%grhormass(1:CP%Nu_mass_eigenstates))
    grhonu=rhomass+State%grhornomass

    om = (State%grhob+State%grhoc)/sqrt(3*(State%grhog+grhonu))
    omtau=om*tau
    Rv=grhonu/(grhonu+State%grhog)

    Rp15=4*Rv+15

    a=tau*State%adotrad*(1+omtau/4)
    a2=a*a

    initv=0

    !  Set adiabatic initial conditions

    chi=1  !Get transfer function for chi
    initv(i_clxg)=-chi*EV%Kf(1)/3*x2*(1-omtau/5)
    initv(i_clxr)= initv(i_clxg)
    initv(i_clxb)=0.75_dl*initv(i_clxg)

    initv(i_qg)=initv(i_clxg)*x/9._dl
    initv(i_qr)=-chi*EV%Kf(1)*(4*Rv+23)/Rp15*x3/27
    initv(i_vb)=0.75_dl*initv(i_qg)
    initv(i_pir)=chi*4._dl/3*x2/Rp15*(1+omtau/4*(4*Rv-5)/(2*Rv+15))
    initv(i_aj3r)=chi*4/21._dl/Rp15*x3
    initv(i_eta)=-chi*2*EV%Kf(1)*(1 - x2/12*(-10._dl/Rp15 + EV%Kf(1)))

    InitVec = initv(:)
    InitVec = -InitVec
    !So we start with chi=-1 as before

    y(ix_etak)= -InitVec(i_eta)*k/2
    !get eta_s*k, where eta_s is synchronous gauge variable

    y(ix_clxc) = -initv(i_clxb)

    !  Dark Field
    y(ix_df) = -chi*EV%Kf(1)/2*x2
    y(ix_df + 1) = 0._dl
    y(ix_df + 2) = 0._dl
    y(ix_df + 3) = 0._dl
    y(ix_df + 4) = 0._dl
    y(ix_df + 5) = 0._dl
    y(ix_df + 6) = 0._dl
    y(ix_df + 7) = 0._dl
    y(ix_df + 8) = 0._dl

    !  Baryons
    y(ix_clxb)=InitVec(i_clxb)
    y(ix_vb)=InitVec(i_vb)

    !  Photons
    y(EV%g_ix)=InitVec(i_clxg)
    y(EV%g_ix+1)=InitVec(i_qg)

    !  Neutrinos
    y(EV%r_ix)=InitVec(i_clxr)
    y(EV%r_ix+1)=InitVec(i_qr)
    y(EV%r_ix+2)=InitVec(i_pir)

    if (EV%lmaxnr>2) then
        y(EV%r_ix+3)=InitVec(i_aj3r)
    endif

    if (CP%Num_Nu_massive == 0) return

    do nu_i = 1, CP%Nu_mass_eigenstates
        EV%MassiveNuApproxTime(nu_i) = Nu_tau_massive(nu_i)
        a_massive =  20000*k/State%nu_masses(nu_i)*CP%Accuracy%AccuracyBoost*CP%Accuracy%lAccuracyBoost
        if (a_massive >=0.99) then
            EV%MassiveNuApproxTime(nu_i)=State%tau0+1
        else if (a_massive > 17.d0/State%nu_masses(nu_i)*CP%Accuracy%AccuracyBoost) then
            EV%MassiveNuApproxTime(nu_i)=max(EV%MassiveNuApproxTime(nu_i),State%DeltaTime(0._dl,a_massive, 0.01_dl))
        end if
        ind = EV%nu_ix(nu_i)
        do  i=1,EV%nq(nu_i)
            y(ind:ind+2)=y(EV%r_ix:EV%r_ix+2)
            if (EV%lmaxnu_tau(nu_i)>2) y(ind+3)=InitVec(i_aj3r)
            ind = ind + EV%lmaxnu_tau(nu_i)+1
        end do
    end do

    end subroutine initial


    subroutine derivs(EV0,n,tau,ay,ayprime)
    !  Evaluate the time derivatives of the scalar perturbations
    use MassiveNu
    use Recombination
    implicit none
    class(EvolutionVars), pointer :: EV
    class(TCambComponent), target :: EV0
    integer n,nu_i
    real(dl), target :: ay(n),ayprime(n)
    real(dl) tau, w
    real(dl) k,k2
    real(dl) opacity
    real(dl) photbar,cs2,pb43,grho,slip,clxgdot, &
        clxcdot,clxbdot,adotdota,gpres,clxrdot,etak
    real(dl) q,aq,v
    real(dl) G11_t,G30_t, wnu_arr(max_nu)

    real(dl) dgq,grhob_t,grhor_t,grhog_t,sigma,polter
    real(dl) qgdot,qrdot,pigdot,pirdot,vbdot,dgrho,adotoa
    real(dl) a,a2,z,clxc,clxb,vb,clxg,qg,pig,clxr,qr,pir
    real(dl) E2, dopacity
    integer l,i,ind, ind2, off_ix, ix
    real(dl) dgs,sigmadot,dz
    real(dl) dgpi, clxnu, gpres_nu
    !non-_flat vars
    real(dl) cothxor !1/tau in _flat case
    real(dl) xe,Trad, Delta_TM, Tmat, Delta_TCMB
    real(dl) delta_p_b, wing_t, wing2_t,winv_t
    real(dl) Delta_source2, polter_line
    real(dl) Delta_xe, Tspin, tau_eps, tau_fac, Tb
    integer lineoff,lineoffpol
    !Variables for source calculation
    real(dl) diff_rhopi, pidot_sum, dgpi_diff, phi
    real(dl) E(2:3), Edot(2:3)
    real(dl) phidot, polterdot, polterddot, octg, octgdot
    real(dl) ddopacity, visibility, dvisibility, ddvisibility, exptau, lenswindow
    real(dl) ISW, quadrupole_source, doppler, monopole_source, tau0, ang_dist
    real(dl) dgrho_de, dgq_de, cs2_de
    real(dl) W_,N_,dB,B2,W2,N2,kN,BW,BN,kB,k2N,k2N2

    select type(EV0)
    class is (EvolutionVars)
        EV => EV0
    end select

    k=EV%k_buf
    k2=EV%k2_buf

    !  Get background scale factor, sound speed and ionisation fraction.
    if (EV%TightCoupling) then
        call EV%ThermoData%Values(tau,a,cs2,opacity,dopacity)
    else
        call EV%ThermoData%Values(tau,a,cs2,opacity)
    end if
    a2=a*a

    etak=ay(ix_etak)

    !  CDM variables
    clxc=ay(ix_clxc)

    !  Baryon variables
    clxb=ay(ix_clxb)
    vb=ay(ix_vb)
    !  Compute expansion rate from: grho 8*pi*rho*a**2

    grhob_t=State%grhob/a
    grhor_t=State%grhornomass/a2
    grhog_t=State%grhog/a2

    !total perturbations: matter terms first, then add massive nu, de and radiation
    !  8*pi*G*(a**2)*SUM[rho_i*clx_i]
    dgrho = grhob_t*clxb + (State%grhoc/a)*clxc
    !  8*pi*G*(a**2)*SUM[(rho_i+P_i)*theta_i]/k
    dgq = 0._dl

    gpres_nu=0

    if (State%CP%Num_Nu_Massive > 0) then
        call MassiveNuVars(EV,ay,a,gpres_nu,dgrho,dgq, wnu_arr)
    end if

    adotoa = State%CP%DarkField%splines%dadtau%GetValue(a)/a

    cothxor=1._dl/tau

    B2 = adotoa*adotoa

    !adotoa=sqrt(grho/3)
    grho = 3*B2

    W_ = State%CP%DarkField%splines%W%GetValue(a)
    N_ = State%CP%DarkField%splines%N%GetValue(a)

    N2 = N_*N_
    BN = adotoa*N_
    W2 = W_*W_
    kN = k*N_
    BW = adotoa*W_
    kB = k*adotoa
    k2N = k2*N_
    k2N2 = k2*N2

    ! 8*pi*G*(a**2)*P
    gpres = gpres_nu + (grhor_t + grhog_t)/3 - (1200*B2*N2 + (26*W_ - 320*BN)*W_)

    dB = -(gpres + B2)/(2._dl + 160*N2)

    ! Recalculate gpres
    gpres = -2*dB - B2

    if (EV%no_nu_multpoles) then
        !RSA approximation of arXiv:1104.2933, dropping opactity terms in the velocity
        !Approximate total density variables with just matter terms
        z=(0.5_dl*dgrho/k + etak)/adotoa
        dz= -adotoa*z - 0.5_dl*dgrho/k
        clxr=-4*dz/k
        qr=-4._dl/3*z
        pir=0
    else
        !  Massless neutrinos
        clxr=ay(EV%r_ix)
        qr  =ay(EV%r_ix+1)
        pir =ay(EV%r_ix+2)
    endif

    pig=0
    if (EV%no_phot_multpoles) then
        if (.not. EV%no_nu_multpoles) then
            z=(0.5_dl*dgrho/k + etak)/adotoa
            dz= -adotoa*z - 0.5_dl*dgrho/k
            clxg=-4*dz/k-4/k*opacity*(vb+z)
            qg=-4._dl/3*z
        else
            clxg=clxr-4/k*opacity*(vb+z)
            qg=qr
        end if
    else
        !  Photons
        clxg=ay(EV%g_ix)
        qg=ay(EV%g_ix+1)
        if (.not. EV%TightCoupling) pig=ay(EV%g_ix+2)
    end if

    !  Photon mass density over baryon mass density
    photbar=grhog_t/grhob_t
    pb43=4._dl/3*photbar

    !  8*pi*G*(a**2)*SUM[rho_i*clx_i] - radiation terms
    dgrho = dgrho + grhog_t*clxg + grhor_t*clxr + State%CP%DarkField%a2X(ix_df,ix_etak,ay,k,k2,W_,N_,B2,W2,N2,kN,BW,BN,k2N,k2N2)

    !  Get sigma (shear) and z from the constraints
    ! have to get z from eta for numerical stability
    z=(0.5_dl*dgrho/k + etak)/adotoa

    ! Recalculate z
    z = z / (1._dl  - State%CP%DarkField%p(79)*N2)

    ayprime(ix_df) = 2*k*z

    ! Recalculate dgrho
    dgrho = 2*k*(adotoa*z - etak)

    !  8*pi*G*(a**2)*SUM[(rho_i+P_i)*theta_i]/k
    dgq = dgq + grhob_t*vb + grhog_t*qg + grhor_t*qr - k*State%CP%DarkField%ia2M(ix_df,ix_etak,ay,ayprime,k,W_,N_,N2,BN,kN)

    ! Recalculate dgq
    dgq = dgq / (1._dl  + 0.5_dl*State%CP%DarkField%p(95)*N2)

    ayprime(ix_etak)=0.5_dl*dgq

    !eta*k equation
    sigma=(z+1.5_dl*dgq/k2)

    clxcdot=-k*z
    ayprime(ix_clxc) = clxcdot

    !  Dark Field equation of motion
    call State%CP%DarkField%dfDeriv(ix_df,ix_etak,ay,ayprime,k,k2,W_,N_,adotoa,dB,B2,BW,BN,kN,kB,k2N)

    !  Baryon equation of motion.
    clxbdot=-k*(z+vb)
    ayprime(ix_clxb)=clxbdot
    !  Photon equation of motion
    clxgdot=-k*(4._dl/3._dl*z+qg)

    !Sources
    Delta_TM = clxg/4
    delta_p_b = cs2*clxb

    Delta_xe = 0

    ! Easy to see instability in k \sim 1e-3 by tracking evolution of vb

    !  Use explicit equation for vb if appropriate

    if (EV%TightCoupling) then
        !  ddota/a
        adotdota=(adotoa*adotoa-gpres)/2

        pig = 32._dl/45/opacity*k*(sigma+vb)

        !  First-order approximation to baryon-photon splip
        slip = - (2*adotoa/(1+pb43) + dopacity/opacity)* (vb-3._dl/4*qg) &
            +(-adotdota*vb-k/2*adotoa*clxg +k*(cs2*clxbdot-clxgdot/4))/(opacity*(1+pb43))

        if (second_order_tightcoupling) then
            ! by Francis-Yan Cyr-Racine simplified (inconsistently) by AL assuming _flat
            !AL: First order slip seems to be fine here to 2e-4

            !  8*pi*G*a*a*SUM[rho_i*sigma_i]
            dgs = grhog_t*pig+grhor_t*pir

            ! Define shear derivative to first order
            sigmadot = -2*adotoa*sigma-dgs/k+etak

            !Once know slip, recompute qgdot, pig, pigdot
            qgdot = k*(clxg/4._dl-pig/2._dl) +opacity*slip

            pig = 32._dl/45/opacity*k*(sigma+3._dl*qg/4._dl)*(1+(dopacity*11._dl/6._dl/opacity**2)) &
                + (32._dl/45._dl/opacity**2)*k*(sigmadot+3._dl*qgdot/4._dl)*(-11._dl/6._dl)

            pigdot = -(32._dl/45._dl)*(dopacity/opacity**2)*k*(sigma+3._dl*qg/4._dl)*(1 + &
                dopacity*11._dl/6._dl/opacity**2 ) &
                + (32._dl/45._dl/opacity)*k*(sigmadot+3._dl*qgdot/4._dl)*(1+(11._dl/6._dl) &
                *(dopacity/opacity**2))

            EV%pigdot = pigdot

        end if

        !  Use tight-coupling approximation for vb
        !  zeroth order approximation to vbdot + the pig term
        vbdot=(-adotoa*vb+cs2*k*clxb + k/4*pb43*(clxg-2*EV%Kf(1)*pig))/(1+pb43)

        vbdot=vbdot+pb43/(1+pb43)*slip
        EV%pig = pig

    else
        vbdot=-adotoa*vb+k*delta_p_b-photbar*opacity*(4._dl/3*vb-qg)
    end if

    ayprime(ix_vb)=vbdot

    if (.not. EV%no_phot_multpoles) then
        !  Photon equations of motion
        ayprime(EV%g_ix)=clxgdot
        qgdot=4._dl/3*(-vbdot-adotoa*vb+delta_p_b*k)/pb43 &
            +EV%denlk(1)*clxg-EV%denlk2(1)*pig
        ayprime(EV%g_ix+1)=qgdot

        !  Use explicit equations for photon moments if appropriate
        if (.not. EV%tightcoupling) then
            E2=ay(EV%polind+2)
            polter = pig/10+9._dl/15*E2 !2/15*(3/4 pig + 9/2 E2)
            ix= EV%g_ix+2
            if (EV%lmaxg>2) then
                pigdot=EV%denlk(2)*qg-EV%denlk2(2)*ay(ix+1)-opacity*(pig - polter) &
                    +8._dl/15._dl*k*sigma
                ayprime(ix)=pigdot
                do  l=3,EV%lmaxg-1
                    ix=ix+1
                    ayprime(ix)=(EV%denlk(l)*ay(ix-1)-EV%denlk2(l)*ay(ix+1))-opacity*ay(ix)
                end do
                ix=ix+1
                !  Truncate the photon moment expansion
                ayprime(ix)=k*ay(ix-1)-(EV%lmaxg+1)*cothxor*ay(ix) -opacity*ay(ix)
            else !_closed case
                pigdot=EV%denlk(2)*qg-opacity*(pig - polter) +8._dl/15._dl*k*sigma
                ayprime(ix)=pigdot
            endif
            !  Polarization
            !l=2
            ix=EV%polind+2
            if (EV%lmaxgpol>2) then
                ayprime(ix) = -opacity*(ay(ix) - polter) - k/3._dl*ay(ix+1)
                do l=3,EV%lmaxgpol-1
                    ix=ix+1
                    ayprime(ix)=-opacity*ay(ix) + (EV%denlk(l)*ay(ix-1)-EV%polfack(l)*ay(ix+1))
                end do
                ix=ix+1
                !truncate
                ayprime(ix)=-opacity*ay(ix) + &
                    k*EV%poltruncfac*ay(ix-1)-(EV%lmaxgpol+3)*cothxor*ay(ix)
            else !_closed case
                ayprime(ix) = -opacity*(ay(ix) - polter)
            endif
        end if
    end if

    if (.not. EV%no_nu_multpoles) then
        !  Massless neutrino equations of motion.
        clxrdot=-k*(4._dl/3._dl*z+qr)
        ayprime(EV%r_ix)=clxrdot
        qrdot=EV%denlk(1)*clxr-EV%denlk2(1)*pir
        ayprime(EV%r_ix+1)=qrdot
        if (EV%high_ktau_neutrino_approx) then
            !ufa approximation for k*tau>>1, more accurate when there are reflections from lmax
            !Method from arXiv:1104.2933
            !                if (.not. EV%TightCoupling) then
            !                 gpres=gpres+ (grhog_t+grhor_t)/3 +grhov_t*w_lam
            !                 adotdota=(adotoa*adotoa-gpres)/2
            !                end if
            !                ddz=(2*adotoa**2 - adotdota)*z  &
            !                  + adotoa/(2*k)*( 6*(grhog_t*clxg+grhor_t*clxr) + 2*((_grhoc/a)*clxc+grhob_t*clxb) ) &
            !                   - 1._dl/(2*k)*( 2*(grhog_t*clxgdot+grhor_t*clxrdot) + (_grhoc/a)*clxcdot + grhob_t*clxbdot )
            !                dz= -adotoa*z - 0.5_dl*dgrho/k
            !                pirdot= -3*pir*cothxor + k*(qr+4._dl/3*z)
            pirdot= -3*pir*cothxor - clxrdot
            ayprime(EV%r_ix+2)=pirdot

            !                pirdot=k*(0.4_dl*qr-0.6_dl*ay(EV%lmaxg+10)+8._dl/15._dl*sigma)
            !                ayprime(EV%lmaxg+9)=pirdot
            !                ayprime(3+EV%lmaxg+7)=k*ay(3+EV%lmaxg+6)- &
            !                                      (3+1)*cothxor*ay(3+EV%lmaxg+7)
            !               ayprime(3+EV%lmaxg+7+1:EV%lmaxnr+EV%lmaxg+7)=0
        else
            ix=EV%r_ix+2
            if (EV%lmaxnr>2) then
                pirdot=EV%denlk(2)*qr- EV%denlk2(2)*ay(ix+1)+8._dl/15._dl*k*sigma
                ayprime(ix)=pirdot
                do l=3,EV%lmaxnr-1
                    ix=ix+1
                    ayprime(ix)=(EV%denlk(l)*ay(ix-1) - EV%denlk2(l)*ay(ix+1))
                end do
                !  Truncate the neutrino expansion
                ix=ix+1
                ayprime(ix)=k*ay(ix-1)- (EV%lmaxnr+1)*cothxor*ay(ix)
            else
                pirdot=EV%denlk(2)*qr +8._dl/15._dl*k*sigma
                ayprime(ix)=pirdot
            end if
        end if
    end if ! no_nu_multpoles

    !  Massive neutrino equations of motion.
    if (State%CP%Num_Nu_massive >0) then
        !DIR$ LOOP COUNT MIN(1), AVG(1)
        do nu_i = 1, State%CP%Nu_mass_eigenstates
            if (EV%MassiveNuApprox(nu_i)) then
                !Now EV%iq0 = clx, EV%iq0+1 = clxp, EV%iq0+2 = G_1, EV%iq0+3=G_2=pinu
                !see astro-ph/0203507
                G11_t=EV%G11(nu_i)/a/a2
                G30_t=EV%G30(nu_i)/a/a2
                off_ix = EV%nu_ix(nu_i)
                w=wnu_arr(nu_i)
                ayprime(off_ix)=-k*z*(w+1) + 3*adotoa*(w*ay(off_ix) - ay(off_ix+1))-k*ay(off_ix+2)
                ayprime(off_ix+1)=(3*w-2)*adotoa*ay(off_ix+1) - 5._dl/3*k*z*w - k/3*G11_t
                ayprime(off_ix+2)=(3*w-1)*adotoa*ay(off_ix+2) - k*(2._dl/3*EV%Kf(1)*ay(off_ix+3)-ay(off_ix+1))
                ayprime(off_ix+3)=(3*w-2)*adotoa*ay(off_ix+3) + 2*w*k*sigma - k/5*(3*EV%Kf(2)*G30_t-2*G11_t)
            else
                ind=EV%nu_ix(nu_i)
                !DIR$ LOOP COUNT MIN(3), AVG(3)
                do i=1,EV%nq(nu_i)
                    q=State%NuPerturbations%nu_q(i)
                    aq=a*State%nu_masses(nu_i)/q
                    v=1._dl/sqrt(1._dl+aq*aq)

                    ayprime(ind)=-k*(4._dl/3._dl*z + v*ay(ind+1))
                    ind=ind+1
                    ayprime(ind)=v*(EV%denlk(1)*ay(ind-1)-EV%denlk2(1)*ay(ind+1))
                    ind=ind+1
                    if (EV%lmaxnu_tau(nu_i)==2) then
                        ayprime(ind)=-ayprime(ind-2) -3*cothxor*ay(ind)
                    else
                        ayprime(ind)=v*(EV%denlk(2)*ay(ind-1)-EV%denlk2(2)*ay(ind+1)) &
                            +k*8._dl/15._dl*sigma
                        do l=3,EV%lmaxnu_tau(nu_i)-1
                            ind=ind+1
                            ayprime(ind)=v*(EV%denlk(l)*ay(ind-1)-EV%denlk2(l)*ay(ind+1))
                        end do
                        !  Truncate moment expansion.
                        ind = ind+1
                        ayprime(ind)=k*v*ay(ind-1)-(EV%lmaxnu_tau(nu_i)+1)*cothxor*ay(ind)
                    end if
                    ind = ind+1
                end do
            end if
        end do

        if (EV%has_nu_relativistic) then
            ind=EV%nu_pert_ix
            ayprime(ind)=+k*a2*qr -k*ay(ind+1)
            ind2= EV%r_ix
            do l=1,EV%lmaxnu_pert-1
                ind=ind+1
                ind2=ind2+1
                ayprime(ind)= -a2*(EV%denlk(l)*ay(ind2-1)-EV%denlk2(l)*ay(ind2+1)) &
                    +   (EV%denlk(l)*ay(ind-1)-EV%denlk2(l)*ay(ind+1))
            end do
            ind=ind+1
            ind2=ind2+1
            ayprime(ind)= k*(ay(ind-1) -a2*ay(ind2-1)) -(EV%lmaxnu_pert+1)*cothxor*ay(ind)
        end if
    end if

    if (associated(EV%OutputSources)) then
        if (EV%TightCoupling .or. EV%no_phot_multpoles) then
            E=0
            Edot=0
        else
            E = ay(EV%polind+2:EV%polind+3)
            Edot = ayprime(EV%polind+2:EV%polind+3)
        end if
        if (EV%no_nu_multpoles) then
            pirdot=0
            qrdot = -4*dz/3
        end if
        if (EV%no_phot_multpoles) then
            pigdot=0
            octg=0
            octgdot=0
            qgdot = -4*dz/3
        else
            if (EV%TightCoupling) then
                if (second_order_tightcoupling) then
                    octg = (3._dl/7._dl)*pig*(EV%k_buf/opacity)
                    E(2) = pig/4 + pigdot*(1._dl/opacity)*(-5._dl/8._dl)
                    E(3) = (3._dl/7._dl)*(EV%k_buf/opacity)*E(2)
                    Edot(2)= (pigdot/4._dl)*(1+(5._dl/2._dl)*(dopacity/opacity**2))
                else
                    pigdot = -dopacity/opacity*pig + 32._dl/45*k/opacity*(-2*adotoa*sigma  &
                        +etak/EV%Kf(1)-  dgpi/k +vbdot )
                    Edot(2) = pigdot/4
                    E(2) = pig/4
                    octg=0
                end if
                octgdot=0
            else
                octg=ay(EV%g_ix+3)
                octgdot=ayprime(EV%g_ix+3)
            end if
        end if
        dgrho_de=0
        dgq_de=0

        dgpi  = grhor_t*pir + grhog_t*pig
        dgpi_diff = 0  !sum (3*p_nu -rho_nu)*pi_nu
        pidot_sum = grhog_t*pigdot + grhor_t*pirdot
        clxnu =0
        if (State%CP%Num_Nu_Massive /= 0) then
            call MassiveNuVarsOut(EV,ay,ayprime,a, adotoa, dgpi=dgpi, clxnu_all=clxnu, &
                dgpi_diff=dgpi_diff, pidot_sum=pidot_sum)
        end if
        diff_rhopi = pidot_sum - (4*dgpi+ dgpi_diff)*adotoa
        phi = -((dgrho +3*dgq*adotoa/k)/EV%Kf(1) + dgpi)/(2*k2)

        if (associated(EV%OutputSources)) then

            EV%OutputSources = 0
            call EV%ThermoData%IonizationFunctionsAtTime(tau, a, opacity, dopacity, ddopacity, &
                visibility, dvisibility, ddvisibility, exptau, lenswindow)

            tau0 = State%tau0
            phidot = (1.0d0/2.0d0)*(adotoa*(-dgpi - 2*k2*phi) + dgq*k - &
                diff_rhopi+ k*sigma*(gpres + grho))/k2
            !time derivative of shear
            sigmadot = -adotoa*sigma - 1.0d0/2.0d0*dgpi/k + k*phi
            !quadrupole source derivatives; polter = pi_g/10 + 3/5 E_2
            polter = pig/10+9._dl/15*E(2)
            polterdot = (1.0d0/10.0d0)*pigdot + (3.0d0/5.0d0)*Edot(2)
            polterddot = -2.0d0/25.0d0*adotoa*dgq/(k*EV%Kf(1)) - 4.0d0/75.0d0*adotoa* &
                k*sigma - 4.0d0/75.0d0*dgpi - 2.0d0/75.0d0*dgrho/EV%Kf(1) - 3.0d0/ &
                50.0d0*k*octgdot*EV%Kf(2) + (1.0d0/25.0d0)*k*qgdot - 1.0d0/5.0d0 &
                *k*EV%Kf(2)*Edot(3) + (-1.0d0/10.0d0*pig + (7.0d0/10.0d0)* &
                polter - 3.0d0/5.0d0*E(2))*dopacity + (-1.0d0/10.0d0*pigdot &
                + (7.0d0/10.0d0)*polterdot - 3.0d0/5.0d0*Edot(2))*opacity
            !Temperature source terms, after integrating by parts in conformal time

            !2phi' term (\phi' + \psi' in Newtonian gauge), phi is the Weyl potential
            ISW = 2*phidot*exptau
            monopole_source =  (-etak/(k*EV%Kf(1)) + 2*phi + clxg/4)*visibility
            doppler = ((sigma + vb)*dvisibility + (sigmadot + vbdot)*visibility)/k
            quadrupole_source = (5.0d0/8.0d0)*(3*polter*ddvisibility + 6*polterdot*dvisibility &
                + (k**2*polter + 3*polterddot)*visibility)/k**2

            EV%OutputSources(1) = ISW + doppler + monopole_source + quadrupole_source
            ang_dist = f_K(tau0-tau)
            if (tau < tau0) then
                !E polarization source
                EV%OutputSources(2)=visibility*polter*(15._dl/8._dl)/(ang_dist**2*k2)
                !factor of four because no 1/16 later
            end if

            if (size(EV%OutputSources) > 2) then
                !Get lensing sources
                if (tau>State%tau_maxvis .and. tau0-tau > 0.1_dl) then
                    EV%OutputSources(3) = -2*phi*f_K(tau-State%tau_maxvis)/(f_K(tau0-State%tau_maxvis)*ang_dist)
                    !We include the lensing factor of two here
                end if
            end if
        end if
    end if

    end subroutine derivs


    end module GaugeInterface
