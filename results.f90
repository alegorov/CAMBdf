    ! Modules used by cmbmain and other routines.

    !     Code for Anisotropies in the Microwave Background
    !     by Antony Lewis (http://cosmologist.info) and Anthony Challinor
    !     See readme.html for documentation.
    !
    !     Based on CMBFAST  by  Uros Seljak and Matias Zaldarriaga, itself based
    !     on Boltzmann code written by Edmund Bertschinger, Chung-Pei Ma and Paul Bode.
    !     Original CMBFAST copyright and disclaimer:
    !
    !     Copyright 1996 by Harvard-Smithsonian Center for Astrophysics and
    !     the Massachusetts Institute of Technology.  All rights reserved.
    !
    !     THIS SOFTWARE IS PROVIDED "AS IS", AND M.I.T. OR C.f.A. MAKE NO
    !     REPRESENTATIONS OR WARRANTIES, EXPRESS OR IMPLIED.
    !     By way of example, but not limitation,
    !     M.I.T. AND C.f.A MAKE NO REPRESENTATIONS OR WARRANTIES OF
    !     MERCHANTABILITY OR FITNESS FOR ANY PARTICULAR PURPOSE OR THAT
    !     THE USE OF THE LICENSED SOFTWARE OR DOCUMENTATION WILL NOT INFRINGE
    !     ANY THIRD PARTY PATENTS, COPYRIGHTS, TRADEMARKS OR OTHER RIGHTS.
    !
    !     portions of this software are based on the COSMICS package of
    !     E. Bertschinger.  See the LICENSE file of the COSMICS distribution
    !     for restrictions on the modification and distribution of this software.

    module results
    use constants, only : const_pi, const_twopi
    use MiscUtils
    use RangeUtils
    use StringUtils
    use MathUtils
    use config
    use model
    use splines
    implicit none
    public

    Type TBackgroundOutputs
        real(dl), allocatable :: H(:), DA(:), rs_by_D_v(:)
    end Type TBackgroundOutputs

    integer, parameter :: derived_age=1, derived_zstar=2, derived_rstar=3, derived_thetastar=4, derived_DAstar = 5, &
        derived_zdrag=6, derived_rdrag=7,derived_kD=8,derived_thetaD=9, derived_zEQ =10, derived_keq =11, &
        derived_thetaEQ=12, derived_theta_rs_EQ = 13
    integer, parameter :: nthermo_derived = 13

    Type lSamples
        integer :: nl = 0
        integer :: lmin = 2
        integer, allocatable :: l(:)
        logical :: use_spline_template = .true.
    contains
    procedure :: Init => lSamples_init
    procedure :: InterpolateClArr
    procedure :: InterpolateClArrTemplated
    end Type lSamples

    logical, parameter :: dowinlens = .false. !not used, test getting CMB lensing using visibility
    integer, parameter :: thermal_history_def_timesteps = 20000

    Type TThermoData
        logical :: HasThermoData = .false. !Has it been computed yet for current parameters?
        !Background thermal history, interpolated from precomputed tables
        integer :: nthermo !Number of table steps
        !baryon temperature, sound speed, ionization fractions, and opacity
        real(dl), dimension(:), allocatable :: tb, cs2, xe, dotmu
        ! e^(-tau) and derivatives
        real(dl), dimension(:), allocatable :: emmu, dcs2,demmu, ddotmu, dddotmu, ddddotmu
        real(dl), dimension(:), allocatable :: ScaleFactor, dScaleFactor, adot, dadot
        real(dl), dimension(:), allocatable :: winlens, dwinlens
        real(dl) tauminn,dlntau
        real(dl) :: tight_tau, actual_opt_depth
        !Times when 1/(opacity*tau) = 0.01, for use switching tight coupling approximation
        real(dl) :: matter_verydom_tau
        real(dl) :: recombination_saha_tau
        !sound horizon and recombination redshifts
        real(dl) :: r_drag0, z_star, z_drag  !!JH for updated BAO likelihood.
        real(dl), dimension(:), allocatable :: step_redshift, rhos_fac, drhos_fac
        real(dl) :: tau_start_redshiftwindows,tau_end_redshiftwindows
        logical :: has_lensing_windows = .false.
        real(dl) recombination_Tgas_tau
        Type(TCubicSpline) :: ScaleFactorAtTime
        !Mapping between redshift and time
        real(dl), private,dimension(:), allocatable :: redshift_time, dredshift_time
        real(dl), private, dimension(:), allocatable :: arhos_fac, darhos_fac, ddarhos_fac
    contains
    procedure :: Init => Thermo_Init
    procedure :: OpacityToTime => Thermo_OpacityToTime
    procedure :: values => Thermo_values
    procedure :: IonizationFunctionsAtTime
    procedure, private :: SetTimeSteps
    end type TThermoData

    !Sources
    Type CalWins
        real(dl), allocatable :: awin_lens(:), dawin_lens(:)
    end Type CalWins

    Type LimberRec
        integer n1,n2 !corresponding time step array indices
        real(dl), dimension(:), allocatable :: k
        real(dl), dimension(:), allocatable :: Source
    end Type LimberRec

    Type ClTransferData
        !Cl transfer function variables
        !values of q for integration over q to get C_ls
        Type (lSamples) :: ls ! l that are computed
        integer :: NumSources
        !Changes -scalars:  2 for just CMB, 3 for lensing
        !- tensors: T and E and phi (for lensing), and T, E, B respectively

        type (TRanges) :: q
        real(dl), dimension(:,:,:), allocatable :: Delta_p_l_k

        !The L index of the lowest L to use for Limber
        integer, dimension(:), allocatable :: Limber_l_min
        !For each l, the set of k in each limber window
        !indices LimberWindow(SourceNum,l)
        Type(LimberRec), dimension(:,:), allocatable :: Limber_windows

        !The maximum L needed for non-Limber
        integer max_index_nonlimber

    end Type ClTransferData

    type TCLdata
        Type(ClTransferData) :: CTransScal, CTransTens, CTransVec

        real(dl), dimension (:,:), allocatable :: Cl_scalar, Cl_tensor, Cl_vector
        !Indices are Cl_xxx( l , Cl_type)
        !where Cl_type is one of the above constants

        real(dl), dimension (:,:,:), allocatable :: Cl_Scalar_Array
        !Indices are Cl_xxx( l , field1,field2)
        !where ordering of fields is T, E, \psi (CMB lensing potential), window_1, window_2...

        !The following are set only if doing lensing
        integer lmax_lensed !Only accurate to rather less than this
        real(dl) , dimension (:,:), allocatable :: Cl_lensed
        !Cl_lensed(l, Cl_type) are the interpolated Cls
    contains
    procedure :: InitCls => TCLdata_InitCls
    end type TCLdata

    Type TTimeSources
        ! values of q to evolve the propagation equations to compute the sources
        type(TRanges) :: Evolve_q
        real(dl), dimension(:,:,:), allocatable :: LinearSrc !Sources and second derivs
        !LinearSrc indices  ( k_index, source_index, time_step_index )
        integer SourceNum, NonCustomSourceNum
        !SourceNum is total number sources (2 or 3 for scalars, 3 for tensors).
    end type TTimeSources

    type, extends(TCAMBdata) :: CAMBdata

        type(CAMBparams) :: CP

        real(dl) ThermoDerivedParams(nthermo_derived)

        !     grhocrit =kappa*a^2*rho_crit(0)
        !     grhornomass=grhor*number of massless neutrino species
        !     taurst,taurend - time at start/end of recombination
        !     dtaurec - dtau during recombination
        !     adotrad - a(tau) in radiation era
        real(dl) grhocrit,dadtau0,grhog,grhor,grhob,grhoc,grhornomass
        real(dl) taurst,dtaurec,taurend,tau_maxvis,adotrad

        real(dl) tau0,chi0 !time today and rofChi(tau0/_curvature_radius)
        real(dl) scale !relative to _flat. e.g. for scaling lSamp%l sampling.

        real(dl) akthom !sigma_T * (number density of protons now)
        real(dl) fHe !n_He_tot / n_H_tot
        real(dl) Nnow
        real(dl) z_eq !defined assuming all neutrinos massless
        !Neutrinos
        real(dl) grhormass(max_nu)
        !     nu_masses=m_nu*c**2/(k_B*T_nu0)
        real(dl) nu_masses(max_nu)

        logical :: get_growth_sigma8 = .true.
        !gets sigma_vdelta, like sigma8 but using velocity-density cross power,
        !in late LCDM f*sigma8 = sigma_vdelta^2/sigma8

        logical :: needs_good_pk_sampling = .false.

        logical ::call_again = .false.
        !if being called again with same parameters to get different thing

        real(dl) :: reion_tau_start, reion_tau_complete
        integer :: reion_n_steps

        Type(TNuPerturbations) :: NuPerturbations

        Type(TBackgroundOutputs) :: BackgroundOutputs

        !Time steps for sampling sources
        Type(TRanges) :: TimeSteps
        !Background interpolation tables for thermal history etc.
        Type(TThermoData) :: ThermoData

        !Matter transfer data
        Type (MatterTransferData):: MT

        !Matter power spectrum for default variable (used for non-linear corrections)
        Type(MatterPowerData), allocatable :: CAMB_PK

        Type(TClData) :: CLdata

        real(dl), dimension(:), allocatable :: optical_depths_for21cm

        Type(TTimeSources), allocatable :: ScalarTimeSources
        integer :: Scalar_C_last = C_PhiE


    contains
    procedure :: DeltaTime => CAMBdata_DeltaTime
    procedure :: TimeOfz => CAMBdata_TimeOfz
    procedure :: SetParams => CAMBdata_SetParams
    procedure :: Free => CAMBdata_Free
    procedure :: grho_no_de
    procedure :: GetReionizationOptDepth
    procedure :: rofChi
    end type CAMBdata

    interface
    FUNCTION state_function(obj, a)
    use precision
    import
    class(CAMBdata) :: obj
    real(dl), intent(in) :: a
    real(dl) :: state_function
    END FUNCTION  state_function
    end interface

    procedure(obj_function), private :: dtauda

    contains


    subroutine derivaWN_optional(this,nvar,tau,ay,ayprime,calculateWN)
    class(CAMBdata) :: this
    integer :: nvar
    real(dl), target :: ay(nvar),ayprime(nvar)
    logical :: calculateWN
    real(dl) :: tau,a,W,N,B

    a = ay(1)
    W = ay(2)
    N = ay(3)
    
    ayprime(1) = sqrt((this%grho_no_de(a) + 6*(W*W)*(a*a)) / (3._dl + 240*(N*N)))

    if (calculateWN) then
        B = ayprime(1)/a

        ayprime(2) = 6*B*W - 40*(B*B)*N
        ayprime(3) = W - 8*B*N
    end if

    end subroutine derivaWN_optional
    

    subroutine derivaWN(this0,n,tau,ay,ayprime)
    !  Evaluate the time derivatives a, W, N
    class(CAMBdata), pointer :: this
    class(TCambComponent), target :: this0
    integer :: n
    real(dl), target :: ay(n), ayprime(n)
    real(dl) :: tau

    select type(this0)
    class is (CAMBdata)
        this => this0
    end select
    
    call derivaWN_optional(this,n,tau,ay,ayprime,.true.)

    end subroutine derivaWN

    function DarkField_splines(this, W0)
    use DarkFieldModule
    use ArrayUtils
    class(CAMBdata) :: this
    real(dl) :: W0
    type(TDarkFieldSplines), allocatable :: DarkField_splines
    procedure(TClassDverk) :: dverk
    real(dl), allocatable :: Xarr(:), dadtau(:), Warr(:), Narr(:)
    integer :: n, bufSize
    real(dl) :: a
    integer :: ind
    real(dl) :: tau,tol1,tauend
    real(dl), parameter :: max_a = 1.1_dl
    real(dl), parameter :: taustart = 0.001_dl
    real(dl), parameter :: step1 = 0.2_dl
    real(dl), parameter :: step2 = 5._dl
    real(dl), parameter :: sigma = (sqrt(65._dl) - 3._dl) / 2
    integer, parameter :: nvar = 3
    real(dl) :: c(24), w(nvar,9), y(nvar), y1(nvar)

    bufSize = 100000
    
    allocate(Xarr(bufSize))
    allocate(dadtau(bufSize))
    allocate(Warr(bufSize))
    allocate(Narr(bufSize))

    ind = 1
    tol1 = 0.0001_dl
    
    a = 0._dl
    Xarr(1) = a
    dadtau(1) = sqrt(this%grho_no_de(a) / 3)
    Warr(1) = 0._dl
    Narr(1) = 0._dl
    
    tau = taustart
    a = dadtau(1) * tau
    y(1) = a
    y(2) = W0*(tau**sigma)
    y(3) = (W0/(sigma + 9._dl))*(tau**(sigma + 1._dl))

    tauend = step1
    n = 2
    
    do while(.true.)
        if (n > bufSize) then
            bufSize = 2 * bufSize
            call realloc_D(Xarr, bufSize)
            call realloc_D(dadtau, bufSize)
            call realloc_D(Warr, bufSize)
            call realloc_D(Narr, bufSize)
        end if
        
        call dverk(this,nvar,derivaWN,tau,y,tauend,tol1,ind,c,nvar,w)
        call derivaWN_optional(this,nvar,tau,y,y1,.false.)
        a = y(1)

        Xarr(n) = a
        dadtau(n) = y1(1)
        Warr(n) = y(2)
        Narr(n) = y(3)

        if (a >= max_a) exit
        tauend = tauend + step1 + a*(step2 - step1)
        n = n + 1
    end do
    
    allocate(DarkField_splines)
    call TCubicSpline_Init(DarkField_splines%dadtau, Xarr,  dadtau, n)
    call TCubicSpline_Init(DarkField_splines%W, Xarr,  Warr, n)
    call TCubicSpline_Init(DarkField_splines%N, Xarr,  Narr, n)
    deallocate(Xarr)
    deallocate(dadtau)
    deallocate(Warr)
    deallocate(Narr)

    end function DarkField_splines
    

    subroutine init_aWN(this)
    class(CAMBdata) :: this
    real(dl), parameter :: initialValue = 1.2e-17_dl
    real(dl), parameter :: maxError = 1e-10_dl
    real(dl) :: max_dy, l, r, m, yl, yr, ym
    type(TDarkFieldSplines), allocatable :: sl, sr, sm

    l = 0._dl
    sl = DarkField_splines(this, l)
    yl = sl%dadtau%GetValue(1._dl)
    
    r = initialValue
    do while (.true.)
        if (allocated(sr)) then
            call sr%Clear()
            deallocate(sr)
        end if
        sr = DarkField_splines(this, r)
        yr = sr%dadtau%GetValue(1._dl)
        if (yr > this%dadtau0) exit
        r = 2 * r
    end do
    
    max_dy = maxError * this%dadtau0
    do while(abs(yr - yl) > max_dy)
        m = (l + r)/2
        if ((m <= l) .or. (m >= r)) exit
        sm = DarkField_splines(this, m)
        ym = sm%dadtau%GetValue(1._dl)
        
        if (ym > this%dadtau0) then
            call sr%Clear()
            deallocate(sr)
            r = m
            sr = sm
            yr = ym
        else
            call sl%Clear()
            deallocate(sl)
            l = m
            sl = sm
            yl = ym
        end if
    end do

    call sr%Clear()
    deallocate(sr)

    this%CP%DarkField%W0 = l
    this%CP%DarkField%splines = sl

    end subroutine init_aWN

    subroutine CAMBdata_SetParams(this, P, error, DoReion, call_again)
    !Initialize background variables; does not yet calculate thermal history
    use constants
    class(CAMBdata), target :: this
    type(CAMBparams), intent(in) :: P
    real(dl) fractional_number, conv
    integer, optional :: error !Zero if OK
    logical, optional :: DoReion
    logical, optional :: call_again
    logical WantReion, calling_again
    integer nu_i,actual_massless
    real(dl) nu_massless_degeneracy, neff_i, eta_k, h2
    real(dl) zpeak, sigma_z, zpeakstart, zpeakend
    !Constants in SI units

    global_error_flag = 0

    if (.not. allocated(P%DarkField)) then
        call GlobalError('DarkField not set', error_darkenergy)
    end if

    if (present(error)) error = global_error_flag
    if (global_error_flag/=0) return

    WantReion = DefaultTrue(DoReion)
    calling_again= DefaultFalse(call_again)

    if (calling_again) then
        this%CP%Accuracy = P%Accuracy
        this%CP%Reion%Reionization = P%Reion%Reionization
        this%CP%WantScalars =P%WantScalars
        this%CP%WantCls = P%WantCls
    else
        this%CP=P
        this%CP%Max_eta_k = max(this%CP%Max_eta_k,this%CP%Max_eta_k_tensor)
    end if

    if (.not. calling_again) then
        this%ThermoData%HasThermoData = .false.
        if (this%CP%Num_Nu_Massive /= sum(this%CP%Nu_mass_numbers(1:this%CP%Nu_mass_eigenstates))) then
            if (sum(this%CP%Nu_mass_numbers(1:this%CP%Nu_mass_eigenstates))/=0) &
                call GlobalError('Num_Nu_Massive is not sum of Nu_mass_numbers', error_unsupported_params)
        end if
10      if (this%CP%Omnuh2 < 1.e-7_dl) this%CP%Omnuh2 = 0
        if (this%CP%Omnuh2==0 .and. this%CP%Num_Nu_Massive /=0) then
            this%CP%Num_Nu_Massless = this%CP%Num_Nu_Massless + this%CP%Num_Nu_Massive
            this%CP%Num_Nu_Massive  = 0
            this%CP%Nu_mass_numbers = 0
        end if

        nu_massless_degeneracy = this%CP%Num_Nu_massless !N_eff for massless neutrinos
        if (this%CP%Num_nu_massive > 0) then
            if (this%CP%Nu_mass_eigenstates==0) &
                call GlobalError('Have Num_nu_massive>0 but no nu_mass_eigenstates', error_unsupported_params)
            if (this%CP%Nu_mass_eigenstates==1 .and. this%CP%Nu_mass_numbers(1)==0) &
                this%CP%Nu_mass_numbers(1) = this%CP%Num_Nu_Massive
            if (all(this%CP%Nu_mass_numbers(1:this%CP%Nu_mass_eigenstates)==0)) this%CP%Nu_mass_numbers=1 !just assume one for all

            !default case of equal heating of all neutrinos
            fractional_number = this%CP%Num_Nu_massless + this%CP%Num_Nu_massive
            actual_massless = int(this%CP%Num_Nu_massless + 1e-6_dl)
            neff_i = fractional_number/(actual_massless + this%CP%Num_Nu_massive)
            nu_massless_degeneracy = neff_i*actual_massless
            this%CP%Nu_mass_degeneracies(1:this%CP%Nu_mass_eigenstates) = &
                this%CP%Nu_mass_numbers(1:this%CP%Nu_mass_eigenstates)*neff_i

            if (abs(sum(this%CP%Nu_mass_fractions(1:this%CP%Nu_mass_eigenstates))-1) > 1e-4) &
                call GlobalError('Nu_mass_fractions do not add up to 1', error_unsupported_params)
        else
            this%CP%Nu_mass_eigenstates = 0
        end if

        if (global_error_flag/=0) then
            if (present(error)) error = global_error_flag
            return
        end if


        call This%ThermoData%ScaleFactorAtTime%Clear()

        !  grho gives the contribution to the expansion rate from: (g) photons,
        !  (r) one flavor of relativistic neutrino (2 degrees of freedom),
        !  (m) nonrelativistic matter (for Omega=1).  grho is actually
        !  8*pi*G*rho/c^2 at a=1, with units of Mpc**(-2).
        !  a=tau(Mpc)*adotrad, with a=1 today, assuming 3 neutrinos.
        !  (Used only to set the initial conformal time.)

        !H0 is in km/s/Mpc

        this%grhocrit = 3*this%CP%h0**2/c**2*1000**2 !3*h0^2/c^2 (=8*pi*G*rho_crit/c^2)
        this%dadtau0 = sqrt(this%grhocrit/3)

        this%grhog = kappa/c**2*4*sigma_boltz/c**3*this%CP%tcmb**4*Mpc**2 !8*pi*G/c^2*4*sigma_B/c^3 T^4
        ! grhog=1.4952d-13*tcmb**4
        this%grhor = 7._dl/8*(4._dl/11)**(4._dl/3)*this%grhog !7/8*(4/11)^(4/3)*grhog (per neutrino species)
        !grhor=3.3957d-14*tcmb**4

        !correction for fractional number of neutrinos, e.g. 3.04 to give slightly higher T_nu hence rhor
        !for massive Nu_mass_degeneracies parameters account for heating from grhor

        this%grhornomass=this%grhor*nu_massless_degeneracy
        this%grhormass=0
        do nu_i = 1, this%CP%Nu_mass_eigenstates
            this%grhormass(nu_i)=this%grhor*this%CP%Nu_mass_degeneracies(nu_i)
        end do
        h2 = (this%CP%H0/100)**2
        this%grhoc=this%grhocrit*this%CP%omch2/h2
        this%grhob=this%grhocrit*this%CP%ombh2/h2

        !  adotrad gives da/dtau in the asymptotic radiation-dominated era:
        this%adotrad = sqrt((this%grhog+this%grhornomass+sum(this%grhormass(1:this%CP%Nu_mass_eigenstates)))/3)

        this%Nnow = this%CP%ombh2/h2*(1-this%CP%yhe)*this%grhocrit*c**2/kappa/m_H/Mpc**2

        this%akthom = sigma_thomson*this%Nnow*Mpc
        !sigma_T * (number density of protons now)

        this%fHe = this%CP%YHe/(mass_ratio_He_H*(1.d0-this%CP%YHe))  !n_He_tot / n_H_tot

        this%z_eq = (this%grhob+this%grhoc)/(this%grhog+this%grhornomass+sum(this%grhormass(1:this%CP%Nu_mass_eigenstates))) -1

        if (this%CP%omnuh2/=0) then
            !Initialize things for massive neutrinos
            call ThermalNuBackground%Init()
            call this%NuPerturbations%Init(P%Accuracy%AccuracyBoost*P%Accuracy%neutrino_q_boost)
            !  nu_masses=m_nu(i)*c**2/(k_B*T_nu0)
            do nu_i=1, this%CP%Nu_mass_eigenstates
                this%nu_masses(nu_i)= ThermalNuBackground%find_nu_mass_for_rho(this%CP%omnuh2/h2*this%CP%Nu_mass_fractions(nu_i)&
                    *this%grhocrit/this%grhormass(nu_i))
            end do
            if (all(this%nu_masses(1:this%CP%Nu_mass_eigenstates)==0)) then
                !All density accounted for by massless, so just use massless
                this%CP%Omnuh2 = 0
                goto 10
            end if
            !Just prevent divide by zero
            this%nu_masses(1:this%CP%Nu_mass_eigenstates) = max(this%nu_masses(1:this%CP%Nu_mass_eigenstates),1e-3_dl)
        else
            this%nu_masses = 0
        end if
        call this%CP%DarkField%Init()
        call init_aWN(this)
        if (global_error_flag==0) this%tau0=this%TimeOfz(0._dl)
        if (global_error_flag==0) then
            this%chi0=this%rofChi(this%tau0)
            this%scale= this%chi0/this%tau0  !e.g. change l sampling depending on approx peak spacing
            if (WantReion) call this%CP%Reion%Init(this)
        end if
    end if

    if (global_error_flag/=0) then
        if (present(error)) error = global_error_flag
        return
    end if

    if (present(error)) then
        error = 0
    else if (FeedbackLevel > 0 .and. .not. calling_again) then
        write(*,'("Om_b h^2             = ",f9.6)') P%ombh2
        write(*,'("Om_nu h^2            = ",f9.6)') P%omnuh2
        if (this%CP%Num_Nu_Massive > 0) then
            write(*,'("N_eff (total)        = ",f9.6)') nu_massless_degeneracy + &
                sum(this%CP%Nu_mass_degeneracies(1:this%CP%Nu_mass_eigenstates))
            do nu_i=1, this%CP%Nu_mass_eigenstates
                conv = k_B*(8*this%grhor/this%grhog/7)**0.25*this%CP%tcmb/eV * &
                    (this%CP%nu_mass_degeneracies(nu_i)/this%CP%nu_mass_numbers(nu_i))**0.25 !approx 1.68e-4
                write(*,'(I2, " nu, g=",f7.4," m_nu*c^2/k_B/T_nu0= ",f9.2," (m_nu= ",f6.3," eV)")') &
                    this%CP%nu_mass_numbers(nu_i), this%CP%nu_mass_degeneracies(nu_i), &
                    this%nu_masses(nu_i),conv*this%nu_masses(nu_i)
            end do
        end if
    end if

    end subroutine CAMBdata_SetParams

    subroutine CAMBdata_Free(this)
    class(CAMBdata) :: this

    call Free_ClTransfer(this%CLdata%CTransScal)
    call Free_ClTransfer(this%ClData%CTransVec)
    call Free_ClTransfer(this%ClData%CTransTens)
    call this%MT%Free()
    if (allocated(this%CAMB_Pk)) deallocate(this%CAMB_PK)

    end subroutine CAMBdata_Free

    function CAMBdata_DeltaTime(this, a1,a2, in_tol)
    class(CAMBdata) :: this
    real(dl) CAMBdata_DeltaTime, atol
    real(dl), intent(IN) :: a1,a2
    real(dl), optional, intent(in) :: in_tol

    atol = PresentDefault(tol/1000/exp(this%CP%Accuracy%AccuracyBoost*this%CP%Accuracy%IntTolBoost-1), in_tol)
    CAMBdata_DeltaTime = Integrate_Romberg(this, dtauda,a1,a2,atol)

    end function CAMBdata_DeltaTime

    function CAMBdata_TimeOfz(this, z, tol)
    class(CAMBdata) :: this
    real(dl) CAMBdata_TimeOfz
    real(dl), intent(in), optional :: tol
    real(dl), intent(IN) :: z

    CAMBdata_TimeOfz= this%DeltaTime(0._dl,1._dl/(z+1._dl), tol)
    end function CAMBdata_TimeOfz

    function rofChi(this,Chi) !sinh(chi) for open.
    class(CAMBdata) :: this
    real(dl) Chi,rofChi

    rofChi=chi
    end function rofChi


    function reion_doptdepth_dz(this,z)
    class(CAMBdata) :: this
    real(dl) :: reion_doptdepth_dz
    real(dl), intent(in) :: z

    reion_doptdepth_dz = this%CP%Reion%x_e(z)*this%akthom*dtauda(this,1._dl/(1._dl+z))

    end function reion_doptdepth_dz

    function grho_no_de(this, a) result(grhoa2)
    !  Return 8*pi*G*rho_no_de*a**4 where rho_no_de includes everything except dark energy.
    class(CAMBdata) :: this
    real(dl), intent(in) :: a
    real(dl) grhoa2, rhonu
    integer nu_i

    grhoa2 = (this%grhoc + this%grhob) * a + this%grhog + this%grhornomass

    if (this%CP%Num_Nu_massive /= 0) then
        !Get massive neutrino density relative to massless
        do nu_i = 1, this%CP%nu_mass_eigenstates
            call ThermalNuBack%rho(a * this%nu_masses(nu_i), rhonu)
            grhoa2 = grhoa2 + rhonu * this%grhormass(nu_i)
        end do
    end if

    end function grho_no_de
    
    function GetReionizationOptDepth(this)
    class(CAMBdata) :: this
    real(dl) GetReionizationOptDepth
    integer n
    real(dl) zstart, zend

    call this%CP%Reion%get_timesteps(n, zstart, zend)
    GetReionizationOptDepth = Integrate_Romberg(this, reion_doptdepth_dz,0.d0,zstart,&
        1d-5/this%CP%Accuracy%AccuracyBoost)

    end function GetReionizationOptDepth

    subroutine lSamples_init(this, State, lmin, max_l)
    ! This subroutines initializes lSet%l arrays. Other values will be interpolated.
    class(lSamples) :: this
    class(CAMBdata), target :: State
    integer, intent(IN) :: lmin,max_l
    integer lind, lvar, step, top, bot, lmin_log
    integer, allocatable :: ls(:)
    real(dl) AScale

    allocate(ls(max_l))
    if (allocated(this%l)) deallocate(this%l)
    this%lmin = lmin
    lmin_log = State%CP%min_l_logl_sampling
    associate(Accuracy => State%CP%Accuracy)
        Ascale=State%scale/Accuracy%lSampleBoost

        if (Accuracy%lSampleBoost >=50) then
            !just do all of them
            lind=0
            do lvar=lmin, max_l
                lind=lind+1
                ls(lind)=lvar
            end do
            this%nl=lind
            allocate(this%l(lind))
            this%l = ls(1:lind)
            return
        end if

        lind=0
        do lvar=lmin, 10
            lind=lind+1
            ls(lind)=lvar
        end do

        do lvar=11, 14
            lind=lind+1
            ls(lind)=lvar
        end do
        if (Accuracy%lSampleBoost > 1) then
            do lvar=15, 37, 1
                lind=lind+1
                ls(lind)=lvar
            end do
        else
            do lvar=15, 37, 2
                lind=lind+1
                ls(lind)=lvar
            end do
        end if

        step = max(nint(5*Ascale),2)
        bot=40
        top=bot + step*10

        do lvar=bot, top, step
            lind=lind+1
            ls(lind)=lvar
        end do

        step=max(nint(20*Ascale),4)
        bot=ls(lind)+step
        top=bot+step*2

        do lvar = bot,top,step
            lind=lind+1
            ls(lind)=lvar
        end do

        if (ls(lind)>=max_l) then
            do lvar=lind,1,-1
                if (ls(lvar)<=max_l) exit
            end do
            lind=lvar
            if (ls(lind)<max_l) then
                lind=lind+1
                ls(lind)=max_l
            end if
        else
            step=max(nint(25*Ascale),4)
            !Get EE right around l=200 by putting extra point at 175
            bot=ls(lind)+step
            top=bot+step

            do lvar = bot,top,step
                lind=lind+1
                ls(lind)=lvar
            end do

            if (ls(lind)>=max_l) then
                do lvar=lind,1,-1
                    if (ls(lvar)<=max_l) exit
                end do
                lind=lvar
                if (ls(lind)<max_l) then
                    lind=lind+1
                    ls(lind)=max_l
                end if
            else
                if (.not. this%use_spline_template) then
                    step=max(nint(42*Ascale),7)
                else
                    step=max(nint(50*Ascale),7)
                end if
                bot=ls(lind)+step
                top=min(lmin_log,max_l)

                do lvar = bot,top,step
                    lind=lind+1
                    ls(lind)=lvar
                end do

                if (max_l > lmin_log) then
                    !Should be pretty smooth or tiny out here
                    step=max(nint(400*Ascale),50)
                    lvar = ls(lind)
                    do
                        lvar = lvar + step
                        if (lvar > max_l) exit
                        lind=lind+1
                        ls(lind)=lvar
                        step = nint(step*1.5) !log spacing
                    end do
                    if (ls(lind) < max_l - 100) then
                        !Try to keep lensed spectra up to specified lmax
                        lind=lind+1
                        ls(lind)=max_l - lensed_convolution_margin
                    else if (ls(lind) - ls(lind-1) > lensed_convolution_margin) then
                        ls(lind)=max_l - lensed_convolution_margin
                    end if
                end if
            end if
            if (ls(lind) /=max_l) then
                lind=lind+1
                ls(lind)=max_l
            end if
            !Not in _flat case so interpolation table is the same when using lower l_max
        end if
    end associate
    this%nl=lind
    allocate(this%l, source=ls(1:lind))

    end subroutine lSamples_init

    subroutine InterpolateClArr(lSet,iCl, all_Cl, max_index)
    class(lSamples), intent(in) :: lSet
    real(dl), intent(in) :: iCl(1:*)
    real(dl), intent(out):: all_Cl(lSet%lmin:*)
    integer, intent(in), optional :: max_index
    integer il,llo,lhi, xi
    real(dl) ddCl(lSet%nl)
    real(dl) xl(lSet%nl)
    real(dl) a0,b0,ho
    integer max_ind

    max_ind = PresentDefault(lSet%nl, max_index)

    if (max_ind > lSet%nl) call MpiStop('Wrong max_ind in InterpolateClArr')

    xl = real(lSet%l(1:lSet%nl),dl)
    call spline_def(xl,iCL,max_ind,ddCl)

    llo=1
    do il=lSet%lmin,lSet%l(max_ind)
        xi=il
        if ((xi > lSet%l(llo+1)).and.(llo < max_ind)) then
            llo=llo+1
        end if
        lhi=llo+1
        ho=lSet%l(lhi)-lSet%l(llo)
        a0=(lSet%l(lhi)-xi)/ho
        b0=(xi-lSet%l(llo))/ho

        all_Cl(il) = a0*iCl(llo)+ b0*iCl(lhi)+((a0**3-a0)* ddCl(llo) &
            +(b0**3-b0)*ddCl(lhi))*ho**2/6
    end do

    end subroutine InterpolateClArr

    subroutine InterpolateClArrTemplated(lSet,iCl, all_Cl, max_ind, template_index)
    class(lSamples), intent(in) :: lSet
    real(dl), intent(in) :: iCl(*)
    real(dl), intent(out):: all_Cl(lSet%lmin:*)
    integer, intent(in) :: max_ind
    integer, intent(in), optional :: template_index
    integer maxdelta, il
    real(dl) DeltaCL(lSet%nl)
    real(dl), allocatable :: tmpall(:)

    if (max_ind > lSet%nl) call MpiStop('Wrong max_ind in InterpolateClArrTemplated')
    if (lSet%use_spline_template .and. present(template_index)) then
        if (template_index<=3) then
            !interpolate only the difference between the C_l and an accurately interpolated template.
            !Using unlensed for template, seems to be good enough
            maxdelta=max_ind
            do while (lSet%l(maxdelta) > lmax_extrap_highl)
                maxdelta=maxdelta-1
            end do
            DeltaCL(1:maxdelta)=iCL(1:maxdelta)- highL_CL_template(lSet%l(1:maxdelta), template_index)

            call lSet%InterpolateClArr(DeltaCl, all_Cl, maxdelta)

            do il=lSet%lmin,lSet%l(maxdelta)
                all_Cl(il) = all_Cl(il) +  highL_CL_template(il,template_index)
            end do

            if (maxdelta < max_ind) then
                !directly interpolate high L where no t  emplate (doesn't effect lensing spectrum much anyway)
                allocate(tmpall(lSet%lmin:lSet%l(max_ind)))
                call InterpolateClArr(lSet,iCl, tmpall, max_ind)
                !overlap to reduce interpolation artefacts
                all_cl(lSet%l(maxdelta-2):lSet%l(max_ind) ) = tmpall(lSet%l(maxdelta-2):lSet%l(max_ind))
                deallocate(tmpall)
            end if
            return
        end if
    end if

    call InterpolateClArr(lSet,iCl, all_Cl, max_ind)

    end subroutine InterpolateClArrTemplated


    subroutine Thermo_values(this,tau, a, cs2b, opacity, dopacity)
    !Compute unperturbed sound speed squared,
    !and ionization fraction by interpolating pre-computed tables.
    !If requested also get time derivative of opacity
    class(TThermoData) :: this
    real(dl), intent(in) :: tau
    real(dl), intent(out) :: a, cs2b, opacity
    real(dl), intent(out), optional :: dopacity
    integer i
    real(dl) d

    d=log(tau/this%tauminn)/this%dlntau+1._dl
    i=int(d)
    d=d-i
    if (i < 1) then
        !Linear interpolation if out of bounds (should not occur).
        write(*,*) 'tau, taumin = ', tau, this%tauminn
        call MpiStop('thermo out of bounds')
    else if (i >= this%nthermo) then
        cs2b=this%cs2(this%nthermo)
        opacity=this%dotmu(this%nthermo)
        a=1
        if (present(dopacity)) then
            dopacity = this%ddotmu(this%nthermo)/(tau*this%dlntau)
        end if
    else
        cs2b=this%cs2(i)+d*(this%dcs2(i)+d*(3*(this%cs2(i+1)-this%cs2(i))  &
            -2*this%dcs2(i)-this%dcs2(i+1)+d*(this%dcs2(i)+this%dcs2(i+1)  &
            +2*(this%cs2(i)-this%cs2(i+1)))))
        opacity=this%dotmu(i)+d*(this%ddotmu(i)+d*(3*(this%dotmu(i+1)-this%dotmu(i)) &
            -2*this%ddotmu(i)-this%ddotmu(i+1)+d*(this%ddotmu(i)+this%ddotmu(i+1) &
            +2*(this%dotmu(i)-this%dotmu(i+1)))))
        a = (this%ScaleFactor(i)+d*(this%dScaleFactor(i)+d*(3*(this%ScaleFactor(i+1)-this%ScaleFactor(i)) &
            -2*this%dScaleFactor(i)-this%dScaleFactor(i+1)+d*(this%dScaleFactor(i)+this%dScaleFactor(i+1) &
            +2*(this%ScaleFactor(i)-this%ScaleFactor(i+1))))))*tau
        if (present(dopacity)) then
            dopacity=(this%ddotmu(i)+d*(this%dddotmu(i)+d*(3*(this%ddotmu(i+1)  &
                -this%ddotmu(i))-2*this%dddotmu(i)-this%dddotmu(i+1)+d*(this%dddotmu(i) &
                +this%dddotmu(i+1)+2*(this%ddotmu(i)-this%ddotmu(i+1))))))/(tau*this%dlntau)
        end if
    end if
    end subroutine Thermo_values

    function Thermo_OpacityToTime(this,opacity)
    class(TThermoData) :: this
    real(dl), intent(in) :: opacity
    integer j
    real(dl) Thermo_OpacityToTime
    !Do this the bad slow way for now..
    !The answer is approximate
    j =1
    do while(this%dotmu(j)> opacity)
        j=j+1
    end do

    Thermo_OpacityToTime = exp((j-1)*this%dlntau)*this%tauminn

    end function Thermo_OpacityToTime

    subroutine Thermo_Init(this, State,taumin)
    !  Compute and save unperturbed baryon temperature and ionization fraction
    !  as a function of time.  With nthermo=10000, xe(tau) has a relative
    ! accuracy (numerical integration precision) better than 1.e-5.
    use constants
    use StringUtils
    class(TThermoData) :: this
    class(CAMBdata), target :: State
    real(dl), intent(in) :: taumin
    integer nthermo
    real(dl) tau01,a0,barssc,dtau
    real(dl) tau,a,a2
    real(dl) adot,fe,thomc0
    real(dl) dtbdla,vfi,cf1,maxvis, vis, z_scale
    integer ncount,i,j1,iv,ns
    real(dl), allocatable :: spline_data(:)
    real(dl) last_dotmu, om
    real(dl) a_verydom
    real(dl) awin_lens1p,awin_lens2p,dwing_lens, rs, DA
    real(dl) a_eq, rs_eq, tau_eq, rstar
    integer noutput
    Type(CalWins), dimension(:), allocatable, target :: RW
    real(dl) Tspin, Trad, rho_fac, window, tau_eps
    integer RW_i, j2
    real(dl) Tb21cm, winamp, z, background_boost
    character(len=:), allocatable :: outstr
    real(dl), allocatable ::  taus(:)
    real(dl), allocatable :: xe_a(:), sdotmu(:), opts(:)
    real(dl), allocatable :: scale_factors(:), times(:), dt(:)
    Type(TCubicSpline) :: dotmuSp
    integer ninverse, nlin
    real(dl) dlna, zstar_min, zstar_max
    real(dl) reion_z_start, reion_z_complete
    Type(CAMBParams), pointer :: CP

    CP => State%CP

    !Allocate memory outside parallel region to keep ifort happy
    background_boost = CP%Accuracy%BackgroundTimeStepBoost*CP%Accuracy%AccuracyBoost
    if (background_boost > 20) then
        write(*,*) 'Warning: very small time steps can give less accurate spline derivatives'
        write(*,*) 'e.g. around reionization if not matched very smoothly'
    end if
    !Higher q starts earlier; scale by log(taumin) so actual step size is not made worse by increasing k_max
    nthermo = nint(thermal_history_def_timesteps*log(1.4e4/taumin)/log(1.4e4/2e-4)*background_boost)
    this%tauminn=0.95d0*taumin
    this%dlntau=log(State%tau0/this%tauminn)/(nthermo-1)

    this%nthermo = nthermo
    allocate(spline_data(nthermo), sdotmu(nthermo))

    if (allocated(this%tb) .and. this%nthermo/=size(this%tb)) then
        deallocate(this%scaleFactor, this%cs2, this%dcs2, this%ddotmu)
        deallocate(this%dscaleFactor, this%adot, this%dadot)
        deallocate(this%tb, this%xe, this%emmu, this%dotmu)
        deallocate(this%demmu, this%dddotmu, this%ddddotmu)
        if (dowinlens .and. allocated(this%winlens)) deallocate(this%winlens, this%dwinlens)
    endif
    if (.not. allocated(this%tb)) then
        allocate(this%scaleFactor(nthermo), this%cs2(nthermo), this%dcs2(nthermo), this%ddotmu(nthermo))
        allocate(this%dscaleFactor(nthermo), this%adot(nthermo), this%dadot(nthermo))
        allocate(this%tb(nthermo), this%xe(nthermo), this%emmu(nthermo),this%dotmu(nthermo))
        allocate(this%demmu(nthermo), this%dddotmu(nthermo), this%ddddotmu(nthermo))
        if (dowinlens) allocate(this%winlens(nthermo), this%dwinlens(nthermo))
    end if

    om = (State%grhob+State%grhoc)/&
        sqrt(3*(State%grhog+sum(State%grhormass(1:CP%Nu_mass_eigenstates))+State%grhornomass))
    a0=this%tauminn*State%adotrad*(1+om*this%tauminn/4)
    ninverse = nint(background_boost*log(1/a0)/log(1/2d-10)*4000)

    nlin = ninverse/2
    allocate(scale_factors(ninverse+nlin))
    allocate(times(ninverse+nlin))
    allocate(dt(ninverse+nlin))
    allocate(taus(nthermo), xe_a(nthermo))

    !$OMP PARALLEL SECTIONS DEFAULT(SHARED)
    !$OMP SECTION
    call CP%Recomb%Init(State,WantTSpin=.false.)    !almost all the time spent here

    !$OMP SECTION

    call splini(spline_data,nthermo)

    this%tight_tau = 0
    this%actual_opt_depth = 0
    ncount=0
    this%z_drag=0.d0
    thomc0= Compton_CT * CP%tcmb**4
    this%r_drag0 = 3.d0/4.d0*State%grhob/State%grhog
    last_dotmu = 0

    this%matter_verydom_tau = 0
    a_verydom = CP%Accuracy%AccuracyBoost*5*(State%grhog+State%grhornomass)/(State%grhoc+State%grhob)
    if (CP%Reion%Reionization) then
        call CP%Reion%get_timesteps(State%reion_n_steps, reion_z_start, reion_z_complete)
        State%reion_tau_start = max(0.05_dl, State%TimeOfZ(reion_z_start, 1d-3))
        !Time when a very small reionization fraction (assuming tanh fitting)
        State%reion_tau_complete = min(State%tau0, &
            State%reion_tau_start+ State%DeltaTime(1/(1+reion_z_start),1/(1.d0+reion_z_complete),1d-3))
    else
        State%reion_tau_start = State%tau0
        State%reion_tau_complete = State%tau0
    end  if
    !  Initial conditions: assume radiation-dominated universe.
    !  Assume that any entropy generation occurs before tauminn.
    !  This gives wrong temperature before pair annihilation, but
    !  the error is harmless.

    !Get scale factor as function of time by inverting tau(a)
    dlna = log(0.2_dl/a0)/(ninverse-1)
    do i=2, ninverse-1
        scale_factors(1+i) = a0*exp((i-1)*dlna)
    end do
    scale_factors(1) = a0
    scale_factors(2) = a0*exp(dlna/3)
    da = 0.8_dl/(nlin-2)
    do i=1, nlin-2
        scale_factors(ninverse+i) = 0.2_dl + (i-1)*da
    end do
    scale_factors(ninverse+nlin-1) = 0.9_dl + 0.1_dl*scale_factors(ninverse+nlin-2)
    scale_factors(ninverse+nlin) = 1
    do i=1, ninverse+nlin
        dt(i) = dtauda(State,scale_factors(i))
    end do
    call this%ScaleFactorAtTime%Init(scale_factors, dt)
    call this%ScaleFactorATTime%IntegralArray(times(2), first_index=2)
    times(1) = this%tauminn
    times(2:) = times(2:) + 2*(sqrt(1 + om*scale_factors(2)/ State%adotrad) -1)/om
    times(ninverse+nlin) = State%tau0
    call This%ScaleFactorAtTime%Init(times, scale_factors)
    taus(1) = this%tauminn
    do i=2,nthermo-1
        taus(i) = this%tauminn*exp((i-1)*this%dlntau)
    end do
    taus(nthermo) = State%tau0
    call this%ScaleFactorAtTime%Array(taus(2:), this%scaleFactor(2:))
    this%scaleFactor(1) = a0
    this%scaleFactor(nthermo) = 1
    this%adot(1) = 1/dtauda(State,a0)

    tau01=this%tauminn
    do i=2,nthermo
        !Get recombination-independent parts of background now as function of conformal time tau
        !Now at the log spaced time steps
        tau=taus(i)
        dtau = tau-tau01
        a = this%scaleFactor(i)
        adot = 1/dtauda(State,a)
        this%adot(i) = adot
        if (this%matter_verydom_tau ==0 .and. a > a_verydom) then
            this%matter_verydom_tau = tau
        end if
        z= 1._dl/a-1._dl
        tau01 =tau
    end do
    !$OMP END PARALLEL SECTIONS

    if (global_error_flag/=0) return

    call CP%Recomb%xe_tm(a0,this%xe(1), this%tb(1))
    barssc=barssc0*(1._dl-0.75d0*CP%yhe+(1._dl-CP%yhe)*this%xe(1))
    this%cs2(1)=4._dl/3._dl*barssc*this%tb(1)
    this%dotmu(1)=this%xe(1)*State%akthom/a0**2


    !$OMP PARALLEL DO DEFAULT(SHARED), SCHEDULE(STATIC,16)
    do i=2,nthermo
        call CP%Recomb%xe_tm(this%scaleFactor(i), xe_a(i), this%tb(i))
    end do

    do i=2,nthermo
        tau =taus(i)
        a = this%scaleFactor(i)
        a2=a*a
        adot=this%adot(i)

        ! If there is re-ionization, smoothly increase xe to the
        ! requested value.
        if (CP%Reion%Reionization .and. tau > State%reion_tau_start) then
            if(ncount == 0) then
                ncount=i-1
            end if
            this%xe(i) = CP%Reion%x_e(1/a-1, tau, this%xe(ncount))
        else
            this%xe(i)=xe_a(i)
        end if

        !  approximate Baryon sound speed squared (over c**2).
        fe=(1._dl-CP%yhe)*this%xe(i)/(1._dl-0.75d0*CP%yhe+(1._dl-CP%yhe)*this%xe(i))
        dtbdla=-2._dl*this%tb(i)
        if (a*this%tb(i)-CP%tcmb < -1e-8) then
            dtbdla= dtbdla -thomc0*fe/adot*(a*this%tb(i)-CP%tcmb)/a**3
        end if
        barssc=barssc0*(1._dl-0.75d0*CP%yhe+(1._dl-CP%yhe)*this%xe(i))
        this%cs2(i)=barssc*this%tb(i)*(1-dtbdla/this%tb(i)/3._dl)

        ! Calculation of the visibility function
        this%dotmu(i)=this%xe(i)*State%akthom/a2

        if (this%tight_tau==0 .and. 1/(tau*this%dotmu(i)) > 0.005) this%tight_tau = tau !0.005
        !Tight coupling switch time when k/opacity is smaller than 1/(tau*opacity)
    end do

    if (CP%Reion%Reionization .and. (this%xe(nthermo) < 0.999d0)) then
        write(*,*)'Warning: xe at redshift zero is < 1'
        write(*,*) 'Check input parameters an Reionization_xe'
        write(*,*) 'function in the Reionization module'
    end if

    !Integrate for optical depth
    call dotmuSp%Init(taus(nthermo:1:-1), this%dotmu(nthermo:1:-1))
    allocate(opts(nthermo))
    call dotmuSp%IntegralArray(opts)
    sdotmu = opts(nthermo:1:-1)
    do j1=1,nthermo
        if (sdotmu(j1)< -69) then
            this%emmu(j1)=1.d-30
        else
            this%emmu(j1)=exp(sdotmu(j1))
        end if
    end do
    z_scale =  COBE_CMBTemp/CP%TCMB
    zstar_min = 700._dl * z_scale
    zstar_max = 2000._dl * z_scale

    iv=0
    vfi=0._dl
    ! Getting the starting and finishing times for decoupling and time of maximum visibility
    if (ncount == 0) then
        cf1=1._dl
        ns=nthermo
    else
        cf1=exp(-sdotmu(ncount))
        ns=ncount
    end if
    maxvis = 0
    do j1=1,ns
        vis = this%emmu(j1)*this%dotmu(j1)
        tau = taus(j1)
        vfi=vfi+vis*cf1*this%dlntau*tau
        if ((iv == 0).and.(vfi > 1.0d-7/CP%Accuracy%AccuracyBoost)) then
            State%taurst=9._dl/10._dl*tau
            iv=1
        elseif (iv == 1) then
            if (vis > maxvis) then
                maxvis=vis
                State%tau_maxvis = tau
            end if
            if (vfi > 0.995) then
                State%taurend=tau
                iv=2
                exit
            end if
        end if
    end do

    if (iv /= 2) then
        call GlobalError('ThemoData Init: failed to find end of recombination',error_reionization)
        return
    end if

    if (dowinlens) then
        vfi=0
        awin_lens1p=0
        awin_lens2p=0
        this%winlens=0
        do j1=1,nthermo-1
            vis = this%emmu(j1)*this%dotmu(j1)
            tau = this%tauminn* taus(j1)
            vfi=vfi+vis*cf1*this%dlntau*tau
            if (vfi < 0.995) then
                dwing_lens =  vis*cf1*this%dlntau*tau / 0.995

                awin_lens1p = awin_lens1p + dwing_lens
                awin_lens2p = awin_lens2p + dwing_lens/(State%tau0-tau)
            end if
            this%winlens(j1)= awin_lens1p/(State%tau0-tau) - awin_lens2p
        end do
    end if

    ! Calculating the timesteps during recombination.

    State%dtaurec=min(State%dtaurec,State%taurst/40)/CP%Accuracy%AccuracyBoost
    if (do_bispectrum .and. hard_bispectrum) State%dtaurec = State%dtaurec / 4

    if (CP%Reion%Reionization) State%taurend=min(State%taurend,State%reion_tau_start)

    if (DebugMsgs) then
        write (*,*) 'taurst, taurend = ', State%taurst, State%taurend
    end if

    !$OMP PARALLEL SECTIONS DEFAULT(SHARED)
    !$OMP SECTION
    call splder(this%dotmu,this%ddotmu,nthermo,spline_data)
    call splder(this%ddotmu,this%dddotmu,nthermo,spline_data)
    call splder(this%dddotmu,this%ddddotmu,nthermo,spline_data)
    !$OMP SECTION
    call splder(this%cs2,this%dcs2,nthermo,spline_data)
    call splder(this%emmu,this%demmu,nthermo,spline_data)
    call splder(this%adot,this%dadot,nthermo,spline_data)
    if (dowinlens) call splder(this%winlens,this%dwinlens,nthermo,spline_data)
    !$OMP SECTION
    this%ScaleFactor(:) = this%scaleFactor/taus !a/tau
    this%dScaleFactor(:) = (this%adot - this%ScaleFactor)*this%dlntau !derivative of a/tau
    call this%SetTimeSteps(State,State%TimeSteps)
    !$OMP END PARALLEL SECTIONS

    this%HasThermoData = .true.
    end subroutine Thermo_Init


    subroutine SetTimeSteps(this,State,TimeSteps)
    !Set time steps to use for sampling the source functions for the CMB power spectra
    class(TThermoData) :: this
    Type(TRanges) :: TimeSteps
    class(CAMBdata) State
    real(dl) dtau0
    integer nri0, nstep
    !Sources
    integer ix,i,nwindow, L_limb
    real(dl) keff, win_end, TimeSampleBoost, delta

    TimeSampleBoost = State%CP%Accuracy%AccuracyBoost*State%CP%Accuracy%TimeStepBoost
    call TimeSteps%Init()

    call TimeSteps%Add_delta(State%taurst, State%taurend, State%dtaurec)

    ! Calculating the timesteps after recombination
    dtau0=State%tau0/500._dl/TimeSampleBoost
    if (do_bispectrum) dtau0 = dtau0/3
    !Don't need this since adding in Limber on small scales

    call TimeSteps%Add_delta(State%taurend, State%tau0, dtau0)

    !Sources

    this%tau_start_redshiftwindows = State%tau0
    this%tau_end_redshiftwindows = 0


    if (State%CP%Reion%Reionization) then
        nri0=int(State%reion_n_steps*State%CP%Accuracy%AccuracyBoost)
        !Steps while reionization going from zero to maximum
        call TimeSteps%Add(State%reion_tau_start,State%reion_tau_complete,nri0)
    end if

    if (global_error_flag/=0) then
        return
    end if

    !Create arrays out of the region information.
    call TimeSteps%GetArray()
    nstep = TimeSteps%npoints

    if (DebugMsgs .and. FeedbackLevel > 0) call WriteFormat('Set %d time steps', nstep)

    end subroutine SetTimeSteps

    subroutine IonizationFunctionsAtTime(this,tau, a, opac, dopac, ddopac, &
        vis, dvis, ddvis, expmmu, lenswin)
    class(TThermoData) :: this
    real(dl), intent(in) :: tau
    real(dl), intent(out):: a, opac, dopac, ddopac, vis, dvis, ddvis, expmmu, lenswin
    real(dl) d, cs2
    integer i

    call this%Values(tau,a,cs2,opac,dopac)

    d=log(tau/this%tauminn)/this%dlntau+1._dl
    i=int(d)
    d=d-i

    if (i < this%nthermo) then
        ddopac=(this%dddotmu(i)+d*(this%ddddotmu(i)+d*(3._dl*(this%dddotmu(i+1) &
            -this%dddotmu(i))-2._dl*this%ddddotmu(i)-this%ddddotmu(i+1)  &
            +d*(this%ddddotmu(i)+this%ddddotmu(i+1)+2._dl*(this%dddotmu(i) &
            -this%dddotmu(i+1)))))-(this%dlntau**2)*tau*dopac) &
            /(tau*this%dlntau)**2
        expmmu=this%emmu(i)+d*(this%demmu(i)+d*(3._dl*(this%emmu(i+1)-this%emmu(i)) &
            -2._dl*this%demmu(i)-this%demmu(i+1)+d*(this%demmu(i)+this%demmu(i+1) &
            +2._dl*(this%emmu(i)-this%emmu(i+1)))))

        if (dowinlens) then
            lenswin=this%winlens(i)+d*(this%dwinlens(i)+d*(3._dl*(this%winlens(i+1)-this%winlens(i)) &
                -2._dl*this%dwinlens(i)-this%dwinlens(i+1)+d*(this%dwinlens(i)+this%dwinlens(i+1) &
                +2._dl*(this%winlens(i)-this%winlens(i+1)))))
        end if
        vis=opac*expmmu
        dvis=expmmu*(opac**2+dopac)
        ddvis=expmmu*(opac**3+3*opac*dopac+ddopac)
    else
        ddopac=this%dddotmu(this%nthermo)
        expmmu=this%emmu(this%nthermo)
        vis=opac*expmmu
        dvis=expmmu*(opac**2+dopac)
        ddvis=expmmu*(opac**3+3._dl*opac*dopac+ddopac)
    end if

    end subroutine IonizationFunctionsAtTime

    subroutine Init_ClTransfer(CTrans)
    !Need to set the Ranges array q before calling this
    Type(ClTransferData) :: CTrans
    integer st

    deallocate(CTrans%Delta_p_l_k, STAT = st)
    call CTrans%q%getArray(.true.)

    allocate(CTrans%Delta_p_l_k(CTrans%NumSources,&
        min(CTrans%max_index_nonlimber,CTrans%ls%nl), CTrans%q%npoints), STAT = st)
    if (st /= 0) call MpiStop('Init_ClTransfer: Error allocating memory for transfer functions')
    CTrans%Delta_p_l_k = 0

    end subroutine Init_ClTransfer

    subroutine Init_Limber(CTrans,State)
    Type(ClTransferData) :: CTrans
    class(CAMBdata) :: State

    call Free_Limber(Ctrans)
    allocate(CTrans%Limber_l_min(CTrans%NumSources))
    CTrans%Limber_l_min = 0

    end subroutine Init_Limber

    subroutine Free_ClTransfer(CTrans)
    Type(ClTransferData) :: CTrans

    if (allocated(CTrans%Delta_p_l_k)) deallocate(CTrans%Delta_p_l_k)
    call CTrans%q%Free()
    call Free_Limber(CTrans)

    end subroutine Free_ClTransfer

    subroutine Free_Limber(CTrans)
    Type(ClTransferData) :: CTrans

    if (allocated(CTrans%Limber_l_min)) deallocate(CTrans%Limber_l_min)
    if (allocated(CTrans%Limber_windows)) deallocate(CTrans%Limber_windows)

    end subroutine Free_Limber


    subroutine TCLdata_InitCls(this, State)
    class(TCLData) :: this
    class(CAMBdata) :: State

    associate(CP=>State%CP)
        call CheckLoadedHighLTemplate
        if (CP%WantScalars) then
            if (allocated(this%Cl_scalar)) deallocate(this%Cl_scalar)
            allocate(this%Cl_scalar(CP%Min_l:CP%Max_l, C_Temp:State%Scalar_C_last), source=0._dl)
        end if
    end associate

    end subroutine TCLdata_InitCls

    end module results
